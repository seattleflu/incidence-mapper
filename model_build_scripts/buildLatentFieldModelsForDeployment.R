# buildModelsForDeployment.R
# this script (and similar others?) controls standardized database queries and model training for web deployment

library(dbViewR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(incidenceMapR)
library(modelServR)
library(modelVisualizeR)


SRC <- 'production'

db <- selectFromDB(queryIn= list(SELECT  =c("*")), source = SRC)
print(mostRecentSample<-max(db$observedData$encountered_date))


fluPathogens <- c('Flu_A_H1','Flu_A_H3','Flu_A_pan','Flu_B_pan','Flu_C_pan')

pathogens <- db$observedData %>% group_by(pathogen) %>% summarize(n = n()) %>% arrange(desc(n))

tmp<-pathogens$pathogen[pathogens$n >= 500 | grepl('Flu',pathogens$pathogen)]
tmp2<-as.list(tmp)
names(tmp2)<-tmp


pathogenKeys <- list(
                     flu=fluPathogens#,
                     # Flu_A_H1 = 'Flu_A_H1', Flu_A_H3 = 'Flu_A_H3', Flu_B_pan = 'Flu_B_pan', Flu_C_pan = 'Flu_C_pan'#,
                     # all='all', other_non_flu = setdiff(pathogens$pathogen,c(fluPathogens,'not_yet_tested','measles','Measles'))#,
                     # rsv = c('RSVA','RSVB'), RSVA='RSVA', RSVB='RSVB', 
                     # AdV='AdV',CoV='CoV',RV='RV'
                     )


factors   <- c('site_type','age_range_coarse_upper','sfs_year')


SOURCE='sfs_domain_geojson'
GEO='residence_regional_name'


siteTypes <- c('childrensHospital','clinic','collegeCampus','childrensClinic','port','retrospective','workplace','publicSpace','self-test')


# 1 week ahead of current week
nowcastWeek <- 1+(isoweek(Sys.time()) %% 52)
nowcastYear <- isoyear(Sys.time())
if (nowcastWeek < isoweek(Sys.time())){
  nowcastYear <- nowcastYear+1
}
nowcastWeek <- paste(nowcastYear ,'-W',sprintf("%02d",nowcastWeek),sep='')


#####################################
###### timeseries latent field models ############
#####################################

# number of subjects with pathogen and factor at residence location 

shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = SOURCE)

for (PATHOGEN in names(pathogenKeys)){

  queryIn <- list(
    SELECT   =list(COLUMN=c('pathogen', factors, GEO,'encountered_week')),
    WHERE    =list(COLUMN='pathogen', IN=pathogenKeys[[PATHOGEN]]),
    WHERE    =list(COLUMN='site_type', IN=siteTypes),
    GROUP_BY =list(COLUMN=c(factors,GEO,"encountered_week")),
    SUMMARIZE=list(COLUMN='pathogen', IN= pathogenKeys[[PATHOGEN]])
  )
  db <- selectFromDB(  queryIn, source=SRC, na.rm=TRUE )
  db <- expandDB(db, shp=shp, currentWeek =nowcastWeek)
  
        
  db <- appendCatchmentModel(db,shp=shp, source=SRC, na.rm=TRUE )

  modelDefinition <- latentFieldModel(db=db, shp=shp)
  model <- modelTrainR(modelDefinition)
  
  print(summary(model$inla))
  
  
  ## rescale latent field model for consistent color code
  reference_scale <- 0.201*.12/.1
  
  
  current_scale <- model$latentField %>% select(encountered_week,residence_regional_name,modeled_intensity_mode) %>%
    group_by(encountered_week) %>% summarize(region_intensity_mode = mean(modeled_intensity_mode))
  current_scale <- max(current_scale$region_intensity_mode)
  
  scale_factor <- reference_scale/current_scale
  
  for (COLUMN in names(model$latentField)[grepl('modeled_',names(model$latentField))]){
    model$latentField[[COLUMN]] <- model$latentField[[COLUMN]] * scale_factor
  }
  
  dir.create('/home/rstudio/seattle_flu/dev_model_diagnostic_plots/', showWarnings = FALSE)
  fname <- paste('/home/rstudio/seattle_flu/dev_model_diagnostic_plots/',paste('inla_latent',PATHOGEN,SOURCE,GEO,'encountered_week',Sys.Date(),sep='-'),'.png',sep='')
  png(filename = fname,width = 6, height = 5, units = "in", res = 300)
  print(
    ggplot(model$latentField) + 
          geom_line(aes_string(x='encountered_week',y="modeled_intensity_median", color=GEO,group =GEO)) + 
          # geom_ribbon(aes_string(x='encountered_week',ymin="modeled_intensity_lower_95_CI", ymax="modeled_intensity_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
          guides(color=FALSE, fill=FALSE) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    )
  dev.off()

  
  dir.create('/home/rstudio/seattle_flu/dev_model_diagnostic_plots/', showWarnings = FALSE)
  fname <- paste('/home/rstudio/seattle_flu/dev_model_diagnostic_plots/',paste('inla_latent',PATHOGEN,SOURCE,GEO,'age',Sys.Date(),sep='-'),'.png',sep='')
  png(filename = fname,width = 6, height = 5, units = "in", res = 300)
  if ('site_age_siteIdx' %in% names(model$inla$summary.random)){
    plotDat <- model$inla$summary.random$site_age_siteIdx
  } else {
    plotDat <- model$inla$summary.random$age_row_iid
  }
  plotDat$group <- unique(model$modeledData$site_type)
  plotDat$color <- unique(model$modeledData$sfs_year)
  plotDat$age <- rep(unique(model$modeledData$age_range_coarse_upper), each = length(unique(model$modeledData$site_type)))
  N<-length(unique(model$modeledData$age_range_coarse_upper))
  plotDat$ageID <- rep(1:N, each = length(unique(model$modeledData$site_type)))
  print(
    ggplot(plotDat) + geom_point(aes(x=ageID, y=mean, color=color), stat='identity') + theme_bw() +
      geom_linerange(aes(x=ageID, ymin=`0.025quant`, ymax=`0.975quant`, color=color)) +
      xlab('age_range_coarse_upper') + ylab('age random effect') + facet_wrap("group")+
      scale_x_continuous(breaks=1:N,labels=as.character(unique(plotDat$age)))
  )
  dev.off()
  
  saveModel(model, storeRDS=FALSE)
  # saveModel(model, storeRDS=TRUE)
  
}

