# buildModelsForDeployment.R
# this script (and similar others?) controls standardized database queries and model training for web deployment

library(dbViewR)
library(incidenceMapR)
library(modelServR)
library(modelVisualizeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

SRC <- 'production'
# SRC <- 'simulated_data'

db <- selectFromDB(queryIn= list(SELECT  =c("*")), source = SRC)

pathogens <- db$observedData %>% group_by(pathogen) %>% summarize(n = n()) %>% arrange(desc(n))

tmp<-pathogens$pathogen[pathogens$n >= 500 | grepl('Flu',pathogens$pathogen)]
tmp2<-as.list(tmp)
names(tmp2)<-tmp


fluPathogens <- c('Flu_A_H1','Flu_A_H3','Flu_A_pan','Flu_B_pan','Flu_C_pan')

# pathogenKeys <- c(list(all='all',
#                   flu=fluPathogens,
#                   rsv=c('RSVA','RSVB'),
#                   other_non_flu = setdiff(pathogens$pathogen,c(fluPathogens,'not_yet_tested','measles'))
#                   ),
#                   tmp2)

pathogenKeys <- list(all='all', flu=fluPathogens, other_non_flu = setdiff(pathogens$pathogen,c(fluPathogens,'not_yet_tested','measles','Measles')))


factors   <- c('site_type','sex','flu_shot')#,'age_range_fine_upper')
# factors   <- c('site','sex','flu_shot')#,'age_range_fine_upper') # site overwhelms ram....


geoLevels <- list(
                   seattle_geojson = c('residence_puma','residence_neighborhood_district_name','residence_cra_name'), #,'residence_census_tract'),
                   wa_geojson = c('residence_puma')
                 )


currentWeek <- paste(isoyear(Sys.time()) ,'-W',isoweek(Sys.time()),sep='')


#####################################
###### timeseries latent field models ############
#####################################

# number of subjects with pathogen and factor at residence location 
for (SOURCE in names(geoLevels)){
  for (GEO in geoLevels[[SOURCE]]){
    
    # SOURCE='seattle_geojson'
    # SOURCE='wa_geojson'
    # GEO='residence_census_tract'
    # GEO='residence_puma'
    # GEO='residence_cra_name'
    # PATHOGEN='flu'
    # PATHOGEN='all'
    # PATHOGEN='rsv'
    # PATHOGEN='other_non_flu'
    
    shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = SOURCE)
    
    for (PATHOGEN in names(pathogenKeys)){

     
      
      queryIn <- list(
        SELECT   =list(COLUMN=c('pathogen', factors, GEO,'encountered_week')),
        WHERE    =list(COLUMN='pathogen', IN=pathogenKeys[[PATHOGEN]]),
        GROUP_BY =list(COLUMN=c(factors,GEO,"encountered_week")),
        SUMMARIZE=list(COLUMN='pathogen', IN= pathogenKeys[[PATHOGEN]])
      )
      
      
      db <- expandDB(selectFromDB(  queryIn, source=SRC, na.rm=TRUE ), shp=shp, currentWeek=currentWeek)
      
      #if you want to add the ILI data to the db
      db <- appendILIDataFc(db, currentWeek)
      
      # training occassionaly segfaults on but it does not appear to be deterministic...
      tries <- 0
      success<-0
      while (success==0 & tries<=2){
        tries <- tries+1
        tryCatch(
          {
            
            db <- appendCatchmentModel(db,shp=shp, source=SRC, na.rm=TRUE )

            modelDefinition <- latentFieldModel(db=db, shp=shp)
            model <- modelTrainR(modelDefinition)
            
            print(summary(model$inla))
            
            dir.create('/home/rstudio/seattle_flu/model_diagnostic_plots/', showWarnings = FALSE)
            fname <- paste('/home/rstudio/seattle_flu/model_diagnostic_plots/',paste('inla_latent',PATHOGEN,SOURCE,GEO,'encountered_week',sep='-'),'.png',sep='')
            png(filename = fname,width = 6, height = 5, units = "in", res = 300)
            print(
              ggplot(model$latentField) + 
                    geom_line(aes_string(x='encountered_week',y="modeled_intensity_median", color=GEO,group =GEO)) + 
                    geom_ribbon(aes_string(x='encountered_week',ymin="modeled_intensity_lower_95_CI", ymax="modeled_intensity_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
                    guides(color=FALSE, fill=FALSE) + 
                    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
              )
            dev.off()
            
            # ggplot(model$modeledData %>% filter(site_type %in% 'retrospective' & flu_shot=='false' & sex=='female')) + 
            #   geom_line(aes_string(x='encountered_week',y="modeled_count_median", color=GEO,group =GEO)) + 
            #   geom_ribbon(aes_string(x='encountered_week',ymin="modeled_count_lower_95_CI", ymax="modeled_count_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
            #   guides(color=FALSE, fill=FALSE) + 
            #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
            # 
            
            saveModel(model, storeRDS=FALSE)
            
            success<-1
            
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
        )
      }
    }
  }
}

