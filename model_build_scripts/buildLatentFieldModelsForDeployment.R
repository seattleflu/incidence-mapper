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
mostRecentSample<-max(db$observedData$encountered_date)

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

pathogenKeys <- list(
                     flu=fluPathogens, Flu_A_H1 = 'Flu_A_H1', Flu_A_H3 = 'Flu_A_H3', Flu_B_pan = 'Flu_B_pan', Flu_C_pan = 'Flu_C_pan' #,
                     # all='all', other_non_flu = setdiff(pathogens$pathogen,c(fluPathogens,'not_yet_tested','measles','Measles')),
                     # rsv = c('RSVA','RSVB'), RSVA='RSVA', RSVB='RSVB', 
                     # AdV='AdV',CoV='CoV',RV='RV'
                     )


factors   <- c('site_type','sex','flu_shot','age_range_coarse_upper')


geoLevels <- list(
                   sfs_domain_geojson = 'residence_regional_name'#,
                   # seattle_geojson = c('residence_puma','residence_neighborhood_district_name','residence_cra_name'), #,'residence_census_tract'),
                   # wa_geojson = c('residence_puma')
                 )

siteTypes <- c('childrensHospital','clinic','collegeCampus','childrensClinic','port','retrospective','workplace','publicSpace','self-test')


# currentWeek <- '2019-W25'
currentWeek <- paste(isoyear(Sys.time()) ,'-W',isoweek(Sys.time()),sep='')


#####################################
###### timeseries latent field models ############
#####################################

# number of subjects with pathogen and factor at residence location 
for (SOURCE in names(geoLevels)){
  for (GEO in geoLevels[[SOURCE]]){
    
    # SOURCE='sfs_domain_geojson'
    # GEO='residence_regional_name'
    # PATHOGEN='flu'
    # PATHOGEN='Flu_C_pan'
    
    # SOURCE='seattle_geojson'
    # SOURCE='wa_geojson'
    # GEO='residence_census_tract'
    # GEO='residence_puma'
    # GEO='residence_cra_name'
    # PATHOGEN='all'
    # PATHOGEN='rsv'
    # PATHOGEN='other_non_flu'
    
    shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = SOURCE)
    
    for (PATHOGEN in names(pathogenKeys)){

      queryIn <- list(
        SELECT   =list(COLUMN=c('pathogen', factors, GEO,'encountered_week')),
        WHERE    =list(COLUMN='pathogen', IN=pathogenKeys[[PATHOGEN]]),
        WHERE    =list(COLUMN='site_type', IN=siteTypes),
        GROUP_BY =list(COLUMN=c(factors,GEO,"encountered_week")),
        SUMMARIZE=list(COLUMN='pathogen', IN= pathogenKeys[[PATHOGEN]])
      )

      db <- expandDB(db<-selectFromDB(  queryIn, source=SRC, na.rm=TRUE ), shp=shp, currentWeek=currentWeek)
      
      
      
      db$observedData <- db$observedData %>% filter(site_type=='retrospective')
      
      db$observedData$positive[db$observedData$encountered_week >= paste('2019-W',isoweek(mostRecentSample), sep='')] <- NaN
      
      
      
      # hack in all flu timeseries
      db2 <- read.csv('all_flu_by_time_query_result_2019-12-11T22_02_07.419Z.csv')
      lineages <- unique(db2$lineage)
      levels(db2$lineage)<-fluPathogens[c(1,2,4,5)]
      names(db2)[1]<-'pathogen'
      names(db2)[2]<-'encountered_week'
      
      db2$encountered_week <- factor(db2$encountered_week, levels=sort(unique(db2$encountered_week)), ordered = TRUE)
      
      
      countData<-list()
      countData$queryList<-list(GROUP_BY=list(COLUMN='encountered_week'))
      countData$observedData <- db2 %>% filter(pathogen %in% pathogenKeys[[PATHOGEN]]) %>% 
        group_by(encountered_week) %>% summarize(positive = sum(case_count), n = sum(case_count))
      
      ## make sure always extrapolates to currentWeek 
      
        mostRecentWeek <- max(countData$observedData$encountered_week)
        
        minYear <- year<-as.numeric(gsub('-W[0-9]{2}','',mostRecentWeek))
        maxYear <- year<-as.numeric(gsub('-W[0-9]{2}','',currentWeek))
        minWeek <- as.numeric(gsub('[0-9]{4}-W','',mostRecentWeek ))
        maxWeek <- as.numeric(gsub('[0-9]{4}-W','',currentWeek )) + (maxYear-minYear)*52 
        
        if(minWeek <maxWeek & minYear<=maxYear){
          weeks <- 1+( (seq(minWeek+1,maxWeek,by=1)-1) %% 52)
          yearBreaks <- c(0,which(diff(weeks)<1), length(weeks))
          years=c()
          for (k in 2:length(yearBreaks)){
            years <- c(years, rep(minYear+(k-2), yearBreaks[k]-yearBreaks[k-1]  ))
          }
          
          forecast_weeks <- paste(years,'-W',sprintf("%02d",weeks),sep='')
        } else {
          forecast_weeks <- NULL
        }
        
        if(!is.null(forecast_weeks)){
          countData$observedData <- countData$observedData %>% tibble::add_row(encountered_week = forecast_weeks, n=0)
        }

      countData$observedData$positive[countData$observedData$encountered_week >= as.character(mostRecentWeek)] <- NaN
      
      
      countData$observedData$time_row<-as.numeric(countData$observedData$encountered_week)
      
      
      countModelDef<-smoothModel(countData)
      countModel <- modelTrainR(countModelDef)
    
      countModel$modeledData$log_all_flu_count_field_effect <- log(countModel$modeledData$modeled_count_mean)

      db$observedData <- db$observedData %>% left_join( countModel$modeledData %>% select(encountered_week,log_all_flu_count_field_effect))
      db$observedData$encountered_week <- factor(db$observedData$encountered_week, levels=sort(unique(db$observedData$encountered_week)), ordered = TRUE)
      
      
      #if you want to add the ILI data to the db
      # latentFieldModel can't handle this correctly right now
      # db <- appendILIDataFc(db, currentWeek)
      
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
                    # geom_ribbon(aes_string(x='encountered_week',ymin="modeled_intensity_lower_95_CI", ymax="modeled_intensity_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
                    guides(color=FALSE, fill=FALSE) + 
                    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
              )
            dev.off()
            
            dir.create('/home/rstudio/seattle_flu/model_diagnostic_plots/', showWarnings = FALSE)
            fname <- paste('/home/rstudio/seattle_flu/model_diagnostic_plots/',paste('inla_latent',PATHOGEN,SOURCE,GEO,'age',sep='-'),'.png',sep='')
            png(filename = fname,width = 6, height = 5, units = "in", res = 300)
            plotDat <- model$inla$summary.random$site_age_siteIdx
            plotDat$group <- unique(model$modeledData$site_type)
            plotDat$age <- rep(unique(model$modeledData$age_range_coarse_upper), each = length(unique(model$modeledData$site_type)))
            N<-length(unique(model$modeledData$age_range_coarse_upper))
            plotDat$ageID <- rep(1:N, each = length(unique(model$modeledData$site_type)))
            print(
              ggplot(plotDat) + geom_point(aes(x=ageID, y=mean), stat='identity') + theme_bw() +
                geom_linerange(aes(x=ageID, ymin=`0.025quant`, ymax=`0.975quant`)) +
                xlab('age_range_coarse_upper') + ylab('age random effect') + facet_wrap("group")+
                scale_x_continuous(breaks=1:N,labels=as.character(unique(plotDat$age)))
            )
            dev.off()
            
              
            # ggplot(model$modeledData %>% filter(site_type %in% 'retrospective' & flu_shot=='false' & sex=='female')) + 
            #   geom_line(aes_string(x='encountered_week',y="modeled_count_median", color=GEO,group =GEO)) + 
            #   geom_ribbon(aes_string(x='encountered_week',ymin="modeled_count_lower_95_CI", ymax="modeled_count_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
            #   guides(color=FALSE, fill=FALSE) + 
            #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
            # 
            
            saveModel(model, storeRDS=FALSE)
            # saveModel(model, storeRDS=TRUE)
            
            success<-1

        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
        )
      }
    }
  }
}

