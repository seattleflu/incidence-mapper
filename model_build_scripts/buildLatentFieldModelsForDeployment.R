# buildModelsForDeployment.R
# this script (and similar others?) controls standardized database queries and model training for web deployment

library(dbViewR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

SRC <- 'production'
# SRC <- 'simulated_data'

db <- selectFromDB(queryIn= list(SELECT  =c("*")), source = SRC)
print(mostRecentSample<-max(db$observedData$encountered_date))
unique(db$observedData$pathogen)

unique(db$observedData$site_type)

as.data.frame(db$observedData %>% group_by(site) %>% summarize(site_type=unique(site_type)),30)



fluPathogens <- c('Flu_A_H1','Flu_A_H3','Flu_A_pan','Flu_B_pan','Flu_C_pan')

pathogens <- db$observedData %>% group_by(pathogen) %>% summarize(n = n()) %>% arrange(desc(n))

tmp<-pathogens$pathogen[pathogens$n >= 500 | grepl('Flu',pathogens$pathogen)]
tmp2<-as.list(tmp)
names(tmp2)<-tmp



# pathogenKeys <- c(list(all='all',
#                   flu=fluPathogens,
#                   rsv=c('RSVA','RSVB'),
#                   other_non_flu = setdiff(pathogens$pathogen,c(fluPathogens,'not_yet_tested','measles'))
#                   ),
#                   tmp2)

pathogenKeys <- list(
                     flu=fluPathogens, Flu_A_H1 = 'Flu_A_H1', Flu_A_H3 = 'Flu_A_H3', Flu_B_pan = 'Flu_B_pan', Flu_C_pan = 'Flu_C_pan',
                     all='all', #other_non_flu = setdiff(pathogens$pathogen,c(fluPathogens,'not_yet_tested','measles','Measles'))#,
                     rsv = c('RSVA','RSVB'), RSVA='RSVA', RSVB='RSVB'#,
                     # AdV='AdV',CoV='CoV',RV='RV'
                     )


# factors   <- c('sex','flu_shot','age_range_coarse_upper','sfs_year')
factors   <- c('age_range_coarse_upper','sfs_year')


geoLevels <- list(
                   sfs_domain_geojson = 'residence_regional_name'#,
                   # seattle_geojson = c('residence_puma','residence_neighborhood_district_name','residence_cra_name'), #,'residence_census_tract'),
                   # wa_geojson = c('residence_puma')
                 )

siteTypes <- c("retrospective","collegeCampus","clinic","home","publicSpace","workplace")

# nowcastWeek <- '2019-W25'
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
for (SOURCE in names(geoLevels)){
  for (GEO in geoLevels[[SOURCE]]){

    SOURCE='sfs_domain_geojson'
    GEO='residence_regional_name'
    PATHOGEN='flu'

    # PATHOGEN='Flu_C_pan'
    # SOURCE='seattle_geojson'
    # SOURCE='wa_geojson'
    # GEO='residence_census_tract'
    # GEO='residence_puma'
    # GEO='residence_cra_name'
    # PATHOGEN='all'
    # PATHOGEN='rsv'
    # PATHOGEN='other_non_flu'
    # PATHOGEN='RSVA'
    
    
    shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = SOURCE)
    
    for (PATHOGEN in names(pathogenKeys)){

      queryIn <- list(
        SELECT   =list(COLUMN=c('pathogen', factors, GEO,'encountered_week','site_type','site')),
        WHERE    =list(COLUMN='pathogen', IN=pathogenKeys[[PATHOGEN]]),
        WHERE    =list(COLUMN='site_type', IN=siteTypes),
        GROUP_BY =list(COLUMN=c(factors,GEO,"encountered_week",'site_type')),
        SUMMARIZE=list(COLUMN='pathogen', IN= pathogenKeys[[PATHOGEN]])
      )
      db <- selectFromDB(  queryIn, source=SRC, na.rm=TRUE )
      db <- expandDB(db, shp=shp, currentWeek =nowcastWeek)
    
      as.data.frame(db$observedData %>% group_by(site,sfs_year) %>% summarize(positive = sum(positive,na.rm=TRUE)))
      as.data.frame(db$observedData %>% group_by(site_type,sfs_year) %>% summarize(positive = sum(positive,na.rm=TRUE)))
      
      # forecast past incomplete data
      # db$observedData$positive[db$observedData$encountered_week >= paste('2019-W',isoweek(mostRecentSample), sep='')] <- NA
      
      
      # db$observedData$positive[db$observedData$encountered_week > '2019-W52']<-NA
      
      # temporarily removing childrens hospital update
      db$observedData$positive[db$observedData$site=='RetrospectiveChildrensHospitalSeattle' &
                                 db$observedData$encountered_week > '2019-W45']<-NA
      db$observedData$positive[db$observedData$site=='retrospective' &
                                 db$observedData$encountered_week > '2019-W45']<-NA
      
      
      
      library(incidenceMapR)
      library(modelServR)
      library(modelVisualizeR)
      
      
      # hack in all flu lab timeseries
      db2 <- read.csv('all_flu_by_time_query_result_2020-01-19T13_15_24.869283-08_00.csv')
      lineages <- unique(db2$lineage)
      levels(db2$lineage)<-fluPathogens[c(1,2,4,5)]
      names(db2)[1]<-'pathogen'
      names(db2)[2]<-'encountered_week'
      
      db2$encountered_week <- factor(db2$encountered_week, levels=sort(unique(db2$encountered_week)), ordered = TRUE)
      
      
      countData<-list()
      countData$queryList<-list(GROUP_BY=list(COLUMN='encountered_week'))
      countData$observedData <- db2 %>% filter(pathogen %in% pathogenKeys[[PATHOGEN]]) %>% 
        group_by(encountered_week) %>% summarize(positive = sum(case_count), n = sum(case_count))
      
      ## make sure always extrapolates to nowcastWeek 
      
        mostRecentWeek <- max(countData$observedData$encountered_week)
        # mostRecentWeek <- '2019-W52'
        
        minYear <- as.numeric(gsub('-W[0-9]{2}','',mostRecentWeek))
        maxYear <- as.numeric(gsub('-W[0-9]{2}','',nowcastWeek))
        minWeek <- as.numeric(gsub('[0-9]{4}-W','',mostRecentWeek ))
        maxWeek <- as.numeric(gsub('[0-9]{4}-W','',nowcastWeek )) + (maxYear-minYear)*52 
        
        if(minWeek <maxWeek & minYear<=maxYear){
          weeks <- 1+( (seq(minWeek+1,maxWeek,by=1)-1) %% 52)
          if(weeks[1]==1 && length(weeks)<52){
            years <- rep(maxYear,length(weeks))
          } else {
            yearBreaks <- c(0,which(diff(weeks)<1), length(weeks))
            years=c()
            for (k in 2:length(yearBreaks)){
              years <- c(years, rep(minYear+(k-2), yearBreaks[k]-yearBreaks[k-1]  ))
            }
          }

          forecast_weeks <- paste(years,'-W',sprintf("%02d",weeks),sep='')
        } else {
          forecast_weeks <- NULL
        }
        
        if(!is.null(forecast_weeks)){
          countData$observedData <- countData$observedData %>% tibble::add_row(encountered_week = forecast_weeks, n=0) %>% distinct(encountered_week,.keep_all = TRUE)
        }
  
      # filter out most recent week, which is likely very data incomplete
      countData$observedData$positive[countData$observedData$encountered_week >= as.character(mostRecentWeek)] <- NA
      tail(countData$observedData)
      
      countData$observedData$time_row<-as.numeric(countData$observedData$encountered_week)
      
      
      countModelDef<-smoothModel(countData)
      countModel <- modelTrainR(countModelDef)
    
      countModel$modeledData$log_all_flu_count_field_effect <- log(countModel$modeledData$modeled_count_median)
      
      summary(countModel$inla)
      plot(countModel$modeledData$modeled_count_median)
      lines(countData$observedData$positive)
      
      db$observedData <- db$observedData %>% left_join( countModel$modeledData %>% select(encountered_week,log_all_flu_count_field_effect))
      db$observedData$encountered_week <- factor(db$observedData$encountered_week, levels=sort(unique(db$observedData$encountered_week)), ordered = TRUE)
      
      
      #if you want to add the ILI data to the db
      # latentFieldModel can't handle this correctly right now
      # db <- appendILIDataFc(db, nowcastWeek)
      
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
            
            
            ## rescale latent field model for consistent color code
            # reference_scale <- 0.201
            reference_scale <- 0.201*.12/.1
            
            
            current_scale <- model$latentField %>% select(encountered_week,residence_regional_name,modeled_intensity_mode) %>%
              group_by(encountered_week) %>% summarize(region_intensity_mode = mean(modeled_intensity_mode))
            current_scale <- max(current_scale$region_intensity_mode)
            
            scale_factor <- reference_scale/current_scale
            
            for (COLUMN in names(model$latentField)[grepl('modeled_',names(model$latentField))]){
              model$latentField[[COLUMN]] <- model$latentField[[COLUMN]] * scale_factor
            }
            
            dir.create('/home/rstudio/seattle_flu/model_diagnostic_plots/', showWarnings = FALSE)
            fname <- paste('/home/rstudio/seattle_flu/model_diagnostic_plots/',paste('inla_latent',PATHOGEN,SOURCE,GEO,'encountered_week',Sys.Date(),sep='-'),'.png',sep='')
            png(filename = fname,width = 6, height = 5, units = "in", res = 300)
            print(
              ggplot(model$latentField) + 
                    geom_line(aes_string(x='encountered_week',y="modeled_intensity_median", color=GEO,group =GEO)) + 
                    # geom_ribbon(aes_string(x='encountered_week',ymin="modeled_intensity_lower_95_CI", ymax="modeled_intensity_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
                    guides(color=FALSE, fill=FALSE) +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
              )
            dev.off()
            
            
            ggplot(model$modeledData %>% group_by_('encountered_week', GEO) %>% summarize(modeled_count_mode=sum(modeled_count_mode, na.rm=TRUE))) + 
              geom_line(aes_string(x='encountered_week',y="modeled_count_mode", color=GEO,group =GEO)) + 
              # facet_wrap('site_type',scales = 'free_y') +
              guides(color=FALSE, fill=FALSE) + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
        
            ggplot(model$modeledData %>% group_by_('encountered_week', GEO) %>% summarize(positive=sum(positive, na.rm=TRUE))) + 
              geom_line(aes_string(x='encountered_week',y="positive", color=GEO,group =GEO)) + 
              # facet_wrap('site_type',scales = 'free_y') +
              guides(color=FALSE, fill=FALSE) + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
            
            
            dir.create('/home/rstudio/seattle_flu/model_diagnostic_plots/', showWarnings = FALSE)
            fname <- paste('/home/rstudio/seattle_flu/model_diagnostic_plots/',paste('inla_latent',PATHOGEN,SOURCE,GEO,'age',Sys.Date(),sep='-'),'.png',sep='')
            png(filename = fname,width = 6, height = 5, units = "in", res = 300)
            if ('site_age_siteIdx' %in% names(model$inla$summary.random)){
              plotDat <- model$inla$summary.random$site_age_siteIdx
            } else {
              plotDat <- model$inla$summary.random$age_row_iid
            }
            plotDat$group <- unique(model$modeledData$site)
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

