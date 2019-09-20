# buildModelsForDeployment.R
# this script (and similar others?) controls standardized database queries and model training for web deployment

library(dbViewR)
library(incidenceMapR)
library(modelServR)
library(modelVisualizeR)
library(dplyr)
library(tidyr)
library(ggplot2)

SRC <- 'production'
# SRC <- 'simulated_data'

fluABtestPathogens <- c('Flu_A_H1','Flu_A_H3','Flu_A_pan','Flu_B_pan')

factors   <- c('site_type')

geoLevels <- list( king_county_geojson = c('residence_puma'),
                   wa_geojson = c('residence_puma'),
                   seattle_geojson = c('residence_neighborhood_district_name')
                 )

#####################################
###### DeDx smooth model ############
#####################################

for (SOURCE in names(geoLevels)){
  for (GEO in geoLevels[[SOURCE]]){

    shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = SOURCE)
    
    queryIn <- list(
      SELECT   =list(COLUMN=c('pathogen',factors, GEO,'encountered_week','age_range_fine_upper')),
      GROUP_BY =list(COLUMN=c(factors,GEO,'encountered_week','age_range_fine_upper')),
      SUMMARIZE=list(COLUMN='pathogen', IN= fluABtestPathogens)
    )
    
    db <- expandDB(selectFromDB(  queryIn, source=SRC, na.rm=TRUE ), shp=shp )
    
    # right now, self-test is missing age and residence data, and publicSpace is missing age data.
    # This is a database issue that will be fixed, so this is hack to fill them in for prediction. 
    
    fill_site_types <- c('publicSpace','self-test')
    for( FILL in setdiff(fill_site_types,unique(db$observedData$site_type)) ){
      if (GEO == 'residence_neighborhood_district_name'){
        fillDat <- expand.grid(encountered_week=sort(unique(db$observedData$encountered_week)),
                               residence_neighborhood_district_name=sort(unique(db$observedData$residence_neighborhood_district_name)),
                               site_type=FILL,
                               age_range_fine_upper = sort(unique(db$observedData$age_range_fine_upper)),
                               pathogen = sort(unique(db$observedData$pathogen)),
                               n=0,
                               positive=NA)
      } else if (GEO=='residence_puma'){
        fillDat <- expand.grid(encountered_week=sort(unique(db$observedData$encountered_week)),
                               residence_puma=sort(unique(db$observedData$residence_puma)),
                               site_type=FILL,
                               age_range_fine_upper = sort(unique(db$observedData$age_range_fine_upper)),
                               pathogen = sort(unique(db$observedData$pathogen)),
                               n=0,
                               positive=NA)
      } else if (GEO=='residence_census_tract'){
        fillDat <- expand.grid(encountered_week=sort(unique(db$observedData$encountered_week)),
                               residence_census_tract=sort(unique(db$observedData$residence_census_tract)),
                               site_type=FILL,
                               age_range_fine_upper = sort(unique(db$observedData$age_range_fine_upper)),
                               pathogen = sort(unique(db$observedData$pathogen)),
                               n=0,
                               positive=NA)
      }
      fillDat$time_row <- match(fillDat$encountered_week, sort(unique(db$observedData$encountered_week)))
      fillDat$age_row <- match(fillDat$age_range_fine_upper, sort(unique(db$observedData$age_range_fine_upper)))
      
      db$observedData <- rbind(db$observedData,fillDat)
    }
      
    # training occassionaly segfaults on but it does not appear to be deterministic...
    tries <- 0
    success<-0
    while (success==0 & tries<=2){
      tries <- tries+1
      tryCatch(
        {
          
          modelDefinition <- smoothModel(db=db, shp=shp)
          model <- modelTrainR(modelDefinition)
          
          print(summary(model$inla))
          
          saveModel(model)
          
          vizSite <- unique(db$observedData$site_type)
          

          for(SITE in vizSite){
          
              
            dir.create('/home/rstudio/seattle_flu/model_diagnostic_plots/', showWarnings = FALSE)
            fname <- paste('/home/rstudio/seattle_flu/model_diagnostic_plots/',paste('DeDx',SOURCE,GEO,SITE,'age_range_fine_upper','encountered_week','mode',sep='-'),'.png',sep='')
            png(filename = fname,width = 6, height = 5, units = "in", res = 300)
            print(
              ggplot(model$modeledData %>% filter(site_type == SITE)) +
                    geom_line(aes_string(x='encountered_week',y="modeled_fraction_median", group=GEO,color=GEO)) +
                    facet_wrap('age_range_fine_upper') +
                    # geom_ribbon(aes_string(x='encountered_week',ymin="modeled_fraction_lower_95_CI", ymax="modeled_fraction_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
                    guides(color=FALSE) +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1))
              )
            dev.off()
            
            fname <- paste('/home/rstudio/seattle_flu/model_diagnostic_plots/',paste('DeDx',SOURCE,GEO,SITE,'age_range_fine_upper','encountered_week','sd',sep='-'),'.png',sep='')
            png(filename = fname,width = 6, height = 5, units = "in", res = 300)
            print(
              ggplot(model$modeledData %>% filter(site_type == SITE)) +
                geom_line(aes_string(x='encountered_week',y="modeled_fraction_sd", group=GEO,color=GEO)) +
                facet_wrap('age_range_fine_upper') +
                guides(color=FALSE) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))
            )
            dev.off()
          }
          
          success<-1
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
      )
    }
  }
}
