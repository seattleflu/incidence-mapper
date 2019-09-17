library(testthat)
library('dbViewR')
library('incidenceMapR')

context("test sample model")


## simulated data with bing query search data, run fixed effects model
queryIn <- list(
  SELECT   =list(COLUMN=c('site_type','residence_census_tract')),
  WHERE    =list(COLUMN='site_type', IN = c('kiosk')),
  GROUP_BY =list(COLUMN=c('site_type','residence_census_tract')),
  SUMMARIZE=list(COLUMN='site_type', IN= c('kiosk'))
)

db <- expandDB( selectFromDB(  queryIn,  source = 'simulated_data' ) )
shp <- masterSpatialDB(shape_level = 'census_tract', source = 'simulated_data', rm_files = TRUE)

modelDefinition <- smoothModel(db=db, shp=shp)
model <- modelTrainR(modelDefinition)
