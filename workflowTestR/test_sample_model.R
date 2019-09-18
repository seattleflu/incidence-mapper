library('dbViewR')
library('incidenceMapR')

## simulated data with bing query search data, run fixed effects model
queryIn <- list(
  SELECT   =list(COLUMN=c('site_type','residence_census_tract')),
  WHERE    =list(COLUMN='site_type', IN = c('kiosk')),
  GROUP_BY =list(COLUMN=c('site_type','residence_census_tract')),
  SUMMARIZE=list(COLUMN='site_type', IN= c('kiosk'))
)

db <- expandDB( selectFromDB(  queryIn,  source = 'simulated_data' ) )
shp <- masterSpatialDB(shape_level = 'census_tract', source = 'simulated_data', rm_files = TRUE)

mycompare <- function(target, current) { all.equal(target, current, tolerance = 1e-3)}


unitizer_sect("modelDefinition Test", compare=testFuns(
  value=mycompare, conditions=mycompare, output=function(x,y) TRUE, 
  message=function(x,y) TRUE, aborted=function(x,y) TRUE), 
              {
                (modelDefinition <- smoothModel(db=db, shp=shp))
              })

unitizer_sect("modelDefinition Test", compare=testFuns(
  value=mycompare, conditions=mycompare, output=function(x,y) TRUE, 
  message=function(x,y) TRUE, aborted=function(x,y) TRUE), 
              {
                (model <- modelTrainR(modelDefinition))
              })


