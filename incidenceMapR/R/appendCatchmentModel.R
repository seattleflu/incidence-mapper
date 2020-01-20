#' appendCatchmentModel: function for adding catchment covariate estimated by smoothModel to dbViewR object for use in latentFieldModel
#' 
#' @param db object from dbViewer with observedData tibble and query
#' @param shp shapefile object from masterShapeDB
#' @param source = 'simulated_data' (default) or 'production
#' @param na.rm = TRUE (default) or FALSE.  Remove NA pathogens
#' @return db with added catchment column to observedData
#' 
#' @import INLA
#' @import dbViewR
#' @import magrittr
#' @import dplyr
#'
#' @export
#'
appendCatchmentModel <- function(db,shp = NULL, source='production', na.rm=TRUE){
  
  validGeoLevels <- c('residence_puma','residence_cra_name','residence_neighborhood_district_name','residence_census_tract','residence_city','residence_regional_name',
                      'work_puma','work_cra_name','work_neighborhood_district_name','work_census_tract','work_city','work_regional_name')
  geo <- validGeoLevels[validGeoLevels %in% names(db$observedData)]
  
  # get pathogen list
  queryIn <- list(
    SELECT   =list(COLUMN=c('pathogen')),
    GROUP_BY =list(COLUMN=c('pathogen')),
    SUMMARIZE=list(COLUMN='pathogen', IN= 'all')
  )
  pathogens <- unique(selectFromDB(  queryIn, source=source, na.rm=na.rm)$observedData$pathogen)
  
  # find catchment maps for each site_type and geoLevel
  # catchment represented by all samples not from target virus. Idea being participation due to other pathogens is assumed
  # to be uncorrelated with pathogen of interest. Social dynamics could violate this assumption.
  
  # making timeseries inference work for pathogen == 'all even if it maybe shouldn't...
    inputPathogen <- unique(db$observedData$pathogen)
    
    if (inputPathogen == 'all'){
      outGroup <- 'all'
    } else if (inputPathogen == 'flu_positive'){
      outGroup <- pathogens[!grepl('flu',pathogens, ignore.case = TRUE)]
    } else if (inputPathogen == 'rsv_positive'){
      outGroup <- pathogens[!grepl('rsv',pathogens, ignore.case = TRUE)]
    } else if (inputPathogen == 'flu_negative'){
      outGroup <- pathogens[grepl('flu',pathogens, ignore.case = TRUE)]
    } else {
      outGroup<-setdiff(pathogens, unlist(strsplit(inputPathogen,'-')))
    }
  
  siteCols = names(db$observedData)[grepl('site',names(db$observedData))]
  queryIn <- list(
    SELECT   =list(COLUMN=c('pathogen',siteCols,geo)),
    WHERE    =list(COLUMN=c('pathogen'), IN = outGroup),
    WHERE    =list(COLUMN=siteCols, IN = unique(db$observedData[[siteCols]])),
    GROUP_BY =list(COLUMN=c(siteCols,geo)),
    SUMMARIZE=list(COLUMN=siteCols, IN= 'all')
  )
  
  catchmentDb <- selectFromDB(  queryIn, source=source, na.rm=na.rm )
  catchmentDb <- expandDB( catchmentDb, shp=shp )
  
  
  # positives as 0 instead of NaN when positive count is total count always (eg catchments) 
  catchmentDb$observedData$positive[is.na(catchmentDb$observedData$positive)]<-0
  
  tmp<-catchmentDb$observedData %>% group_by(!!as.name(siteCols)) %>% summarize(positive = sum(positive,na.rm=TRUE))
  dropList <- tmp %>% filter(positive<10)
  for (k in 1:nrow(dropList)){
    catchmentDb$observedData <- catchmentDb$observedData %>%  filter(!((!!as.name(siteCols) %in% dropList$site[k])))
  }
  
  catchmentModelDefinition <- smoothModel(db=catchmentDb, shp=shp)
  catchmentModel <- modelTrainR(catchmentModelDefinition)
  
  # drop sites without enough data to estimate catchment
  # edge case for new sites (like publicSpace)
  for (COLUMN in siteCols){
    db$observedData <- db$observedData %>% filter((!!as.name(COLUMN)) %in% unique(catchmentDb$observedData[[COLUMN]]))
  }
  
  # append catchment as intercept covariate
  db$observedData <- db$observedData %>% left_join(catchmentModel$modeledData %>% select(!!as.name(siteCols), geo, modeled_count_median))
  names(db$observedData)[names(db$observedData) %in% 'modeled_count_median'] <- 'catchment'
  db$observedData$catchment <- log(db$observedData$catchment) - mean(log(db$observedData$catchment))
  
  return(db)
}
