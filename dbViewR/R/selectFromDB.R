#' selectFromDB: function for fetching data from research DB
#' (Currently pulls from simulated data at https://github.com/seattleflu/simulated-data)
#'
#' STANDARD DB QUERIES WILL ALL LIKELY MIGRATE TO THE HUTCH RESARCH DB BEFORE PRODUCTION.
#'
#' @param queryList  list or json specifying query  (See example)
#' @param source source database, one of: 'simulated_data' (default) or 'production'
#' @param credentials_path path to your pg_service and pgpass file for production database
#' @param na.rm = FALSE (default) Drop rows with NA from dataset as incidenceMapR will ignore them anyway
#' @param shp = dbViewR::masterSpatialDB(shape_level = 'census_tract', source = 'king_county_geojson') (Needed hack until higher-level shape labels are in database)
#' @return dbViewR list with query and observedData table that has been prepared for defineModels.R
#'
#' @import jsonlite
#' @import dplyr
#' @import lubridate
#' @import DBI
#' @import RPostgres
#' @importFrom RCurl getURL
#' @importFrom magrittr %>%
#' @importFrom lazyeval interp
#' @import tidyr 
#'
#' @export
#' @examples
#' return h1n1pdm summary by time and location
#' queryJSON <- jsonlite::toJSON(
#'   list(
#'       SELECT   =list(COLUMN=c('pathogen','encountered_week','residence_puma','residence_census_tract')),
#'       GROUP_BY =list(COLUMN=c('encountered_week','residence_puma','residence_census_tract')),
#'       SUMMARIZE=list(COLUMN='pathogen', IN= c('h1n1pdm'))
#'       )
#'    )
#'    db <- selectFromDB( queryJSON )
#'
selectFromDB <- function( queryIn = jsonlite::toJSON(
                            list(
                              SELECT   =list(COLUMN=c('pathogen','encountered_date','residence_puma','residence_census_tract')),
                              GROUP_BY =list(COLUMN=c('encountered_week','residence_puma','residence_census_tract')),
                              SUMMARIZE=list(COLUMN='pathogen', IN= c('h1n1pdm'))
                            )
                          ), source = 'production', 
                          credentials_path = '/home/rstudio/seattle_flu',
                          na.rm = FALSE
                          ){

  if(class(queryIn) == "json"){
    queryList <- jsonlite::fromJSON(queryIn)
  } else if(class(queryIn) == "list"){
    queryList<-queryIn
  }
  

  # connect to database
  if(source == 'simulated_data'){

    rawData <- RCurl::getURL("https://raw.githubusercontent.com/seattleflu/simulated-data/master/simulated_subject_database.csv")
    db <- read.table(text = rawData, header=TRUE, sep=",", stringsAsFactors = FALSE)

  } else if(source == 'production'){

    # Define standard Pg environment variables for our connection files.
    #
    # These are set here instead of passed in via the Dockerfile or
    # docker-compose.rstudio.yaml file, because this typically runs as an
    # rstudio-managed user with its own shell environment separate from the
    # base Docker user.
    Sys.setenv(PGSERVICEFILE = file.path(credentials_path, ".pg_service.conf"),
               PGPASSFILE    = file.path(credentials_path, ".pgpass"))

    # Connect to database using service definition and credentials in files
    # defined by the environment variables above.
    rawData <- DBI::dbConnect(RPostgres::Postgres(), service = "seattleflu-production")

    db <- DBI::dbGetQuery(rawData, "select distinct * from shipping.incidence_model_observation_v2;")
    
    # db <- DBI::dbGetQuery(rawData, paste('select distinct * from shipping.incidence_model_observation_v1 encounter',
    #                                      'left join shipping.presence_absence_result_v1 taq',
    #                                      'on encounter.sample = taq.sample',
    #                                      ';'),sep=' ') 
    
    # this logic should be substantially rethought, as I'm mixing sql and dplyr in confusing ways, but it will have to do for now!
    
    # get all samples and nest
    db2 <- DBI::dbGetQuery(rawData, paste('select distinct * from shipping.presence_absence_result_v1',
                                          ';'),sep=' ') 

    names(db2)[names(db2)=='target'] <- 'pathogen'
    
    # count pathogens found and tests performed
    db2 <- db2 %>% group_by(sample) %>%
      mutate(number_pathogens_found = sum(present), number_pathogens_tested = n())
    
    # add in "undetected" pathogen for samples that were tested but had no detections
    db3 <- db2 %>% group_by(sample,number_pathogens_found,number_pathogens_tested) %>% filter(all(present == FALSE) &  all(number_pathogens_tested>0)) %>%
      summarize(pathogen = 'undetected') %>% mutate(present=TRUE)
    
    # join undetecteds with positives
    db4 <- bind_rows(db2 %>% filter(present == TRUE), db3)
    
    # join with encounter list, using nice formatting
    db <- db %>% left_join(db4) 
    
    # put in "not_yet_tested" for samples with no test results
    idx<-is.na(db$number_pathogens_tested)
    db$number_pathogens_tested[idx] <- 0
    db$pathogen[idx] <- 'not_yet_tested'
    db$present[idx] <- TRUE
    
    # drop present column if all TRUE
    if(all(db$present==TRUE)){
      db <- db %>% select(-present)
    }
    
    DBI::dbDisconnect(rawData)


  } else {
     print('unknown source database!')
  }

 
  # combine PCR targets that describe one pathogen family
  db$pathogen[db$pathogen %in% c('Adenovirus_pan_2','Adenovirus_pan_1','AdV_1of2','AdV_2of2') ] <- 'AdV'
  db$pathogen[db$pathogen %in% c('12 Rhinovirus_pan_2','11 Rhinovirus_pan_1','RV_1of2','RV_2of2') ] <- 'RV'
  db$pathogen[db$pathogen %in% c('Influenza_B','Flu_b_pan') ] <- 'Flu_B_pan'
  db$pathogen[db$pathogen %in% c('flu_A_pan','Flu_a_pan') ] <- 'Flu_A_pan'
  db$pathogen[db$pathogen %in% c('AP324NU') ] <- 'Flu_C_pan'
  db$pathogen[db$pathogen %in% c('hPIV1','hPIV2','hPIV1_hPIV2','hPIV3','hPIV4','hPIV3_hPIV4') ] <- 'hPIV'
  db$pathogen[db$pathogen %in% c('CoV_229E_CoV_OC43','CoV_HKU1_CoV_NL63','CoV_OC43','CoV_229E','CoV_HKU1','CoV_NL63') ] <- 'CoV'
  db$pathogen[db$pathogen %in% c('APZTD4A','S. pneumoniae_APZTD4A') ] <- 'S.pneumoniae'
  db$pathogen[db$pathogen %in% c('AI5IRK5','M. pneumoniae_AI5IRK5') ] <- 'M.pneumoniae'
  db$pathogen[db$pathogen %in% c('C. pneumoniae_AI1RW2H') ] <- 'C.pneumoniae'
  db$pathogen[db$pathogen %in% c('AP7DPVF','EnterovirusA_B 1_AP7DPVF') ] <- 'EV_pan'
  db$pathogen[db$pathogen %in% c('Enterovirus-D_APFVK4U','enterovirus-D_APFVK4U') ] <- 'EV_D68'
  db$pathogen[db$pathogen %in% c('AI20U8U') ] <- 'B.pertussis'
  db$pathogen[db$pathogen %in% c('APKA3DE') ] <- 'Mumps'
  
  # filter out nested PCR targets to retain high-level target only
  # Flu A
  keepTargetList <- unique(db$sample[db$pathogen %in% c("Flu_A_H1","Flu_A_H3")])
  dropTargetList <- unique(db$sample[db$pathogen %in% c("Flu_A_pan")])
  
  dropSampleList <- intersect(dropTargetList,keepTargetList)
  
  db <- db %>% filter( !(sample %in% dropSampleList & db$pathogen %in% c("Flu_A_pan")))
  
  # enterovirus
  keepTargetList <- unique(db$sample[db$pathogen %in% c("EV_D68")])
  dropTargetList <- unique(db$sample[db$pathogen %in% c("EV_pan")])
  
  dropSampleList <- intersect(dropTargetList,keepTargetList)
  
  db <- db %>% filter( !(sample %in% dropSampleList & db$pathogen %in% c("EV_pan")))
  
  # filter out controls
  db <- db %>% filter( !(pathogen %in% c('Hs04930436_g1','Ac00010014_a1')))
  
  # format flu_shot NA
  db$flu_shot[is.na(db$flu_shot)] <- 'unknown'
  db$flu_shot <- tolower(db$flu_shot)
  
  
  # run query
  # this logic will probably move to sql queries in the database instead of dplyr after....
    if(queryList$SELECT !="*"){
      
 
      #(Needed hack until higher-level shape labels are in database)
        if ( any( grepl('residence',queryList$SELECT$COLUMN) | grepl('work',queryList$SELECT$COLUMN) ) ){
          if (! any( grepl('regional_name',queryList$SELECT$COLUMN) | grepl('cra_name',queryList$SELECT$COLUMN) | grepl('neighbo',queryList$SELECT$COLUMN) ) ){
            shp = dbViewR::masterSpatialDB(shape_level = 'census_tract', source = 'wa_geojson')
          } else if (any( grepl('regional_name',queryList$SELECT$COLUMN) )){
            shp = dbViewR::masterSpatialDB(shape_level = 'census_tract', source = 'sfs_domain_geojson')
          } else {
            shp = dbViewR::masterSpatialDB(shape_level = 'census_tract', source = 'king_county_geojson')
          }
        
          # append higher-level spatial labels
          # this feature will eventually be in the database, but it's needed for now to index to pumas, cra_name, etc
          nestedVariables <- c('regional_name','cra_name','neighborhood_district_name','puma','city')
          
          for( COLUMN in nestedVariables){
            COLNAME <- paste0('residence_',COLUMN)
            if( ('residence_census_tract' %in% names(db))  & !(COLNAME %in% names(db)) & (COLNAME %in% names(shp))){
              db[[COLNAME]] <- as.character(shp[[COLNAME]][match(db$residence_census_tract,shp$residence_census_tract)])
            }
            COLNAME <- paste0('work_',COLUMN)
            if( ('work_census_tract' %in% names(db)) & !(COLNAME %in% names(db)) & (COLNAME %in% names(shp))){
              db[[COLNAME]] <- as.character(shp[[COLNAME]][match(db$work_census_tract,shp$work_census_tract)])
            }
          }
        }
      
      ## real flow starts here
        
      db <- db %>% dplyr::select(dplyr::one_of(queryList$SELECT$COLUMN))

      for(FILTER in which(grepl('WHERE',names(queryList)))){

        if( any(grepl('IN',names(queryList[[FILTER]])))){

          if(any(queryList[[FILTER]]$IN != 'all')){
            filter_criteria <- lazyeval::interp(~y %in% x, .values=list(y = as.name(queryList[[FILTER]]$COLUMN), x = queryList[[FILTER]]$IN))
            db <- db %>% dplyr::filter_(filter_criteria)
          }
        
        } else if( any(grepl('BETWEEN',names(queryList[[FILTER]])))){

          filter_criteria_low <- lazyeval::interp(~y >= x, .values=list(y = as.name(queryList[[FILTER]]$COLUMN), x = queryList[[FILTER]]$BETWEEN[1]))
          filter_criteria_high <- lazyeval::interp(~y <= x, .values=list(y = as.name(queryList[[FILTER]]$COLUMN), x = queryList[[FILTER]]$BETWEEN[2]))

          db <- db %>% dplyr::filter_(filter_criteria_low)  %>% dplyr::filter_(filter_criteria_high)

        }
      }
    
      if('GROUP_BY' %in% names(queryList)){
        db<- db %>% dplyr::group_by_(.dots=queryList$GROUP_BY$COLUMN)
      }
      
      if('SUMMARIZE' %in% names(queryList)){
  
        if (queryList$SUMMARIZE$IN != 'all'){
          summary_criteria <- lazyeval::interp(~sum(y %in% x), .values=list(y = as.name(queryList$SUMMARIZE$COLUMN), x = queryList$SUMMARIZE$IN))
        } else {
          summary_criteria <- lazyeval::interp(~n())  # must always output n and positive for downstream interpretation!
        }
          
        db <- db %>% dplyr::summarise_(n = lazyeval::interp(~n()), positive = summary_criteria) 
      }
  
      # "pathogen" column is required for incidenceMapR model definitions
      # this logic should use a key-value config file instead of being hard-coded
        if(!('pathogen' %in% queryList$GROUP_BY$COLUMN)){
          if( 'pathogen' %in% queryList$WHERE$COLUMN){
            
            if (all(grepl('flu',queryList$WHERE$IN,ignore.case = TRUE)) & length(queryList$WHERE$IN)>1){
              db$pathogen <- 'flu_positive'
            } else if (all(grepl('rsv',queryList$WHERE$IN,ignore.case = TRUE)) & length(queryList$WHERE$IN)>1) {
              db$pathogen <- 'rsv_positive'
            } else if (!any(grepl('flu',queryList$WHERE$IN,ignore.case = TRUE) |
                            grepl('not_yet_tested',queryList$WHERE$IN,ignore.case = TRUE) |
                            grepl('measles',queryList$WHERE$IN,ignore.case = TRUE) ) & length(queryList$WHERE$IN)>1
                       ){
              db$pathogen <- 'flu_negative'
            } else {
              db$pathogen <- paste(queryList$WHERE$IN['pathogen' %in% queryList$WHERE$COLUMN],collapse='-')
            }
            
            
          } else if ('pathogen' %in% queryList$SUMMARIZE$COLUMN) {
            db$pathogen <- paste(queryList$SUMMARIZE$IN['pathogen' %in% queryList$SUMMARIZE$COLUMN],collapse='-')
          } else {
            db$pathogen<-'all' 
          }
        }
    
    }
  
  
  # type harmonization
    for( COLUMN in names(db)[names(db) %in% c('residence_census_tract','residence_cra_name','residence_puma','residence_neighborhood_district_name','residence_city','residence_regional_name',
                                              'work_census_tract','work_cra_name','work_puma','work_neighborhood_district_name','work_city','work_regional_name')]){
      db[[COLUMN]] <- as.character(db[[COLUMN]])
    }

    if('encountered_week' %in% names(db)){
      db$encountered_week <- factor(db$encountered_week, levels=sort(unique(db$encountered_week)), ordered = TRUE)
    }
  
  # drop rows with NA since incidenceMapR (INLA) will ignore them anyway
    if(na.rm){
      
      # fixes #Error in charToDate(x): character string is not in a standard unambiguous format
      dateIdx<- (sapply(db,class)=="Date")
      db[,dateIdx] <- as.character(db[,dateIdx])
      
      db <- db %>% replace(.=='NA', NA) %>% tidyr::drop_na()
      
      # restore type
      if (any(dateIdx)){
        db[,dateIdx] <- as.Date(db[,dateIdx])
      }
    
    }
    
  summarizedData <- list(observedData = db,queryList = c(queryList))

  return(summarizedData)
}



