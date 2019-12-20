#' expandDB: function to expand database to fill in unobserved elements for prediction/interpolation
#'
#' @param db tibble with valid column names for INLA model
#' @param shp sf shapefile object with reference geodata
#' @return db tibble with filled in NaNs or Zeros as required
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom lazyeval interp
#'
#' @export
#' @examples
#'
expandDB <- function( db = dbViewR::selectFromDB(),
                      shp = dbViewR::masterSpatialDB(shape_level = 'census_tract', source = 'king_county_geojson'),
                      currentWeek = NULL, 
                      minWeek = NULL){
  
  
  # list of valid column data for expanding and joining
    
    validColumnData <- list()
    
    # encountered_week
    if ("encountered_week" %in% names(db$observedData)){
      
      weeks <- unique(sort(db$observedData$encountered_week))
      # mostRecentWeekWithData <- as.character(weeks[length(weeks)])
      
      if(is.null(minWeek)){
        minYear <- as.numeric(gsub('-W[0-9]{2}','',weeks[1]))
        minWeek <- as.numeric(gsub('[0-9]{4}-W','',weeks[1] ))
      } else {
        minYear <- as.numeric(gsub('-W[0-9]{2}','',minWeek))
      }

      if(is.null(currentWeek)){
        maxYear <- as.numeric(gsub('-W[0-9]{2}','',weeks[length(weeks)]))
        maxWeek <- as.numeric(gsub('[0-9]{4}-W','',weeks[length(weeks)] )) + (maxYear-minYear)*52  + 4 # 4 week look-ahead
      } else {
        maxYear <- as.numeric(gsub('-W[0-9]{2}','',currentWeek)) 
        maxWeek <- as.numeric(gsub('[0-9]{4}-W','',currentWeek )) + (maxYear-minYear)*52 #+ 4 # 4 week look-ahead
      }
      
      weeks <- 1+( (seq(minWeek,maxWeek,by=1)-1) %% 52)
      yearBreaks <- c(0,which(diff(weeks)<1), length(weeks))
      years=c()
      for (k in 2:length(yearBreaks)){
          years <- c(years, rep(minYear+(k-2), yearBreaks[k]-yearBreaks[k-1]  ))
      }    
      validColumnData$encountered_week <- paste(years,'-W',sprintf("%02d",weeks),sep='')
      validColumnData$time_row <- 1:length(validColumnData$encountered_week)
    }

    # age
      validColumnData$age = seq(0,90,by=1)
      validColumnData$age_bin_fine_lower = c(0,0.08,0.5,1,seq(5,90,by=5))
      validColumnData$age_bin_fine_upper = c(0.08,0.5,1,seq(5,90,by=5))
      validColumnData$age_bin_coarse_lower = c(0,0.5,5,18,65)
      validColumnData$age_bin_coarse_upper = c(0.5,5,18,65,90)
      
    
    # age bin
    # if(any(grepl('age',names(db$observedData)))) {
    #   validColumnData$age_bin <- seq(0,90,by=5)
    #   validColumnData$age_row <- 1:length(validColumnData$age_bin)
    # }

    # geography
      validColumnData$residence_census_tract = shp$residence_census_tract
      validColumnData$residence_cra_name = sort(unique(shp$residence_cra_name)) 
      validColumnData$residence_neighborhood_district_name = sort(unique(shp$residence_neighborhood_district_name))
      validColumnData$residence_puma = sort(unique(shp$residence_puma))
      validColumnData$residence_city = sort(unique(shp$residence_city))
      validColumnData$residence_regional_name = sort(unique(shp$residence_regional_name))
      validColumnData$work_census_tract = shp$work_census_tract
      validColumnData$work_cra_name = sort(unique(shp$work_cra_name))
      validColumnData$work_neighborhood_district_name = sort(unique(shp$work_neighborhood_district_name))
      validColumnData$work_puma = sort(unique(shp$work_puma))
      validColumnData$work_city = sort(unique(shp$work_city))
      validColumnData$work_regional_name = sort(unique(shp$work_regional_name))
      
      
      # NA handling
      validColumnData$residence_cra_name = validColumnData$residence_cra_name[validColumnData$residence_cra_name!='NA']
      validColumnData$residence_neighborhood_district_name = validColumnData$residence_neighborhood_district_name[validColumnData$residence_neighborhood_district_name!='NA']
      validColumnData$residence_city = validColumnData$residence_city[validColumnData$residence_city!='NA']
      validColumnData$residence_regional_name = validColumnData$residence_regional_name[validColumnData$residence_regional_name!='NA']
      validColumnData$work_cra_name = validColumnData$work_cra_name[validColumnData$work_cra_name!='NA']
      validColumnData$work_neighborhood_district_name = validColumnData$work_neighborhood_district_name[validColumnData$work_neighborhood_district_name!='NA']
      validColumnData$work_city = validColumnData$work_city[validColumnData$work_city!='NA']
      validColumnData$work_regional_name = validColumnData$work_regional_name[validColumnData$work_regional_name!='NA']
      
      
    # factors (these don't get interpolated by the models, so we only want the valid levels for the dataset at hand)
      factorNames <- names(db$observedData)[ !( (names(db$observedData) %in% c('age','n','positive')) | 
                                                grepl('residence_',names(db$observedData)) | 
                                                grepl('work_',names(db$observedData)) |
                                                grepl('encounter',names(db$observedData))  )]
      for ( COLUMN in factorNames){ 
        validColumnData[[COLUMN]] <- sort(unique(db$observedData[[COLUMN]]))
      }

  
  # don't expand on nested shape variables
    if (any(grepl('census_tract',names(db$observedData)))){
      nestedVariables <- c('residence_cra_name','residence_neighborhood_district_name','residence_puma','residence_city', 'residence_regional_name',
                           'work_cra_name','work_name','work_puma','work_city', 'work_regional_name')
    } else {
      nestedVariables <- c()
    }

  # expand.grid for non-nested variables
    colIdx <- ( names(validColumnData) %in% names(db$observedData) ) &  !( names(validColumnData) %in% nestedVariables) 
    tmp<-expand.grid(validColumnData[colIdx],stringsAsFactors = FALSE)
    
  # join
  if('encountered_week' %in% names(db$observedData)){
    db$observedData$encountered_week <- as.character(db$observedData$encountered_week)  
  }
  db$observedData <- dplyr::left_join(tmp,db$observedData, by=names(validColumnData)[colIdx])

  # order weeks  
  if('encountered_week' %in% names(db$observedData)){
    db$observedData$encountered_week <- factor(db$observedData$encountered_week, levels=sort(unique(db$observedData$encountered_week)), ordered = TRUE)
  }
  
  # sample size as zero instead of NaN
  if ("n" %in% names(db$observedData)){
    db$observedData$n[is.na(db$observedData$n)]<-0

    # deal with NAs based on when sites turned on
    if (any(grepl('site', names(db$observedData), ignore.case = TRUE))){
      COLUMN <- names(db$observedData)[grepl('site', names(db$observedData), ignore.case = TRUE)]
      if (length(COLUMN)>1){
        if( any(COLUMN=='site')){
          COLUMN <- COLUMN[COLUMN=='site']
        } else {
          COLUMN <- COLUMN[COLUMN=='site_type']
        }
      }
      siteValues <- unique(db$observedData[[COLUMN]])
      for (SITE in siteValues){
        siteIdx <- (db$observedData[[COLUMN]]==SITE) 
        # when time is a variable, only set positive to 0 when n==0 if we're sure the site was operational and providing data
        # ideally, this would pull from a data source with gold-standard dates
        if('encountered_week' %in% names(db$observedData)){
          minWeek=as.character(min(db$observedData$encountered_week[siteIdx & !is.na(db$observedData$positive)]))
          maxWeek=as.character(max(db$observedData$encountered_week[siteIdx & !is.na(db$observedData$positive)]))
          db$observedData$positive[siteIdx & is.na(db$observedData$positive) &
                                   (db$observedData$encountered_week>=minWeek) &
                                   (db$observedData$encountered_week <= maxWeek)]<-0
        } else {
          idx <- is.na(db$observedData$positive)
          if(all(db$observedData$positive[siteIdx & !idx]==db$observedData$n[siteIdx & !idx])){
            db$observedData$positive[siteIdx & idx]<-0
          }
        }
      }
    }
  }
  
  
  
  # nested variables
    colIdx <- which(( names(validColumnData) %in% names(db$observedData) ) & ( names(validColumnData) %in% nestedVariables ) )
    for( k in colIdx){
      colName <- names(validColumnData)[k]
      db$observedData[[colName]] <- shp[[colName]][match(db$observedData[[nestedVariables[nestedVariables == colName]]],as.character(shp[[nestedVariables[nestedVariables == colName]]]))]
    }
  
  # row indices for INLA
  if(any(grepl('encountered_week',names(db$observedData)))){
    db$observedData$time_row <- validColumnData$time_row[match(db$observedData$encountered_week,validColumnData$encountered_week)]
  }

  # if(any(grepl('age_bin',names(db$observedData)))){
  #   db$observedData$age_row <- validColumnData$age_row[match(db$observedData$age_bin,validColumnData$age_bin)]
  # }
  if(any(grepl('age_range_fine_lower',names(db$observedData)))){
    db$observedData$age_row <- match(db$observedData$age_range_fine_lower,validColumnData$age_range_fine_lower)
  } else if (any(grepl('age_range_fine_upper',names(db$observedData)))){
    db$observedData$age_row <- match(db$observedData$age_range_fine_upper,validColumnData$age_range_fine_upper)
  } else if (any(grepl('age_range_coarse_lower',names(db$observedData)))){
    db$observedData$age_row <- match(db$observedData$age_range_coarse_lower,validColumnData$age_range_coarse_lower)
  } else if (any(grepl('age_range_coarse_upper',names(db$observedData)))){
    db$observedData$age_row <- match(db$observedData$age_range_coarse_upper,validColumnData$age_range_coarse_upper)
  } 
  
  return(db)
}
