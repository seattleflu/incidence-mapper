#' appendILIDataFc: function for adding ILI Data and short fc to db$observedData 
#'
#' @param db object from dbViewer with observedData tibble and query
#' @param currentWeek current week
#' 
#' @import dplyr
#' 
#' @return db with added ILI Data + fc tibble
#' 
#' @export

appendILIDataFc <- function(db, currentWeek){
  
  ILI_forecast_dat <- forecast_ILI(currentWeek = currentWeek)
  
  #ILI_forecast_dat includes the actual past ILI values and the forecasted data for the past year up to currentWeek
  
  #reduce the ILI_forecast_dat so that it returns actual ILI values if available, and only the forecast data for the recent weeks where actual ILI unavailable
  ILI_forecast_dat$ILI_all <- dplyr::if_else(is.na(ILI_forecast_dat$ILI), ILI_forecast_dat$fc_ILI, ILI_forecast_dat$ILI)
  ILI_forecast_dat <- ILI_forecast_dat %>% dplyr::select(encountered_week, ILI_all) %>% dplyr::rename(ILI = ILI_all)
  
  #add to the db
  db$observedData <- dplyr::left_join(db$observedData, ILI_forecast_dat, by = c("encountered_week"="encountered_week"))
  
  return(db)
}
