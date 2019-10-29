#' @title Forecast ILI using ARIMA  
#' September 2019
#' Helper functions first, actual function is forecast_ILI below


#' plotforecasts: internal function for saving plot of ILI and ILI forecast with 95%CI
#'
#' @param dat ILI data to plot
#' @param fitted_param ILI
#' @param ggtitlestr title of plot
#' @param outpath where to save
#'  
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import gridExtra
#' @import tibble
#'
#' @export
#' @examples
#' plot ILI data and forecast
#'   plotforecasts(dat_WAILI, "ILI", "ARIMA forecast ILI", "ARIMA_fc_ILI")
 
plotforecasts<- function(dat, fitted_param , ggtitlestr = NULL , outpath = NULL){
  
  fc <- dat %>% select(matches("fc_|upper_|lower_")) %>% colnames() %>% c(fitted_param)
  bounds_only <- dat %>% select(matches("upper_|lower_")) %>% colnames()
  data_toplot<- select(dat, "encountered_week", fc) %>% tidyr::gather(variable, value, -encountered_week, -lower_95, -upper_95)
  
  #everything in its own facet except the ILI which is overlaid on all facets
  all_otherdata <- data_toplot[-grep(paste("^", fitted_param, "$", sep=""), data_toplot$variable), ]
  p<-ggplot(aes(x=encountered_week, y=value), data=all_otherdata)+
    geom_point(aes(color="fc")) + geom_line(aes(group=variable,  color="fc")) + facet_wrap(~variable, scales="free_y") + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    geom_point(aes(group=1,  color="n"),data=transform(data_toplot[grep(paste("^", fitted_param, "$", sep=""), data_toplot$variable), ], variable=NULL)) +
    geom_line(aes(group=1,  color="n"),data=transform(data_toplot[grep(paste("^", fitted_param, "$", sep=""), data_toplot$variable), ], variable=NULL)) +
    scale_color_manual(labels = c("fc", fitted_param), values = c("red", "gray")) + 
    ylab(fitted_param) + ggtitle(ggtitlestr) + 
    geom_ribbon(aes(ymin=lower_95, ymax=upper_95, fill="95% CI", group=variable),  alpha=0.2)
 
  
  if(!is.null(outpath)){
    ggsave(paste(outpath, ".png", sep=""), plot=p, width=7, height=5, units="in")
  } else{
    print(p) 
  }
}


#' calcfc: internal function for calculating the ILI fit with no regressors
#' takes logit of ILI data for modeling, inverse logit is done before returning ILI estimates
#' uses num_prior_weeks - number of points to fit for each data point in the 2018-2019 season
#' forward estimates for num_future_weeks - number of weeks to make future ILI estimates
#' returns ILI forward estimates for num_future_weeks- point estimate and the 95% CI upper and lower
#'
#' @param dat dataframe with week, and ILI as columns
#' @param num_prior_data number of previous data points (weeks) to use to fit the ARIMA model
#' @param num_future_data number of weeks forward to make ILI estimates
#' 
#' @import lubridate
#' @import reshape2
#' @import dplyr
#' @import lubridate
#' @import forecast
#' @import tibble
#' @import boot
#' 
#' @return ILI_estimate - dataframe with point estimate and 95% CI upper and lower for num_future_weeks
calcfc<- function(dat, num_prior_data=52, currentWeek=currentWeek){
  
  
  #take logit of ILI data
  dat$ILI<-boot::logit(dat$ILI)
  
  # regularize weeks with 0 reported ILI
  idx <- is.finite(dat$ILI)
  dILI <- median(abs(diff(dat$ILI[idx])))
  for(k in which(!idx)){
    if (k>1){
      dat$ILI[k] <- dat$ILI[k-1]-dILI
    } else {
      dat$ILI[k] <- dat$ILI[k+1]-dILI
    }
  }
  
  
  ## calculate forecast date range
  currentWeek <- paste(lubridate::isoyear(Sys.time()) ,'-W',lubridate::isoweek(Sys.time()),sep='')
  last_ILI_week <- max(dat$encountered_week)
  
  minYear <- year<-as.numeric(gsub('-W[0-9]{2}','',last_ILI_week))
  maxYear <- year<-as.numeric(gsub('-W[0-9]{2}','',currentWeek))
  minWeek <- as.numeric(gsub('[0-9]{4}-W','',last_ILI_week ))
  maxWeek <- as.numeric(gsub('[0-9]{4}-W','',currentWeek )) + (maxYear-minYear)*52 
  
  weeks <- 1+( (seq(minWeek+1,maxWeek,by=1)-1) %% 52)
  yearBreaks <- c(0,which(diff(weeks)<1), length(weeks))
  years=c()
  for (k in 2:length(yearBreaks)){
    years <- c(years, rep(minYear+(k-2), yearBreaks[k]-yearBreaks[k-1]  ))
  }
  
  forecast_weeks <- paste(years,'-W',sprintf("%02d",weeks),sep='')

  #set up data frame to store forecasts
  fc_master <- slice(dat, 1:n()) %>% select(encountered_week) 
  end_fc_time <- nrow(fc_master)
  
  if(!is.null(forecast_weeks)){
    fc_master <- fc_master %>% tibble::add_row(encountered_week = forecast_weeks)
  }
  fc_master <- fc_master %>% tibble::add_column(fc = NA, upper = NA, lower=NA)
  
  num_future_data <- nrow(fc_master) - end_fc_time

  
  #run forecasting for each week + future forecast
  for (current_time in (num_prior_data+1):end_fc_time){
    
    #start forecasting with 1 year of data 
    dat_s <- slice(dat, (current_time-num_prior_data):current_time)
    
    #fit model and forecast 1 wk ahead for all weeks
    fit <- auto.arima(dat_s$ILI) 
    fc <- forecast(fit, h=1) #forecast 1 wk ahead
    if (current_time == nrow(dat) && !is.null(forecast_weeks)){
      fc <- forecast(fit, h=num_future_data) #forecast more weeks ahead for last time point
    }
    
    
    #merge fitted value with original data for plotting 
    fc_master_timerange <- (current_time+1):(current_time+1)
    if(current_time == nrow(dat) && !is.null(forecast_weeks)){
      fc_master_timerange <- (current_time+1):(current_time+num_future_data)
    }
    
    #store the forecast values in the master list
    fc_master$fc[fc_master_timerange] <-data.matrix(fc$mean)
    fc_master$upper[fc_master_timerange] <- data.matrix(fc$upper[,2])
    fc_master$lower[fc_master_timerange] <- data.matrix(fc$lower[,2])
    
    #inverse logit the forecast and the range 
    fc_master$fc[fc_master_timerange] <-boot::inv.logit(fc_master$fc[fc_master_timerange])
    fc_master$upper[fc_master_timerange] <-boot::inv.logit(fc_master$upper[fc_master_timerange])
    fc_master$lower[fc_master_timerange] <-boot::inv.logit(fc_master$lower[fc_master_timerange])
    
    
  }
  #rename columns
  colnames(fc_master) <- c("encountered_week", "fc_ILI", "upper_95", "lower_95")
  
  return(fc_master)
}




#' forecast_ILI: function to extend CDC ILI data to desired week via ARIMA forecast
#' 
#' @param num_prior_data number of previous data points to use to fit the ARIMA model
#' @param num_future_data number of weeks forward from the latest observation in the CDC database to make ILI estimates
#' @return ILI_estimate - dataframe with point estimate and 95% CI upper and lower for num_future_data weeks
#' 
#' @import cdcfluview
#' @import boot
#' @import magrittr
#' @import dplyr
#' @import lubridate
#'
#' @export
#' @examples
#' return ILI estimate 2 weeks ahead
#'    ILI_forecast <- forecast_ILI(num_prior_data=52, num_future_weeks = 2)
#'    
#'    
forecast_ILI <- function(num_prior_data=52, currentWeek = NULL){
    
  # current week being run
  if( is.null(currentWeek)){
    currentWeek <- paste(lubridate::isoyear(Sys.time()) ,'-W',lubridate::isoweek(Sys.time()),sep='')
  }
  
  #load the ILI Data
  ILI<- cdcfluview::ilinet(region=c("state"), years = 2017:(lubridate::isoyear(Sys.time())+1)) #all states through this season

  #get only WA rows and unweighted ILI columns 
  WA_ILI<- ILI %>% filter(grepl('Washington', region)) %>% select(unweighted_ili, week_start) %>% dplyr::rename(ILI = unweighted_ili)
  
  # add isoweek as ordered factor
  WA_ILI$encountered_week <- paste(lubridate::isoyear(WA_ILI$week_start) ,'-W',sprintf("%02d",lubridate::isoweek(WA_ILI$week_start)),sep='')
  WA_ILI$encountered_week <- factor(WA_ILI$encountered_week, levels=sort(unique(WA_ILI$encountered_week)), ordered = TRUE)
  
  #convert ILI to fraction btwn 0 and 1
  WA_ILI$ILI <- WA_ILI$ILI/100
  

  #forecast WA ILI 1 wk ahead for weeks where ILI data is known using ARIMA
  #combine ILI data with forecast ILI data
  dat_WAILI<-left_join(calcfc(WA_ILI, num_prior_data=num_prior_data, currentWeek=currentWeek),
                       WA_ILI,by = 'encountered_week') 

  return(dat_WAILI)
}