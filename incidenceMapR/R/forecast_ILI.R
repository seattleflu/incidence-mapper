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
#' @import reshape2
#' @import plyr
#' @import dplyr
#' @import gridExtra
#' @import tibble
#'
#' @export
#' @examples
#' plot ILI data and forecast
#'   plotforecasts(dat_WAILI, "ILI", "ARIMA forecast ILI", "ARIMA_fc_ILI")
 
plotforecasts<- function(dat, fitted_param, ggtitlestr, outpath){
  
  fc <- dat %>% select(matches("fc_|upper_|lower_")) %>% colnames() %>% c(fitted_param)
  bounds_only <- dat %>% select(matches("upper_|lower_")) %>% colnames()
  data_toplot<- select(dat, "encountered_week", fc) %>% melt(id.vars=c("encountered_week", bounds_only))
  
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
  
  ggsave(paste(outpath, ".pdf", sep=""), plot=p, width=7, height=5, units="in")
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
#' @import MMWRweek
#' @import reshape2
#' @import plyr
#' @import dplyr
#' @import lubridate
#' @import forecast
#' @import tibble
#' @import boot
#' 
#' @return ILI_estimate - dataframe with point estimate and 95% CI upper and lower for num_future_weeks
calcfc<- function(dat, num_prior_data=52, num_future_data = 1){
  
  
  #take logit of ILI data
  dat$ILI<-logit(dat$ILI)
  
  #calculate the extra_dates for the forecast_weeks (TODO fix when it crosses to year 2020)
  last_week<-dat$week_start[nrow(dat)]
  extra_dates_day<-last_week + weeks(1:num_future_data)
  extra_dates <- paste(MMWRweek(extra_dates_day)[,1],"-W",
                                sprintf("%02d",MMWRweek(extra_dates_day)[,2]), sep="") %>% data.frame()#add CDC week number
  
  colnames(extra_dates)<-c("encountered_week")  
  
  
  #set up data frame to store forecasts
  fc_master <- slice(dat, 1:n()) %>% select(encountered_week)
  end_fc_time<-nrow(dat)-1
  if(!is.null(extra_dates))
  {
    fc_master<-rbind(fc_master, extra_dates[1])
    end_fc_time <- end_fc_time + 1
  }
  
  #column names to add 
  fc_master<-add_column(fc_master, fc = NA, upper = NA, lower=NA)
  
  #run forecasting for each week + future forecast
  for (current_time in (num_prior_data+1):end_fc_time){
    
    #start forecasting with 1 year of data 
    dat_s <- slice(dat, (current_time-num_prior_data):current_time)
    
    #fit model and forecast 1 wk ahead for all weeks
    fit <- auto.arima(dat_s$ILI) 
    fc <- forecast(fit, h=1) #forecast 1 wk ahead
    if (current_time == nrow(dat) && !is.null(extra_dates)){
      fc <- forecast(fit, h=num_future_data) #forecast more weeks ahead for last time point
    }
    
    
    #merge fitted value with original data for plotting 
    fc_master_timerange <- (current_time+1):(current_time+1)
    if(current_time == nrow(dat) && !is.null(extra_dates)){
      fc_master_timerange <- (current_time+1):(current_time+num_future_data)
    }
    
    #store the forecast values in the master list
    fc_master$fc[fc_master_timerange] <-data.matrix(fc$mean)
    fc_master$upper[fc_master_timerange] <- data.matrix(fc$upper[,2])
    fc_master$lower[fc_master_timerange] <- data.matrix(fc$lower[,2])
    
    #inverse logit the forecast and the range 
    fc_master$fc[fc_master_timerange] <-inv.logit(fc_master$fc[fc_master_timerange])
    fc_master$upper[fc_master_timerange] <-inv.logit(fc_master$upper[fc_master_timerange])
    fc_master$lower[fc_master_timerange] <-inv.logit(fc_master$lower[fc_master_timerange])
    
    
  }
  #rename columns
  colnames(fc_master) <- c("week", "fc_ILI", "upper_95", "lower_95")
  
  return(fc_master)
}




#' 
#' 
#' @param num_prior_data number of previous data points to use to fit the ARIMA model
#' @param num_future_data number of weeks forward from the latest observation in the CDC database to make ILI estimates
#' @return ILI_estimate - dataframe with point estimate and 95% CI upper and lower for num_future_data weeks
#' 
#' @import cdcfluview
#' @import MMWRweek
#' @import boot
#'
#' @export
#' @examples
#' return ILI estimate 2 weeks ahead
#'    ILI_forecast <- forecast_ILI(num_prior_data=52, num_future_weeks = 2)
#'    
#'    
forecast_ILI <- function(num_prior_data=52, num_future_weeks = 2){
    
  #load the ILI Data
  ILI<- ilinet(region=c("state"), years = c(2017, 2018, 2019)) #all states for 2018-19 flu season

  #get only WA rows and unweighted ILI columns 
  WA_ILI<- ILI %>% filter(grepl('Washington', region)) %>% select(unweighted_ili, week_start) %>% dplyr::rename(ILI = unweighted_ili) 
  #add CDC week number
  WA_ILI$encountered_week <- paste(MMWRweek(WA_ILI$week_start)[,1],"-W",
                                sprintf("%02d",MMWRweek(WA_ILI$week_start)[,2]), sep="")#add CDC week number 
  
  #convert ILI to fraction btwn 0 and 1
  WA_ILI$ILI <- WA_ILI$ILI/100
  
  #forecast WA ILI 1 wk ahead for weeks where ILI data is known using ARIMA
  #combine ILI data with forecast ILI data
  dat_WAILI<-left_join(calcfc(WA_ILI, num_prior_data=num_prior_data, num_future_data=num_future_weeks),
                       WA_ILI, 
                      by = c("week" = "encountered_week" )) 
  dat_WAILI<-dplyr::rename(dat_WAILI, encountered_week=week)
  
  #optional save a plot of the ILI fit and  forecast
  plotforecasts(dat_WAILI, "ILI", "ARIMA forecast ILI", "ARIMA_fc_ILI")
  
  return(dat_WAILI)
}