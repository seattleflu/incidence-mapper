# explore forecast_ILI.R

currentWeek <-  paste(lubridate::isoyear(Sys.time()) ,'-W',lubridate::isoweek(Sys.time()),sep='')

forecast_dat <- forecast_ILI(currentWeek = currentWeek)

#optional save a plot of the ILI fit and  forecast
plotforecasts(forecast_dat, fitted_param = "ILI", ggtitlestr = "ARIMA forecast ILI", outpath = '/home/rstudio/seattle_flu/model_diagnostic_plots/forecast_WA_ILI')
