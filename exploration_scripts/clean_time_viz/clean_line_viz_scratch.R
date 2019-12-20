# SFS pretty time series plot

# script assumes a latentField model named "model" is loaded in your workspace!

library(gglot2)
library(dplyr)
library(ggthemes)

GEO='residence_regional_name'
PATHOGEN='flu'
shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = 'sfs_domain_geojson')

plotDat <- model$latentField %>% select(encountered_week,residence_regional_name,modeled_intensity_median,modeled_intensity_sd) %>%
  group_by(encountered_week) %>% summarize(region_intensity_median = mean(modeled_intensity_median))

weeks <- paste('2019-W',sprintf('%02d',isoweek(c('2019-01-01','2019-02-01','2019-03-01',
                   '2019-04-01','2019-05-01','2019-06-01',
                   '2019-07-01','2019-08-01','2019-09-01',
                   '2019-10-01','2019-11-01','2019-12-01'))), sep='')
weekLabels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',
             'Sep','Oct','Nov','Dec')

plotLabel <- data.frame(encountered_week='2019-W40',
                        region_intensity_median=0.21,
                        label='Seattle Flu Intensity')

p<-ggplot(plotDat) + 
  geom_line(aes(x=encountered_week,y=region_intensity_median, color=GEO,group =GEO), 
            size=1.25, color='#095da8') + 
  guides(color=FALSE, fill=FALSE) + # geom_rangeframe() +  
  theme_hc(base_family='sans', base_size = 12) +
  scale_x_discrete(breaks=weeks[seq(1,12,by=2)], labels=weekLabels[seq(1,12,by=2)]) +
  xlab('2019')+
  ylab('') +
  geom_text(data=plotLabel,aes(x=encountered_week,y=region_intensity_median, label=label)) +
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2), 
                     labels=c('none','low','moderate','high','peak'))
p

fname <- paste('/home/rstudio/seattle_flu/exploration_scripts/clean_time_viz/',paste('inla_latent',PATHOGEN,SOURCE,GEO,'encountered_week',sep='-'),'.pdf',sep='')
pdf(file = fname,width = 5, height = 3.5)
print(p)
dev.off()
