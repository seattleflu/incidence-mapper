# plot 3 maps

library(dbViewR)
library(incidenceMapR)
library(lubridate)
library(magrittr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

# get legend: https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

fluPathogens <- c('Flu_A_H1','Flu_A_H3','Flu_A_pan','Flu_B_pan','Flu_C_pan')

# setting
SOURCE='sfs_domain_geojson'
GEO = 'residence_regional_name'


# SOURCE='seattle_geojson'
# GEO = 'residence_neighborhood_district_name'

factors   <- c('site_type')#,'age_range_fine_upper')
shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = SOURCE)

# SOURCE='seattle_geojson'
# GEO = 'residence_neighborhood_district_name'

# PATHOGEN='Flu_A_H1'
# PATHOGEN='Flu_A_H3'
# PATHOGEN <- fluPathogens
# PATHOGEN='RSVA'
# PATHOGEN='RSVB'
# PATHOGEN=c('RSVA','RSVB')
PATHOGEN='all'
# PATHOGEN='CoV'

# db <- selectFromDB(queryIn= list(SELECT  =c("*")), source = SRC)
# pathogens <- db$observedData %>% group_by(pathogen) %>% summarize(n = n()) %>% arrange(desc(n))
# PATHOGEN = setdiff(pathogens$pathogen,c(fluPathogens,'not_yet_tested','measles','Measles'))



# build model
queryIn <- list(
  SELECT   =list(COLUMN=c('pathogen', factors, GEO,'encountered_week')),
  WHERE    =list(COLUMN='pathogen', IN=PATHOGEN),
  GROUP_BY =list(COLUMN=c(factors,GEO,"encountered_week")),
  SUMMARIZE=list(COLUMN='pathogen', IN= PATHOGEN)
)

db <- expandDB(selectFromDB(  queryIn, source='production', na.rm=TRUE ), shp=shp)
db <- appendCatchmentModel(db,shp=shp, source='production', na.rm=TRUE )

modelDefinition <- latentFieldModel(db=db, shp=shp)
model <- modelTrainR(modelDefinition)
summary(model$inla)

## line plot
ggplot(model$latentField) + 
  geom_line(aes_string(x='encountered_week',y="modeled_intensity_median", color=GEO,group =GEO)) + 
  # geom_ribbon(aes_string(x='encountered_week',ymin="modeled_intensity_lower_95_CI", ymax="modeled_intensity_upper_95_CI", fill=GEO,group =GEO),alpha=0.1) +
  guides(color=FALSE, fill=FALSE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



## set up map environment
plotDat <- right_join(model$latentField %>% group_by_(.dots =c(GEO,'encountered_week')) %>% summarise(modeled_intensity_median = median(modeled_intensity_median)),shp, by=GEO)
bbox<-sf::st_bbox(shp$geometry)

# #h1
# plotDat$modeled_intensity_median <- pmin(0.10,plotDat$modeled_intensity_median)/0.10
# 
# #h3
# plotDat$modeled_intensity_median <- pmin(0.12,plotDat$modeled_intensity_median)/0.12

#all
plotDat$modeled_intensity_median <- pmin(1,plotDat$modeled_intensity_median)/1

mapSettings <- ggplot() + xlim(-122.45,-122.05) + ylim(47.45,47.8) +
  theme_bw() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.grid.major=element_line(colour="transparent"), panel.border = element_blank()) +
  geom_sf(data=shp,size=0.1,aes(fill=NaN))

colorLimits=c(0,1)
colorBreaks=seq(0,5,by=1)/5
# colorLimits<-c(min(plotDat$modeled_intensity_median,na.rm=TRUE),max(plotDat$modeled_intensity_median,na.rm=TRUE))
# colorBreaks<-unique(c(min(colorLimits),round(1e4*seq(sqrt(min(colorLimits)),sqrt(max(colorLimits)), length.out = 5)^2)/1e4, max(colorLimits)))

## choose weeks to show
week <- c('2019-W02','2019-W05','2019-W08', '2019-W11','2019-W14')
p <- list()
for(k in 1:length(week)){
  p[[k]] <- mapSettings + geom_sf(data=plotDat %>% filter(encountered_week == week[k]),size=0, aes(fill=modeled_intensity_median))  +
    guides(fill=guide_legend(title="modeled_intensity")) +
    viridis::scale_fill_viridis(na.value="transparent",breaks=colorBreaks,limits=colorLimits)+
    ggtitle(week[k])
  if(k<length(week)){
    p[[k]] <- p[[k]] + theme(legend.position='none')
  }
}
legend <- g_legend(p[[k]]) 
p[[k]] <- p[[k]] + theme(legend.position='none')

png(filename=paste('Year1',paste(PATHOGEN,collapse='-',sep=''),'fourMaps.png',sep='-'),res=600, units='in', width=6.5,height=4)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]], p[[5]], legend,nrow=1)
dev.off()
