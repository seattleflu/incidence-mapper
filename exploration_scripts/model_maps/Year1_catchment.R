# summarize catchments from year 1

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

# setting
SOURCE='seattle_geojson'
GEO = 'residence_neighborhood_district_name'


# build model
shp <- masterSpatialDB(shape_level = gsub('residence_','',GEO), source = SOURCE)
queryIn <- list(
  SELECT   =list(COLUMN=c('pathogen', GEO,'site_type')),
  GROUP_BY =list(COLUMN=c(GEO,'site_type')),
  SUMMARIZE=list(COLUMN='pathogen', IN= 'all')
)

db <- expandDB(selectFromDB(  queryIn, source='production', na.rm=TRUE ), shp=shp)


## set up map environment
plotDat <- right_join(db$observedData,shp, by=GEO)
plotDat <- plotDat %>% group_by(site_type) %>% mutate(fraction = n/max(n))

bbox<-sf::st_bbox(shp$geometry)

mapSettings <- ggplot() + #xlim(c(min(122.5, -bbox[1]),max(121.7,-bbox[3]))) + ylim(c(max(47.17,bbox[2]),min(47.76,bbox[4]))) +
  theme_bw() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.grid.major=element_line(colour="transparent"), panel.border = element_blank()) +
  geom_sf(data=shp,size=0.1,aes(fill=NaN)) 


# colorLimits<-c(min(plotDat$n,na.rm=TRUE),max(plotDat$n,na.rm=TRUE))
# colorBreaks<-unique(c(min(colorLimits),round(1e0*seq(sqrt(min(colorLimits)),sqrt(max(colorLimits)), length.out = 5)^2)/1e0, max(colorLimits)))

p <- mapSettings + geom_sf(data=plotDat ,size=0, aes(fill=fraction, group='site_type'))  +
  guides(fill=guide_legend(title="Year 1\n relative recruitment")) +
  viridis::scale_fill_viridis(na.value="transparent") + #,breaks=colorBreaks,limits=colorLimits)+ 
  facet_wrap('site_type')
p

png(filename=paste('Year1_catchment_map_by_site_type.png',sep='-'),res=600, units='in', width=10,height=6)
p
dev.off()


plotDat2 <- plotDat %>% group_by(!!rlang::parse_expr(GEO)) %>% summarize(n=sum(n)) %>% mutate(fraction = n/max(n)) %>% right_join(shp,by=GEO)
p2 <- mapSettings + geom_sf(data=plotDat2 ,size=0, aes(fill=n))  +
  guides(fill=guide_legend(title="Year 1 count")) +
  viridis::scale_fill_viridis(na.value="transparent",breaks=seq(350,1600,by=200))
p2

png(filename=paste('Year1_catchment_map_all.png',sep='-'),res=600, units='in', width=6.5,height=4)
p2
dev.off()
