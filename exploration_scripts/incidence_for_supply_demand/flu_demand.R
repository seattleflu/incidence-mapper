# script to transform unnormalized incidence into example infectious demand curves
library(dbViewR)
library(tidycensus)
library(magrittr)
library(dplyr)
library(ggplot2)

# pull in incidence map intensity model
dat <- list()
dat$observedData <- read.csv('../../test_model_store/8723e4ec4636fed9d78b85aa2b394074.csv')
dat$observedData$encountered_week <- factor(dat$observedData$encountered_week,
                                            levels=sort(unique(dat$observedData$encountered_week)),
                                            ordered=TRUE)

# filter to year 1
dat$observedData <- dat$observedData %>% filter(encountered_week <= '2019-W26')

# pull in census tract to regional name lookup
shp_all_tracts <- masterSpatialDB(shape_level = 'census_tract', source = 'sfs_domain_geojson')
shp_sfs_tracts <- shp_all_tracts %>% filter(residence_regional_name != "NA")

# join tracts       
dat$observedData <- dat$observedData %>% right_join(as.data.frame(shp_sfs_tracts) %>% 
                                                      mutate(residence_regional_name = regional_name) %>%
                                                      select(residence_regional_name,residence_census_tract) 
                                                    )
# join total population census data
populationSizeAdjustmentFactor <- (724745/610333)^(9/7) 
popVariables <- 'P001001' # total population
dat <- dat %>% addCensusData(geography = "tract", variables=popVariables,
                                       state = "WA", county = NULL,
                                       source = "decennial", year=2010,
                                       credentials_path = '/home/rstudio/performance-metrics') 
names(dat$observedData)[names(dat$observedData)=='P001001.2010'] = "population_2010"
names(dat$observedData)
dat$observedData$population_2019 <- round(dat$observedData$population_2010 * populationSizeAdjustmentFactor)

# scale intensity to cumulative infected is 10% of total regional pop
dat$observedData <- dat$observedData %>%
  mutate(new_infections = modeled_intensity_median * population_2019 )
sum(dat$observedData$new_infections)

incidenceScaleFactor <- 0.1 * sum(dat$observedData$population_2019[dat$observedData$encountered_week=='2019-W01']) /
  sum(dat$observedData$new_infections)
dat$observedData <- dat$observedData %>%
  mutate(new_infections = round(new_infections * incidenceScaleFactor ))

sum(dat$observedData$new_infections)


# plot
png(filename = 'modeled_all_flu_infections_2018-2019_census_tract_13Jan2020.png',width = 15, height = 10, units = "in", res = 300)
print(
  ggplot(dat$observedData) + geom_line(aes(x=encountered_week,y=new_infections, group=residence_census_tract, color=residence_regional_name)) +
    facet_wrap('residence_regional_name') + guides(color=FALSE) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
)
dev.off()

# clean for export
dat$observedData <- dat$observedData  %>% select(encountered_week, residence_regional_name,
                                                 residence_census_tract,population_2019,new_infections)

write.csv(dat$observedData,'modeled_all_flu_infections_2018-2019_census_tract_13Jan2020.csv',row.names = FALSE)

  