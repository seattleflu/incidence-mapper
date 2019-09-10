# DeDx output format script

library(modelServR)
library(jsonlite)
library(ggplot2)

dat<- loadModelFileById('da4335e1dff3f2763dbe52e441fe9739', model_store_dir='./test_model_store')

# reproduce query that generated file
fluABtestPathogens <- c('Flu_A_H1','Flu_A_H3','Flu_A_pan','Flu_B_pan')

queryIn <- list(
  SELECT   =list(COLUMN=c('pathogen','site_type','residence_puma','encountered_week','age_range_fine_upper')),
  GROUP_BY =list(COLUMN=c('site_type','residence_puma','encountered_week','age_range_fine_upper')),
  SUMMARIZE=list(COLUMN='pathogen', IN= fluABtestPathogens)
)

# filter to a few weeks in january
dat <- dat %>% filter(encountered_week %in% c('2019-W01','2019-W02','2019-W03','2019-W04'))

# drop time and age bin columns and empirical fraction
dat <- dat[,-c(5,8,9,10)]


# since self-test demographic interpolation isn't hooked up yet, what is most similar to self-test?
# demogQuery <- list(
#   SELECT   =list(COLUMN=c('site_type','age_range_fine_upper')),
#   GROUP_BY =list(COLUMN=c('site_type','age_range_fine_upper')),
#   SUMMARIZE=list(COLUMN='site_type', IN= 'all')
# )
# 
# #VPN down...  
# demog <- selectFromDB(demogQuery, credentials_path = './')

# cached
demog <- read.csv('C:/Users/mfamulare/Dropbox (IDM)/seattle-flu-study/cached_data/SFS-08052019-1649-distinct_from_shipping.incidence_model_observation_v1_join_shipping.presence_absence_result_v1.csv')
levels(demog$site_type) <- c(levels(demog$site_type),'self-test')
demog$site_type[is.na(demog$site_type)] <- 'self-test'
demog <- demog %>% distinct(sample,.keep_all=TRUE) %>% group_by(site_type,age_range_fine_upper) %>% 
  summarize(n=n()) %>% mutate(fraction = n/sum(n))

ggplot(demog) + geom_line(aes(x=age_range_fine_upper, y=fraction, group=site_type,color=site_type)) + facet_wrap('site_type')

# no ages for self test!
# best_fit <- data.frame(site_type=unique(demog$site_type))
# for (k in 1:nrow(best_fit)){
#   print(demog[demog$site_type==best_fit$site_type[k],]$fraction)
#   best_fit$score = sum((demog[demog$site_type==best_fit$site_type[k],]$fraction - demog[demog$site_type=='self-test',]$fraction)^2)
# }
# best_fit

# just assume self-test looks like "collegeCampus" for now...
levels(dat$site_type) <- c(levels(dat$site_type),'self-test')
dat$site_type[dat$site_type == 'collegeCampus'] <- 'self-test'

# drop all but self-test
dat <- dat[dat$site_type=='self-test',]

dat$residence_puma <- factor(dat$residence_puma)
ggplot(dat) + geom_line(aes(x=age_range_fine_upper,y=modeled_fraction_mode,group=residence_puma,color=residence_puma)) + facet_wrap('encountered_week')


# example output
saveList = list(database_query=queryIn, 
                model = 'inla_observed', 
                spatial_domain ="king_county_geojson_puma",
                created = '2019-09-10 13:50:18Z',
                model_output=dat)

write_json(saveList,'./test_model_store/da4335e1dff3f2763dbe52e441fe9739.json',pretty=TRUE)



