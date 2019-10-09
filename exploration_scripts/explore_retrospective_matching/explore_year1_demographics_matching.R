# explore retrospective matching

library(dbViewR)
library(dplyr)
library(magrittr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)


### select retrospective data at puma scale
  shp <- masterSpatialDB(shape_level = 'census_tract', source = 'wa_geojson')
  
  queryIn <- list(
    SELECT = list(COLUMN = '*')
  )
  
  db <- selectFromDB(queryIn)
  
  db$observedData$residence_puma <- as.character(shp$residence_puma[match(db$observedData$residence_census_tract,shp$residence_census_tract)])
  
  
  fluABtestPathogens <- c('Flu_A_H1','Flu_A_H3','Flu_A_pan','Flu_B_pan')
  
  names(db$observedData)
  
  dat <- db$observedData %>% mutate(flu_positive =(pathogen %in% fluABtestPathogens)) %>%
    filter(site_type %in% c('retrospective')) %>% filter(site != 'Unknown') %>% distinct(sample, flu_positive,.keep_all = TRUE) %>%
    filter(!is.na(residence_puma))
  dups <- intersect(dat$sample[dat$flu_positive],dat$sample[!dat$flu_positive])  
  dat <- dat[dat$flu_positive | !(dat$sample %in% dups), ]
    
##############################################################################################
#### data exploration of demographic marginal distributions for flu-positives vs flu-negatives
############################################################################################## 
# we are looking to see if the demographic details of flu-negs and flu-pos are similar, 
# so that we can use simple sampling.   If any are different, we need to figure out why, 
# and if they should be different (like age-dependence of different diseases) or should not
# like sex. 
    
# marginal distributions each week.
  
factors = c('age_range_fine_upper','sex','race','residence_puma') # dropped: hispanic_or_latino, flu_shot

datPlotHolder<-list()
for(FACTOR in factors){
  
  # FACTOR='sex'
  # FACTOR='age_range_fine_upper'
  # FACTOR='residence_puma'
  # FACTOR='race'
  
  plotDat <- dat %>%
    group_by_('site','flu_positive','encountered_week',FACTOR) %>%  
    summarize(n=n()) %>% group_by_('site','flu_positive',FACTOR) %>% mutate(demog_frac = n/sum(n)) %>%
    arrange_('encountered_week','flu_positive','site',FACTOR)
  
  plotDat <- plotDat %>% group_by_('site','encountered_week','flu_positive') %>%
    mutate(mean_factor = n/sum(n))
  
  
  datPlotHolder[[FACTOR]]<-ggplot(plotDat) + geom_line(aes_string(x='encountered_week',y='demog_frac',group='flu_positive',color='flu_positive')) +
      facet_grid(as.formula(paste('site ~',FACTOR))) +
      theme(axis.text.x = element_text(angle = 90))+ ggtitle(FACTOR) +
    theme_bw() + scale_y_continuous(expand=c(0,0))
  
}

for(k in 1:length(factor)){ print(datPlotHolder[k])}
  
# visually, flu- and flu+ look a bit different. how significant is this, given noise, and what might real differences mean?
  
  
###########################################################################################
## Manhattan plots: p values of two-sided tests comparing distirbutions.
## how frequently are the flu- and flu+ each week significantly different?
## expected answer under null of no differences in p<0.05 happens less than 5% of the time
  
pvals <- expand.grid(site=unique(dat$site), 
                     encountered_week = sort(unique(dat$encountered_week)),
                     factor=factors)
  
plotHolder=list()

for(FACTOR in factors){
    
  if(FACTOR == 'sex'){ # sex fisher test
    
    plotDat <- dat %>%
      group_by_('site','flu_positive','encountered_week',FACTOR) %>%  
      summarize(n=n()) %>% group_by(site,flu_positive) %>% mutate(demog_frac = n/sum(n)) %>%
      arrange_('encountered_week','flu_positive','site',FACTOR)
    
    for(k in which(pvals$factor==FACTOR)){
      vec <- c(
        plotDat$n[plotDat$site == pvals$site[k] &
                    plotDat$encountered_week == pvals$encountered_week[k] & 
                    plotDat$flu_positive == TRUE & 
                    plotDat$sex ==  'male' ],
        plotDat$n[plotDat$site == pvals$site[k] &
                    plotDat$encountered_week == pvals$encountered_week[k] & 
                    plotDat$flu_positive == FALSE & 
                    plotDat$sex ==  'male'],
        plotDat$n[plotDat$site == pvals$site[k] &
                    plotDat$encountered_week == pvals$encountered_week[k] & 
                    plotDat$flu_positive == TRUE & 
                    plotDat$sex ==  'female' ],
        plotDat$n[plotDat$site == pvals$site[k] &
                    plotDat$encountered_week == pvals$encountered_week[k] & 
                    plotDat$flu_positive == FALSE & 
                    plotDat$sex ==  'female']
      )
      
      if(length(vec)==4){
        tbl <- matrix(vec, ncol=2, byrow=TRUE)
        pvals$p[k]<-fisher.test(tbl,alternative="two.sided")$p.value
      } else {
        pvals$p[k]<-NaN
      }
      
    }
    
  } else { #ks-test for everything else

    for(k in which(pvals$factor==FACTOR)){
      
      x<-as.numeric(factor(dat[[FACTOR]][dat$site == pvals$site[k] &
                                        dat$encountered_week == pvals$encountered_week[k] & 
                                        dat$flu_positive == TRUE],levels=sort(unique(dat[[FACTOR]]))))
      y<-as.numeric(factor(dat[[FACTOR]][dat$site == pvals$site[k] &
                                        dat$encountered_week == pvals$encountered_week[k] & 
                                        dat$flu_positive == FALSE],levels=sort(unique(dat[[FACTOR]]))))
      x<-x[!is.na(x)]
      y<-y[!is.na(y)]
      
      if(length(x)>0 & length(y)>0 ){
        pvals$p[k]<-ks.test(x,y,alternative="two.sided")$p.value
      } else {
        pvals$p[k]<-NaN
      }
      
    }
    
  }
  plotHolder[[FACTOR]] <-
    ggplot(pvals %>% filter(factor==FACTOR)) + 
    geom_point(aes(x=encountered_week,y=p)) + facet_wrap('site') +
    geom_hline(aes(yintercept=0.05),linetype='dashed') +
    geom_smooth(aes(x=as.numeric(as.factor(encountered_week)),y=p), se=FALSE) +
    theme(axis.text.x = element_text(angle = 90)) + ggtitle(FACTOR) +
    theme_bw() + scale_y_continuous(expand=c(0,0))
  
}


# summarize pvalues

excludeIdx <- pvals$factor =='age_range_fine_upper' |  # age should vary between flu and non-flu
  (pvals$site=='RetrospectiveChildrensHospitalSeattle' & pvals$factor=='race') # exclude seattlechildrens because race is all NULL

# histogram if null is true has asymptotic uniform distribution.
hist(pvals$p[!excludeIdx ],20)
# pretty much!  Heavy concetration at p=1 comes from cells with very small counts.

# if null is true, about 5% of p values should be less than 0.05
mean(pvals$p[!excludeIdx]<0.05,na.rm=T)
# yup!

# plots of pvalues 
grid.arrange(grobs=plotHolder[2:4],nrow=1)


# overall punchline is this data is consistent with the hypothesis that there are no weekly
# demographic differences by location, sex, or race for flu+ and flu- retrospectives. But location shows the 
# biggest differences, so we explore that below.


# AGE should vary between flu+ and flu-.

hist(pvals$p[pvals$factor =='age_range_fine_upper' ],20)

mean(pvals$p[pvals$factor =='age_range_fine_upper' ]<0.05,na.rm=T)

plotHolder[1]

# age is  not compatible witht the null, especially for seattle childrens.  Makes sense!



  
#########################################################  
# why is location different between flu- and flu+ ?
  
  shp <- masterSpatialDB(shape_level = 'puma', source = 'wa_geojson')
  
  bbox<-sf::st_bbox(shp$geometry)
  mapSettings <- ggplot() + theme_bw() +
    theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.grid.major=element_line(colour="transparent"), panel.border = element_blank())
  p<-mapSettings + geom_sf(data=shp,size=0.1,aes(fill=NaN))


  for(SITE in unique(dat$site)){
  
    # flu pos
      queryIn <- list(
        SELECT   =list(COLUMN=c('pathogen', 'residence_puma','site')),
        WHERE    =list(COLUMN='pathogen', IN= fluABtestPathogens),
        WHERE    =list(COLUMN='site', IN= SITE),
        GROUP_BY =list(COLUMN=c('residence_puma')),
        SUMMARIZE=list(COLUMN='pathogen', IN= 'all')
      )
      db <- expandDB(selectFromDB(  queryIn, source='production', na.rm=TRUE ), shp=shp )
      
      plotDat <- right_join(db$observedData,shp, by='residence_puma')
      plotDat$positive[plotDat$n==0]<-NaN
      
      p1 <- p + geom_sf(data=plotDat,size=0, aes(fill=positive))  +
        guides(fill=guide_legend(title="flu-pos")) +
        scale_fill_viridis(na.value="transparent",trans = "sqrt")
    
    # flu neg
      queryIn <- list(
        SELECT   =list(COLUMN=c('pathogen', 'residence_puma','site')),
        WHERE    =list(COLUMN='pathogen', IN= setdiff(unique(dat$pathogen),fluABtestPathogens)),
        WHERE    =list(COLUMN='site', IN= SITE),
        GROUP_BY =list(COLUMN=c('residence_puma')),
        SUMMARIZE=list(COLUMN='pathogen', IN= 'all')
      )
      db <- expandDB(selectFromDB(  queryIn, source='production', na.rm=TRUE ), shp=shp )
      
      plotDat <- right_join(db$observedData,shp, by='residence_puma')
      plotDat$positive[plotDat$n==0]<-NaN
      
      p2 <- p + geom_sf(data=plotDat,size=0, aes(fill=positive))  +
        guides(fill=guide_legend(title="flu-neg")) +
        scale_fill_viridis(na.value="transparent",trans = "sqrt")
      

    grid.arrange(p1,p2,nrow=1, top=textGrob(SITE))
  }
# plausible reason in flu-negative catchments are just much bigger. This could be due to volume of data 
# and too-small statistical coverage, or real effects in care seeking. Either way, let's look within King County next

  
# location in Seattle-metro only
  corePumas=as.character(c(11601:11605, # seattle
                           11606,11701:11704,  # everett, south snohomish
                           11607:11611, # king county eastside clockwise around lake washington
                           11612:11614 # south king
                        ))
  df <- data.frame(residence_puma = shp$residence_puma, in_metro = as.numeric(shp$residence_puma %in% corePumas))
  df$in_metro[df$in_metro==0]<- NA
  plotDat <- right_join(df,shp, by='residence_puma')
  p + geom_sf(data=plotDat,size=0, aes(fill=in_metro))  + guides(fill=FALSE) + scale_fill_continuous(na.value="transparent")
  
  dat2 <- db$observedData %>% filter(residence_puma %in% corePumas) %>% 
    mutate(flu_positive =(pathogen %in% fluABtestPathogens)) %>%
    filter(site_type %in% c('retrospective')) %>% filter(site != 'Unknown') %>% distinct(sample,.keep_all = TRUE) %>%
    filter(!is.na(residence_puma))
  
  FACTOR = 'residence_puma'
  
  pvals <- expand.grid(site=unique(dat2$site))
  
  for(k in 1:nrow(pvals)){
    
    x<-as.numeric(factor(dat2$residence_puma[dat2$site == pvals$site[k] &
                                              dat2$flu_positive == TRUE],levels=sort(unique(dat2$residence_puma))))
    y<-as.numeric(factor(dat2$residence_puma[dat2$site == pvals$site[k] &
                                              dat2$flu_positive == FALSE],levels=sort(unique(dat2$residence_puma))))
    x<-x[!is.na(x)]
    y<-y[!is.na(y)]
    
    if(length(x)>0 & length(y)>0 ){
      pvals$p[k]<-ks.test(x,y,alternative="two.sided")$p.value
    } else {
      pvals$p[k]<-NaN
    }
    
  }
  
  ggplot(pvals) + geom_point(aes(x=site,y=p)) +
    geom_hline(aes(yintercept=0.05),linetype='dashed') +
    theme(axis.text.x = element_text(angle = 90)) 
  
  pvals
  # Seattle and north to everett has similar distributions between flu+ and flu-, but
  # as we add in the east-side and south king county, the differences at the larger hospitals
  # manifest.  Eyeball test says the distributions are still mostly similar. 
  # The question that remains for me is this:
  #   TO WHAT EXTENT DO WE THINK LOCATION DIFFERENCES AFFECT CARE-SEEKING BEHAVIOR vs REGIONAL DIFFERENCES 
  #   INCIDENCE?
  # 
  # To figure that out, I think we have to match for symptom severity. 
  

  
###################################################
## how similar are demographics given matched severity?
###################################################
  
shp <- masterSpatialDB(shape_level = 'census_tract', source = 'wa_geojson')
db <- selectFromDB(list(SELECT = list(COLUMN = '*')))
db$observedData$residence_puma <- as.character(shp$residence_puma[match(db$observedData$residence_census_tract,shp$residence_census_tract)])

# there are no symptoms for retrospectives!
dat <- db$observedData %>% filter(!is.na(symptoms) & (site_type == 'retrospective'))
  

