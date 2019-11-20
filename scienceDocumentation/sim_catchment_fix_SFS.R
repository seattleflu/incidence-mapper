########################################################################################
# Author:   Roy Burstein
# Date:     August 2019
# Purpose:  Simulate a simple version of Mike's catchment model to answer the question:
#           is the remaining 'latent field' proportional to incidence?
########################################################################################

library(data.table)
library(magrittr)
library(ggplot2)
library(INLA)
library(gridExtra)
library(grid)

set.seed(123457)
########################################################################################
## Simulate data (admin, site, pathogen cases by week)

### set params

# number of admin areas
n_admins <- 5

# number of pathogens tracked
n_pathogens <- 10

# number of site types
n_sitetypes <- 5

# TODO, need N as well (site and admin specific population) - draw from a binomial.. 

# set up param table
params <- expand.grid(admin    = paste0('AD_',1:n_admins),
                      pathogen = paste0('PATH_',letters[1:n_pathogens])) %>%
              data.table()


# for now, we assume cases are independent (the current method of using total other pathogens depends on this assumption, I think)
# later we can change this
# make them zero inflated but if have cases some random number of cases
params[, totalcases := rbinom(.N, 1, 0.9) * rpois(.N, exp(rnorm(.N, 5, 1)))]


# imagine an under-reporting matrix by path, admin, site. (suppress week for now). parameter is called rho_paw
# assume these are randomly drawn for now (not likely true obvioulsy.. )
ur <- expand.grid(admin    = paste0('AD_',1:n_admins),
                  pathogen = paste0('PATH_',letters[1:n_pathogens]), # The assumption is that it is equal across pathogens
                  sitetype = paste0('SITE_',1:n_sitetypes)) %>% data.table()

# if want a as rho
#ur[, rho := rbeta(.N,2,2)]
ur[, rho := ((rbinom(.N,1,.5)*rbeta(.N,5,1)+rbinom(.N,1,.5)*rbeta(.N,5,1))/2)]

# if want seperate a and s rhos
#ura <- data.table(admin    = paste0('AD_',1:n_admins),      rho = rbeta(n_admins,2,2) )
#ur <- merge(ur, ura, by = 'admin') #, urs, by ='sitetype') #urs <- data.table( sitetype = paste0('SITE_',1:n_sitetypes),rhos = rbeta(n_sitetypes,2,2) )

#ur[, rho := rhoa*rhos]


# Simulate overall incidence over a 20 week period for each admin area and pathogen (again independently)
d <- data.table()
for(r in 1:nrow(params)){
  
  Nc  <- params$totalcases[r] 
    
  if(Nc != 0){
  
    # an epi curve ( just do a simple gaussian distribution, moving mean and sd a bit for each admin-pathogen combo)
    peak <- abs(round(rnorm(1, 14, 3)))
    stdv <- abs(round(rnorm(1, 4,  1)))
    week <- ceiling(abs(round(rnorm(Nc, mean = peak, sd = stdv))))
  #  true_intensity <- data.table(true_intensity=dnorm(1:max(week),peak,stdv),week=1:max(week))
    
    # for each case, draw a sitetype it was detected at (for now assume equal probability detection at each site)
    # this is pretty unrealistic in the real data so we'll need to think aboiut the effect of this assumption (what is most cases are from one site one admin (ie UW))
    #styp <- paste0('SITE_',ceiling(runif(Nc,0,n_sitetypes)))
    styp <- paste0('SITE_',ceiling(((rbinom(Nc,1,.5)*rbeta(Nc,5,1)+rbinom(Nc,1,.5)*rbeta(Nc,5,1))/2)*n_sitetypes+.001))
    # append to d
    out <- data.table(admin = params$admin[r], pathogen = params$pathogen[r], week= week, sitetype = styp)
   # out <- merge(out, true_intensity, by = 'week', all.x=T)
    d   <- rbind(d, out)
  } 
}


# now, randomly choose if each case was recorded or not
dtrue <- merge(d, ur, by = c('admin','sitetype'), all.x = TRUE)
dtrue[, captured := rbinom(.N,1,rho)]

dobs <- dtrue[captured == 1]

dtrue <- dtrue[order(admin,pathogen,week)]
dobs  <- dobs[order(admin,pathogen,week)]


# plot to see whats happening
#ggplot(dtrue, aes(week,fill=sitetype)) + geom_histogram() + facet_grid(admin~pathogen)
g1<-ggplot(dobs, aes(week,fill=pathogen)) + geom_histogram() + ylab('Observed Count') +
  facet_grid(admin~sitetype) + theme_bw() + theme(legend.position = 'none')
g2<-ggplot(dtrue, aes(week,fill=pathogen)) + geom_histogram() + ylab('True Count') +
  facet_grid(admin~sitetype) + theme_bw() + theme(legend.position = 'none')
g1
g2


# get catchment as currently defined (all pathogens by Admin, site other than the pathogen in question)
# late can do a smooth model like mike did for this and then wont get any zeroes 
catch <- data.table()
for(p in unique(dobs$pathogen)){
  tmp <- dobs[pathogen != p]
  if(nrow(tmp)>0){
    tmp <- tmp[, .(catchment = .N+0.001), by = .(admin, sitetype)]
    tmp$pathogen <- p
    catch <- rbind(catch, tmp)
  } 
}

# expand out and fill in zeroes
# This could lead to problems, since avg over S/A could miss a whole outbreak if its in one S/A?
expp <-  expand.grid(admin    = paste0('AD_',1:n_admins),
                     pathogen = paste0('PATH_',letters[1:n_pathogens]),
                     sitetype = paste0('SITE_',1:n_sitetypes)) %>% data.table()
catch <- merge(expp, catch, by = c('admin','pathogen','sitetype'), all.x=T)
catch[, id := 1:.N]
catch[is.na(catchment),     catchment := 0]
m <- inla(formula = catchment ~ f(id, model = 'iid'), 
          data = catch, family = 'poisson',control.predictor = list(compute=TRUE,link=1))
catch$catchment <- exp(m$summary.linear.predictor$mean) ## NOTE TO MIKE: THIS WAY OF IMPUTING ONLY WORKS WITHOUT FE. OTHERWISE OVERESTIMATING CATCHMENT SIZE
#catch[is.na(catchment), catchment := 0.1]

# get actual catchment/propensity as well

#dobs <- merge(dobs, catch, by = c('admin','sitetype','pathogen'), all.x = T)
#dobs[is.na(catchment), catchment := 0.1] # in cases with no other pathogens present

# aggregate date by week admin site pathogen for the mode
dagg <- dobs[, .(cases = .N), by = .(admin,sitetype,pathogen,week)]

# exapnd the aggregated data for all possible combos
expanded <-  expand.grid(admin    = paste0('AD_',1:n_admins),
                         pathogen = paste0('PATH_',letters[1:n_pathogens]),
                         sitetype = paste0('SITE_',1:n_sitetypes),
                         week     = 1:max(d$week)) %>% data.table()
dagg <- merge(expanded, dagg, by = c('admin', 'pathogen', 'sitetype', 'week'), all.x = TRUE)

dagg <- merge(dagg, catch, by = c('admin','sitetype','pathogen'), all.x = T)


dagg[is.na(catchment) , catchment := 0.1]
dagg[is.na(cases),     cases := 0]

dagg[, catchment := log(catchment)-mean(log(catchment))]

# MERGE CATCHMENT ERROR TODO LOOK INTO IT


########################################################################################
## Run model, try to do similar to mike's iid model as possible

## model one pathogen at a time.
path_to_model <- 'PATH_a'

# subset
inputData <- dagg[pathogen == path_to_model]

# priors
hyper <- list()
#hyper$global <- list(prec = list( prior = "pc.prec", param = 100, alpha = 0.01))
hyper$local  <- list(prec = list( prior = "pc.prec", param = 1/100, alpha = 0.01)) # these priors make the RWs all look the same
#hyper$age    <- list(prec = list( prior = "pc.prec", param = 1, alpha = 0.01))
#hyper$time   <- list(prec = list( prior = "pc.prec", param = 1/100, alpha = 0.01))

family <- 'poisson'

outcome <- inputData$cases

# initialize formula  
formula <- as.formula('cases ~ 1  + sitetype + catchment') # catchment

# time as a latent fied
inputData$time_row_rw2 <- inputData$week
inputData$time_row_IID <- inputData$week

inputData$admin_row <- match(inputData$admin,unique(inputData$admin))
inputData$time_row_admin <- inputData$week
  
# DIAGONAL: https://groups.google.com/forum/#!topic/r-inla-discussion-group/WqnLwPP6Ges
formula <- update(formula,  ~ . + f(admin_row, model = 'iid', graph = NULL , diagonal=5e-09, 
                                     group = time_row_admin, control.group=list(model="rw2")))
formula <- update(formula,  ~ . + f(time_row_rw2, model='rw2', diagonal=5e-09)) # + hyper= hyper$time
                 # f(time_row_IID, model='iid', hyper=hyper$local,  constr = TRUE) )
# Need a closer look at priors here, the defaults are doing much better

#fit
model <- INLA::inla(formula           = formula,
                    family            = family, 
                    data              = inputData, 
                    control.predictor = list(compute=TRUE,link=1),
                    control.compute   = list(config=TRUE,dic=TRUE),
                    verbose           = TRUE,
                    keep              = FALSE,
                    control.inla      = list(int.strategy="auto", strategy = "gaussian", cmin=0))


# admin_row will be n_admin time n_weeks long. time_row_rw2 will be n_weeks long   
out <- data.table(intensity = model$summary.random$admin_row$mean +rep(model$summary.random$time_row_rw2$mean, each=n_admins),
                  admin     = rep( paste0('AD_',1:n_admins), max(dagg$week)),
                  week      = rep(1:max(dagg$week), each=n_admins))
out <- merge(out,  dobs[pathogen == path_to_model, .(observed_cases = .N), by = .(admin, week)], by = c('admin','week'), all.x = TRUE)
out <- merge(out, dtrue[pathogen == path_to_model, .(true_cases     = .N, rho = mean(rho)), by = .(admin, week)], by = c('admin','week'), all.x = TRUE)

out[is.na(observed_cases), observed_cases := 0]
out[is.na(true_cases),     true_cases := 0]

# do the other fixed effects variables equal rho*N?

# rank admins at each time point?
rank <- copy(out)[order(week,intensity)]
rank[, intensity_rank := 1:.N, by = week]
rank <- rank[order(week,observed_cases)]
rank[, observed_cases_rank := 1:.N, by = week]
rank <- rank[order(week,true_cases)]
rank[, true_cases_rank := 1:.N, by = week]


# plot smoothed time series
g1=ggplot(out, aes(y=exp(intensity), x=week, color=admin, group=admin)) + geom_line() + theme_bw() + theme(legend.position = 'none')
g2=ggplot(out, aes(y=true_cases, x=week, color=admin, group=admin)) + geom_line() + theme_bw() + theme(legend.position = 'none')
g3=ggplot(out, aes(y=observed_cases, x=week, color=admin, group=admin)) + geom_line() + theme_bw() + theme(legend.position = 'none')
g4=ggplot(out, aes(x=log(true_cases), y=intensity, color=admin)) + geom_point() + theme_bw() + theme(legend.position = 'none')
g5=ggplot(out, aes(x=log(observed_cases), y=intensity, color=admin)) + geom_point() + theme_bw() + theme(legend.position = 'none')

g6 = ggplot(rank[week < 25], aes(week,intensity_rank,group=admin,color=admin),alpha=.4) + 
  geom_line() + theme_bw() + # admin %in% paste0('AD_',1:5) & 
  theme(legend.position = 'none') + ylab('ADMIN RANK') +
  geom_line(aes(week,true_cases_rank+.1,group=admin,color=admin), lty = 'dashed')
grid.arrange(g1,g6,g2,g3,g4,g5, layout_matrix = rbind(c(1,2),c(3,4),c(5,6)))


summary(lm(exp(intensity)~true_cases,data=out))$r.squared 
summary(lm(exp(intensity)~observed_cases,data=out))$r.squared


# explore other model compopnents, formalize whats happening.. 

# how does it scale with more pathogens? 
# how robust is it to changing under-reporting patterns?



#########################################
# what if we fit a simpler model on catchment only and looked at residual?
m <- glm(cases ~ catchment + sitetype , family = 'poisson', data = inputData) # admin all pathogens needed in here to resolve site-admin specific biases, current data is subset to one pathogen
m3 <- mgcv::gam(cases ~ catchment + sitetype + s(week, by = admin), family = 'poisson', data = inputData) #

m2 <- glm(cases ~ catchment + pathogen + sitetype*admin, family = 'poisson', data = dagg) # admin all pathogens needed in here to resolve site-admin specific biases, current data is subset to one pathogen
dd <- data.table(inputData, pred = exp(predict(m2))[dagg$pathogen==path_to_model])
dd <- data.table(inputData, pred = exp(predict(m)))
dd <- data.table(inputData, pred = exp(predict(m3))) #, terms = c('catchment', 'sitetype', '(Intercept)'))))

dd <- dd[, .(cases=sum(cases),pred=sum(pred),resid=sum(cases)-sum(pred)), by = .(week,admin)]

g7<-ggplot(dd, aes(x=week,y=resid,color=admin,group=admin)) + geom_line() + theme_bw() + theme(legend.position = 'none')
grid.arrange(g1,g7,g2,g3)






# old lincomb stuff
if(TRUE==FALSE){
  validLatentFieldColumns <- c('admin_row','time_row_admin','time_row_rw2') #,'time_row_IID')
  
  # linear combination of pathogen and latent fields
  
  # find unique rows after discarding factors that are being averaged over
  lc.data <- data.frame(inputData[,names(inputData) %in% validLatentFieldColumns, with = FALSE])
  lc.rowIdx <- !duplicated(lc.data)
  lc.data <- lc.data[lc.rowIdx,]
  
  # generate list of desired linear combinations # https://groups.google.com/forum/#!topic/r-inla-discussion-group/_e2C2L7Wc30
  lcIdx=c()
  spentColumn<-rep(FALSE,length(validLatentFieldColumns))
  for(COLUMN in validLatentFieldColumns){
    if(!spentColumn[validLatentFieldColumns %in% COLUMN]) {
      lcIdx[[COLUMN]] <- inla.idx(lc.data[[COLUMN]])          
    }
    
    spentColumn[validLatentFieldColumns %in% COLUMN]<-TRUE
  }
  
  
  # generate list of desired linear combinations # https://groups.google.com/forum/#!topic/r-inla-discussion-group/_e2C2L7Wc30
  lc.latentField <- vector("list", nrow(lc.data))
  
  w<-vector("list", length(names(lcIdx))+1)
  w[[length(names(lcIdx))+1]]<-1 #pathogen
  
  for(k in 1:nrow(lc.data)){
    
    for(n in 1:length(names(lcIdx))){
      w[[n]]<-rep(0,nrow(lc.data))
      w[[n]][lcIdx[[n]][k]]<-1
    }
    names(w) <- c(names(lcIdx),'(Intercept)')
    
    lc <- inla.make.lincomb(w)
    names(lc)<- paste0('latent_field',k)
    lc.latentField[k]<-lc
    lc.data$latentField[k]<-names(lc)
    
  }
  
  #lc.latentField <- inla.make.lincombs(lc.data)
  
  #lc.latentField <- inla.make.lincombs(time_row_rw2 = diag(lc.data$time_row_rw2), 
  #                          admin_row = diag(lc.data$admin_row), 
  #                         time_row_admin = diag(lc.data$time_row_admin))
  
  # get original values for linear combination categories
  lc.colIdx <- (names(inputData) %in% c('admin','week'))
  lc.data <-inputData[which(lc.rowIdx),c('admin','week'),with=FALSE]
  
}

