#####################################################
##
##   BSTS - FINDING DEFAULT STATE SPACE
##
##    Run script for simulations of
##     - ...
##
#####################################################


rm(list=ls())
##
library(plyr);library(dplyr)
library(tidyr)
library(CausalImpact)
library(bsts)
library(tibble)
library(did)
library(sarima)
library(ggpubr)
library(Boom)
library(BoomSpikeSlab)
library(cowplot)
library(forecast)
library(coda)
library(e1071)
library(qualV)

## Directories
dir_proj <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'
dir_ext <- 'D:\\BSTS_external'
dir_plot <- file.path(dir_proj, 'plots')
dir_r <- file.path(dir_proj,'R')

## Set working directory
setwd(dir_proj)


##==============================
##  file prefix for saving images, writing outputs, etc.
##-----------------------------
prefix <- '__DEBUG__bsts-illus_default-st-sp_'

##==============================
## Load simulation functions
##------------------------------
## source(file.path(dir_r,'internal_intervention_sim.R')) ## previous version

## Actor index vectorized simulation
source(file.path(dir_r,'single_intervention_sim_vec.R')) 
## Setting up and adding state space components to state.space list
source(file.path(dir_r,'bsts_helper_functions.R')) 
## BSTS vs. DiD comparison and sensitivity analysis
source(file.path(dir_r,'bsts_did_comparison_functions.R')) 



# ## DEBUG ###
# effect.types=c('constant','quadratic','geometric')
# sim.id=round(10*as.numeric(Sys.time()))
# plot.show=F
# plot.save=F
# save.items.dir=NA
# bsts.niter=1000
# ##----------


################################################################################################################
################################################################################################################


##=======================================
##
##
##  1.  FIND OPTIONAL STATE SPACE CONFIG 
##      TO SERVE AS DEFAULT 
##
##
##=======================================
## Static defaults
noise.level <- 1.3
b4 <- 1.0  ## seasonal component weight
b5 <- 0.04 ## linear growth trend
dgp.nseasons= 52
dgp.freq= 1

## Variables Ranges to grid search (optional)
ns <- list(200) ## 400
sim.lengths <- list(520)
treat.rules <- list('random')  ## 'below.benchmark'
seasonalities <- list(TRUE)   ## c(TRUE,  FALSE )
prior.sd.scenarios <- list('sd.low') ## list('sd.low','sd.high')  ## sd.low
expect.mod.sizes <- list(2)
## FOCAL CONSTRUCT
dgp.ars <- list(0)  ## 0.6  ## .1,.2,.4
## STATE SPACE CONFIGURATIONS
# st.sp.lists <- list(
#   ## ## LEVEL
#   # `1`=c('AddLocalLevel'),
#   ## ## TREND
#   # `2`=c('AddLocalLinearTrend'),
#   # # `3`=c('AddStudentLocalLinearTrend'),
#   # ## ## LEVEL + SLOPE ( + AR SLOPE DRIFT)
#   # `4`=c('AddSemilocalLinearTrend'),
#   # ## ## SEASONAL
#   `5`= c('AddTrig', 'AddRegression'),
#   `5b`=c('AddSeasonal', 'AddRegression'),
#   # ## ## AUTOCORRELATION
#   # # list('AddAr'),
#   # `6`=c('AddAr'),
#   # ##------ COMBINATIONS --------------
#   # ## AR & SEASONALITY
#   # `7`= c('AddAr','AddTrig'),
#   # `7b`=c('AddAr','AddSeasonal'),
#   # # `7a`=c('AddAutoAr','AddTrig')#,
#   # ## LEVEL + ...
#   `8`= c('AddLocalLevel','AddTrig', 'AddRegression'),
#   `8b`=c('AddLocalLevel','AddSeasonal', 'AddRegression'),
#   # # `9`=c('AddLocalLevel','AddAr'),
#   # # `10`=c('AddLocalLevel','AddTrig','AddAutoAr'),
#   # # ## (LEVEL + SLOPE) + ...
#   `11`= c('AddLocalLinearTrend','AddTrig', 'AddRegression'),
#   `11b`=c('AddLocalLinearTrend','AddSeasonal', 'AddRegression'),
#   # `12`=c('AddLocalLinearTrend','AddAr'),
#   # `12`= c('AddLocalLinearTrend','AddAr','AddTrig'),
#   # `12b`=c('AddLocalLinearTrend','AddAr','AddSeasonal'),
#   # `13`= c('AddStudentLocalLinearTrend','AddTrig'),
#   # `13b`=c('AddStudentLocalLinearTrend','AddSeasonal')#,
#   # `14`=c('AddStudentLocalLinearTrend','AddAr'),
#   # # `14`c('AddStudentLocalLinearTrend','AddTrig','AddAr'),
#   ## (LEVEL + SLOPE ( + AR1 SLOPE DRIFT)) + ...
#   `15`= c('AddSemilocalLinearTrend','AddTrig', 'AddRegression'),
#   `15b`=c('AddSemilocalLinearTrend','AddSeasonal', 'AddRegression')#,
# )

##**DEBUG**
st.sp.lists <- list(
  `11b`=c('AddLocalLinearTrend','AddSeasonal')  ## 'AddSeasonal'
)

##
simlist <- list()
## NUMBER OF ACTORS (number of timeseries)
for (d in 1:length(ns)) {
  n <- ns[[ d ]]
  ## SIMULATION LENGTHS - NUMBER OF PERIODS
  for (f in 1:length(sim.lengths)) {
    npds <- sim.lengths[[ f ]]
    intpd <-  round( npds * 5/6 ) 
    ## AUTOCORRELATION VALUES
    for (g in 1:length(dgp.ars)) {
      dgp.ar <- dgp.ars[[ g ]]
      
      ## SEASONALITY
      for (h in 1:length(seasonalities)) {
        seasonality <- seasonalities[[ h ]]
        
        ## ENDOGENEITY (SELF-SELECTION)
        for (i in 1:length(treat.rules)) {
          treat.rule <- treat.rules[[ i ]]
          
          ##--------------------------------------------
          ## Setup state space configurations
          ##--------------------------------------------
          bsts.state.specs <- list()
          ## PRIOR NOISE / UNCERTAINTY (PRIOR STDEV)
          for (j in 1:length(prior.sd.scenarios)) {
            prior.sd.scenario <- prior.sd.scenarios[[ j ]]
            
            ## STATE SPACE COMPONENTS CONFIGURATION
            for (k in 1:length(st.sp.lists)) {
              st.sp.vec <- st.sp.lists[[ k ]]
              
              bsts.state.config <- list()
              for (kk in 1:length(st.sp.vec)) {
                ## SKIP AddSharedLocalLevel() for now...
                .id <- length(bsts.state.config)+1
                bsts.state.config[[ .id ]] <- getStateSpaceConfBySimScenario(st.sp.vec[ kk ], prior.sd.scenario)#, ## c('sd.high','sd.low')
                # x.spikslab.prior=x.spikslab.prior, ## X  values for Boom::SpikeSlabPrior()
                # y.spikslab.prior=y.spikslab.prior)
                # .id <- .id + 1
              }
              bsts.state.specs[[ paste(st.sp.vec, collapse='|') ]] <- bsts.state.config
            }
            
            ## EXPECTED MODEL SIZES
            for (l in 1:length(expect.mod.sizes)) {
              expect.mod.size <- expect.mod.sizes[[ l ]]
              
              ##--------------------------------------------
              ## Append simulation configuration to simlist
              ##--------------------------------------------
              key <- sprintf('d%s|f%s|g%s|h%s|i%s|j%s|l%s', d,f,g,h,i,j,l)
              cat(sprintf('\n%s\n',key))
              # .idx <- sprintf('ar%s',dgp.ar)
              simlist[[ key ]] <- list(
                n = n,    ## Number of firms
                npds = npds,  ## Number of periods
                intpd = intpd, ## #intervention period
                noise.level = noise.level, ## stdev of simulated noise terms
                prior.sd.scenario = prior.sd.scenario, ## BSTS Prior SD scenario (high vs. low uncertainty in priors
                treat.rule = treat.rule, 
                treat.prob = ifelse(treat.rule=='random', 0.5, 1), 
                treat.threshold = ifelse(treat.rule=='random', 1, 0.5),
                ## Dynamic treatment effect  (quadratic polynomial)
                w0 = 1.5,  ## constant
                w1 = 0.13, ## linear
                w2 = -.05 / sqrt(npds), ## quadratic
                w2.shift = -round( sqrt(npds)*.7 ), ## quadratic shift rightward (make all of U-shape after intervention)
                ##
                b4 = b4,   ## seasonal component weight
                b5 = b5, ##
                b9 = dgp.ar  , ## autocorrelation
                seasonality = seasonality,
                dgp.nseasons= ifelse(seasonality, dgp.nseasons, NA), 
                dgp.freq= ifelse(seasonality, dgp.freq, NA),
                bsts.state.specs=bsts.state.specs,
                expect.mod.size=expect.mod.size,
                # bsts.state.specs=list(list(AddSemilocalLinearTrend),list(AddSemilocalLinearTrend,AddStudentLocalLinearTrend)),
                rand.seed = 13579
              )
              
            }
            
          } ## // end prior.sd.scenarios loop
          
        } ## // end treat.rules loop
        
      } ## // end seasonalities loop
      
    } ## // end ARs loop
    
  } ## // end nps loop  (sim lengths)
  
} ## // end ns loop (number of actors)




##
effect.types = c('constant') #c('geometric','quadratic')
##
bsts.niter.start <- 200 ## 5000
bsts.niter.max   <- 200 ## 8e4
## ID for the simulation (to search/filter all simulation figures, RDS files, etc.)
sim.id <- round(10*as.numeric(Sys.time()))


## RUN SIMULATION -  SIMULATE TIME SERIES
simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               sim.id = sim.id,
                               plot.show = F, plot.save = F )

## RUN BSTS and compare to DID
simlist.files <- runSimCompareBstsDiD(simlist, 
                                      effect.types = effect.types,
                                      sim.id = sim.id,
                                      save.items.dir= dir_ext,
                                      bsts.niter = bsts.niter.start,
                                      bsts.max.iter= bsts.niter.max
)  ## D:\\BSTS_external






# ########################## END ##########################################





##=======================================
##
##
##  2.  FIND ALTERNATIVE STATE SPACE CONFIG 
##      FOR COMPARISON 
##
##
##=======================================
## Static defaults
noise.level <- 1.3
b4 <- 1.0  ## seasonal component weight
b5 <- 0.04 ## linear growth trend
dgp.nseasons= 52
dgp.freq= 1

## Variables Ranges to grid search (optional)
ns <- list(200) ## 400
sim.lengths <- list(520)
treat.rules <- list('random')  ## 'below.benchmark'
seasonalities <- list(TRUE)   ## c(TRUE,  FALSE )
prior.sd.scenarios <- list('sd.low') ## list('sd.low','sd.high')  ## sd.low
expect.mod.sizes <- list(2)
## FOCAL CONSTRUCT
dgp.ars <- list(0)  ## 0.6  ## .1,.2,.4
## STATE SPACE CONFIGURATIONS
st.sp.lists <- list(
  `1`=c('AddLocalLevel', 'AddRegression'),
  `2`=c('AddLocalLinearTrend', 'AddRegression'),
  `4`=c('AddSemilocalLinearTrend', 'AddRegression'),
  `6`= c('AddAr', 'AddRegression'),
  `7`= c('AddAr','AddTrig', 'AddRegression'),
  `7b`=c('AddAr','AddSeasonal', 'AddRegression')
)

##
simlist <- list()
## NUMBER OF ACTORS (number of timeseries)
for (d in 1:length(ns)) {
  n <- ns[[ d ]]
  ## SIMULATION LENGTHS - NUMBER OF PERIODS
  for (f in 1:length(sim.lengths)) {
    npds <- sim.lengths[[ f ]]
    intpd <-  round( npds * 5/6 ) 
    ## AUTOCORRELATION VALUES
    for (g in 1:length(dgp.ars)) {
      dgp.ar <- dgp.ars[[ g ]]
      
      ## SEASONALITY
      for (h in 1:length(seasonalities)) {
        seasonality <- seasonalities[[ h ]]
        
        ## ENDOGENEITY (SELF-SELECTION)
        for (i in 1:length(treat.rules)) {
          treat.rule <- treat.rules[[ i ]]
          
          ##--------------------------------------------
          ## Setup state space configurations
          ##--------------------------------------------
          bsts.state.specs <- list()
          ## PRIOR NOISE / UNCERTAINTY (PRIOR STDEV)
          for (j in 1:length(prior.sd.scenarios)) {
            prior.sd.scenario <- prior.sd.scenarios[[ j ]]
            
            ## STATE SPACE COMPONENTS CONFIGURATION
            for (k in 1:length(st.sp.lists)) {
              st.sp.vec <- st.sp.lists[[ k ]]
              
              bsts.state.config <- list()
              for (kk in 1:length(st.sp.vec)) {
                ## SKIP AddSharedLocalLevel() for now...
                .id <- length(bsts.state.config)+1
                bsts.state.config[[ .id ]] <- getStateSpaceConfBySimScenario(st.sp.vec[ kk ], prior.sd.scenario)#, ## c('sd.high','sd.low')
                # x.spikslab.prior=x.spikslab.prior, ## X  values for Boom::SpikeSlabPrior()
                # y.spikslab.prior=y.spikslab.prior)
                # .id <- .id + 1
              }
              bsts.state.specs[[ paste(st.sp.vec, collapse='|') ]] <- bsts.state.config
            }
            
            ## EXPECTED MODEL SIZES
            for (l in 1:length(expect.mod.sizes)) {
              expect.mod.size <- expect.mod.sizes[[ l ]]
              
              ##--------------------------------------------
              ## Append simulation configuration to simlist
              ##--------------------------------------------
              key <- sprintf('d%s|f%s|g%s|h%s|i%s|j%s|l%s', d,f,g,h,i,j,l)
              cat(sprintf('\n%s\n',key))
              # .idx <- sprintf('ar%s',dgp.ar)
              simlist[[ key ]] <- list(
                n = n,    ## Number of firms
                npds = npds,  ## Number of periods
                intpd = intpd, ## #intervention period
                noise.level = noise.level, ## stdev of simulated noise terms
                prior.sd.scenario = prior.sd.scenario, ## BSTS Prior SD scenario (high vs. low uncertainty in priors
                treat.rule = treat.rule, 
                treat.prob = ifelse(treat.rule=='random', 0.5, 1), 
                treat.threshold = ifelse(treat.rule=='random', 1, 0.5),
                ## Dynamic treatment effect  (quadratic polynomial)
                w0 = 1.5,  ## constant
                w1 = 0.13, ## linear
                w2 = -.05 / sqrt(npds), ## quadratic
                w2.shift = -round( sqrt(npds)*.7 ), ## quadratic shift rightward (make all of U-shape after intervention)
                ##
                b4 = b4,   ## seasonal component weight
                b5 = b5, ##
                b9 = dgp.ar  , ## autocorrelation
                seasonality = seasonality,
                dgp.nseasons= ifelse(seasonality, dgp.nseasons, NA), 
                dgp.freq= ifelse(seasonality, dgp.freq, NA),
                bsts.state.specs=bsts.state.specs,
                expect.mod.size=expect.mod.size,
                # bsts.state.specs=list(list(AddSemilocalLinearTrend),list(AddSemilocalLinearTrend,AddStudentLocalLinearTrend)),
                rand.seed = 13579
              )
              
            }
            
          } ## // end prior.sd.scenarios loop
          
        } ## // end treat.rules loop
        
      } ## // end seasonalities loop
      
    } ## // end ARs loop
    
  } ## // end nps loop  (sim lengths)
  
} ## // end ns loop (number of actors)




##
# effect.types = c('constant') #c('geometric','quadratic')

# # ## BSTS MCMC iterations
# # bsts.niter <- 2e4
# ## ID for the simulation (to search/filter all simulation figures, RDS files, etc.)
# if ( ! 'sim.id' %in% ls())
#   sim.id <- round(10*as.numeric(Sys.time()))


## RUN SIMULATION -  SIMULATE TIME SERIES
simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               sim.id = sim.id,
                               plot.show = F, plot.save = F )

## RUN BSTS and compare to DID
simlist.files <- runSimCompareBstsDiD(simlist, 
                                      effect.types = effect.types,
                                      sim.id = sprintf('_%s-ALT_', sim.id),
                                      save.items.dir= dir_ext,
                                      bsts.niter = bsts.niter.start,
                                      bsts.max.iter= bsts.niter.max
)  ## D:\\BSTS_external









###############################################################
###############################################################





# 'D:\\BSTS_external'

# filename <- '__GRIDSEARCH_output__n200_pd520_niter50000_16724204820_d1f1g1h1i1j1l1.rds'
# filename <- '__GRIDSEARCH_output__n200_pd520_niter1e+05_16726555451_d1f1g1h1i1j1l1.rds'
filename <- '__GRIDSEARCH_output__n200_pd520_niter80000_16727165757_d1f1g1h1i1j1l1.rds'
simlist <- readRDS( file.path( dir_ext, 'bsts_default_state_space_20230103', filename ) )


res.tbl <- simlist$compare$res.tbl$constant

intpd <- which( res.tbl[[1]]$event.time == 0)
npds <- nrow(res.tbl[[1]])

x <- x.orig <- simlist$compare$bsts$constant
sslist <- lapply(1:length(x), function(i) x[[i]]$CausalImpact$model$bsts.model )
CompareBstsModels( sslist )

##
rm(simlist)

# z <- rowSums(abs(sslist[[1]]$one.step.prediction.errors[,1:(intpd-1)]))

res <- t( sapply(1:length(x), function(i) summary(rowSums(abs(sslist[[i]]$one.step.prediction.errors[,1:(intpd-1)]))) ) )
print(round(res, 3))

sscomp <- abs(as.data.frame(res))
sscomp$ss.id <- 1:nrow(sscomp)
sscomp$ss.comps <- sapply(1:nrow(sscomp), function(i) paste(st.sp.lists[[i]], collapse = '|') )



## PRE ----
# x <- simlist$compare$bsts$constant
## ***SLOW*** [10 - 15 minutes]
ss.fits <- lapply(1:length(st.sp.lists), function(i){
  cat(sprintf(' %s ',i))
  idx <- 1:(intpd-1)
  mod <- x[[i]]$CausalImpact$model$bsts.model
  newdata <- cbind(treatment_y_outcome=mod$original.series[idx], mod$predictors[idx,])
  burn <- round( .2 * mod$niter )
  fit <- predict(mod, newdata=newdata, horizon = 1, burn = burn) ## one-step ahead prediction
  yhat.t <- fit$mean
  y.t <- fit$original.series[idx]
  return(list(yhat.t=yhat.t, 
              y.t=y.t#, 
              # fit=fit, 
              # mod=mod
  ))
})

## Save comparison list (slow to run)
saveRDS(ss.fits, file=file.path(dir_ext, 'bsts_default_state_space_comparison_1step-pred-err_list.rds'))


## alias for constant MCC curve shape state space list
# x <- simlist$compare$bsts$constant
## sMAPE - symmetric MAPE
sscomp$bsts.pre.onestep.smape <- laply(ss.fits, function(l){
  yhat.t <- l$yhat.t
  y.t <- l$y.t
  100 * mean( abs(yhat.t - y.t) / (abs(y.t) - abs(yhat.t))  , na.rm=T)
})
sscomp$bsts.pre.onestep.asmape <- abs( sscomp$bsts.pre.onestep.smape )
sscomp$bsts.pre.onestep.mape <- laply(ss.fits, function(l){
  yhat.t <- l$yhat.t
  y.t <- l$y.t
  100 * mean( abs( (yhat.t - y.t) / y.t ), na.rm=T )
})
sscomp$bsts.pre.onestep.mae <- laply(ss.fits, function(l){
  err <- l$yhat.t - l$y.t
  mean( abs( err ), na.rm=T ) 
})
sscomp$bsts.pre.onestep.mse <- laply(ss.fits, function(l){
  err <- l$yhat.t - l$y.t
  mean( err^2, na.rm=T )
})
sscomp$bsts.pre.onestep.rmse <- laply(ss.fits, function(l){
  err <- l$yhat.t - l$y.t
  sqrt( mean( err^2, na.rm=T ) )
})
sscomp$bsts.pre.cae <- laply(ss.fits, function(l){
  err <- l$yhat.t - l$y.t
  sum( abs(err), na.rm=T ) 
})

## POST --
# res.tbl <- simlist$compare$res.tbl$constant
## 
sscomp$bsts.post.smape <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  dgp <- res.tbl[[i]]$b3.att[idx]
  bsts <- res.tbl[[i]]$bsts.point.effect[idx]
  100 * mean( abs(bsts - dgp) / (abs(dgp) - abs(bsts))  , na.rm=T)
})
sscomp$bsts.post.asmape <- abs( sscomp$bsts.post.smape )
sscomp$bsts.post.mape <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  dgp <- res.tbl[[i]]$b3.att[idx]
  bsts <- res.tbl[[i]]$bsts.point.effect[idx]
  eps <- bsts - dgp
  100 * mean( abs( eps / dgp), na.rm=T)
})
## RMSE - Root Mean Square Error
sscomp$bsts.post.rmse <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  dgp <- res.tbl[[i]]$b3.att[idx]
  bsts <- res.tbl[[i]]$bsts.point.effect[idx]
  eps <- bsts - dgp
  sqrt( mean( eps^2, na.rm=T) )
})
sscomp$bsts.post.mse <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  dgp <- res.tbl[[i]]$b3.att[idx]
  bsts <- res.tbl[[i]]$bsts.point.effect[idx]
  eps <- bsts - dgp
  mean( eps^2, na.rm=T)
})
sscomp$bsts.post.mae <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  dgp <- res.tbl[[i]]$b3.att[idx]
  bsts <- res.tbl[[i]]$bsts.point.effect[idx]
  eps <- bsts - dgp
  mean( abs( eps ), na.rm=T)
})
sscomp$bsts.post.cae <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  dgp <- res.tbl[[i]]$b3.att[idx]
  bsts <- res.tbl[[i]]$bsts.point.effect[idx]
  err <- bsts - dgp
  sum( abs(err), na.rm=T ) 
})



sscomp$has.season <- sapply(1:nrow(sscomp),function(i) 1*grepl(pattern = '(AddTrig|AddSeasonal)', x = sscomp$ss.comps[i], ignore.case = T, perl = T))


# 
# sscomp <- sscomp[order(sscomp$bsts.pre.onestep.mape, decreasing = F),]
# sscomp$bsts.pre.onestep.mape.rank <- 1:nrow(sscomp)
# print(sscomp)
# 
# sscomp <- sscomp[order(sscomp$bsts.post.mape, decreasing = F),]
# sscomp$bsts.post.mape.rank <- 1:nrow(sscomp)
# print(sscomp)
# 
# sscomp <- sscomp[order(sscomp$bsts.pre.onestep.rmse, decreasing = F),]
# sscomp$bsts.pre.onestep.rmse.rank <- 1:nrow(sscomp)
# print(sscomp)
# 
# sscomp <- sscomp[order(sscomp$bsts.post.rmse, decreasing = F),]
# sscomp$bsts.post.rmse.rank <- 1:nrow(sscomp)
# print(sscomp)


# sscomp$rank.diff <- sscomp$pre.rank - sscomp$post.rank
# print(sscomp)

write.csv(sscomp, file = file.path(dir_ext,'bsts_default_state_space_14config_summary_pre-post.csv'))
# ########################## END ##########################################



