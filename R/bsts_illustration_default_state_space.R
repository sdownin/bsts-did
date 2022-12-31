#####################################################
##
##   Dynamic Causal Inference Simulations
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
prefix <- 'bsts-illus_default-st-sp_'

##==============================
## Load simulation functions
##------------------------------
# source(file.path(dir_r,'internal_intervention_sim.R'))
source(file.path(dir_r,'single_intervention_sim_vec.R'))  ## Actor index vectorized simulation
##
source(file.path(dir_r,'bsts_helper_functions.R')) ## Setting up and adding state space components to state.space list


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
##  0.  FIND OPTIONAL STATE SPACE CONFIG 
##      TO SERVE AS DEFAULT 
##
##
##=======================================
## Static defaults
noise.level <- 1.5
b4 <- 1  ## seasonal component weight
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
  ## ## LEVEL
  `1`=c('AddLocalLevel'),
  ## ## TREND
  `2`=c('AddLocalLinearTrend'),
  # # `3`=c('AddStudentLocalLinearTrend'),
  # ## ## LEVEL + SLOPE ( + AR SLOPE DRIFT)
  `4`=c('AddSemilocalLinearTrend'),
  # ## ## SEASONAL
  `5`= c('AddTrig'),
  `5b`=c('AddSeasonal'),
  # ## ## AUTOCORRELATION
  # # list('AddAr'),
  `6`=c('AddAr'),
  # ##------ COMBINATIONS --------------
  # ## AR & SEASONALITY
  `7`= c('AddAr','AddTrig'),
  `7b`=c('AddAr','AddSeasonal'),
  # # `7a`=c('AddAutoAr','AddTrig')#,
  # ## LEVEL + ...
  `8`= c('AddLocalLevel','AddTrig'),
  `8b`=c('AddLocalLevel','AddSeasonal'),
  # # `9`=c('AddLocalLevel','AddAr'),
  # # `10`=c('AddLocalLevel','AddTrig','AddAutoAr'),
  # # ## (LEVEL + SLOPE) + ...
  `11`= c('AddLocalLinearTrend','AddTrig'),
  `11b`=c('AddLocalLinearTrend','AddSeasonal'),
  # `12`=c('AddLocalLinearTrend','AddAr'),
  # `12`= c('AddLocalLinearTrend','AddAr','AddTrig'),
  # `12b`=c('AddLocalLinearTrend','AddAr','AddSeasonal'),
  # `13`= c('AddStudentLocalLinearTrend','AddTrig'),
  # `13b`=c('AddStudentLocalLinearTrend','AddSeasonal')#,
  # `14`=c('AddStudentLocalLinearTrend','AddAr'),
  # # `14`c('AddStudentLocalLinearTrend','AddTrig','AddAr'),
  ## (LEVEL + SLOPE ( + AR1 SLOPE DRIFT)) + ...
  `15`= c('AddSemilocalLinearTrend','AddTrig'),
  `15b`=c('AddSemilocalLinearTrend','AddSeasonal')#,
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
                intpd = intpd, ## 60% pre-intervention training / 40% post-interventihttp://127.0.0.1:37489/graphics/plot_zoom_png?width=1200&height=900on
                noise.level = noise.level, ## stdev of simulated noise terms
                prior.sd.scenario = prior.sd.scenario, ## BSTS Prior SD scenario (high vs. low uncertainty in priors
                treat.rule = treat.rule, 
                treat.prob = ifelse(treat.rule=='random', 0.5, 1), 
                treat.threshold = ifelse(treat.rule=='random', 1, 0.5),
                b4 = b4,   ## seasonal component weight
                b5 = b5, ##
                b9 = dgp.ar  , ## autocorrelation
                seasonality = seasonality,
                dgp.nseasons= ifelse(seasonality, dgp.nseasons, NA), 
                dgp.freq= ifelse(seasonality, dgp.freq, NA),
                bsts.state.specs=bsts.state.specs,
                expect.mod.size=expect.mod.size,
                # bsts.state.specs=list(list(AddSemilocalLinearTrend),list(AddSemilocalLinearTrend,AddStudentLocalLinearTrend)),
                rand.seed = 7531
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
bsts.niter <- 1e4
## ID for the simulation (to search/filter all simulation figures, RDS files, etc.)
sim.id <- round(10*as.numeric(Sys.time()))


## RUN SIMULATION -  SIMULATE TIME SERIES
simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               sim.id = sim.id,
                               plot.show = T, plot.save = T )

## RUN BSTS and compare to DID
simlist.files <- runSimCompareBstsDiD(simlist, 
                                      effect.types = effect.types,
                                      sim.id = sim.id,
                                      save.items.dir= dir_ext,
                                      bsts.niter=bsts.niter * 5
)  ## D:\\BSTS_external







# 'D:\\BSTS_external'

filename <- '__GRIDSEARCH_output__n200_pd520_niter50000_16724204820_d1f1g1h1i1j1l1.rds'
simlist <- readRDS( file.path( dir_ext, filename ) )


x <- simlist$compare$bsts$constant
sslist <- lapply(1:length(x), function(i) x[[i]]$CausalImpact$model$bsts.model )
CompareBstsModels( sslist )

z <- rowSums(abs(sslist[[1]]$one.step.prediction.errors[,1:(intpd-1)]))

res <- t( sapply(1:length(x), function(i) summary(rowSums(abs(sslist[[i]]$one.step.prediction.errors[,1:(intpd-1)]))) ) )
print(round(res, 3))

sscomp <- abs(as.data.frame(res))
sscomp$ss.id <- 1:nrow(sscomp)
sscomp$ss.comps <- sapply(1:nrow(sscomp), function(i) paste(st.sp.lists[[i]], collapse = '|') )

res.tbl <- simlist$compare$res.tbl$constant

## PRE ----
x <- simlist$compare$bsts$constant
ss.fits <- lapply(1:length(st.sp.lists), function(i){
  cat(sprintf(' %s ',i))
  idx <- 1:(intpd-1)
  mod <- x[[i]]$CausalImpact$model$bsts.model
  newdata <- cbind(treatment_y_outcome=mod$original.series[idx], mod$predictors[idx,])
  burn <- round( .2 * mod$niter )
  fit <- predict(mod, newdata=newdata, horizon = 1, burn = burn) ## one-step ahead prediction
  yhat.t <- fit$mean
  y.t <- fit$original.series[idx]
  return(list(yhat.t=yhat.t, y.t=y.t, fit=fit, mod=mod))
})
saveRDS(ss.fits, file=file.path(dir_ext, 'bsts_default_state_space_comparison_sMAPE.rds'))
## alias for constant MCC curve shape state space list
x <- simlist$compare$bsts$constant
## sMAPE - symmetric MAPE
sscomp$bsts.pre.onestep.smape <- lapply(1:length(ss.fits), function(l){
  yhat.t <- l['yhat.t']
  y.t <- l['y.t']
  100 * mean( abs(yhat.t - y.t) / (abs(y.t) - abs(yhat.t))  , na.rm=T)
})
sscomp$bsts.pre.onestep.mape <- lapply(1:length(ss.fits), function(l){
  yhat.t <- l['yhat.t']
  y.t <- l['y.t']
  100 * mean( abs( (yhat.t - y.t) / y.t ), na.rm=T )
})
sscomp$bsts.pre.onestep.mae <- lapply(1:length(ss.fits), function(l){
  err <- l['yhat.t'] - l['y.t']
  mean( abs( err, na.rm=T ) )
})
sscomp$bsts.pre.onestep.mse <- lapply(1:length(ss.fits), function(l){
  err <- l['yhat.t'] - l['y.t']
  mean( err^2, na.rm=T )
})
sscomp$bsts.pre.onestep.rmse <- lapply(1:length(ss.fits), function(l){
  err <- l['yhat.t'] - l['y.t']
  sqrt( mean( err^2, na.rm=T ) )
})
# ## MAPE - Mean Absolute Percentage Error
# sscomp$bsts.pre.onestep.mape <- sapply(1:nrow(sscomp), function(i){
#   # idx <- 1:(intpd-1)
#   # # eps.mat <- sslist[[i]]$one.step.prediction.errors[,idx]
#   # ## pointwise mean of pre-intervention 1-step ahead prediction error
#   # eps.t <- colMeans( sslist[[i]]$one.step.prediction.errors[,idx], na.rm=T )
#   # ## pointwise observations pre-intervention
#   # y.t <-  sslist[[i]]$original.series[idx] 
#   # ## mean absolute percentage error
#   # 100 * mean( abs(eps.t / y.t), na.rm=T)
# })
## RMSE - Root Mean Square Error
sscomp$bsts.pre.onestep.rmse <- sapply(1:nrow(sscomp), function(i){
  idx <- 1:(intpd-1)
  att <- res.tbl[[i]]$b3.att[idx]
  eps <- res.tbl[[i]]$bsts.point.effect[idx] - att
  sqrt( mean( eps^2, na.rm=T) )
})

## POST --

## MAPE - Mean Absolute Percentage Error
sscomp$bsts.post.mape <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  att <- res.tbl[[i]]$b3.att[idx]
  eps <- res.tbl[[i]]$bsts.point.effect[idx] - att
  mean( abs(eps / att), na.rm=T)
})
## RMSE - Root Mean Square Error
sscomp$bsts.post.rmse <- sapply(1:nrow(sscomp), function(i){
  idx <- intpd:npds
  att <- res.tbl[[i]]$b3.att[idx]
  eps <- res.tbl[[i]]$bsts.point.effect[idx] - att
  sqrt( mean( eps^2, na.rm=T) )
})



sscomp <- sscomp[order(sscomp$bsts.pre.onestep.mape, decreasing = F),]
sscomp$bsts.pre.onestep.mape.rank <- 1:nrow(sscomp)
print(sscomp)

sscomp <- sscomp[order(sscomp$bsts.post.mape, decreasing = F),]
sscomp$bsts.post.mape.rank <- 1:nrow(sscomp)
print(sscomp)

sscomp <- sscomp[order(sscomp$bsts.pre.onestep.rmse, decreasing = F),]
sscomp$bsts.pre.onestep.rmse.rank <- 1:nrow(sscomp)
print(sscomp)

sscomp <- sscomp[order(sscomp$bsts.post.rmse, decreasing = F),]
sscomp$bsts.post.rmse.rank <- 1:nrow(sscomp)
print(sscomp)


# sscomp$rank.diff <- sscomp$pre.rank - sscomp$post.rank
# print(sscomp)

write.csv(sscomp, file = file.path(dir_ext,'bsts_default_state_space_14config_summary_pre-post.csv'))
# ########################## END ##########################################
