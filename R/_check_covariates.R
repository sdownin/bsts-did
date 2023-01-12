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
noise.level <- 1.5
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
expect.mod.sizes <- list(5)
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
  `8b`=c('AddLocalLevel','AddSeasonal', 'AddRegression'),
  `11b`=c('AddLocalLinearTrend','AddSeasonal', 'AddRegression')
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
effect.types = c('quadratic')  ## constant','geometric
##
bsts.niter.start <- 500 ## 5000
bsts.niter.max   <- 500 ## 8e4
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
                                      bsts.max.iter= bsts.niter.max,
                                      bsts.n.cov.cats = 3
)  ## D:\\BSTS_external


bsts.n.cov.cats## GET RESULT FROM STORAGE
. <- readRDS(simlist.files[[1]]$file[1])
bsts.model <- getBstsModelFromSimlist(list(.), key = 1, effect.type = 'quadratic')
par(mar=c(2,2.5,2,1))
plot(bsts.model, 'components')
plot(bsts.model, 'predictors')
plot(bsts.model, 'coefficients')


l <- readRDS(file.path(dir_ext, '__GRIDSEARCH_output__n100_pd520_niter30000_16734236256_d1f1g1h1i1j1l1.rds'))
# compare <- l$compare
# l$compare <- NULL
# l$compare <- compare
# compare <- NULL

# l$bsts.state.specs

CompareBstsModels(list(getBstsModelFromSimlist(list(l), key=1, effect.type = 'constant'),
                       getBstsModelFromSimlist(list(l), key=1, effect.type = 'geometric'),
                       getBstsModelFromSimlist(list(l), key=1, effect.type = 'quadratic')
                       ))

PlotBstsComponents(
  getBstsModelFromSimlist(list(l), key=1, effect.type = 'constant'),
)





# ## Get the  CausalImpact object from the BSTS vs. DiD comparison simlist object
# ##
# getCausalImpactFromSimlist <- function(simlist, key=NA, 
#                                        effect.type='quadratic',
#                                        state.space.list.id=1
# ) {
#   if (is.na(key)) {
#     if (length(names(simlist))==0) {
#       names(simlist) <- as.character( 1:length(simlist) )
#     }
#     key <- names(simlist)[1]
#   }
#   bstslist <- simlist[[key]]$compare$bsts[[ effect.type ]]
#   return( bstslist[[ state.space.list.id ]]$CausalImpact )
# }
# 
# 
# ##
# ##  Get the  CausalImpact object from the BSTS vs. DiD comparison simlist object
# ##
# getBstsModelFromSimlist <- function(simlist, key=NA, 
#                                     effect.type='quadratic',
#                                     state.space.list.id=1
















