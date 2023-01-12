#####################################################
##
##   PERIOD AGGREGATION - GRANULARITY 
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
library(stringdist)
##
library(RColorBrewer)
library(scico)
library(ggplot2)

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
prefix <- 'bsts-illus_pd-agg_'

##==============================
## Load simulation functions
##------------------------------
# source(file.path(dir_r,'internal_intervention_sim.R'))  ## previous version

## Actor index vectorized simulation
source(file.path(dir_r,'single_intervention_sim_vec.R')) 
## Setting up and adding state space components to state.space list
source(file.path(dir_r,'bsts_helper_functions.R')) 
## BSTS vs. DiD comparison and sensitivity analysis
source(file.path(dir_r,'bsts_did_comparison_functions.R')) 




################################################################################################################
################################################################################################################







##=======================================
##
##
##  1.  PERIOD AGGREGATION
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
ns <- list(200, 400) ## 200 ## list( 100, 200 ) ## ** INTERACT SIZE (N) with GRANULARITY
sim.lengths <- list(520)
treat.rules <- list('random')  ## 'below.benchmark'
seasonalities <- list(TRUE)   ## c(TRUE,  FALSE )
prior.sd.scenarios <- list('sd.low') ## list('sd.low','sd.high')  ## sd.low
expect.mod.sizes <- list(2)
## FOCAL CONSTRUCT
dgp.ars <- list(0)  ## 0.6  ## .1,.2,.4
## STATE SPACE CONFIGURATIONS
st.sp.lists <- list(
  `8b`=c('AddLocalLevel','AddSeasonal')
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
effect.types <- c('constant','geometric','quadratic') ##  c('constant','geometric','quadratic') ## c('quadratic')  
## Scale bsts iterations in increments (doubling) starting at 10k (more iterations for models that converge more slowly)
bsts.niter.start <- 5000  ## 10k
bsts.niter.max   <- 8e4  ## 80k
##
sim.id <- round(10*as.numeric(Sys.time()))




## RUN SIMULATION -  SIMULATE TIME SERIES
simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               sim.id = sim.id,
                               plot.show = F, plot.save = FALSE )

## RUN BSTS and compare to DID
simlist.files <- runSimCompareBstsDiD(simlist, 
                                      effect.types = effect.types,
                                      sim.id = sim.id,
                                      save.items.dir= dir_ext,
                                      bsts.niter = bsts.niter.max, ## **START at MAX niter for large npds **
                                      bsts.max.iter= bsts.niter.max
                                      )  ## D:\\BSTS_external





# ## LOAD INDIVIDUAL SIMULATION COMPARISON LIST
# key <- sprintf('d%sf%sg%sh%si%sj%s', 1,1,1,1,1,1)
# simx <- readRDS(file.path(dir_ext,sprintf('__GRIDSEARCH_output__%s_%s.rds',simlist[[1]]$sim$id, key )))


####################################################
##---------------------
## ORIGINAL DATA 520 pds  (Weekly data,  10 years)
##---------------------
##----------------------
## Aggregate every 52 periods (yearly) == 10
##---------------------
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 52 )
simlistx.files <- runSimCompareBstsDiD(simlistx, 
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max) 
##----------------------
## Aggregate every 40 periods (every 40 weeks) == 13
##---------------------
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 40 )
simlistx.files <- runSimCompareBstsDiD(simlistx, 
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max ) 
##----------------------
## Aggregate every 26 periods (half-yearly or semi-yearly) == 20
##---------------------
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 26 )
simlistx.files <- runSimCompareBstsDiD(simlistx, 
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max ) 
##----------------------
## Aggregate every 20 periods (every 20 weeks) == 26
##---------------------
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 20 )
simlistx.files <- runSimCompareBstsDiD(simlistx, 
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max) 
##----------------------
## Aggregate every 10 periods (every 10 weeks) == 52
##---------------------
# simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
#                                plot.show = F, plot.save = F )
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 10 )
simlistx.files <- runSimCompareBstsDiD(simlistx, 
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max) 
##----------------------
## Aggregate every 8 periods (bi-monthly, every 2 months) == 65
##---------------------
# simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
#                                plot.show = F, plot.save = F )
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 8 )
simlistx.files <- runSimCompareBstsDiD(simlistx,
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max)
##----------------------
## Aggregate every 5 periods (every 5 weeks) == 104
##---------------------
# simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
#                                plot.show = F, plot.save = F )
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 5 )
simlistx.files <- runSimCompareBstsDiD(simlistx,
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max)
##----------------------
## Aggregate every 4 periods (monthly, every 4 weeks) == 130
##---------------------
# simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
#                                plot.show = F, plot.save = F )
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 4 )
simlistx.files <- runSimCompareBstsDiD(simlistx, 
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max) 
##----------------------
## Aggregate every 2 periods (biweekly, every 2 weeks) == 260
##---------------------
# simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
#                                plot.show = F, plot.save = F )
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 2 )
simlistx.files <- runSimCompareBstsDiD(simlistx,
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.max,  ## **START at MAX niter for large npds **
                                       bsts.max.iter= bsts.niter.max)

## *** CAUTION *** CAN THROW ERRORS(?) - RUN LAST
##----------------------
## Aggregate every 104 periods (biannually) == 5
##---------------------
simlistx <- simlist
simlistx <- runSimUpdateSimlist(simlistx, effect.types=effect.types, sim.id=sim.id, plot.show=F, plot.save=F)
simlistx <- updateSimlistAggregateSimDfPd(simlistx, pd.agg = 104 )
simlistx.files <- runSimCompareBstsDiD(simlistx, 
                                       effect.types = effect.types,
                                       save.items.dir= dir_ext,
                                       sim.id = sim.id,
                                       bsts.niter = bsts.niter.start,
                                       bsts.max.iter= bsts.niter.max )












##############################################
##
##
##
## AGGREGATION x MCC CURVE INTERACTION
##
##
##
##############################################
# gsdir  <-  'D:\\BSTS_external\\MCC_curve_interact_aggregation_comparison'  ## 'D:\\BSTS_external\\bsts_did_comparison_summary_gridsearch'
gsdir <- file.path(dir_ext, 'bsts_period_aggregation_MCC_curve_interaction_20230103')
# gsdir <- 'D:\\BSTS_external\\MCC_curve_interact_aggreagation_comparision_n200-400'

files <- dir(gsdir, pattern = '\\.rds$')
suml <- list()
sdf <- data.frame(stringsAsFactors = F)
effect.types <- c('constant','quadratic','geometric')
for (i in 1:length(files)) {
  
  x <- readRDS(file.path(gsdir, files[i]))
  
  intpd <- x$intpd
  npds <- x$npds
  n <- x$n
  npds.post <- npds - intpd + 1
  
  for (k in 1:length(effect.types)) {
    cat(sprintf('i=%s %s k=%s\n',i,files[i],k))
    
    effect.type <- effect.types[k]    
    
    parallel.p <- x$compare$did[[effect.type]]$attgt$Wpval
    parallel.W <- x$compare$did[[effect.type]]$attgt$W
    
    nss <- length(x$compare$bsts[[effect.type]])
    
    for (h in 1:length(nss)) {
      
      # x$compare$bsts[[effect.type]][[h]]$CausalImpact
      
      # x$compare$res.tbl
      
      res.tbl <- x$compare$res.tbl[[effect.type]][[h]]
      
      
      npds.int <- length(intpd:npds)
      
      dsig <- 1 * (res.tbl$did.point.conf.low[intpd:npds] * res.tbl$did.point.conf.high[intpd:npds] > 0 ) ## positive product means both same sign, not containing zero --> significant
      bsig <- 1 * (res.tbl$bsts.point.effect.lower[intpd:npds] * res.tbl$bsts.point.effect.upper[intpd:npds] > 0)
      
      did.bsts.sigvec.dist.cos <- sum(dsig * bsig) / (sqrt(sum(dsig^2))*sqrt(sum(bsig^2)))
      did.bsts.sigvec.dist.ham <- hamming.distance(dsig, bsig)
      
      ## time until first significance
      ##   0 = significant in 1st treatment period (at t_event = 0)
      tsig.1.did <- ifelse(sum(dsig==1)>0, which(dsig == 1)[1]  , NA)
      tsig.1.bsts <- ifelse(sum(bsig==1)>0 , which(bsig == 1)[1]  , NA)
      
      tsig.1.tnotsig.1.did <- NA
      tsig.1.tnotsig.1.bsts <- NA
      
      ## time from 1st significance to 1st insignificance 
      if (!is.na(tsig.1.did)) {
        id.1notsig.after.1sig.did <- which( dsig[(tsig.1.did+1):npds.int] == 0 )
        tsig.1.tnotsig.1.did <- ifelse(length(id.1notsig.after.1sig.did)>0, id.1notsig.after.1sig.did[1] + tsig.1.did , NA)
      }
      if (!is.na(tsig.1.bsts)) {
        id.1notsig.after.1sig.bsts <- which( bsig[(tsig.1.bsts+1):npds.int] == 0 )
        tsig.1.tnotsig.1.bsts <- ifelse(length(id.1notsig.after.1sig.bsts)>0, id.1notsig.after.1sig.bsts[1] + tsig.1.bsts , NA)
      }
      
      
      ## longest pointwise significance span
      tsig.long.did <- stringdist(paste(dsig,collapse=''), 
                                  paste(rep(1, length(dsig)),collapse=''), 
                                  method='lcs')
      tsig.long.bsts <- stringdist(paste(bsig,collapse=''), 
                                   paste(rep(1, length(bsig)),collapse=''), 
                                   method='lcs')
      
      ## Period of last pointwise signifance
      tsig.last.did <- max(which(dsig == 1))
      tsig.last.bsts<- max(which(bsig == 1))
      
      # x$compare$bsts[[effect.type]][[h]]$CausalImpact$summary
      # x$compare$bsts[[effect.type]][[h]]$CausalImpact$series
      
      st.sp.l <- x$compare$bsts[[effect.type]][[h]]$CausalImpact$model$bsts.model$state.specification
      st.sp.str <- paste(sapply(st.sp.l, function(x)x$name), collapse = '|')
      
      
      if ( is.null(parallel.p) ) {
        # cat(sprintf(' skiping i=%s, parallel trends test p-val is missing (NULL)\n',i))
        # next
        parallel.p <- NA
      }
      if (is.null(parallel.W)) {
        parallel.W <- NA
      }
      
      att.dgp <- mean( res.tbl$b3.att[intpd:npds], na.rm = T)
      att.dgp.tot <- sum( res.tbl$b3.att[intpd:npds], na.rm=T)
      
      causimp <- getCausalImpactFromSimlist(list(x), key=1, effect.type=effect.type) 
      
      att.bsts <- causimp$summary$AbsEffect[1]
      att.bsts.tot <-  causimp$summary$AbsEffect[2]
      
      att.bsts.lo <- mean(causimp$series$point.effect.lower[intpd:npds], na.rm=T)
      att.bsts.lo.tot <- causimp$series$cum.effect.lower[npds]
      att.bsts.hi <- mean(causimp$series$point.effect.upper[intpd:npds], na.rm=T)
      att.bsts.hi.tot <- causimp$series$cum.effect.upper[npds]
      
      # b3.bsts.tot <- sum( res.tbl$bsts.point.effect[intpd:npds], na.rm=T)
      
      # b3.did <- x$compare$did[[effect.type]]$agg.es$overall.att
      
      attgt <- x$compare$did[[effect.type]]$attgt
      agg.group <- aggte(attgt, type = 'group', bstrap = T)
      
      att.did <- agg.group$overall.att
      att.did.tot <- sum(res.tbl$did.estimate[intpd:npds], na.rm=T)
      
      att.did.lo <- mean(res.tbl$did.point.conf.low[intpd:npds], na.rm=T)
      att.did.lo.tot <- sum(res.tbl$did.point.conf.low[intpd:npds], na.rm=T)
      
      att.did.hi <- mean(res.tbl$did.point.conf.high[intpd:npds], na.rm=T)
      att.did.hi.tot <- sum(res.tbl$did.point.conf.high[intpd:npds], na.rm=T)
      
      att.did.pre.rmse <- sqrt(mean( (res.tbl$did.estimate[1:(intpd-1)] - res.tbl$b3.att[1:(intpd-1)])^2, na.rm=T ))
      att.did.pre.abs.err.tot <- sum( abs(res.tbl$did.estimate[1:(intpd-1)] - res.tbl$b3.att[1:(intpd-1)]), na.rm=T)
      
      att.did.post.rmse <- sqrt(mean( (res.tbl$did.estimate[intpd:npds] - res.tbl$b3.att[intpd:npds])^2, na.rm=T ))
      att.did.post.abs.err.tot <- sum( abs(res.tbl$did.estimate[intpd:npds] - res.tbl$b3.att[intpd:npds]), na.rm=T)
      
      att.did.pre.mape <- mean( abs( (res.tbl$did.estimate[1:(intpd-1)] - res.tbl$b3.att[1:(intpd-1)]) / res.tbl$did.estimate[1:(intpd-1)] ), na.rm=T )
      
      bsts.posterior.p <- causimp$summary$p[1]
      # bsts.mape <- 
      bsts.post.rmse <- sqrt(mean( (causimp$series$response[intpd:npds] - causimp$series$point.pred[intpd:npds])^2 )) 
      bsts.post.abs.err.tot <- sum( abs(causimp$series$response[intpd:npds] - causimp$series$point.pred[intpd:npds]), na.rm=T)
      
      bsts.pre.rmse <- sqrt(mean( (causimp$series$response[1:(intpd-1)] - causimp$series$point.pred[1:(intpd-1)])^2 ))
      bsts.pre.abs.err.tot <- sum( abs(causimp$series$response[1:(intpd-1)] - causimp$series$point.pred[1:(intpd-1)]) , na.rm=T)
      
      bsts.pre.1step <- causimp$model$bsts.model$one.step.prediction.errors[,1:(intpd-1)]
      bsts.pre.1step.tot <- sum(bsts.pre.1step, na.rm=T)
      bsts.pre.1step.rmse <- sqrt(mean( colMeans(bsts.pre.1step)^2 ))
      # plot(rowMeans(bsts.pre.1step),type='l')
      # colMeans(bsts.pre.1step)
      
      bsts.pre.mape <- mean( abs( (causimp$series$response[1:(intpd-1)] - causimp$series$point.pred[1:(intpd-1)]) / causimp$series$response[1:(intpd-1)] ), na.rm=T )
      
      bsts.niter <- causimp$model$bsts.model$niter
      
      bsts.cnt.prior.inclu <- sum( causimp$model$bsts.model$prior$prior.inclusion.probabilities > .5, na.rm=T)
      
      
      ## Period of first cumulative significance (is there an effect at all? If yes, when? )
      agg.es <- x$compare$did[[effect.type]]$agg.es
      did.tbl <- tidy(agg.es)
      did.tbl <- did.tbl[which(did.tbl$event.time >= 0),]
      idx.did.cumu.sig <- which( cumsum(did.tbl$conf.low) > 0 ) ## lower confidence interval > 0 == significant 
      tsig.cumu.1.did <- ifelse( length(idx.did.cumu.sig)>0, idx.did.cumu.sig[1], NA)
      
      bsts.tbl <- causimp$series[intpd:npds, ]
      idx.bsts.cumu.sig <- which( bsts.tbl$cum.effect.lower > 0  )
      tsig.cumu.1.bsts <- ifelse( length(idx.bsts.cumu.sig)>0, idx.bsts.cumu.sig[1], NA)
      
      # (x$compare$did[[effect.type]]$agg.es$crit.val.egt * x$compare$did[[effect.type]]$agg.es$overall.se)
      
      idf <- data.frame(file=files[i], intpd=intpd, npds=npds, n=n, npds.per.n=npds/n,
                        npds.postint=npds.int,
                        effect.type=effect.type, st.sp=st.sp.str, bsts.niter=bsts.niter,
                        parallel.p=parallel.p, parallel.W=parallel.W, 
                        att.bsts=att.bsts,  
                        att.bsts.lo=att.bsts.lo,  
                        att.bsts.hi=att.bsts.hi, 
                        att.bsts.tot=att.bsts.tot, 
                        att.bsts.lo.tot=att.bsts.lo.tot, 
                        att.bsts.hi.tot=att.bsts.hi.tot,
                        #
                        bsts.post.rmse=bsts.post.rmse,  
                        bsts.post.abs.err.tot=bsts.post.abs.err.tot,
                        bsts.pre.rmse=bsts.pre.rmse,
                        bsts.pre.abs.err.tot=bsts.pre.abs.err.tot,
                        bsts.pre.1step.tot=bsts.pre.1step.tot,
                        bsts.pre.1step.rmse=bsts.pre.1step.rmse,
                        #
                        bsts.pre.mape=bsts.pre.mape,
                        ##
                        bsts.posterior.p=bsts.posterior.p, 
                        bsts.niter=bsts.niter,
                        bsts.cnt.prior.inclu=bsts.cnt.prior.inclu,
                        ##
                        att.did=att.did, 
                        att.did.tot=att.did.tot,
                        att.did.lo=att.did.lo, 
                        att.did.lo.tot=att.did.lo.tot,
                        att.did.hi=att.did.hi, 
                        att.did.hi.tot=att.did.hi.tot,
                        ##
                        att.did.pre.rmse=att.did.pre.rmse,
                        att.did.pre.abs.err.tot=att.did.pre.abs.err.tot,
                        att.did.post.rmse=att.did.post.rmse,
                        att.did.post.abs.err.tot=att.did.post.abs.err.tot,
                        #
                        att.did.pre.mape=att.did.pre.mape,
                        ##
                        did.bsts.sigvec.dist.cos=did.bsts.sigvec.dist.cos,
                        did.bsts.sigvec.dist.ham=did.bsts.sigvec.dist.ham,
                        did.bsts.sigvec.dist.cos.scaled= did.bsts.sigvec.dist.cos / length(intpd:npds),
                        did.bsts.sigvec.dist.ham.scaled= did.bsts.sigvec.dist.ham / length(intpd:npds),
                        ## period of 1st pointwise significance
                        tsig.1.did=tsig.1.did,
                        tsig.1.bsts=tsig.1.bsts,
                        tsig.1.did.lag = tsig.1.did - 1,  ## delay before 1st pointwise significance
                        tsig.1.bsts.lag = tsig.1.bsts - 1, ## delay before 1st pointwise significance
                        ## period of last pointwise significance
                        tsig.last.did = tsig.last.did, 
                        tsig.last.bsts = tsig.last.bsts,
                        ## longest pointwise significance period
                        tsig.long.did=tsig.long.did,
                        tsig.long.bsts=tsig.long.bsts,
                        ## 1st not significant period after 1st significant period - For Use in Geometric effect precision illustration
                        tsig.1.tnotsig.1.did=tsig.1.tnotsig.1.did,
                        tsig.1.tnotsig.1.bsts=tsig.1.tnotsig.1.bsts,
                        ##
                        tsig.cumu.1.did = tsig.cumu.1.did,
                        tsig.cumu.1.bsts = tsig.cumu.1.bsts,
                        ## NO ENTROPY? Check if still pointwise significant at end
                        did.no.entropy = ( tsig.last.did == npds.int ),
                        bsts.no.entropy = ( tsig.last.bsts == npds.int ),
                        ##
                        bias.bsts=(att.bsts - att.dgp), 
                        bias.did=(att.did - att.dgp),
                        bsts.advantage=(abs(att.did - att.dgp) / abs(att.bsts - att.dgp))
      )
      sdf <- rbind(sdf, idf)
      
    }
    
  }
  
}
##
sdf$agg.span <- max(sdf$npds) / sdf$npds
sdf$agg.granularity <- sdf$npds / max(sdf$npds) 
sdf$agg.granularity.log <- log(sdf$agg.granularity)
##
write.csv(sdf, file=file.path(gsdir,'bsts_period_aggregation_MCC_curve_interaction_n200-400_SUMMARY.csv'))
# cdf <- read.csv(file.path(gsdir,'bsts_did_period_aggregation_summary.csv'))
# cdf <- read.csv(file.path(gsdir,'_prev','bsts_did_period_aggregation_summary.csv'))
# # att.did <- agg.es$overall.att
# # att.bsts <- impact_amount$summary$AbsEffect[1]

# cdf$perf.adv.model <- 'DID'
# cdf$perf.adv.model[cdf$bsts.advantage > 1] <- "BSTS"
# cdf$total.obs <- cdf$n * cdf$npds
# cdf$actors <- as.factor(cdf$n)

ns <- c(200, 400)
aggcor <- data.frame()
for (i in 1:length(effect.types)) {
  effect.type <- effect.types[i]
  for (j in 1:length(ns)) {
    n <- ns[j]
    idx <- which( sdf$effect.type==effect.type & sdf$n == n )
    sdfx <- sdf[idx, ]
    ##
    aggcor <- rbind(aggcor, data.frame(
      nobs.cor = nrow(sdfx),
      nactors=n,
      effect.type=effect.type,
      gran.mean = mean( sdfx$agg.granularity, na.rm=T),
      gran.sd = sd( sdfx$agg.granularity, na.rm=T),
      
      gran.log.mean = mean( sdfx$agg.granularity.log, na.rm=T),
      gran.log.sd = sd( sdfx$agg.granularity.log, na.rm=T),
      
      bsts.niter.conv.mean = mean(sdfx$bsts.niter, na.rm=T),
      bsts.niter.conv.sd = sd(sdfx$bsts.niter, na.rm=T),
      bsts.adv.mean = mean( abs(sdfx$bsts.advantage) , na.rm=T),
      bsts.adv.sd = sd( abs(sdfx$bsts.advantage) , na.rm=T),
      
      bias.bsts.mean = mean( sdfx$bias.bsts, na.rm=T),
      bias.bsts.sd = sd( sdfx$bias.bsts, na.rm=T),
      bias.did.mean = mean( sdfx$bias.did, na.rm=T),
      bias.did.sd = sd( sdfx$bias.did, na.rm=T),
      
      cor.bsts.adv.gran= cor(abs(sdfx$bsts.advantage), sdfx$agg.granularity),
      cor.bsts.adv.gran.testp = cor.test(abs(sdfx$bsts.advantage), sdfx$agg.granularity)$p.value,
      
      cor.bsts.bias.gran= cor(sdfx$bias.bsts, sdfx$agg.granularity), 
      cor.bsts.bias.gran.testp= cor.test(sdfx$bias.bsts, sdfx$agg.granularity)$p.value, 
      
      cor.did.bias.gran= cor(sdfx$bias.did, sdfx$agg.granularity),
      cor.did.bias.gran.testp =   cor.test(sdfx$bias.did, sdfx$agg.granularity)$p.value, 
      
      ##---
      
      cor.bsts.adv.gran.log= cor(abs(sdfx$bsts.advantage), sdfx$agg.granularity.log),
      cor.bsts.adv.gran.log.testp = cor.test(abs(sdfx$bsts.advantage), sdfx$agg.granularity.log)$p.value,
      
      cor.bsts.bias.gran.log= cor(sdfx$bias.bsts, sdfx$agg.granularity.log), 
      cor.bsts.bias.gran.log.testp= cor.test(sdfx$bias.bsts, sdfx$agg.granularity.log)$p.value, 
      
      cor.did.bias.gran.log= cor(sdfx$bias.did, sdfx$agg.granularity.log),
      cor.did.bias.gran.log.testp =   cor.test(sdfx$bias.did, sdfx$agg.granularity.log)$p.value
    ))
  }
}
write.csv(aggcor, file=file.path(gsdir, 'bsts_did_period_aggregation_correlation_MCC_curve_interaction.csv'))




########################################
########################################



# ggplot(data=sdf, aes(x=npds, y=parallel.p+1, colour=effect.type)) + 
#   geom_point() + facet_wrap( . ~ factor(n))  + 
#   geom_smooth(method='lm') +
#   scale_y_log10()




















########################################
##  PERIOD AGGREGATION - SUMMARIZE
##   N=200,400
##   T=[5] 10, 13, 20, 26, 52, 65, 104, 130, 260, 520
########################################
# gsdir  <-  'D:\\BSTS_external\\bsts_did_aggregation_comparison'  ## 'D:\\BSTS_external\\bsts_did_comparison_summary_gridsearch'
gsdir <- file.path(dir_ext, 'bsts_period_aggregation_MCC_curve_interaction_20230103')
files <- dir(gsdir, pattern = '\\.rds$')
suml <- list()
cdf <- data.frame(stringsAsFactors = F)
effect.type <- 'quadratic'
for (i in 1:length(files)) {
  cat(sprintf('i=%s %s\n',i,files[i]))
  x <- readRDS(file.path(gsdir, files[i]))
  intpd <- x$intpd
  npds <- x$npds
  n <- x$n
  p <- x$compare$did[[ effect.type ]]$attgt$Wpval
  W <- x$compare$did[[ effect.type ]]$attgt$W
  
  bsts.model <- getBstsModelFromSimlist(list(x), key=1, effect.type = effect.type)
  
  causimp <- getCausalImpactFromSimlist(list(x), effect.type = effect.type)
  
  
  if ( is.null(p) ) {
    cat(sprintf(' skiping i=%s, parallel trends test p-val is missing (NULL)\n',i))
    next
  }
  
  b3.bsts<-  causimp$summary$AbsEffect[1]
  b3.tot.bsts <-  causimp$summary$AbsEffect[2]
  
  b3.did <- x$compare$did[[ effect.type ]]$agg.es$overall.att
  
  b3.dgp <- mean( x$compare$res.tbl[[ effect.type ]][[1]]$b3.att[intpd:npds], na.rm = T)
  
  idf <- data.frame(file=files[i], intpd=intpd, npds=npds, n=n, bsts.niter=bsts.model$niter,
                    parallel.p=p, parallel.W=W, b3.bsts=b3.bsts,
                    b3.tot.bsts=b3.tot.bsts, b3.did=b3.did, 
                    bias.bsts=(b3.bsts - b3.dgp), bias.did=(b3.did - b3.dgp),
                    bsts.advantage=(abs(b3.did - b3.dgp) / abs(b3.bsts - b3.dgp))
                    )
  cdf <- rbind(cdf, idf)
}
write.csv(cdf, file=file.path(gsdir,'bsts_did_period_aggregation_summary.csv'))
# cdf <- read.csv(file.path(gsdir,'bsts_did_period_aggregation_summary.csv'))
# cdf <- read.csv(file.path(gsdir,'_prev','bsts_did_period_aggregation_summary.csv'))
# # att.did <- agg.es$overall.att
# # att.bsts <- impact_amount$summary$AbsEffect[1]


read.csv(file.path(gsdir,'bsts_did_period_aggregation_summary.csv'))


cdf$perf.adv.model <- 'DID'
cdf$perf.adv.model[cdf$bsts.advantage > 1] <- "BSTS"
cdf$total.obs <- cdf$n * cdf$npds
cdf$actors <- as.factor(cdf$n)

########################################
########################################


ggplot(data=cdf, aes(x=npds/n, y=parallel.W, colour=perf.adv.model)) + 
  geom_point(aes(size=npds, shape=actors)) + 
  scale_fill_brewer() +
  scale_size(range = c(3,14)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  geom_smooth(method = "lm", alpha=.1) +
  scale_y_log10() + scale_x_log10() +
  theme_bw() + labs(colour='Performance\nAdvantage', size='Periods', shape='Actors')
ggsave(filename='did_vs_bsts_compare_pd_aggregation_wald-stat_npds-n_colour_shape_size.png')

ggplot(data=cdf, aes(x=npds/n, y=parallel.W)) + 
  geom_point(aes(size=total.obs, colour=perf.adv.model)) + 
  scale_fill_brewer() +
  scale_size(range = c(3,14)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  geom_smooth(method = "lm", alpha=.1) +
  scale_y_log10() + scale_x_log10() +
  theme_bw() + labs(colour='Performance\nAdvantage', size='Observations (N*T)')
ggsave(filename='did_vs_bsts_compare_pd_aggregation_wald-stat_npds-n.png')

ggplot(data=cdf, aes(x=npds/n, y=parallel.W)) + 
  geom_point(aes(size=total.obs, colour=perf.adv.model)) + 
  scale_fill_brewer() +
  scale_size(range = c(3,14)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  # geom_smooth(method = "lm", alpha=.1) +
  ylab('Parallel Trends Pretest - Wald Chisq') + xlab('T/N (Num. Periods Per Actor)') +
  scale_y_log10(limits=c(2,4000)) + scale_x_log10() +
  theme_bw() + labs(colour='Performance\nAdvantage', size='Observations (N*T)')
ggsave(filename='did_vs_bsts_compare_pd_aggregation_wald-stat_npds-n_NOLINE.png')


w.qs <- quantile(cdf$npds, probs = c(.333, .67, 1))
cdf$n.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<w.qs[1],'1Low',ifelse(cdf$npds[i]<w.qs[2],'2Mid','3High')))
ggplot(data=cdf, aes(x=npds/n, y=parallel.W, colour=perf.adv.model)) + 
  geom_point(aes(size=n)) + 
  scale_fill_brewer() + facet_wrap( . ~ n.cat) +
  scale_size(range = c(3,14)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  # geom_smooth(method = "lm", alpha=.1) +
  ylab('Parallel Trends Pretest - Wald Chisq') + xlab('N (Num. Actor)') +
  scale_y_log10(limits=c(2,10000)) + scale_x_log10() +
  theme_bw() #+ labs(colour='Performance\nAdvantage', size='Observations (N*T)')




# p.qs <- quantile(cdf$parallel.p, probs = c(1/4, 3/4, 1))
# npds.qs <- quantile(cdf$npds, probs = c(1/2, 1))
# cdf$W.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i]<w.qs[1],'1Low',ifelse(cdf$parallel.W[i]<w.qs[2],'2Mid','3High')))
# cdf$W.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i]<30,'W<30','W>=30'))
cdf$p.cat    <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.p[i]<.001,'1 Not Parallel p<0.001',ifelse(cdf$parallel.p[i]<0.5,'2 Maybe Parallel .001<p<.5','3 Parallel p>.5')))
cdf$n.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$n[i] <= 100, 'N<=200','N>200'))
cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<20,'1 Short Series T<20','2 Long Series T>=20'))
# cdf$w.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i]<20,'3Low W',ifelse(cdf$parallel.W[i]<60,'2Mid W','1High W')))
# cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<npds.qs[1],'1Low T',ifelse(cdf$npds[i]<npds.qs[2],'2Mid T','3High T')))
# cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<npds.qs[1],'1Low T',ifelse(cdf$npds[i]<npds.qs[2],'2Mid T','3High T')))
cdf.sub <- cdf[cdf$n>50, ]
ggplot(data=cdf.sub, aes(x=npds/n, y=parallel.W)) + 
  geom_point(aes(size=n, colour=perf.adv.model), alpha=.7) + 
  scale_fill_brewer() + 
  facet_grid( p.cat ~ npds.cat  ) +
  scale_size(range = c(3,10)) +
  # scale_color_manual(values=c('#00B050','#0070C0','#619CFF')) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  # geom_smooth(method = "lm", alpha=.1) +
  # ylab('Parallel Trends Pretest - Wald Chisq') + xlab('N (Num. Actor)') +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() #+ labs(colour='Performance\nAdvantage', size='Observations (N*T)')




# cdf$p.cat    <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.p[i]<.001,'1 Not Parallel p<0.001',ifelse(cdf$parallel.p[i]<0.5,'2 Maybe Parallel .001<p<.5','3 Parallel p>.5')))
cdf$p.cat    <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.p[i]<.01,'1 Reject Parallel p<.01','2 NOT Reject Parallel p>=.01'))
cdf$n.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$n[i] <= 100, 'N<=200','N>200'))
cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<20,'1 Short Series T<20','2 Long Series T>=20'))
cdf.sub <- cdf[cdf$n>50, ]
ggplot(data=cdf.sub, aes(x=npds/n, y=parallel.W)) + 
  geom_point(aes(size=n, colour=perf.adv.model), alpha=.7) + 
  scale_fill_brewer() + 
  facet_grid( p.cat ~ npds.cat  ) +
  scale_size(range = c(3,8)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() #+ labs(colour='Performance\nAdvantage', size='Observations (N*T)')
ggsave(filename='did_vs_bsts_compare_pd_aggregation_wald-stat_npds-n_NOLINE_FACET-Pcat-Tcat.png')

# cdf$p.cat    <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.p[i]<.001,'1 Not Parallel p<0.001',ifelse(cdf$parallel.p[i]<0.5,'2 Maybe Parallel .001<p<.5','3 Parallel p>.5')))
cdf$w.cat    <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i]>100,'1 Reject Parallel W>100','2 NOT Reject Parallel W<=100'))
cdf$n.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$n[i] <= 100, 'N<=200','N>200'))
cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<20,'1 Short Series T<20','2 Long Series T>=20'))
cdf.sub <- cdf[cdf$n>50, ]
ggplot(data=cdf.sub, aes(x=npds/n, y=parallel.W)) + 
  geom_point(aes(size=n, colour=perf.adv.model), alpha=.7) + 
  scale_fill_brewer() + 
  facet_grid( w.cat ~ npds.cat  ) +
  scale_size(range = c(3,8)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() #+ labs(colour='Performance\nAdvantage', size='Observations (N*T)')
ggsave(filename='did_vs_bsts_compare_pd_aggregation_wald-stat_npds-n_NOLINE_FACET-Wcat-Tcat.png')

# cdf$p.cat    <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.p[i]<.01,'1 Not Parallel p<0.01','2 Possibly Parallel p>=0.01'))
# cdf$n.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$n[i] <= 100, 'N<=200','N>200'))
# cdf$W.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i] > 100, '2 Not Parallel W>100','1 Possibly Parallel W<100'))
# cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<20,'1 Short Series T<20','2 Long Series T>=20'))
cdf.sub <- cdf[cdf$n>50, ]
ggplot(data=cdf.sub, aes(x=npds/n, y=parallel.W)) + 
  geom_point(aes(size=n, colour=perf.adv.model), alpha=.7) + 
  scale_fill_brewer() + 
  # facet_grid( W.cat ~ npds.cat  ) +
  scale_size(range = c(3,10)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() #+ labs(colour='Performance\nAdvantage', size='Observations (N*T)')
ggsave(filename='did_vs_bsts_compare_pd_aggregation_wald-stat_npds-n_NOLINE_NOFACET.png')












cdf$p.cat    <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.p[i]<.001,'p<0.001','p>=0.001'))
cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<20,'T<20','T>=20'))
# cdf$w.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i]<20,'3Low W',ifelse(cdf$parallel.W[i]<60,'2Mid W','1High W')))
# cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<npds.qs[1],'1Low T',ifelse(cdf$npds[i]<npds.qs[2],'2Mid T','3High T')))
# cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<npds.qs[1],'1Low T',ifelse(cdf$npds[i]<npds.qs[2],'2Mid T','3High T')))
# 
ggplot(data=cdf, aes(x=n, y=parallel.W, colour=perf.adv.model)) + 
  geom_point(aes(size=n)) + 
  scale_fill_brewer() + 
  facet_grid( p.cat ~ npds.cat  ) +
  scale_size(range = c(3,10)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  # geom_smooth(method = "lm", alpha=.1) +
  # ylab('Parallel Trends Pretest - Wald Chisq') + xlab('N (Num. Actor)') +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() #+ labs(colour='Performance\nAdvantage', size='Observations (N*T)')





n.qs <- quantile(cdf$n, probs = c(.5, 1))
npds.qs <- quantile(cdf$npds, probs = c(1/3, 2/3, 1))
# cdf$W.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i]<w.qs[1],'1Low',ifelse(cdf$parallel.W[i]<w.qs[2],'2Mid','3High')))
# cdf$W.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$parallel.W[i]<30,'W<30','W>=30'))
cdf$n.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$n[i]<n.qs[1],'2Low N','1High N'))
# cdf$n.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$n[i]<n.qs[1],'3Low N',ifelse(cdf$n[i]<n.qs[2],'2Mid N','1High N')))
cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<npds.qs[1],'1Low T',ifelse(cdf$npds[i]<npds.qs[2],'2Mid T','3High T')))
# cdf$npds.cat <- sapply(1:nrow(cdf), function(i)ifelse(cdf$npds[i]<npds.qs[1],'1Low T','2High T', T')))
ggplot(data=cdf, aes(x=parallel.p, y=parallel.W, colour=perf.adv.model)) + 
  geom_point(aes(size=npds*n)) + 
  scale_fill_brewer() + 
  facet_grid( n.cat ~ npds.cat ) +
  scale_size(range = c(3,14)) +
  scale_color_manual(values=c('#00B050','#0070C0')) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  # geom_smooth(method = "lm", alpha=.1) +
  # ylab('Parallel Trends Pretest - Wald Chisq') + xlab('N (Num. Actor)') +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() #+ labs(colour='Performance\nAdvantage', size='Observations (N*T)')



ggplot(data=cdf, aes(x=npds, y=n, fill=parallel.W)) + 
  geom_tile() + 
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  # geom_smooth(method = "lm") +
  # scale_fill_viridis_c() +
  scale_y_log10() + scale_x_log10() +
  # theme(label='Performance\nAdvantage') +
  theme_bw()


ggplot(data=cdf, aes(x=parallel.p, y=parallel.W)) + 
  geom_point() + 
  # geom_text(data = cdf, position = ) +
  # scale_colour_gradient2(low='red',mid='white',high='blue',midpoint=1) +
  # scale_color_brewer(palette = "RdYlBu") +
  # geom_smooth(method = "lm") +
  # scale_fill_viridis_c() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()




# cdf$go.for.bsts <- 0
# cdf$go.for.bsts[cdf$parallel.W > 10] <- 1











########################################
##
##
##  3.2.4. SENSITIVITY ANALYSIS
##         
##
##  PRIOR DISTRIBUTIONS COMPARISON
##
##
##
########################################

## Static Defaults
# n <- 100  ## NUmber of actors (i.e., number of timeseries)
# npds <- 60 #100
# intpd <- round( npds * 2/3 )
noise.level <- 1.5
# treat.rule <- 'random' # 'below.benchmark' ## 'random'
# treat.prob <-  0.5
# treat.threshold <- NA
b4 <- 1
b5 <- 0.04
dgp.nseasons= 12
dgp.freq= 1
## Scenarios
# lags = list(c(1),c(2),c(3)) 
# lags <- list( c(1) ) ##list(NULL)
## 
ns <- list(200) ## 400
sim.lengths <- list(240)
treat.rules <- list('random', 'below.benchmark')  ## 'below.benchmark'
seasonalities <- list(TRUE, FALSE)   ## c(TRUE,  FALSE )
prior.sd.scenarios <- list('sd.low') ## list('sd.low','sd.high')  ## sd.low
expect.mod.sizes <- list(1,3,5,7)
## FOCAL CONSTRUCT
dgp.ars <- list(0, 0.1, 0.2)  ## 0.6  ## .1,.2,.4
## STATE SPACE CONFIGURATIONS
st.sp.lists <- list(
  ## ## LEVEL
  # `1`=c('AddLocalLevel'),
  ## ## TREND
  # `2`=c('AddLocalLinearTrend'),
  # `3`=c('AddStudentLocalLinearTrend'),
  ## ## LEVEL + SLOPE ( + AR SLOPE DRIFT)
  # `4`=c('AddSemilocalLinearTrend'),
  ## ## SEASONAL
  # `5`=c('AddTrig')#,
  ## ## AUTOCORRELATION
  # list('AddAr'),
  `6`=c('AddAr'),
  ##------ COMBINATIONS --------------
  ## AR & SEASONALITY
  `7`=c('AddAr','AddTrig')#,
  # `7a`=c('AddAutoAr','AddTrig')#,
  ## LEVEL + ...
  # `8`=c('AddLocalLevel','AddTrig')#,
  # `9`=c('AddLocalLevel','AddAr'),
  # `10`=c('AddLocalLevel','AddTrig','AddAutoAr'),
  # ## (LEVEL + SLOPE) + ...
  # `11`=c('AddLocalLinearTrend','AddTrig')#,
  # `12`=c('AddLocalLinearTrend','AddAr'),
  # # `12`=c('AddLocalLinearTrend','AddTrig','AddAr'),
  # `13`=c('AddStudentLocalLinearTrend','AddTrig'),
  # `14`=c('AddStudentLocalLinearTrend','AddAr'),
  # # `14`c('AddStudentLocalLinearTrend','AddTrig','AddAr'),
  ## (LEVEL + SLOPE ( + AR1 SLOPE DRIFT)) + ...
  # `15`=c('AddTrig','AddSemilocalLinearTrend')#,
)



##
simlist <- list()
## NUMBER OF ACTORS (number of timeseries)
for (d in 1:length(ns)) {
  n <- ns[[ d ]]
  ## SIMULATION LENGTHS - NUMBER OF PERIODS
  for (f in 1:length(sim.lengths)) {
    npds <- sim.lengths[[ f ]]
    intpd <-  round( npds * 2/3 )
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
                intpd = intpd, ## 60% pre-intervention training / 40% post-intervention
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
effect.types = c('constant','geometric','quadratic') #c('quadratic')
##
bsts.niter <- 1e4
##
sim.id <- round(10*as.numeric(Sys.time()))


## RUN SIMULATION -  SIMULATE TIME SERIES
simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               sim.id = sim.id,
                               plot.show = FALSE, plot.save = FALSE )

## RUN BSTS and compare to DID
simlist.files <- runSimCompareBstsDiD(simlist, 
                                      effect.types = effect.types,
                                      sim.id = sim.id,
                                      save.items.dir= dir_ext,
                                      bsts.niter=bsts.niter * 8
)  ## D:\\BSTS_external


































































































simdf <- simx$sim$df

head(simdf)
actors <- sort(unique(simdf$actor[!is.na(simdf$match_id)]))
# pd.aggs <- c(4, 12)
# for (i in 1:length(actors)) {
#   actor <- actors[i]
#   for (j in 1:length(pd.aggs)) {
#     pd.agg <- pd.aggs[j]
#     
#   }
# }

# ## FUNCTION FOR AGGREGATING DATA TO DIFFERENT PERIOD
# tmpdf <- simdf
# pd.agg <- 4
# 
# ###
# ## Aggregate Aggregated Timeseries Simulation Panel Dataframe
# ###
# getAggregatedSimPanelDf <- function(tmpdf, pd.agg, 
#                                     na.rm=TRUE) {
#   if (is.na(pd.agg) | pd.agg <= 1) {
#     cat(sprintf('\npd.agg=%s is either NA or <= 1. Try setting pd.agg >= 2.\n',pd.agg))
#   }
#   ##
#   ts <- sort(unique(tmpdf$t))
#   npds.orig <- length(ts)
#   npds.new <- npds.orig / pd.agg
#   actors <- sort(unique(simdf$actor[!is.na(simdf$match_id)]))
#   
#   aggmap <- data.frame(t.old=1:npds.orig,
#                        t.new=rep(1:npds.new, each=pd.agg))
#   
#   intpd.old <- unique(simdf$t[which(simdf$t.post.intpd==1)])[1]
#   intpd.new <- aggmap$t.new[which(aggmap$t.old == intpd.old)]
#   
#   
#   tmpdf$t.agg <- NA
#   for (i in 1:nrow(aggmap)) {
#     idx <- which( tmpdf$t==aggmap$t.old[i] )
#     tmpdf$t.agg[idx] <- aggmap$t.new[i]
#   }
#   
#   if (na.rm) {
#     tmpdf <- tmpdf[ !is.na(tmpdf$match_id), ]
#   }
#   
#   aggdf <- tmpdf %>% 
#     group_by(effect.type, t.agg, actor) %>% 
#     summarize(
#       n_actors = n(),
#       group = paste(unique(group), collapse = '|'), 
#       match_id = paste(unique(match_id), collapse = '|'), 
#       match_pd = paste(unique(match_pd), collapse = '|'), 
#       ##
#       y = mean(y, na.rm=T),
#       ##
#       x1 = mean(x1, na.rm=T),
#       x2 = mean(x2, na.rm=T),
#       x3 = mean(x3, na.rm=T),
#       #
#       b1 = mean(b1, na.rm=T),
#       b2 = mean(b2, na.rm=T),
#       b3 = mean(b3, na.rm=T),
#       ##
#       c1 = mean(c1, na.rm=T),
#       c2 = mean(c2, na.rm=T),
#       c3 = mean(c3, na.rm=T),
#       ##
#       season.val = mean(season.val, na.rm=T),
#       u = mean(u, na.rm=T),
#       v = mean(v, na.rm=T)
#     )
#   
#   t.aggs <- sort(unique(aggdf$t.agg))
#   aggdf$t.agg.post.intpd <- NA
#   for (i in 1:length(t.aggs)) {
#     t.new <- t.aggs[ i ]
#     aggdf$t.agg.post.intpd[which(aggdf$t.agg==t.new)] <- ( t.new - (intpd.new - 1) )
#   }
#   
#   return(aggdf)
# }
# 
# ###
# ## Update Simlist configurations and simulated panel dataframes for aggregated periods
# ###
# updateSimlistAggregateSimDfPd <- function(simlist, pd.agg, na.rm=TRUE) {
#   
#   for (i in 1:length(simlist)) 
#   {
#     simx <- simlist[[ i ]]
#     
#     if (is.null(simx$sim)) {
#       cat(sprintf('\nsimlist item i=%s is missing "sim" object. First call runSimUpdateSimlist().\n',i))
#     }
#     
#     aggdf <- getAggregatedSimPanelDf(simx$sim$df, pd.agg=pd.agg, na.rm=na.rm)
#     simlist[[i]]$sim <- aggdf
#     simlist[[i]]$npds <- length(unique(aggdf$t.agg))
#     simlist[[i]]$intpd <- unique(aggdf$t.agg[which(aggdf$t.agg.post.intpd==1)])[1]
#   }
#     
#   return(simlist)
# }





out <- getAggregatedSimPanelDf(simdf, pd.agg = 4)
dim(out)
unique(out$t.agg)
out <- getAggregatedSimPanelDf(simdf, pd.agg = 12)
dim(out)
unique(out$t.agg)

simx$compare$did$quadratic$attgt$W[1]
simx$compare$did$quadratic$attgt$Wpval[1]




## CODA DIAGNOSTICS

geweke.diag(rowMeans(err5 ^ 2), .1, .5)
geweke.diag(coefmc5, .1, .5)

heidel.diag(rowMeans(err5 ^ 2))
heidel.diag(coefmc5)


pp <- postPredChecks(simx$compare$bsts$quadratic[[1]]$CausalImpact)
std.res <- pp$std.res

qqnorm(rowMeans(std.res), main = "Std.Residual QQ-plot, Y")
qqline(rowMeans(std.res))
Acf(rowMeans(std.res), main = "");title(main='Std.Residual ACF, Y')


causimp <- simx$compare$bsts$quadratic[[1]]$CausalImpact
p.posterior <- causimp$summary$p[1]

hist(causimp$model$bsts.model$sigma.obs)

causimp$model$bsts.model$one.step.prediction.errors

predictors <- causimp$model$bsts.model$predictors

response <- causimp$series$response
y <- response

niter <- nrow(causimp$model$bsts.model$coefficients)
burn <- round( niter * .2 )

## Newdata (post-intervention data) to predict via BSTS
newdata <- causimp$model$bsts.model$predictors[1:(intpd-1), ]
newdata <- cbind(response=response[1:(intpd-1)], newdata)

# predict.mbsts()
post.pred <- predict.bsts(causimp$model$bsts.model, newdata = newdata , burn = burn)
post.pred.dist <- post.pred$distribution
post.pred.mean <- post.pred$mean ##apply(post.pred$distribution, c(1), mean)

## before intervention period bool dummy
.ind <- (1:npds) < intpd

par(mfrow=c(2,2), mar=c(2,2,2.5,1))

## DENSITIES
plot(density(y[.ind]),  main = "Density comparison, Y")
lines(density(post.pred.mean), col='blue')

## HISTOGRAMS & BAYESIAN P-VALUES
# max.distrib <- apply(post.pred, c(2, 3), max)
max.distrib <- apply(post.pred$distribution, 1, max)
pvalue <- sum(max.distrib >= max(y[.ind]))/length(max.distrib)
hist(max.distrib, 30, col = "lightblue", border = "grey", 
     main = paste0("Bayesian p-val (Max Y) = ", round(pvalue, 2)),
     xlab = "Max. in-sample forecasts")
abline(v = max(y[.ind]), col = "darkblue", lwd = 3)


# Residual plots
y.rep <- matrix(y[.ind], length(y[.ind]), (niter - burn),  byrow = FALSE)
# res <- y.rep - (post.pred.dist - CausalMBSTS$mcmc$eps.samples[, i, ])
res <-  y.rep - t(post.pred.dist)
std.res <- res / apply(res,1,sd)   ## [residual i] / [ stdev of residual i]
# std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
qqnorm(rowMeans(std.res), main = "Std.Residual QQ-plot, Y")
qqline(rowMeans(std.res))
Acf(rowMeans(std.res), main = "");title(main='Std.Residual ACF, Y')

# par(mfrow=c(2,2))
# for (i in 1:4) hist(pred$distribution[,i])


#########################################
## Posterior Predictive Checks
########################################


# 
# 
# plotChecks <- function(CausalMBSTS, int.date) {
#   dates <- CausalMBSTS$dates
#   ind <- dates < int.date
#   y <- CausalMBSTS$y
#   post.pred <- CausalMBSTS$predict$post.pred.0
#   post.pred.mean <- apply(post.pred, c(1, 2), mean)
#   
#   for (i in 1:dim(y)[2]) {
#     # Density of posterior mean vs density of the data before intervention
#     plot(density(y[ind, i]), xlab = "", ylab = "",
#          main = paste0("Density comparison, Y", i))
#     lines(density(post.pred.mean[, i]), col = "blue")
#     
#     # Histograms & Bayesian p-value
#     max.distrib <- apply(post.pred, c(2, 3), max)
#     pvalue <- sum(max.distrib[i, ] >= max(y[ind, i]))/ncol(max.distrib)
#     hist(max.distrib[i, ], 30, col = "lightblue", border = "grey", main = paste0("Bayesian p-value = ",
#                                                                                  round(pvalue, 2), ", Y", i), xlab = "Max. in-sample forecasts")
#     abline(v = max(y[ind, i]), col = "darkblue", lwd = 3)
#     
#     # Residual plots
#     y.rep <- matrix(y[ind, i], nrow(y[ind, ]), (CausalMBSTS$mcmc$niter - CausalMBSTS$mcmc$burn),
#                     byrow = FALSE)
#     res <- (y.rep - (post.pred[, i, ] - CausalMBSTS$mcmc$eps.samples[, i, ]))
#     std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
#     qqnorm(rowMeans(std.res), main = paste0("Residual QQ-plot, Y", i))
#     qqline(rowMeans(std.res))
#     Acf(rowMeans(std.res), main = paste0("Residual ACF, Y", i))
#   }
# }




## Skip if 'random' treat.rule and treat.threshold index ii > 0
##   (random treat.rule already in simlist; no need to replicate bc treat.threshold not used)
# if (treat.rule=='random' & ii>1) {
#   cat(sprintf(' skipping redundant config for treat.rule=random, ii=%s\n',ii))
#   next
# }
# bsts.state.config <- list(
# AddTrig=getStateSpaceConf('AddTrig',
#                           period=dgp.nseasons,
#                           frequencies=dgp.freq,
#                           sigma.prior = SdPrior(sigma.guess =.01, sample.size =.1, initial.value =.1, upper.limit = Inf),
#                           # initial.state.prior=NormalPrior(mu=.1, sigma=1, initial.value = .1),
#                           # initial.state.prior = MvnDiagonalPrior(mean.vector = c(.1, .1), 
#                           #                                        sd.vector = c(1, 1)),
#                           method='harmonic'
#                           ),
# AddLocalLevel=getStateSpaceConf('AddLocalLevel',
#                   sigma.prior = SdPrior(sigma.guess =.01, sample.size =.1, initial.value = .1, upper.limit = Inf),
#                   initial.state.prior = NormalPrior(mu=.01,sigma = .1,initial.value = .1))#,
# AddAr=getStateSpaceConf('AddAr',
#                  lags=lag,
#                  sigma.prior=sigma.prior,
#                  initial.state.prior=initial.state.prior)
# AddAutoAr=getStateSpaceConf('AddAutoAr',
#                             lags=lag,
#                             prior=SpikeSlabArPrior(
#                               lag,
#                               prior.inclusion.probabilities = GeometricSequence( lag, initial.value = .9, discount.factor = .3),
#                               prior.mean = rep(0, lag),
#                               prior.sd = GeometricSequence(lag, initial.value = .9, discount.factor = .3),
#                               sdy=.5,
#                               prior.df = 1,
#                               expected.r2 = .5,
#                               sigma.upper.limit = Inf,
#                               truncate = TRUE)
#                             )#,
# AddLocalLinearTrend=getStateSpaceConf('AddLocalLinearTrend',
#                                       level.sigma.prior= SdPrior(sigma.guess = .01, sample.size = .01, initial.value = .01, upper.limit = Inf),
#                                       slope.sigma.prior = SdPrior(sigma.guess = .01, sample.size = .01, initial.value = .01, upper.limit = Inf),
#                                       initial.level.prior = NormalPrior(mu = .01, sigma = .1, initial.value = .01),
#                                       initial.slope.prior = NormalPrior(mu = .01, sigma = .1, initial.value = .01)
#                                       )#,
# AddStudentLocalLinearTrend=getStateSpaceConf('AddStudentLocalLinearTrend', save.weights=TRUE,
#                               level.sigma.prior= SdPrior(sigma.guess =.01, sample.size = .01, initial.value = .01, upper.limit = Inf),
#                               # level.nu.prior= GammaPrior(a = 1.5, b = 4, prior.mean = 3/4, initial.value = 3/4), ### Inherits: DoubleModel(),
#                               level.nu.prior= LognormalPrior(mu = .01, sigma = .01, initial.value = NULL), ### Inherits: DoubleModel(),
#                               slope.sigma.prior = SdPrior(sigma.guess =.01, sample.size = .01, initial.value = .01, upper.limit = Inf),
#                               # slope.nu.prior =GammaPrior(a = 2.5, b = 2, prior.mean = 5/4, initial.value =5/4),
#                               slope.nu.prior= LognormalPrior(mu = .01, sigma = .01, initial.value = NULL), ### Inherits: DoubleModel(),
#                               initial.level.prior = NormalPrior(mu = .01, sigma = .01, initial.value = .01),
#                               initial.slope.prior = NormalPrior(mu = .01, sigma = .01, initial.value = .01)
#                             )#,
# AddSemilocalLinearTrend=getStateSpaceConf('AddSemilocalLinearTrend',
#                                           level.sigma.prior= SdPrior(sigma.guess=.1, sample.size = .01, initial.value = .01, upper.limit = Inf),
#                                           slope.mean.prior = NormalPrior(mu=.01, sigma = .1, initial.value = .01),
#                                           slope.ar1.prior = Ar1CoefficientPrior(mu = .01, sigma = .1, force.stationary = F, force.positive = F, initial.value = .6),
#                                           slope.sigma.prior = SdPrior(sigma.guess =.1, sample.size = .01, initial.value = .01, upper.limit = Inf),
#                                           # slope.nu.prior =GammaPrior(a = 2.5, b = 2, prior.mean = 5/4, initial.value =5/4),
#                                           # level.nu.prior= GammaPrior(a = 1.5, b = 4, prior.mean = 3/4, initial.value = 3/4), ### Inherits: DoubleModel(),
#                                           initial.level.prior = NormalPrior(mu = .01, sigma = .1, initial.value = .01),
#                                           initial.slope.prior = NormalPrior(mu = .01, sigma = .1, initial.value = .01)
#                                           )
# ) ## end list bsts.state.config



######################################################
######################################################
######################################################
######################################################
################# TODO: META SIM #####################
######################################################
######################################################
######################################################
######################################################

### META SIMULATION TO COMPUTE EXPECTATION OF BIAS
###  AS AVG OF SIMULATION RUN BIAS [ATT - DGP]
runMetaSim(simlist, metric='att.error', nruns=10)


######################################################
######################################################
######################################################
##################### END TODO #######################
######################################################
######################################################
######################################################





######################################################### DEBUG #####
simi <- readRDS(file.path(dir_ext, '__GRIDSEARCH_output__16689319176_1.rds'))

bstss <- lapply(simi$compare$bsts,function(z)z[[1]]$CausalImpact$model$bsts.model)
CompareBstsModels(bstss)




##########################################
##########################################
##
##  Summarize Simulation Scenarios simlist
##
##########################################
#########################################


# graphics.off()

effect.types <- c('constant','quadratic','geometric')
compdf <- data.frame(stringsAsFactors = F)
for (i in 1:length(simlist)) 
{
  simi <- simlist[[i]]
  intpd <- simi$intpd
  
  # TODO: ADD 
  # bsts.res.list <- simi$compare$bsts
  # for (j in 1:length(bsts.res.list)) {
  #   bsts.res.j <- bsts.res.list[[ j ]]
  
    for (k in 1:length(effect.types)) 
    {   
      effect.type <- effect.types[k]
      res.tbl <- simi$compare$res.tbl[[ effect.type ]]
      
      agg.es <- simi$compare$did[[ effect.type ]]$agg.es
      impact_amount <- simi$compare$bsts[[ effect.type ]]$CausalImpact
      
      ## Treatment Effect
      gatt.b3 <- mean( as.numeric( res.tbl$b3.att[intpd:nrow(res.tbl)] ) )
      gatt.did <- agg.es$overall.att
      gatt.bsts <- impact_amount$summary$AbsEffect[1]
      
      ## TODO: Change to loop j over BSTS configs list
      ## BSTS State Specifications
      # bsts.state.spec.j <- simi$compare$bsts[[ j ]][[ effect.type ]]$CausalImpact$model$bsts.model$state.specification  ## with all bsts state specs
      bsts.state.spec.j <- simi$compare$bsts[[ effect.type ]]$CausalImpact$model$bsts.model$state.specification
      bsts.state.comps <- sapply(bsts.state.spec.j,function(z)class(z)[1])
      
      ## Simulation scenario (row) dataframe
      df.ijk <- data.frame(
        n = simi$n,
        npds = simi$npds,
        intpd = simi$npds,
        effect.type = effect.type,
        noise.level = simi$noise.level,
        treat.rule = simi$treat.rule,
        treat.prob = simi$treat.prob,
        treat.threshold = simi$treat.threshold,
        b4 = simi$b4,
        b5 = simi$b5,
        b9 = simi$b9,
        dgp.nseasons = simi$dgp.nseasons,
        dgp.freq = simi$dgp.freq,
        ##
        bsts.state.comps = paste(bsts.state.comps, collapse = '|'),
        bsts.state.comps.len = length(bsts.state.comps),
        ## Performance Measures (OVERALL GENERAL ATT)
        gatt.b3 = gatt.b3,
        gatt.did = gatt.did,
        gatt.bsts = gatt.bsts,
        ##
        bsts.cumu.pred.err.mean =  mean( simi$compare$bsts[[effect.type]]$cumu.pred.error, na.rm = T),
        bsts.cumu.pred.err.sum =  sum( simi$compare$bsts[[effect.type]]$cumu.pred.error, na.rm = T),
        # bsts.cumu.pred.err.mean =  mean( simi$compare$bsts[[ j ]][[effect.type]]$cumu.pred.error, na.rm = T), ## with all bsts state specs
        # bsts.cumu.pred.err.sum =  sum( simi$compare$bsts[[ j ]][[effect.type]]$cumu.pred.error, na.rm = T), ## with all bsts state specs
        ###
        did.b3.diff = (gatt.did - gatt.b3),
        bsts.b3.diff = (gatt.bsts - gatt.b3),
        ##
        rand.seed = simi$rand.seed
      )
      
      # ## Add BSTS state components in loop
      # # ]$compare$bsts[[ h ]][[effect.type]]$CausalImpact
      # bsts.res.list <- simi$compare$bsts[[ h ]][[effect.type]]$
      # for (l in 1:length(bests.res.list)) {
      #   cih <- bsts.res.list[[ h ]][[effect.type]]$CausalImpact
      #   gatt.bsts.h <- impact_amount$summary$AbsEffect[1]
      #   df.ijk <- cbind()
      # }
      
      ## Append row(s)
      compdf <- rbind(compdf, df.ijk)
      
    }
    
  # }
  
}

compdf$dgp.freq[is.na(compdf$dgp.freq)] <- 0
compdf$dgp.nseasons[is.na(compdf$dgp.nseasons)] <- 0

View(compdf)


##===========================
##  PLOTTING HEATMAPS
##--------------------------
# # Library
# library(ggplot2)
# # Dummy data
# x <- LETTERS[1:20]
# y <- paste0("var", seq(1,20))
# data <- expand.grid(X=x, Y=y)
# data$Z <- runif(400, 0, 5)
library(stringr)
heat.x.cols <- c('b9','b5','effect.type','bsts.state.comps')
heat.y.cols <- c('dgp.nseasons','treat.rule','treat.threshold','noise.level')
compdf$heat.x <- apply(compdf[,heat.x.cols],1,function(x)paste(x,collapse = '|'))
compdf$heat.y <- apply(compdf[,heat.y.cols],1,function(x)paste(x,collapse = '|')) 

.base <- compdf[, ! names(compdf) %in% c('did.b3.diff','bsts.b3.diff')]
.a <- .base
.a$stats.type <- 'bsts'
.a$b3.diff <- compdf$bsts.b3.diff
.b <- .base
.b$stats.type <- 'did'
.b$b3.diff <- compdf$did.b3.diff
compdf.stack <- rbind(.a,.b)
col.scale.vals <- c(
  min(compdf.stack$b3.diff[compdf.stack$b3.diff < 0], na.rm=T),
  median(compdf.stack$b3.diff[compdf.stack$b3.diff < 0], na.rm=T),
  0,
  median(compdf.stack$b3.diff[compdf.stack$b3.diff > 0], na.rm=T),
  max(compdf.stack$b3.diff[compdf.stack$b3.diff > 0], na.rm=T)
)
## PLOT FACET HEATMAP
ggplot(compdf.stack, aes(factor(heat.x), factor(heat.y), fill= b3.diff)) +
  geom_tile() + facet_wrap( . ~ factor(stats.type)) +
  xlab(paste(heat.x.cols,collapse = ' | ')) + 
  ylab(paste(heat.y.cols,collapse=' | ')) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  scale_fill_gradient2() #+
  # scale_fill_gradientn(colours=c('red','yellow','white','cyan','blue'), 
  #                      # values=rescale(c(-1,0-.Machine$double.eps,0,0+.Machine$double.eps,1)),
  #                      values=rescale(col.scale.vals)
  #                      )
  # scale_fill_gradient(low="white", high="blue") 

#
par(mfrow=c(1,2))
# BSTS ACCURACY Heatmap 
ggplot(compdf, aes(factor(heat.x), factor(heat.y), fill= bsts.b3.diff)) +
  geom_tile() + 
  xlab(paste(heat.x.cols,collapse = ' | ')) + 
  ylab(paste(heat.y.cols,collapse=' | ')) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
## DID ACCURACY HEATMAP
ggplot(compdf, aes(factor(heat.x), factor(heat.y), fill= did.b3.diff)) +
  geom_tile() + 
  xlab(paste(heat.x.cols,collapse = ' | ')) + 
  ylab(paste(heat.y.cols,collapse=' | ')) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
# Heatmap 
# ggplot(compdf, aes(factor(effect.type), factor(treat.rule), fill= bsts.b3.diff)) + 
#   geom_tile() + facet_wrap(factor(dgp.neasons) ~ factor(b5))







##=====================================
##  Structural Changes
##-------------------------------------
# ## devtools::install_github("KevinKotze/tsm")
# ## devtools::install_github("cran/fArma")
# library(tsm)
# library(fArma)
library(strucchange)
i <- 3
simi <- simlist[[ i ]]
simidfs <- simi$sim$df.summary
par(mfrow=c(3,2))
for (.effect.type in effect.types) {
  for (.group in c('treatment','control')) {
    y <- simidfs$med[simidfs$effect.type==.effect.type & simidfs$group==.group]
    dat <- data.frame(cbind(y[-1], y[-(length(y))]))
    colnames(dat) <- c("ylag0", "ylag1")
    fs <- Fstats(ylag0 ~ ylag1, data = dat)
    print( breakpoints(fs) ) # where the breakpoint is
    print( sctest(fs, type = "supF") )  # the significance of this breakpoint
    plot(fs, main=sprintf('%s %s',.group,.effect.type), ylim=c(0,100)); abline(v=(simi$intpd - 1)/simi$npds)
  }
}






# ## Results comparison table
# res.tbl <- cbind(bsts.res[ ,c('point.effect','point.effect.lower','point.effect.upper')],
#                  did.res[ ,c('term','event.time','estimate','point.conf.low','point.conf.high')],
#                  b3.treat=b3diff$treat,
#                  b3.ctrl=b3diff$ctrl,
#                  b3.att=b3diff$diff 
# )
# ## ROUND RESULTS TABLE (round numeric columns)
# # res.tbl <- res.tbl[ ! is.na(res.tbl$estimate), ]
# num.cols <- names(res.tbl)[ ! names(res.tbl) %in% c('term','event.time') ]
# res.tbl[ , num.cols] <- round( as.numeric(res.tbl[ , num.cols]), 4)
# # for (i in 1:length(num.cols)) {
# #   res.tbl[ , num.cols[i] ] <- round( as.numeric( res.tbl[ , num.cols[i] ] ), 4)
# # }
# ## MOVE ID COLUMNS TO FRONT
# .col.idx <- which(names(res.tbl) %in% c('term','event.time'))
# res.tbl4 <- cbind(res.tbl[, .col.idx],  res.tbl[, -.col.idx] )
# # View(res.tbl4)

##PLOT INCLUSION PROBABILITIES
# png(filename = sprintf('%s_BSTS_inclusion_probs_n%s_pd%s_ss%s_%s_%s_%s.png',
#                        prefix,n,npds,h,key.strip,effect.type,sim.id))
# plot(impact_amount$model$bsts.model,'coefficients', main=sprintf('%s %s', key,effect.type))
# dev.off()

## PLOT DYNAMIC EFFECTS COMPARISON - DID vs. BSTS vs. DGP
png(filename = sprintf('%s_BSTS_dynamic_treatment_effect_comparison_n%s_pd%s_ss%s_%s_%s_%s.png',
                       prefix,n,npds,h,key.strip,effect.type,sim.id))
res.tbl.filename <- sprintf('%s: DGP = %.3f; DiD = %.3f;  BSTS = %.3f',
                            key.strip, mean(res.tbl$b3.att[intpd:nrow(b3.att)]), agg.es$overall.att, impact_amount$summary$AbsEffect[1])
# ################################
# ### PRE DEBUG ###
# ################################
# getSimlistDebug <- function(x){
#   return(list(
#     n = actor.sizes[x$g],
#     npds = 100,
#     noise.level = noise.levels[x$h],
#     treat.rule = treat.rules[x$i],
#     treat.threshold = treat.thresholds[x$ii],
#     treat.prob = ifelse(treat.rules[x$i]=='random', 0.5, 1),
#     intpd = intpds[x$j],
#     autocor = autocors[x$k],
#     seasonality = seasonalities[x$m],
#     linear.trend = linear.trends[x$n],
#     state.conf = state.configs[x$p]
#   ))
# }
# config1 <- list(g=1,h=1,i=1,ii=1,j=1,k=1,m=1,n=1,p=1)
# config2 <- list(g=2,h=2,i=2,ii=2,j=2,k=2,m=2,n=2,p=2)
# simlist.debug <- list(
#   getSimlistDebug(config1),
#   getSimlistDebug(config2)
# )
# simlist <- simlist.debug
# ##
# .debug.sim <- runSimBstsDiDComparison(simlist.debug, plot.show = T, plot.save = T)
# 
# 
# # simlist <- .debug.sim
# sim <- simlist[[1]]
# effect.types = effect.types
# n = sim$n ## NUMBER OF FIRMS
# npds = sim$npds## NUMBER OF PERIODS
# intpd = sim$intpd ## intervention after first section
# ystart = 0
# treat.rule = ifelse(is.null(sim$treat.rule), NA, sim$treat.rule)
# treat.prob =  ifelse(is.null(sim$treat.prob), NA, sim$treat.prob) #0.95,  ## 0.6
# treat.threshold = ifelse(is.null(sim$treat.threshold), NA, sim$treat.threshold)  # 0.015
# sim.id = sim.id ## defaults to timestamp
# ##
# noise.level = ifelse(is.null(sim$noise.level), 0, sim$noise.level)
# ##
# b4 = ifelse(is.null(sim$b4), 0, sim$b4) ## past performance
# b5 = ifelse(is.null(sim$b5), 0, sim$b5) ## growth (linear function of time t)
# b9 = ifelse(is.null(sim$b9), 0, sim$b9) ## Autocorrelation
# ## # PEFORMANCE [Y] FUNCTION PARAMETERS
# b0 = .001 ## intercept
# b1 = .001 ## treatment dummy
# b2 = .001 ## post intervention dummy
# # b3 = .001, ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
# b4 = 0 ## spillover of past performance on current performance (how much of treatment effect persists across periods)
# b5 = .01 ## growth rate (linear effect of time: proportion of time t added to linear combination in yFunc() performance )
# ## Covariates
# b6 = 1 ## Age [1,2,3,...]
# b7 = 0 ## type [0,1]
# b8 = 1 ## level [0,1,2]
# b9 = 0
# ## Dynamic treatment effect function parameters
# w0 = 1.7  ## constant
# w1 = 0.18  ## linear
# w2 = -0.009 ## quadratic  ## -0.005, ## ***** made steeper curve= -0.008 *****
# ## TODO: CHECK w2.shift SENSITIVITY - this shifts quadratic curve several periods to the right so that treatment effect increases slowly
# w2.shift = -round(.07*sim$npds) ## optimal value here is likely a function of the combination of treatment effect function parameters
# ##
# nseasons  = ifelse(is.null(sim$dgp.nseasons), NA, sim$dgp.nseasons)
# season.frequency = ifelse(is.null(sim$freq), NA, sim$freq)
# ## Plotting
# plot.show = T ## TRUE
# plot.save = T  ## TRUE
# rand.seed = 54321












# ########################################
# ##            DEBUG
# ##
# #######################################
# # key <- '3a.rand.trendAr.ssTrendAr'
# key <- '1b.rand.trend.ssNone'
# xsim <- simlist[[key]]$sim
# xres <- simlist[[key]]$results
# 
# 
# xbstsmod <- xres$quadratic$bsts.CausalImpact$model$bsts.model
# plot(xres$quadratic$bsts.CausalImpact)
# plot(xbstsmod)







# 
# 
# ###########################################################
# ##   DEBUG
# ###########################################################
# 
# key <- '2b.rand.trendAr.ssAutoAr'
# effect.type <- 'quadratic'
# sim.id <- 'DEBUGSIMID'
# 
# simdf <- simlist[[key]]$sim$df
# simdf <- simdf[simdf$effect.type == effect.type, ]
# # simdf <- simdf[simdf$effect.type==effect.type, ]
# ## Set group name 'gname' field, where 0 = control, # = period of treatment
# simdf$gname <- 0
# simdf$gname[simdf$group=='treatment'] <- simdf$match_pd[simdf$group=='treatment']
# 
# ## Compute Multiperiod DiD for Avg Treatment Effect on Treated
# ccattgt <- att_gt(yname = "y", ## "Y",
#                   tname = "t",
#                   idname = "actor",
#                   gname = "gname",
#                   xformla = ~c1 + c2 + c3,
#                   data = simdf #,
#                   # panel = F
# )
# # ccattgt
# 
# ## Get first treatment group actor
# tr.actor.1 <- simdf$actor[which(simdf$group=='treatment')[1]]
# 
# ## SIMPLE AGGREGATION (OVERALL EFFECT) ATT
# agg.simple <- aggte(ccattgt, type='simple')
# # summary(agg.simple)
# 
# ## DYNAMIC EFFECTS AND EVENT STUDIES
# agg.es <- aggte(ccattgt, type = "dynamic")
# # summary(agg.es)
# # tidy(agg.es)
# ggdid(agg.es)
# # ggsave(filename = sprintf('%s_did_dynamic_effect_%s_%s_%s.png',
# #                           prefix,key,effect.type,sim.id))
# 
# 
# ##-----------------------------
# 
# ## Correlation of simulated to inferred
# cormat <- cor(cbind(simdf$b3[simdf$actor==tr.actor.1][-1], agg.es$att.egt))
# # ## MATPLOT of 
# # png(filename = sprintf('single_intervention_DiD_BSTS_DGP_comparison_%s.png',sim.id),
# #     width = 6,height = 6,units = 'in',res=300)
# # matplot(cbind(simdf$b3[simdf$actor==tr.actor.1][-1], agg.es$att.egt), type='o',pch=1:2)
# # dev.off()
# ## 
# sigmamat <- cor(cbind(simdf$x1[simdf$t==intpd], ## treatment dummy at intervention period
#                       simdf$y[simdf$t==(intpd-1)] ## performance variable at period before intervention
# ))
# 
# tmp.cordf <- data.frame(cor.type=c('dgp.did','x1t.ytm1'),
#                         effect.type=rep(effect.type,2),
#                         cor=c(cormat[2,1],sigmamat[2,1]))
# 
# ## SAVE TO OUTPUT LIST
# simlist[[key]]$cordf <- rbind(simlist[[key]]$cordf, tmp.cordf)  ## off-diagonal element of symmetric correlation matrix
# 
# # endog <- data.frame(threshold=c(1/2, 1/3,1/4,1/5,1/6),
# #                     cor=c(-0.2523,-.2030,-.2086,-.2107,-.2214))
# 
# # ## Density plot 
# # ggplot(simdf,aes(x=y, colour=effect.type)) + geom_density(alpha=0.1, size=1.4)
# 
# 
# ##====================
# ## BSTS Timseries Setup
# ##--------------------
# ## Aggregate into timeseries dataframe
# tsdf <- simdf %>%
#   filter( ! is.na(match_id)) %>%
#   group_by(t, t.post.intpd, effect.type, match_pd, gname, group, group.color) %>%
#   dplyr::summarize(
#     n_in_pd = n(),
#     actors = paste(unique(actor), collapse = '|'),
#     y_mean = mean(y, na.rm=T),
#     y_sum = sum(y, na.rm=T),
#     y_sd = sd(y, na.rm=T),
#     y_min = min(y, na.rm=T),
#     y_max = max(y, na.rm=T),
#     x1_mean = mean(x1, na.rm=T),
#     x2_mean = mean(x2, na.rm=T),
#     x3_mean = mean(x3, na.rm=T),
#     ##
#     c1_mean = mean(c1, na.rm=T),
#     c2_mean = mean(c2, na.rm=T),
#     c3_mean = mean(c3, na.rm=T),
#     #
#     b1_mean = mean(b1, na.rm=T),
#     b2_mean = mean(b2, na.rm=T),
#     b3_mean = mean(b3, na.rm=T),
#     #
#     u_mean = mean(u, na.rm=T),
#     v_mean = mean(v, na.rm=T)
#   )
# tsdf$.id <- 1:nrow(tsdf)
# 
# ## MAKE WIDE TIMESERIES FOR treatment,control groups in n periods
# val.cols <- c('y_mean','y_sum','y_min','y_max','y_sd',
#               'x1_mean','x2_mean','x3_mean',
#               'c1_mean','c2_mean','c3_mean',
#               'b1_mean','b2_mean','b3_mean',
#               'u_mean','v_mean')
# ts <- unique(tsdf$t)
# groups <- unique(tsdf$group)
# ## init timeseries dataframe - wide
# tsdfw <- data.frame(t=ts,  stringsAsFactors = F)
# for (j in 1:length(groups)) {
#   id.j <- which( tsdf$group == groups[j] ) 
#   for (k in 1:length(val.cols)) {
#     df.col <- data.frame( tsdf[ id.j , val.cols[k] ] )
#     names(df.col) <- sprintf('%s_%s',groups[j],val.cols[k])
#     tsdfw <- cbind(tsdfw,  df.col)
#   }
# }
# 
# 
# # Set up pre- and post-treatment period
# # pre.period <- as.Date(c("2013-01-01","2016-01-25"))
# pre.period <- c(1, intpd-1)  
# # post.period <- as.Date(c("2016-01-26","2018-01-01"))
# post.period <- c(intpd, npds) 
# 
# # # BSTS causal effect analysis using CausalImpact package
# # # CausalImpact option: 
# # # niter: Number of MCMC samples to draw. More samples lead to more accurate inferences. Defaults to 1000.
# # # standardize.data: Whether to standardize all columns of the data before fitting the model (i.e., setting Bayes priors), Defaults to TRUE.
# # # prior.level.sd: Prior standard deviation of the Gaussian random walk of the local level. Defaults to 0.01.
# # # nseasons: Period of the seasonal components. Default to 1.
# # # dynamic.regression: Whether to include time-varying regression coefficients. Defaults to FALSE.
# # impact_amount <- CausalImpact(amount.impact,pre.period,post.period,alpha=0.1, model.args = list(niter = 5000))
# # summary(impact_amount)
# # plot(impact_amount)
# dat <- tsdfw[,c('treatment_y_mean','control_y_mean','control_y_sd',#'control_y_sum', 'control_y_min','control_y_max',
#                 'control_c1_mean','control_c2_mean','control_c3_mean'#,
#                 # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
#                 # 'control_u_mean','control_v_mean'
# )]
# ## Train on y pre-treatment but NA's post-treatment
# y.pre.treat.NAs.post.treat <- c(dat$treatment_y_mean[1:(intpd-1)], rep(NA,npds-intpd+1))
# ## Then use the post-treatment response for causal impact estimation
# post.period.response <- dat$treatment_y_mean[intpd:npds]
# ## Covariates (predictors) - Matrix for "formula = y ~ predictors" argument
# predictors <- dat[, ! names(dat) %in% 'treatment_y_mean'] ## remove response; convert to matrix
# # ## Covariates (predictors) - Dataframe for "data" argument
# # predictors <- as.matrix(predictors)
# 
# ##----------------------------
# ## State Space Configuration
# ##----------------------------
# ## State Space Config list
# st.sp <- list()
# # st.sp <- AddAutoAr(st.sp, y.pre.treat.NAs.post.treat)
# # st.sp <- AddStudentLocalLinearTrend(st.sp, y.pre.treat.NAs.post.treat)
# bsts.state.specs <- list(AddAutoAr, AddStudentLocalLinearTrend)
# nss <- length(bsts.state.specs)
# if (nss > 0)
# {
#   for (j in 1:nss) {
#     addStateComponentFunc <- bsts.state.specs[[j]]
#     st.sp <- addStateComponentFunc(st.sp, y.pre.treat.NAs.post.treat)
#   }
# } else {
#   st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
# }
# print(st.sp)
# 
# ## BSTS model
# bsts.model <- bsts(y.pre.treat.NAs.post.treat ~ . ,
#                    state.specification = st.sp,
#                    data = predictors,
#                    niter = 5000)
# # ## BSTS model for Dynamic Regression
# # bsts.model <- bsts(y.pre.treat.NAs.post.treat,
# #                    state.specification = st.sp,
# #                    niter = 5000)
# 
# ## Use BSTS prediction of counterfactual to estimate CausalImpact
# impact_amount <- CausalImpact(bsts.model=bsts.model,
#                               post.period.response = post.period.response,
#                               alpha=0.05, model.args = list(niter = 5000))
# # ##



























































# ###=========================================
# ##  Aggregate actor series into 1 total
# ##------------------------------------------
# # actors <- sort(unique(simdf$actor))
# # for (i in 1:length(actors)) {
# #   actor <- actors[i]
# #   
# # }
# 
# ## Aggregate into timeseries dataframe
# tsdf <- simdf %>%
#   filter( ! is.na(match_id)) %>%
#   group_by(t, t.post.intpd, effect.type, match_pd, gname, group, group.color) %>%
#   dplyr::summarize(
#     n_in_pd = n(),
#     actors = paste(unique(actor), collapse = '|'),
#     y_mean = mean(y, na.rm=T),
#     y_sum = sum(y, na.rm=T),
#     y_sd = sd(y, na.rm=T),
#     y_min = min(y, na.rm=T),
#     y_max = max(y, na.rm=T),
#     x1_mean = mean(x1, na.rm=T),
#     x2_mean = mean(x2, na.rm=T),
#     x3_mean = mean(x3, na.rm=T),
#     ##
#     c1_mean = mean(c1, na.rm=T),
#     c2_mean = mean(c2, na.rm=T),
#     c3_mean = mean(c3, na.rm=T),
#     #
#     b1_mean = mean(b1, na.rm=T),
#     b2_mean = mean(b2, na.rm=T),
#     b3_mean = mean(b3, na.rm=T),
#     #
#     u_mean = mean(u, na.rm=T),
#     v_mean = mean(v, na.rm=T)
#   )
# 
# tsdf$.id <- 1:nrow(tsdf)
#   
# ##
# ##
# ##  ***** START HERE *****
# ##
# ##
# 
# ## Timeseries Dataframe - Wide
# tsdfw <- tsdf %>% pivot_wider(# id_cols = .(t,group),
#                      # id_expand = FALSE,
#                      names_from = c(group),
#                      # names_prefix = "",
#                      # names_sep = "_",
#                      # names_glue = NULL,
#                      # names_sort = FALSE,
#                      # names_vary = "fastest",
#                      # names_expand = FALSE,
#                      # names_repair = "check_unique",
#                      values_from = c('y_mean','y_max','y_min'),   #,
#                      # values_fill = NULL,
#                      values_fn = mean,
#                      # unused_fn = NULL
#                      )
# 
# # id.cols <- c('t','group')
# val.cols <- c('y_mean','y_sum','y_min','y_max','y_sd',
#               'x1_mean','x2_mean','x3_mean',
#               'c1_mean','c2_mean','c3_mean',
#               'b1_mean','b2_mean','b3_mean',
#               'u_mean','v_mean'
#               )
# ts <- unique(tsdf$t)
# groups <- unique(tsdf$group)
# tsdfw <- data.frame(t=ts,  stringsAsFactors = F)
# for (j in 1:length(groups)) {
#   id.j <- which( tsdf$group == groups[j] ) 
#   # sprintf('%s_%s_%s',ts[i],groups[j],)
#   # z <- tsdf[id.ij, ]
#   # tsdfw <- cbind(tsdfw, group=rep(groups[j],length(ts)))
#   for (k in 1:length(val.cols)) {
#     df.col <- data.frame( tsdf[ id.j , val.cols[k] ] )
#     names(df.col) <- sprintf('%s_%s',groups[j],val.cols[k])
#     tsdfw <- cbind(tsdfw,  df.col)
#   }
# }
# 
# 
# # Set up pre- and post-treatment period
# # pre.period <- as.Date(c("2013-01-01","2016-01-25"))
# pre.period <- c(1, intpd-1)  
# # post.period <- as.Date(c("2016-01-26","2018-01-01"))
# post.period <- c(intpd, npds) 
# 
# # # BSTS causal effect analysis using CausalImpact package
# # # CausalImpact option: 
# # # niter: Number of MCMC samples to draw. More samples lead to more accurate inferences. Defaults to 1000.
# # # standardize.data: Whether to standardize all columns of the data before fitting the model (i.e., setting Bayes priors), Defaults to TRUE.
# # # prior.level.sd: Prior standard deviation of the Gaussian random walk of the local level. Defaults to 0.01.
# # # nseasons: Period of the seasonal components. Default to 1.
# # # dynamic.regression: Whether to include time-varying regression coefficients. Defaults to FALSE.
# # impact_amount <- CausalImpact(amount.impact,pre.period,post.period,alpha=0.1, model.args = list(niter = 5000))
# # summary(impact_amount)
# # plot(impact_amount)
# dat <- tsdfw[,c('treatment_y_mean','control_y_mean','control_y_sd','control_y_sum', 'control_y_min','control_y_max',
#                 'control_c1_mean','control_c2_mean','control_c3_mean',
#                 # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
#                 'control_u_mean','control_v_mean'
#                 )]
# impact_amount <- CausalImpact(dat, pre.period,post.period,
#                               alpha=0.05, model.args = list(niter = 5000))
# summary(impact_amount)
# plot(impact_amount)
# bsts.res <- impact_amount$series
# 
# did.res <- tidy(agg.es)
# # plot(did.res)
# 
# 
# 
# 
# 
# # ## AVERAGE TREATMENT EFFECT USED IN SIMULATION
# # simdf
# 
# tr.actors <- unique(simdf$actor[which(simdf$group=='treatment')])  
# co.actors <- unique(simdf$actor[which(simdf$group=='control')])
# 
# b3diff <- data.frame(
#   treat=simdf %>% filter(group=='treatment' & actor==tr.actors[1]) %>% mutate(treat=b3) %>% dplyr::select(treat),
#   ctrl=simdf %>% filter(group=='control' & actor==co.actors[1]) %>% mutate(ctrl=b3) %>% dplyr::select(ctrl),
#   diff=NA
# )
# b3diff$diff <- b3diff$treat - b3diff$ctrl
# # simdf %>% 
# #   filter(group=='control' & actor==tr.actors[1]) %>% 
# #   dplyr::mutate(b3.diff=b3.treat-b3.ctrl) %>% 
# #   dplyr::select(b3.diff) 
# 
# 
# res.tbl <- cbind(bsts.res[ ,c('point.effect','point.effect.lower','point.effect.upper')],
#                  did.res[ ,c('term','event.time','estimate','point.conf.low','point.conf.high')],
#                  b3.treat=b3diff$treat,
#                  b3.ctrl=b3diff$ctrl,
#                  b3.att=b3diff$diff 
#                  )
# # res.tbl <- res.tbl[ ! is.na(res.tbl$estimate), ]
# num.cols <- names(res.tbl)[ ! names(res.tbl) %in% c('term','event.time') ]
# for (i in 1:length(num.cols)) {
#   res.tbl[ , num.cols[i] ] <- round( as.numeric( res.tbl[ , num.cols[i] ] ), 4)
# }
# res.tbl4 <- cbind(res.tbl[,4:5],  res.tbl[,-(4:5)] )
# View(res.tbl4)
# 
# plot(impact_amount$model$bsts.model,'coefficients')
# 
# matplot(res.tbl4[,c('point.effect','estimate','b3.att')],type='l',lty=c(1,2,4),lwd=c(1,1,2),col=c('black','red','blue'),
#         main=sprintf('ATT[DGP] = %.3f;  ATT[DiD] = %.3f;  ATT[BSTS] = %.3f',
#                      mean(b3diff$diff[intpd:nrow(b3diff)]), agg.es$overall.att, impact_amount$summary$AbsEffect[1]))
# legend('topright',legend=c('BSTS','DiD','DGP'),col=c('black','red','blue'),lty=c(1,2,4),lwd=c(1,1,2)) 



# group_by(order_pd_num, cohort_pd_num, group, age_cat, sex, married) %>%    ## order_pd_t0
# summarise(
#   n_in_pd = n(),
#   id_covar = paste(c(unique(age_cat),unique(sex),unique(married)),collapse = '_'),
#   order_cnt_percap_mean = sum(order_cnt,na.rm=T)/n(),
#   order_sum_percap_mean = sum(order_sum,na.rm=T)/n(),
#   ##
#   order_cnt_mean = mean(order_cnt, na.rm=T),
#   order_cnt_max = max(order_cnt, na.rm=T),
#   order_cnt_min = min(order_cnt, na.rm=T),
#   order_cnt_sd = sd(order_cnt, na.rm=T),
#   order_cnt_tot = sum(order_cnt, na.rm=T),
#   order_sum_mean = mean(order_sum, na.rm=T),
#   order_sum_max = max(order_sum, na.rm=T),
#   order_sum_min = min(order_sum, na.rm=T),
#   order_sum_sd = sd(order_sum, na.rm=T),
#   order_sum_tot = sum(order_sum, na.rm=T),
#   age_mean = mean(age, na.rm=T),
#   age_cat_mode = mode(age_cat),
#   married_y_prop = sum(married=='Y',na.rm=T)/n(),
#   sex_f_prop = sum(sex=='F',na.rm=T)/n()
# ) %>%
# filter(
#   order_pd_num >  quantile(ms$order_pd_num, .025),
#   order_pd_num <= quantile(ms$order_pd_num, .975)
# )

  
# 
# 
# 
# ##===========================================
# ## CHECK BSTS FOR HOW TO HANDLE 
# ##------------------------------------------
# 
# library(bsts)
# data(iclaims)
# 
# ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
# model <- bsts(initial.claims$iclaimsNSA,
#                state.specification = ss,
#                niter = 1000)
# 
# ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
# ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
# model1 <- bsts(initial.claims$iclaimsNSA,
#                state.specification = ss,
#                niter = 1000)
# plot(model)
# plot(model1, "components")  # plot(model1, "comp") works too!
# plot(model1, "help")
# 
# pred1 <- predict(model1, horizon = 12)
# plot(pred1, plot.original = 156)
# 
# # Fit a bsts model with expected model size 1, the default.
# model2 <- bsts(iclaimsNSA ~ .,
#                state.specification = ss,
#                niter = 1000,
#                data = initial.claims)
# 
# 
# # Fit a bsts model with expected model size 5, to include more coefficients.
# model3 <- bsts(iclaimsNSA ~ .,
#                state.specification = ss,
#                niter = 1000,
#                data = initial.claims,
#                expected.model.size = 5)  # Passed to SpikeSlabPrior.
# 
# plot(model2, 'comp')
# plot(model3, 'comp')
# 
# plot(model2, 'coef')
# plot(model3, 'coef')
# 
# 
# #####------------------------------
# library(tidyverse, quietly = TRUE)
# library(bsts, quietly = TRUE)    
# data(iclaims)
# .data <- initial.claims
# claims <- .data$iclaimsNSA
# plot(claims, ylab = "")
# 
# (model_components <- list())
# 
# summary(model_components <- AddLocalLinearTrend(model_components, 
#                                                 y = claims))
# summary(model_components <- AddSeasonal(model_components, y = claims, 
#                                         nseasons  = 52))
# fit <- bsts(claims, model_components, niter = 2000)
# 
# burnin <- 500 # Throw away first 500 
# tibble(
#   date = as.Date(time(claims)),
#   trend = colMeans(fit$state.contributions[-(1:burnin),"trend",]),
#   seasonality = colMeans(fit$state.contributions[-(1:burnin),"seasonal.52.1",])) %>%
#   gather("component", "value", trend, seasonality) %>%
#   ggplot(aes(x = date, y= value)) + 
#   geom_line() + theme_bw() + 
#   theme(legend.title = element_blank()) + ylab("") + xlab("") +
#   facet_grid(component ~ ., scales="free") + guides(colour=FALSE) +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0))
# 
# pred <- predict(fit, horizon = 100, burn = burnin, quantiles = c(.05, .95))
# plot(pred)
# 
# errors <- bsts.prediction.errors(fit, burn = 1000)
# PlotDynamicDistribution(errors$in.sample)
# 
# fit2 <- bsts(iclaimsNSA ~ ., state.specification = model_components, 
#              data = initial.claims, niter = 1000)
# 
# # ## 
# # agg.ca <- aggte(ccattgt, type = "calendar")
# # summary(agg.ca)
# # ggdid(agg.ca)
# # 
# # ## 
# # agg.gr <- aggte(ccattgt, type = "group")
# # summary(agg.gr)
# # ggdid(agg.gr)
# 
# ##---------------
# ## // end debug
# ##======================
# 
# 
# 
# 
# ##--------------------------------------
# ##  COMBINED PLOT FACET GRID
# ##   - Simulation Scenario by dynamic treatment effect shape
# ##--------------------------------------
# dfx.t0.summary <- data.frame(stringsAsFactors = F)
# dfx.att <- data.frame(stringsAsFactors = F)
# for (i in 1:length(simlist)) {
#   dfx.sim.i <- simlist[[i]]$sim$df.t0.summary
#   dfx.sim.i$scenario <- names(simlist)[i]
#   dfx.t0.summary <- rbind(dfx.t0.summary,  dfx.sim.i)
#   ##
#   dfx.att.i <- simlist[[i]]$sim$df.att
#   dfx.att.i$scenario <- names(simlist)[i]
#   dfx.att <- rbind(dfx.att,  dfx.att.i)
# }
# 
# pall <- ggplot(dfx.t0.summary, aes(x=t0, y=med, color=group)) +
#   geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
#   geom_line(size=1.2) +
#   # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
#   geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
#   ylab('Y') +
#   xlim(c(-intpd+2, npds-intpd-2)) +
#   facet_grid( scenario ~ effect.type ) +
#   theme_bw() + theme(legend.position='top') 
# print(pall)
# pall.file <- sprintf('single_intervention_staggered_DiD_COMBINED_%s.png',sim.id)   ### ******* ??? *************
# ggsave(filename=file.path(dir_plot, pall.file), plot=pall,
#        width=10,heigh=10,dpi=300,units='in')
# 
# ## Average Treatment Effect on the Treated (empirical difference between treated and controlled by day from simulation)
# patt <- ggplot(dfx.att, aes(x=t0, y=y.att)) +
#   # geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
#   geom_line(size=1.2) +
#   # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
#   geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
#   ylab('ATT') +
#   xlim(c(-intpd+2, npds-intpd-2)) +
#   facet_grid( scenario ~ effect.type ) +
#   theme_bw() + theme(legend.position='top') 
# print(patt)
# 
# ########################## END ##########################################
