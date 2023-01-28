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
noise.level <- 1.0   ### This 
# b4 <- 1.0  ## seasonal component weight
# b5 <- 1.0 ## Weight of localLevel component in outcome
dgp.nseasons= 52
dgp.freq= 1
post.intpd.portion <- 0.6  ## after 6 of ten years, the intervention happens

## Variables Ranges to grid search (optional)
ns <- list(200) ## 400
sim.lengths <- list(520)
treat.rules <- list('random')  ## 'below.benchmark'
seasonalities <- list(TRUE)   ## c(TRUE,  FALSE )
prior.sd.scenarios <- list('sd.low') ## list('sd.low','sd.high')  ## sd.low
# cov.scenarios <- list(low =list(c1=.01, c2=.05, c3=.1),
#                       high=list(c1=.1, c2=.2, c3=.3))   ### try decreasing variance of random walk 0.1 --> 0.01, while increasing their weights to 1.0
# cov.scenarios <- list(mixed=list(c1=.2, c2=.35, c3=.5))
# cov.scenarios <- list(high=list(c1=.5, c2=.6, c3=.7), 
#                       mid=list(c1=.2, c2=.3, c3=.4),
#                       low=list(c1=.05, c2=.1, c3=.15))
covariates.types <- list('control')  ## random, control
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
  # `8b`=c('AddLocalLevel','AddSeasonal', 'AddRegression'),
  # `8b`=c('AddSemilocalLinearTrend','AddSeasonal', 'AddRegression'),
  # `8b`=c('AddLocalLevel','AddSeasonal', 'AddRegression','AddAr')#,
  # `11b`=c('AddLocalLinearTrend','AddSeasonal', 'AddRegression'),
  `15b`=c('AddSemilocalLinearTrend','AddSeasonal', 'AddRegression')#,
  # `16`=c('AddAr','AddSeasonal', 'AddRegression')
)

## AVERGAING RESULTS ACROSS SIMULATION RUNS STARTING AT DIFFERENT RNG SEEDS
# sample(1:1e7, size = 30, replace = F)
rand.seeds <- list(
  5326233, 7989884, 5084888, 9598677, 1805267, 
  312661, 7783005, 2783942, 5969818, 3065768, 
  8923248, 3389679, 1027026, 203809, 485520, 
  971698, 689247, 6497856, 8462323, 5669916, 
  2663673, 79909, 6082183, 632897, 9153304, 
  86896, 134245, 2263407, 8802814, 465415
)
# rand.seeds <- list(
#   # 2601460,  583397, 9156861, 1413506, 8349641,
#   # 5928146, 2060318, 6458401, 2882649, 2367099,
#   # 9473581, 2961190, 3064348, 6267096, 6636439,
#   # 623676, 4718279, 9386470, 3895028, 4637202,
#   # 2115445, 3522413, 8239197,  6755386, 183790
# ) ## 9473581, 2961190, 3064348, 6267096, 6636439
# # rand.seeds <- list( 623676 )

##
simlist <- list()
## NUMBER OF ACTORS (number of timeseries)
for (d in 1:length(ns)) {
  n <- ns[[ d ]]
  ## SIMULATION LENGTHS - NUMBER OF PERIODS
  for (f in 1:length(sim.lengths)) {
    npds <- sim.lengths[[ f ]]
    intpd <-  round( npds * post.intpd.portion ) 
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
            
            ## COVARIATE SCENARIOS
            for (l in 1:length(covariates.types)) {
              cov.type <- covariates.types[[ l ]]
              
              for (r in 1:length(rand.seeds)) {
                run.rand.seed <- rand.seeds[[ r ]]
                
                
                # for (m in 1:length(bsts.n.cov.cats.list)) {
                #   bsts.n.cov.cats <- bsts.n.cov.cats.list[[ m ]]
                
                ##--------------------------------------------
                ## Append simulation configuration to simlist
                ##--------------------------------------------
                key <- sprintf('d%s|f%s|g%s|h%s|i%s|j%s|l%s|r%s', d,f,g,h,i,j,l,r)
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
                  # cov.scenario = cov.scenario, 
                  # cov.scenario.key = cov.scenario.key,
                  ## Dynamic treatment effect  (quadratic polynomial)
                  w0 = 1.0,  ## constant
                  w1 = 0.03, ## linear
                  w2 =  -0.00333 /sqrt(npds), ## ## quadratic
                  # w2.shift = -round( sqrt(npds)*.85 ), ## quadratic shift rightward (make all of U-shape after intervention)
                  ##
                  # b4 = b4,   ## seasonal component weight  (default 1)
                  # b5 = b5,   ## local level component weight (default 1)
                  b9 = dgp.ar  , ## autocorrelation
                  seasonality = seasonality,
                  dgp.nseasons= ifelse(seasonality, dgp.nseasons, NA), 
                  dgp.freq= ifelse(seasonality, dgp.freq, NA),
                  bsts.state.specs=bsts.state.specs,
                  covariates.type = cov.type,
                  # bsts.state.specs=list(list(AddSemilocalLinearTrend),list(AddSemilocalLinearTrend,AddStudentLocalLinearTrend)),
                  rand.seed = run.rand.seed
                )
                
                # }  // end m list n.cov.cats (?)
                
              }
              
            } ## // end l loop over cov.scenarios
            
          } ## // end prior.sd.scenarios loop
          
        } ## // end treat.rules loop
        
      } ## // end seasonality loop
      
    } ## // end ARs loop
    
  } ## // end nps loop  (sim lengths)
  
} ## // end ns loop (number of actors)


##---------------- RUN SIM GENERATE DATA SET -----------------------

##
effect.types = c('quadratic')  ## c('quadratic','geometric','constant')  ##cov.cols.need.fill.bool  ## constant','geometric
## ID for the simulation (to search/filter all simulation figures, RDS files, etc.)
sim.id <- round(10*as.numeric(Sys.time()))


## RUN SIMULATION -  SIMULATE TIME SERIES
simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               sim.id = sim.id,
                               plot.show = F, plot.save = F,
                               dgp.prior.sd.weight=.01 ) ## c('random','conditional','control')

##------------------- RUN BSTS; COMPARE vs. DID ---------------------
## `bsts.ctrl.cats` control group category definitions:
##  [NA]  = NO control group -- just covariate time series
##  [0,1] = 1 control group  with all untreated actors   +  covariate time series    
##  [1+ ] = Synthetic control series created by binning each covariate to make control group (N categories per covariates; 
##          num.groups = Cartesian product of all binned covariates, 
##          (e.g., bsts.ctrl.cats=3 for 3 covs (c1,c2,c3) --> 3*3*3=27 series to choose from for synthetic controls)
bsts.ctrl.cats.list <- list(NA,1)  #list(1, NA) ## NA=no control;
## BSTS expected model size (for spike-and-slab priors)
bsts.expect.mod.sizes <- list(3)  #list(1, 3, 5, 7) # 1  ## list(7, 4, 1)
## MCMC Iterations
bsts.niter.start <- 3000 ## 5000
bsts.niter.max   <- 3000 ## 8e4
## LOOP OVER BSTS MODEL COMPARISONS ON SAME SIMULATED DATA (SAME DGP SCENARIO)
for (m in 1:length(bsts.ctrl.cats.list)) {
  for (r in 1:length(bsts.expect.mod.sizes)) {
    ## RUN BSTS and compare to DID
    simlist.files <- runSimCompareBstsDiD(
      simlist, 
      effect.types = effect.types,
      sim.id = sim.id,
      save.items.dir= dir_ext,
      bsts.niter = bsts.niter.start,
      bsts.max.iter= bsts.niter.max,
      bsts.ctrl.cats = bsts.ctrl.cats.list[[ m ]],
      bsts.expect.mod.size = bsts.expect.mod.sizes[[ r ]]
    )  ## D:\\BSTS_external
    
    # break   ##**DEBUG**
  }
}

 




##----------------------- END ---------------------------------------





# ##================ DEBUG ====================================
# # x <- readRDS(file.path(dir_ext, '____DEBUG__bsts-illus_default-st-sp__GRIDSEARCH_output__n200_pd520_niter1000_covCats1_msize3_cov-control_16742929058_d1f1g1h1i1j1l1r1.rds'))
# x <- readRDS(file.path(dir_ext, '____DEBUG__bsts-illus_default-st-sp__GRIDSEARCH_output__n200_pd520_niter2000_covCatsNA_msize3_cov-control_16742962441_d1f1g1h1i1j1l1r4.rds'))
# simdf <- x$sim$df
# View(simdf)
# dfx <-  simdf[ (simdf$t > 2 & simdf$t < intpd ), names(simdf) %in% c('y','c1','c2','c3')]
# cor(dfx)
# ##===========================================================


sim.id <- 16743474366


runfiles <- dir(dir_ext, pattern=sprintf('%s.+\\.rds$', sim.id))
df <- data.frame(stringsAsFactors = F)
for (i in 1:length(runfiles)) {
  cat(sprintf('\n%s',runfiles[i]))
  l <- readRDS(file.path(dir_ext, runfiles[i]))
  simdf <- l$sim$df
  intpd <- l$intpd
  npds <- l$npds
  idx.pre <- 1:(intpd-1)
  idx.post <- intpd:npds
  df.pre <- simdf[idx.pre, ]
  df.post <- simdf[ -idx.pre, ]
  ##
  bsts.model <- l$compare$bsts$quadratic[[1]]$CausalImpact$model$bsts.model
  ##
  y <- bsts.model$original.series
  y.pre <- y[idx.pre]
  y.post <- y[ -idx.pre ]
  df.bsts <- cbind(y, bsts.model$predictors)
  df.bsts.pre <- df.bsts[idx.pre, ]
  df.bsts.post <- df.bsts[idx.post, ]
  
  niter <- bsts.model$niter
  burn <- round(niter * 0.2)
  idx.noburn <- (burn+1):niter
  err.pre <- bsts.model$one.step.prediction.errors[idx.noburn, idx.pre]
  
  res <- l$compare$res.tbl$quadratic[[1]][idx.post, ]
  bias <- res$bsts.point.effect - res$b3.att
  
  did.bias <- res$did.estimate - res$b3.att

  tmpdf <- data.frame(
    sim.id = l$sim$id,
    rand.seed = l$rand.seed,
    conv.check.all = l$compare$bsts$quadratic[[1]]$convcheck$converged.all,
    conv.check.prop = l$compare$bsts$quadratic[[1]]$convcheck$converged.prop,
    prior.sd.scenario = l$prior.sd.scenario,
    cov.has.control = ifelse('y_control' %in% dimnames(bsts.model$predictors)[[2]], 1, 0),
    covariates.types = l$covariates.type,
    effect.type = names(l$compare$bsts)[1],
    cor.simdf.y.c1= cor(df.pre[,c('y','c1')])[2,1],
    cor.simdf.y.c2= cor(df.pre[,c('y','c2')])[2,1],
    cor.simdf.y.c3= cor(df.pre[,c('y','c3')])[2,1],
    cor.bstsdf.y.c1= cor(df.bsts.pre[,c('y','c1_mean')])[2,1],
    cor.bstsdf.y.c2= cor(df.bsts.pre[,c('y','c2_mean')])[2,1],
    cor.bstsdf.y.c3= cor(df.bsts.pre[,c('y','c3_mean')])[2,1],
    err.pre.me= mean(err.pre, na.rm=T),
    err.pre.mae= mean( abs(err.pre), na.rm=T),
    err.pre.rmse= sqrt(mean(err.pre^2, na.rm=T)),
    err.pre.cumabs = sum(abs(err.pre), na.rm = T),
    bias.post.me = mean(bias, na.rm=T),
    bias.post.mae=mean( abs(bias), na.rm=T),
    bias.post.rmse=sqrt( mean(bias^2, na.rm=T)),
    bias.post.mape=mean(abs(bias)/res$b3.att, na.rm=T),
    ## ADD DID
    did.post.me = mean(did.bias, na.rm = T),
    did.post.mae = mean( abs(did.bias), na.rm = T),
    did.post.rmse = sqrt( mean(did.bias^2, na.rm = T))
  )
  df <- rbind(df, tmpdf)
}
dim(df)
View(df)
write.csv(df, file=file.path(dir_ext, sprintf('bsts_runs_%s.csv',sim.id)), row.names = F)

ggplot(df, aes(x=err.pre.me, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  xlim(c(-.1,.1)) + 
  geom_density(alpha=.15) +
  geom_vline(xintercept=0,lty=3) + theme_bw()
ggplot(df, aes(x=err.pre.me, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  xlim(c(-.1,.1)) + 
  geom_histogram(alpha=.15, binwidth = 0.01, position = 'identity') +
  geom_vline(xintercept=0,lty=3) + theme_bw()
###
ggplot(df, aes(x=bias.post.me, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  xlim(c(-1,1)) + 
  geom_density(alpha=.15) + 
  geom_vline(xintercept = 0, lty=3) + theme_bw()
ggplot(df, aes(x=bias.post.me, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  xlim(c(-1,1)) + 
  geom_histogram(alpha=.15, binwidth = 0.1, position = 'identity') + 
  geom_vline(xintercept = 0, lty=3) + theme_bw()



#####
# df.bsts.did <- df %>% group_by(rand.seed, cov.has.control)
cov.type.ctrl <- 'Has covariate with same DGP as Y (has a "control")'
cov.type.no <- 'NO covariate with same DGP as Y (No "control")'
df.bsts.did <- rbind(
  data.frame(model='bsts',
             err.post.me=df$bias.post.me[df$cov.has.control==1], 
             cov.has.ctrl=1, cov.type=cov.type.ctrl),
  data.frame(model='bsts',
             err.post.me=df$bias.post.me[df$cov.has.control==0], 
             cov.has.ctrl=0, cov.type=cov.type.no),
  data.frame(model='did',
             err.post.me=df$did.post.me[df$cov.has.control==1], 
             cov.has.ctrl=1, cov.type=cov.type.ctrl),
  data.frame(model='did',
             err.post.me=NA, 
             cov.has.ctrl=0, cov.type=cov.type.no)
)
df.bsts.did.s <- df.bsts.did %>% 
  group_by(cov.has.ctrl,model) %>%
  summarize(
    mean=mean(err.post.me, na.rm=T),
    sd=sd(err.post.me, na.rm=T)#,
    # skew=skewness(err.post.me, na.rm=T, type=2)
  )
# # A tibble: 4 x 3
# # Groups:   cov.has.ctrl [2]
#   cov.has.ctrl model     mean
#          <dbl> <chr>    <dbl>
# 1            0 bsts   -0.0337
# 2            0 did   NaN     
# 3            1 bsts    0.0183
# 4            1 did     0.0189
################
# df.bsts.did$covariate.type <- 'NO Control (No Covariate like Pre-Intervention Y)'
# df.bsts.did$covariate.type[df.bsts.did$cov.has.ctrl==1] <- 'HAS Control (Covariate like Pre-Intervention Y)'
# df.bsts.did$err.post.me[which(df.bsts.did$cov.has.ctrl==1 & df.bsts.did$model=='did')] <- NA
p.bsts.did.err <- ggplot(df.bsts.did, aes(x=err.post.me, fill=model, colour=model)) + 
  xlim(c(-1,1)) + 
  geom_density(alpha=.15) +
  geom_vline(xintercept=0,lty=3) + 
  facet_wrap(. ~ cov.type) + theme_bw()
ggsave(plot=p.bsts.did.err, 
       filename = file.path(dir_ext,'bsts_did_compare_post_err_with_control_facet_plot.png'), 
       units = 'in', height = 6, width=8, dpi=300)

######
ggplot(df, aes(x=bias.post.me, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  xlim(c(-1,1)) + 
  geom_histogram(alpha=.15, bins=11) + 
  geom_vline(xintercept = 0, lty=3) + theme_bw()

ggplot(df, aes(x=bias.post.me, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()

ggplot(df, aes(x=bias.post.rmse, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()


ggplot(df, aes(x=err.pre.cumabs, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()

ggplot(df, aes(x=err.pre.me, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()
ggplot(df, aes(x=err.pre.mae, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()
ggplot(df, aes(x=bias.post.rmse, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()


ggplot(df, aes(x=cor.bstsdf.y.c1, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()
ggplot(df, aes(x=cor.bstsdf.y.c2, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()
ggplot(df, aes(x=cor.bstsdf.y.c3, fill=factor(cov.has.control), colour=factor(cov.has.control))) + 
  geom_histogram(alpha=.15, bins=15) + geom_density(alpha=.15) + theme_bw()






cor.cols <- c(
'cor.bstsdf.y.c1',
'cor.bstsdf.y.c2',
'cor.bstsdf.y.c3',
'bias.post.me',
'cov.has.control'
)
round(cor(df[,cor.cols]), 3)











##_------------------------------------------------------------

# n  <- 100
# c1 <- rnorm(n, 1, .5)
# c2 <- rnorm(n, 1, .5)
# c3 <- rnorm(n, 1, .5)
# X  <- cbind(c1, c2, c3)
# sig.mat <- matrix(c( 1,.3,.1,
#                     .3, 1,.2,
#                     .1,.2, 1), ncol=3, byrow = T)
# mu.vec <- c(.4, .6, .5)
# rmv.mat <- mvtnorm::rmvnorm(3, 
#                             mean = mu.vec, 
#                             sigma = sig.mat) 
# # corplot
# 
# Xmvt <- cbind(x1,x2,x3) %*% sig

# 
# #######################################################################
# #######################################################################
# 
# 
# 
# 
# ##========================= DEBUG ===================================
# ## GET RESULT FROM STORAGE
# . <- readRDS(simlist.files[[1]]$file[1])
# bsts.model <- .$quadratic[[1]]$CausalImpact$model$bsts.model
# # bsts.model <- getBstsModelFromSimlist(list(.), key = 1, effect.type = 'quadratic')
# par(mar=c(2,2.5,2,1))
# plot(bsts.model, 'components')
# plot(bsts.model, 'predictors')
# plot(bsts.model, 'coefficients')
# 
# 
# l <- readRDS(file.path(dir_ext, '__GRIDSEARCH_output__n100_pd520_niter30000_16734236256_d1f1g1h1i1j1l1.rds'))
# # compare <- l$compare
# # l$compare <- NULL
# # l$compare <- compare
# # compare <- NULL
# 
# # l$bsts.state.specs
# 
# CompareBstsModels(list(getBstsModelFromSimlist(list(l), key=1, effect.type = 'constant'),
#                        getBstsModelFromSimlist(list(l), key=1, effect.type = 'geometric'),
#                        getBstsModelFromSimlist(list(l), key=1, effect.type = 'quadratic')
#                        ))
# 
# PlotBstsComponents(
#   getBstsModelFromSimlist(list(l), key=1, effect.type = 'constant'),
# )
# 
# 
# 
# ################
# a <- readRDS(file.path(dir_ext, '__GRIDSEARCH_output__n400_pd520_niter10000_covCats3_16735656074_d1f1g1h1i1j1l1.rds'))
# ################

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
















