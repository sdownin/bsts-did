


library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(CausalImpact)
library(bsts)
library(did)
library(ggpubr)
library(cowplot)
library(coda)
library(Boom)
library(BoomSpikeSlab)
library(e1071)
# library(sarima)
library(forecast)
# library(qualV)


## Directories
dir_proj <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'
dir_ext <- 'D:\\BSTS_external'
dir_plot <- file.path(dir_proj, 'plots')
dir_r <- file.path(dir_proj,'R')

##==============================
##  file prefix for saving images, writing outputs, etc.
##-----------------------------
prefix <- 'bsts-illus_vignette_'
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



## SIMULATION LIST
simlist <- list()

## State space configurations list to illustrate
st.sp.lists <- list(
  # `0_null` = c('AddStaticIntercept'),
  `1_level` = c('AddLocalLevel'),
  `2_season` = c('AddLocalLevel','AddSeasonal'),
  `3_regression`= c('AddLocalLevel','AddSeasonal','AddRegression'),
  `4a_trend` = c('AddSemilocalLinearTrend','AddSeasonal','AddRegression'),
  `4b_ar` = c('AddLocalLevel','AddAr','AddSeasonal','AddRegression'),
  `5a_trend_noreg` = c('AddLocalLinearTrend'),
  `5b_trend_season_noreg` = c('AddLocalLinearTrend','AddSeasonal')
)

## Static defaults
noise.level <- 1.3
b4 <- 1.0  ## seasonal component weight
b5 <- 0.04 ## yearly growth rate
dgp.nseasons= 52
dgp.freq= 1
## Variables Ranges to grid search (optional)
ns <- list(200) ## 200 ## list( 100, 200 ) ## ** INTERACT SIZE (N) with GRANULARITY
sim.lengths <- list(520)
treat.rules <- list('random')  ## 'below.benchmark'
seasonalities <- list(TRUE)   ## c(TRUE,  FALSE )
prior.sd.scenarios <- list('sd.low') ## list('sd.low','sd.high')  ## sd.low
expect.mod.sizes <- list(5)
## FOCAL CONSTRUCT
dgp.ars <- list(0)
n <- 200
npds <- 520
intpd <- round( 520 * (5/6) )


##-------------------
##       1 Null
##-------------------
# key <- '1_null'
key <- '5a_trend_noreg'


##--------------------------------------------
## Setup state space configurations
##--------------------------------------------
bsts.state.specs <- list()
## STATE SPACE COMPONENTS CONFIGURATION
st.sp.vec <- st.sp.lists[[ key ]]
bsts.state.config <- list()
for (kk in 1:length(st.sp.vec)) {
  .id <- length(bsts.state.config)+1
  bsts.state.config[[ .id ]] <- getStateSpaceConfBySimScenario(st.sp.vec[ kk ], 'sd.low')#, ## c('sd.high','sd.low')
}
bsts.state.specs[[ paste(st.sp.vec, collapse='|') ]] <- bsts.state.config


##--------------------------------------------
## Setup state space configurations
##--------------------------------------------
simlist[[ key ]] <- list(
  ##--------Simulation settings--------------
  n = n,    ## Number of firms
  npds = npds,  ## Number of periods
  intpd = intpd, ## #intervention period = (5/6)'ths of the total periods 
  noise.level = 1.3, ## stdev of simulated noise terms
  prior.sd.scenario = 'sd.low', ## BSTS Prior SD scenario (high vs. low uncertainty in priors
  treat.rule = 'random', 
  treat.prob =  0.5,  ## ifelse(treat.rule=='random', 0.5, 1), 
  treat.threshold = 1, ## ifelse(treat.rule=='random', 1, 0.5),
  seasonality = TRUE,
  dgp.nseasons= 52,  ## ifelse(seasonality, dgp.nseasons, NA), 
  dgp.freq=  1, ##ifelse(seasonality, dgp.freq, NA),
  rand.seed = 13579,
  ## Dynamic treatment effect  (quadratic polynomial)
  w0 = 1.5,  ## constant
  w1 = 0.13, ## linear
  w2 = -.05 / sqrt(npds), ## quadratic
  w2.shift = -round( sqrt(npds)*.7 ), ## quadratic shift rightward (make all of U-shape after intervention)
  ## focal parameters of outcome function (manipulated in sensitivity analysis)
  b4 = 1.0,   ## seasonal component weight
  b5 = 0.04,  ## yearly growth rate
  b9 = 0,     ## autocorrelation
  ##-----------BSTS state space----------------
  ##  BSTS state space parameters
  bsts.state.specs = bsts.state.specs, 
  expect.mod.size = 5  # number of covariates included in BSTS model of counterfactual
)



##
effect.types <- c('quadratic') ##  c('constant','geometric','quadratic') ## c('quadratic')  

##
sim.id <- round(10*as.numeric(Sys.time()))


## RUN SIMULATION -  SIMULATE TIME SERIES
simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               sim.id = sim.id,
                               plot.show = F, plot.save = F )





source(file.path(dir_r,'bsts_did_comparison_functions.R')) 
## Scale bsts iterations in increments (doubling) starting at 10k (more iterations for models that converge more slowly)
bsts.niter.start <- 2e2  ## 
bsts.niter.max   <- 2e3  ## 

# ## RUN BSTS and compare to DID
# simlist <- fitBstsUpdateSimlist(simlist,
#                                 effect.types = effect.types,
#                                 sim.id = sim.id,
#                                 save.items.dir= dir_ext,
#                                 bsts.niter = bsts.niter.start, ## **START at MAX niter for large npds **
#                                 bsts.max.iter= bsts.niter.max
# )  ## D:\\BSTS_external







## cache original bsts.niter for dynamic niter updating if MCMC convergence failed
bsts.niter <- 200
bsts.niter.orig <- bsts.niter

# print("runSimBstsDiDComparison()::SIMLIST INPUT:")
# print(simlist)
if (length(simlist) > 0 & length(names(simlist))==0) {
  names(simlist) <- 1:length(simlist)
}

## Simulation ID
if (is.na(sim.id)) {
  sim.id <- simlist[[1]]$sim$id
  if (is.null(sim.id) | is.na(sim.id)) {
    sim.id <- round(10*as.numeric(Sys.time()))
  } 
} 

## IF save simlist items is NA, then save images to work_dir
## else save images to save.items.fir
save.items.dir <- getwd()
save.img.dir <- ifelse(is.na(save.items.dir), getwd(), save.items.dir)

##===============================
##  BSTS State Specification Comparison 
##------------------------------
# for (i in 1:length(simlist))
# {
i <- 1


  key <- names(simlist)[i]
  key.strip <- gsub('[|]','',key,ignore.case = F, perl = T)
  cat(sprintf('\n%s, %s\n',i, key))
  
  simlist[[key]]$cordf <- data.frame()
  simlist[[key]]$compare <- list(bsts=list())
  
  ## simulation output from simulation scenario = simlist[[key]]
  npds <- simlist[[key]]$npds
  intpd <- simlist[[key]]$intpd
  n <- simlist[[key]]$n
  
  ## Simulation object (containing the simulated timeseries)
  sim <- simlist[[key]]$sim
  # sim.id <- simlist[[key]]$sim$id
  
  
  ## list of BSTS State component lists
  bsts.state.specs <- simlist[[key]]$bsts.state.specs
  if (length(names(bsts.state.specs))==0) {
    names(bsts.state.specs) <- 1:length(bsts.state.specs)
  }
  
  # ## BSTS expected model size for spike-and-slab priors
  expect.mod.size <- ifelse(is.null(simlist[[key]]$expect.mod.size), NA, simlist[[key]]$expect.mod.size)
  
  
  # ## Dynamic Treatment Effect Type shapes
  # for (k in 1:length(effect.types)) 
  # {
    
    k <- 1
    
    effect.type <- effect.types[k]
    simdf <- sim$df[sim$df$effect.type == effect.type, ]
    
    ## Set group name 'gname' field, where 0 = control, # = period of treatment
    simdf$match_pd <- as.numeric(simdf$match_pd)
    simdf$gname <- 0
    simdf$gname[simdf$group=='treatment'] <- simdf$match_pd[simdf$group=='treatment']
    ## Remove NAs
    simdf <- simdf[!is.na(simdf$match_id), ]
    
    
    ## Init output list wihtin simulations list
    simlist[[key]]$compare$bsts[[effect.type]] <- list()
    
    
    ##====================
    ## BSTS Timseries Setup
    ##--------------------
    ## Aggregate into timeseries dataframe
    tsdf <- simdf %>%
      dplyr::filter( ! is.na(match_id)) %>%
      group_by(t, t.post.intpd, effect.type, match_pd, gname, group, group.color) %>%
      dplyr::summarize(
        n_in_pd = n(),
        actors = paste(unique(actor), collapse = '|'),
        y_outcome = mean(y, na.rm=T),
        y_sum = sum(y, na.rm=T),
        y_sd = sd(y, na.rm=T),
        y_min = min(y, na.rm=T),
        y_max = max(y, na.rm=T),
        y_skew = skewness(y, na.rm = T, type = 2), ## moment-based distribution
        y_kurt = kurtosis(y, na.rm = T, type = 2), ## moment-based distribution
        ##
        x1_mean = mean(x1, na.rm=T),
        x2_mean = mean(x2, na.rm=T),
        x3_mean = mean(x3, na.rm=T),
        ##
        c1_mean = mean(c1, na.rm=T),
        c2_mean = mean(c2, na.rm=T),
        c3_mean = mean(c3, na.rm=T),
        #
        c1_sd = sd(c1, na.rm=T),
        c2_sd = sd(c2, na.rm=T),
        c3_sd = sd(c3, na.rm=T),
        ##
        c1_skew = skewness(c1, na.rm=T, type = 2),
        c2_skew = skewness(c2, na.rm=T, type = 2),
        c3_skew = skewness(c3, na.rm=T, type = 2),
        #
        c1_kurt = skewness(c1, na.rm=T, type = 2),
        c2_kurt = skewness(c2, na.rm=T, type = 2),
        c3_kurt = skewness(c3, na.rm=T, type = 2),
        #
        b1_mean = mean(b1, na.rm=T),
        b2_mean = mean(b2, na.rm=T),
        b3_mean = mean(b3, na.rm=T),
        #
        u_mean = mean(u, na.rm=T),
        v_mean = mean(v, na.rm=T)
      )
    tsdf$.id <- 1:nrow(tsdf)
    
    ## MAKE WIDE TIMESERIES FOR treatment,control groups in n periods
    val.cols <- c('y_outcome','y_sum','y_min','y_max','y_sd',
                  'y_skew','y_kurt',
                  'x1_mean','x2_mean','x3_mean',
                  'c1_mean','c2_mean','c3_mean',
                  'c1_sd','c2_sd','c3_sd',
                  'c1_skew','c2_skew','c3_skew',
                  'c1_kurt','c2_kurt','c3_kurt',
                  'b1_mean','b2_mean','b3_mean',
                  'u_mean','v_mean')
    ts <- unique(tsdf$t)
    groups <- unique(tsdf$group)
    ## init timeseries dataframe - wide
    tsdfw <- data.frame(t=ts,  stringsAsFactors = F)
    for (jj in 1:length(groups)) {
      id.j <- which( tsdf$group == groups[ jj ] ) 
      for (kk in 1:length(val.cols)) {
        df.col <- data.frame( tsdf[ id.j , val.cols[ kk ] ] )
        names(df.col) <- sprintf('%s_%s',groups[ jj ],val.cols[ kk ])
        tsdfw <- cbind(tsdfw,  df.col)
      }
    }
    
    
    # Set up pre- and post-treatment period
    # pre.period <- as.Date(c("2013-01-01","2016-01-25"))
    pre.period <- c(1, intpd-1)  
    # post.period <- as.Date(c("2016-01-26","2018-01-01"))
    post.period <- c(intpd, npds) 
    
    # # BSTS causal effect analysis using CausalImpact package
    # # CausalImpact option: 
    # # niter: Number of MCMC samples to draw. More samples lead to more accurate inferences. Defaults to 1000.
    # # standardize.data: Whether to standardize all columns of the data before fitting the model (i.e., setting Bayes priors), Defaults to TRUE.
    # # prior.level.sd: Prior standard deviation of the Gaussian random walk of the local level. Defaults to 0.01.
    # # nseasons: Period of the seasonal components. Default to 1.
    # # dynamic.regression: Whether to include time-varying regression coefficients. Defaults to FALSE.
    # impact_amount <- CausalImpact(amount.impact,pre.period,post.period,alpha=0.1, model.args = list(niter = 5000))
    # summary(impact_amount)
    # plot(impact_amount)
    dat <- tsdfw[,c('treatment_y_outcome','control_y_outcome','control_y_sd',#'control_y_sum', 'control_y_min','control_y_max',
                    'control_y_max','control_y_skew', 'control_y_kurt',
                    'control_c1_mean','control_c2_mean','control_c3_mean',
                    'control_c1_sd','control_c2_sd','control_c3_sd',
                    # 'control_c1_skew','control_c2_skew',
                    'control_c3_skew',
                    # 'control_c1_kurt','control_c2_kurt',
                    'control_c3_kurt'#,
                    # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
                    # 'control_u_mean','control_v_mean'
    )]
    ## Train on y pre-treatment but NA's post-treatment
    y.pre.treat.NAs.post.treat <- c(dat$treatment_y_outcome[1:(intpd-1)], rep(NA,npds-intpd+1))
    ## Then use the post-treatment response for causal impact estimation
    post.period.response <- dat$treatment_y_outcome[intpd:npds]
    ## Covariates (predictors) - Matrix for "formula = y ~ predictors" argument
    predictors <- dat[, ! names(dat) %in% 'treatment_y_outcome'] ## remove response; convert to matrix
    # ## Covariates (predictors) - Dataframe for "data" argument
    # predictors <- as.matrix(predictors) 
    
    ## ADD temporal trend to covariates
    predictors$covar_temp_trend <- (1:npds) + rnorm(npds, 0, noise.level)   ## * simlist$`0.rand.base`$b5
    
    ## ***DEBUG***
    print("bsts.state.specs:")
    print(bsts.state.specs)
    
    ##----------------------------
    ## State Space Configuration
    ##----------------------------
    ## LOOP OVER BSTS STATE SPECIFICATIONS FOR THE SAME SIMULATION RUN
    # for (h in 1:length(bsts.state.specs)) 
    # {
      h <- 1
      
      simlist[[key]]$compare$bsts[[effect.type]][[ h ]] <- list()
      
      ## h'th BSTS state space configuration (state component list)
      state.conf <- bsts.state.specs[[ h ]]
      
      ## names of 
      state.comps <- unname(sapply(state.conf, function(x){
        ifelse(class(x)=='list' & 'name' %in% names(x), x$name, 'ERROR x$name not set in state.conf[x]')
      }, simplify = T))
      
      ## Loop over state space components
      nss <- length(state.conf)
      ## State Space Config list
      st.sp <- list()
      if (nss > 0) {
        for (jj in 1:nss) {
          state.conf.item <- state.conf[[ jj ]]
          if (class(state.conf.item)=='list') {
            
            if(state.conf.item$name=='AddSharedLocalLevel') {
              state.conf$coefficient.prior <- SpikeSlabPrior(
                x = predictors,
                y = dat$treatment_y_mean, ##**NOTE** USING ALL Y VALS (not just y.pre.treat.NAs.post.treat)
                expected.r2 = .5, ## [.5]
                prior.df = .01, ##[.01]
                expected.model.size = 2, ## [1]
                prior.information.weight = .01, ## [.01]
                diagonal.shrinkage = 0.5,  ## [.5] setting=0 --> Zellner's G-prior
                # optional.coefficient.estimate = NULL,
                max.flips = -1,  ## <= 0 means all indicators will be sampled
                # prior.inclusion.probabilities = NULL,
                sigma.upper.limit = Inf
              )
              
            }
            ## [[ IF NOT AddSharedLocalLevel(), NOT INCLUDING SPIKESLABPRIOR ]]
            st.sp <- updateStateSpaceAddComponentFromConfList(st.sp,  
                                                              y.pre.treat.NAs.post.treat,  
                                                              state.conf.item)
            cat(sprintf('add to state.space: %s\n',state.conf.item$name))
          }
        }
      } else {
        ## Default in CausalImpact package
        st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
      }
      # print(st.sp)
      
      cat(sprintf('\nRunning BSTS model estimation for state.conf h=%s\n',h))
      
      ##-------------------------------
      ## Regression component of model
      if ('AddRegression' %in% state.comps) {
        bsts.input.form <- y.pre.treat.NAs.post.treat ~ . ## with regression
      } else {
        bsts.input.form <- y.pre.treat.NAs.post.treat    ## without regression
        # predictors <- data.frame(intercept=rep(1, npds)) ##predictors[1] ## keep 1st column as dataframe
      }
      
      ##---------------------------------------------------------
      ## RUN BSTS WITH DYNAMIC niter BASED ON CONVERGENCE 
      isConverged <- FALSE
      isMaxIter <- FALSE
      hasBstsError <- FALSE
      bsts.niter <- bsts.niter.orig  ## reset to original bsts.niter input value
      while ( !isConverged  &  !isMaxIter & !hasBstsError  ) {
        
        # ##**DEBUG**
        # print('st.sp')
        # print(st.sp)
        # ##
        
        ## BSTS model
        bsts.model <- tryCatch(expr = {
          bsts(formula = bsts.input.form,
               state.specification = st.sp,
               data = predictors,
               expected.model.size = expect.mod.size,
               niter = 1000)
        },
        error=function(e) {
          message(sprintf('bsts() error: %s', as.character(e)))
        },
        warning=function(w) {
          message(sprintf('bsts() warning: %s', as.character(w)))
        },
        finally={
          ##PASS
        })
        
        
        ## skip if bsts() function threw an error (not return 'bsts' object)
        if ( class(bsts.model) == 'bsts' ) {
          ## Use BSTS prediction of counterfactual to estimate CausalImpact
          impact_amount <- CausalImpact(bsts.model=bsts.model,
                                        post.period.response = post.period.response,
                                        alpha=0.05, model.args = list(niter = bsts.niter))
          causimp <- impact_amount
          ## POSTERIOR PREDICTIVE CHECKS
          ppcheck.filename <- file.path(save.img.dir,
                                        sprintf('%s_bsts_post_pred_checks_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
                                                prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id))
          convcheck <- postPredChecks(impact_amount, filename=ppcheck.filename, return.val = T)
          # convcheck <- postPredChecksPreIntervention(impact_amount, filename=ppcheck.filename, return.val = T)
          ##
          # print(convcheck)
          ##
          cat(sprintf('\nBSTS niter = %s',bsts.niter))
          cat(convcheck$summary)
          ## UPDATE CONVERGENCE CHECK FLAG - ELSE RERUN WITH INCREASED bsts.niter
          # isConverged <- convcheck$converged.all
          conv.tol <- 0.8
          conv.min.iter.thresh <- 4e4 ## 40k
          # isConverged <- convcheck$converged.prop >= conv.tol
          isConverged <- convcheck$converged.all | (convcheck$converged.prop >= conv.tol & bsts.niter >= conv.min.iter.thresh) ## don't allow incomplete check below minimum threshold of bsts.niter = 10k 
          print(convcheck$converged)
          cat(sprintf('Converged proportion = %.3f (tol = %.3f) (min.iter.converg.thresh=%s)\nConverged status = %s\n\n',
                      convcheck$converged.prop, conv.tol, conv.min.iter.thresh,  isConverged))
        } else {
          hasBstsError <- TRUE
        }
        
        
        ## UPDATE niter --> increase by factor of 2
        if ( !isConverged ) {
          # .niter <- round( bsts.niter * (2 - convcheck$converged.prop) ) ## scale niter increment by proportion of failed checks
          bsts.niter.new <- round( (bsts.niter * 2) / 10) * 10  ## double niter and make divisible by 10 (for bsts ping=10 to divide evenly into niter)
          if ( bsts.niter.new <= bsts.max.iter ) {
            bsts.niter <- bsts.niter.new  
          } else {
            isMaxIter <- TRUE
          }
        }
        
        
      }
      
      if (hasBstsError) {
        next
      }
      ##-------------------------------------------------------------
      
      
      ##
      plot(bsts.model, main=sprintf('BSTS Plot: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
      PlotBstsComponents(bsts.model) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
      PlotBstsSize(bsts.model, main=sprintf('BSTS Size: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
      ##
      
      #                         key,effect.type,sim.id))
      p.bsts.impact.all <- plot(impact_amount, c('original','pointwise','cumulative')) # pointwise','cumulative
      ggsave(filename = file.path(save.img.dir,
                                  sprintf('%s_bsts_CausalImpact_plot_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
                                          prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id))
      )
      # dev.off()
      
      ##-------------------
      ## OVERALL CONFIDENCE INTERVALS AND PVALUES
      ## OVERALL CAUSAL P-VALUE
      pval.bsts.general <- impact_amount$summary$p[2]
      # pval.did.general <- ccattgt
      ci.bsts.general <- c(impact_amount$summary$AbsEffect.lower[1], impact_amount$summary$AbsEffect.upper[1]) 
      ## MUST USE GROUP AGGREGATION TO GET critical.val.egt
      ## (we only have one group [i.e., one treatment time])
      # ci.did.general <- c(agg.group$overall.att - (agg.group$overall.se * agg.group$crit.val.egt),
      #                     agg.group$overall.att + (agg.group$overall.se * agg.group$crit.val.egt))
      ##-------------------
      
      # ## DID
      # did.res <- tidy(agg.es)
      ## BSTS
      bsts.res <- impact_amount$series
      # .nidx <- which(names(bsts.res) %in% c('point.effect','point.effect.lower','point.effect.upper'))
      bsts.res <- bsts.res[ , which(names(bsts.res) %in% c('point.effect','point.effect.lower','point.effect.upper')) ]
      names(bsts.res) <- c('bsts.point.effect','bsts.point.effect.lower','bsts.point.effect.upper')
      
      # plot(did.res)
      
      # ## AVERAGE TREATMENT EFFECT USED IN SIMULATION
      # simdf
      
      tr.actors <- unique(simdf$actor[which(simdf$group=='treatment')])  
      co.actors <- unique(simdf$actor[which(simdf$group=='control')])
      
      ##
      # att.b3 <- mean(b3diff$diff[intpd:npds])
      # att.did <- agg.es$overall.att
      att.bsts <- impact_amount$summary$AbsEffect[1]