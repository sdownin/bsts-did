#####################################################
##
##   Dynamic Causal Inference Simulations
##
##    Run script for simulations of
##     - internal intervention (self-selection) endogeneity
##     - ...
##
#####################################################

library(dplyr)
library(tidyr)
library(CausalImpact)
library(bsts)
library(tibble)
library(did)
library(sarima)

## Directories
dir_proj <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'   
dir_plot <- file.path(dir_proj, 'plots')
dir_r <- file.path(dir_proj,'R')

setwd(dir_proj)


##==============================
##  file prefix for saving images, writing outputs, etc.
##-----------------------------
prefix <- 'single_intervention-statespace-config_'

##==============================
## Load simulation functions
##------------------------------
# source(file.path(dir_r,'internal_intervention_sim.R'))
source(file.path(dir_r,'single_intervention_sim_vec.R'))  ## vectorized -- in development






####################################
##  MAIN COMPARISON FUNCTION
##  - 1. runSimSingleInterventionEffectComparison() on simlist
##  - 2. DiD & BSTS results
##  - 3. comparison of DiD & BSTS performance
##   returns full simlist
######################################
runSimBstsDiDComparison <- function(simlist,     ## n, npds, intpd moved into simlist elements
                                    effect.types=c('constant','quadratic','geometric'), 
                                    sim.id=round(10*as.numeric(Sys.time())),
                                    plot.show=F, plot.save=F)
{
  ##----------------------------
  ## Run Simulation List
  ##----------------------------
  # sim.id <- ifelse( is.na(sim.id), round(as.numeric(Sys.time())), sim.id)
  for (i in 1:length(simlist)) {
    key <- names(simlist)[i]
    sim <- simlist[[key]]
    cat(sprintf('\nScenario label: %s\n\n', key))
    ##  
    set.seed(sim$rand.seed)
    ##
    simlist[[key]]$sim <- runSimSingleInterventionEffectComparison(
      effect.types = effect.types,
      n = sim$n, ## NUMBER OF FIRMS
      npds = sim$npds, ## NUMBER OF PERIODS
      intpd = sim$intpd, ## intervention after first section
      ystart = 0,
      treat.rule = ifelse(is.null(sim$treat.rule), NA, sim$treat.rule),
      treat.prob =  ifelse(is.null(sim$treat.prob), NA, sim$treat.prob), #0.95,  ## 0.6
      treat.threshold = ifelse(is.null(sim$treat.threshold), NA, sim$treat.threshold),  # 0.015
      sim.id = sim.id, ## defaults to timestamp
      ##
      noise.level = ifelse(is.null(sim$noise.level), 0, sim$noise.level),
      ##
      b4 = ifelse(is.null(sim$b4), 0, sim$b4), ## past performance
      b5 = ifelse(is.null(sim$b5), 0, sim$b5), ## growth (linear function of time t)
      b9 = ifelse(is.null(sim$b9), 0, sim$b9), ## Autocorrelation
      ## Dynamic treatment effect function parameters
      w0 = 1.7,  ## constant
      w1 = 0.18,  ## linear
      w2 = -0.009,  ## quadratic  ## -0.005, ## ***** made steeper curve= -0.008 *****
      ## TODO: CHECK w2.shift SENSITIVITY - this shifts quadratic curve several periods to the right so that treatment effect increases slowly  
      w2.shift = -round(.07*sim$npds),  ## optimal value here is likely a function of the combination of treatment effect function parameters
      ##
      nseasons  = ifelse(is.null(sim$nseasons), NA, sim$nseasons), 
      season.frequency = ifelse(is.null(sim$season.frequency), NA, sim$season.frequency),
      ## Plotting
      plot.show = plot.show, ## TRUE
      plot.save = plot.save  ## TRUE
    )
  }
  
  # ## Save simulation list as serialized data file
  # simlist.file <- sprintf('single_intervention_SIMLIST_selection_endog_%s_%s.rds',
  #                         length(simlist), sim.id)
  # saveRDS(simlist, file = file.path(dir_plot, simlist.file))
  
  
  
  ##===============================
  ##  BSTS State Specification Comparison 
  ##------------------------------
  for (i in 1:length(simlist))
  {
    key <- names(simlist)[i]
    cat(sprintf('\n%s, %s\n',i, key))
    
    simlist[[key]]$cordf <- data.frame()
    simlist[[key]]$results <- list()
    
    for (k in 1:length(effect.types)) 
    {
      effect.type <- effect.types[k]
      simdf <- simlist[[key]]$sim$df
      simdf <- simdf[simdf$effect.type == effect.type, ]
      # simdf <- simdf[simdf$effect.type==effect.type, ]
      ## Set group name 'gname' field, where 0 = control, # = period of treatment
      simdf$gname <- 0
      simdf$gname[simdf$group=='treatment'] <- simdf$match_pd[simdf$group=='treatment']
      
      ## Compute Multiperiod DiD for Avg Treatment Effect on Treated
      ccattgt <- att_gt(yname = "y", ## "Y",
                        tname = "t",
                        idname = "actor",
                        gname = "gname",
                        xformla = ~c1 + c2 + c3,
                        data = simdf #,
                        # panel = F
      )
      # ccattgt
      
      ## Get first treatment group actor
      tr.actor.1 <- simdf$actor[which(simdf$group=='treatment')[1]]
      
      ## SIMPLE AGGREGATION (OVERALL EFFECT) ATT
      agg.simple <- aggte(ccattgt, type='simple')
      # summary(agg.simple)
      
      ## DYNAMIC EFFECTS AND EVENT STUDIES
      agg.es <- aggte(ccattgt, type = "dynamic")
      # summary(agg.es)
      # tidy(agg.es)
      ggdid(agg.es)
      ggsave(filename = sprintf('%s_did_dynamic_effect_%s_%s_%s.png',
                                prefix,key,effect.type,sim.id))
      
      
      ##-----------------------------
      
      ## Correlation of simulated to inferred
      cormat <- cor(cbind(simdf$b3[simdf$actor==tr.actor.1][-1], agg.es$att.egt))
      # ## MATPLOT of 
      # png(filename = sprintf('single_intervention_DiD_BSTS_DGP_comparison_%s.png',sim.id),
      #     width = 6,height = 6,units = 'in',res=300)
      # matplot(cbind(simdf$b3[simdf$actor==tr.actor.1][-1], agg.es$att.egt), type='o',pch=1:2)
      # dev.off()
      ## 
      sigmamat <- cor(cbind(simdf$x1[simdf$t==intpd], ## treatment dummy at intervention period
                            simdf$y[simdf$t==(intpd-1)] ## performance variable at period before intervention
      ))
      
      tmp.cordf <- data.frame(cor.type=c('dgp.did','x1t.ytm1'),
                              effect.type=rep(effect.type,2),
                              cor=c(cormat[2,1],sigmamat[2,1]))
      
      ## SAVE TO OUTPUT LIST
      simlist[[key]]$cordf <- rbind(simlist[[key]]$cordf, tmp.cordf)  ## off-diagonal element of symmetric correlation matrix
      
      # endog <- data.frame(threshold=c(1/2, 1/3,1/4,1/5,1/6),
      #                     cor=c(-0.2523,-.2030,-.2086,-.2107,-.2214))
      
      # ## Density plot 
      # ggplot(simdf,aes(x=y, colour=effect.type)) + geom_density(alpha=0.1, size=1.4)
      
      
      ##====================
      ## BSTS Timseries Setup
      ##--------------------
      ## Aggregate into timeseries dataframe
      tsdf <- simdf %>%
        filter( ! is.na(match_id)) %>%
        group_by(t, t.post.intpd, effect.type, match_pd, gname, group, group.color) %>%
        dplyr::summarize(
          n_in_pd = n(),
          actors = paste(unique(actor), collapse = '|'),
          y_mean = mean(y, na.rm=T),
          y_sum = sum(y, na.rm=T),
          y_sd = sd(y, na.rm=T),
          y_min = min(y, na.rm=T),
          y_max = max(y, na.rm=T),
          x1_mean = mean(x1, na.rm=T),
          x2_mean = mean(x2, na.rm=T),
          x3_mean = mean(x3, na.rm=T),
          ##
          c1_mean = mean(c1, na.rm=T),
          c2_mean = mean(c2, na.rm=T),
          c3_mean = mean(c3, na.rm=T),
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
      val.cols <- c('y_mean','y_sum','y_min','y_max','y_sd',
                    'x1_mean','x2_mean','x3_mean',
                    'c1_mean','c2_mean','c3_mean',
                    'b1_mean','b2_mean','b3_mean',
                    'u_mean','v_mean')
      ts <- unique(tsdf$t)
      groups <- unique(tsdf$group)
      ## init timeseries dataframe - wide
      tsdfw <- data.frame(t=ts,  stringsAsFactors = F)
      for (j in 1:length(groups)) {
        id.j <- which( tsdf$group == groups[j] ) 
        for (k in 1:length(val.cols)) {
          df.col <- data.frame( tsdf[ id.j , val.cols[k] ] )
          names(df.col) <- sprintf('%s_%s',groups[j],val.cols[k])
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
      dat <- tsdfw[,c('treatment_y_mean','control_y_mean','control_y_sd',#'control_y_sum', 'control_y_min','control_y_max',
                      'control_c1_mean','control_c2_mean','control_c3_mean'#,
                      # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
                      # 'control_u_mean','control_v_mean'
      )]
      ## Train on y pre-treatment but NA's post-treatment
      y.pre.treat.NAs.post.treat <- c(dat$treatment_y_mean[1:(intpd-1)], rep(NA,npds-intpd+1))
      ## Then use the post-treatment response for causal impact estimation
      post.period.response <- dat$treatment_y_mean[intpd:npds]
      ## Covariates (predictors) - Matrix for "formula = y ~ predictors" argument
      predictors <- dat[, ! names(dat) %in% 'treatment_y_mean'] ## remove response; convert to matrix
      # ## Covariates (predictors) - Dataframe for "data" argument
      # predictors <- as.matrix(predictors) 
      
      ## ADD temporal trend to covariates
      predictors$covar_temp_trend <- (1:npds) + rnorm(npds, 0, .001)   ## * simlist$`0.rand.base`$b5
      
      ##----------------------------
      ## State Space Configuration
      ##----------------------------
      nss <- length(simlist[[key]]$bsts.state.spec)
      ## State Space Config list
      st.sp <- list()
      if (nss > 0) {
        for (j in 1:nss) {
          addStateComp <- simlist[[key]]$bsts.state.spec[[j]]
          if (class(addStateComp)=='function') {
            ## STATE COMPONENT FUNCTIONS (excluding seasonality)
            st.sp <- addStateComp(st.sp, y.pre.treat.NAs.post.treat) 
          } else if ('name' %in% names(addStateComp) & addStateComp$name=='AddTrig') {
            ## SEASONALITY
            st.sp <- AddTrig(st.sp, y.pre.treat.NAs.post.treat, 
                             period = addStateComp$nseasons, 
                             frequencies = addStateComp$freq,
                             method = 'harmonic')
          } else {
            cat('\nWarning state component excluded (not a state function or seasonal effect list\n')
          }
        }
      } else {
        ## Default in CausalImpact package
        st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
      }
      # print(st.sp)
      
      ## BSTS model
      bsts.model <- bsts(y.pre.treat.NAs.post.treat ~ . ,
                         state.specification = st.sp,
                         data = predictors,
                         niter = 5000)
      # ## BSTS model for Dynamic Regression
      # bsts.model <- bsts(y.pre.treat.NAs.post.treat,
      #                    state.specification = st.sp,
      #                    niter = 5000)
      
      ## Use BSTS prediction of counterfactual to estimate CausalImpact
      impact_amount <- CausalImpact(bsts.model=bsts.model,
                                    post.period.response = post.period.response,
                                    alpha=0.05, model.args = list(niter = 5000))
      # ##
      # summary(impact_amount)
      # summary(impact_amount$model$bsts.model)
      # plot(impact_amount)
      
      # summary(impact_amount)
      # png(filename=sprintf('single_intervention_BSTS_CausalImpact_plot_%s_%s_%s.png',
      #                         key,effect.type,sim.id))
      plot(impact_amount, c('original','pointwise','cumulative'))
      ggsave(filename = sprintf('%s_bsts_CausalImpact_plot_%s_%s_%s.png',
                                prefix,key,effect.type,sim.id))
      # dev.off()
      
      
      bsts.res <- impact_amount$series
      did.res <- tidy(agg.es)
      # plot(did.res)
      
      # ## AVERAGE TREATMENT EFFECT USED IN SIMULATION
      # simdf
      
      tr.actors <- unique(simdf$actor[which(simdf$group=='treatment')])  
      co.actors <- unique(simdf$actor[which(simdf$group=='control')])
      
      b3diff <- data.frame(
        treat=simdf %>% filter(group=='treatment' & actor==tr.actors[1]) %>% mutate(treat=b3) %>% dplyr::select(treat),
        ctrl=simdf %>% filter(group=='control' & actor==co.actors[1]) %>% mutate(ctrl=b3) %>% dplyr::select(ctrl),
        diff=NA
      )
      b3diff$diff <- b3diff$treat - b3diff$ctrl
      # simdf %>% 
      #   filter(group=='control' & actor==tr.actors[1]) %>% 
      #   dplyr::mutate(b3.diff=b3.treat-b3.ctrl) %>% 
      #   dplyr::select(b3.diff) 
      
      ## Results comparison table
      res.tbl <- cbind(bsts.res[ ,c('point.effect','point.effect.lower','point.effect.upper')],
                       did.res[ ,c('term','event.time','estimate','point.conf.low','point.conf.high')],
                       b3.treat=b3diff$treat,
                       b3.ctrl=b3diff$ctrl,
                       b3.att=b3diff$diff 
      )
      ## ROUND RESULTS TABLE (round numeric columns)
      # res.tbl <- res.tbl[ ! is.na(res.tbl$estimate), ]
      num.cols <- names(res.tbl)[ ! names(res.tbl) %in% c('term','event.time') ]
      res.tbl[ , num.cols] <- round( as.numeric(res.tbl[ , num.cols]), 4)
      # for (i in 1:length(num.cols)) {
      #   res.tbl[ , num.cols[i] ] <- round( as.numeric( res.tbl[ , num.cols[i] ] ), 4)
      # }
      ## MOVE ID COLUMNS TO FRONT
      .col.idx <- which(names(res.tbl) %in% c('term','event.time'))
      res.tbl4 <- cbind(res.tbl[, .col.idx],  res.tbl[, -.col.idx] )
      # View(res.tbl4)
      
      ##PLOT INCLUSION PROBABILITIES
      png(filename = sprintf('%s_BSTS_inclusion_probs_%s_%s_%s.png',
                             prefix,key,effect.type,sim.id))
      plot(impact_amount$model$bsts.model,'coefficients', main=sprintf('%s %s', key,effect.type))
      dev.off()
      
      ## PLOT DYNAMIC EFFECTS COMPARISON - DID vs. BSTS vs. DGP
      png(filename = sprintf('%s_BSTS_dynamic_treatment_effect_comparison_%s_%s_%s.png',
                             prefix,key,effect.type,sim.id))
      res.tbl.filename <- sprintf('%s: DGP = %.3f; DiD = %.3f;  BSTS = %.3f',
                                  key, mean(b3diff$diff[intpd:nrow(b3diff)]), agg.es$overall.att, impact_amount$summary$AbsEffect[1])
      matplot(x = res.tbl4$event.time, y=res.tbl4[,c('point.effect','estimate','b3.att')],
              type='l',lty=c(1,2,4),lwd=c(1,1,2),col=c('black','red','blue'),
              main=res.tbl.filename, ylab='ATT',xlab='t')
      legend('topright',legend=c('BSTS','DiD','DGP'),col=c('black','red','blue'),lty=c(1,2,4),lwd=c(1,1,2)) 
      dev.off()
      
      ##===============================================================
      ## 1-step ahead prediction error
      bsts.pred.er <- bsts.prediction.errors(impact_amount$model$bsts.model)$in.sample[,(1:(intpd-1))]
      
      ## Append results to output list
      simlist[[key]]$results[[effect.type]]$did.attgt <- ccattgt
      simlist[[key]]$results[[effect.type]]$did.agg.simple <- agg.simple
      simlist[[key]]$results[[effect.type]]$did.agg.es <- agg.es
      simlist[[key]]$results[[effect.type]]$bsts.CausalImpact <- impact_amount
      simlist[[key]]$results[[effect.type]]$bsts.cumu.pred.error <-  cumsum(colSums(abs(bsts.pred.er)))
      simlist[[key]]$results[[effect.type]]$res.tbl <- res.tbl4
      
    }
    
  } #; dev.off()
  
  return(simlist)
}



################################################################################################################
################################################################################################################





#######################################################
##
## SIMULATION SETUP
##
#######################################################
n <- 300    ## Number of firms
npds <- 100 ## 100  ## Number of periods
intpd <- round( npds * 0.6 )   ## 60% pre-intervention training / 40% post-intervention

## dynamic treatment effect types (shapes)
effect.types <- c('quadratic') ## 'constant', 'geometric'

# Specify simulation parameters
simlist <- list(
  # `0.rand.base`=list(
  #   treat.rule='random',  ## random or below.benchmark (endogenous self-selection)
  #   treat.prob=0.5,       ## probability of being in treatment group
  #   treat.threshold=NA,   ## if treat.rule=below.benchmark, treat.threshold is the quantile below which firms evaluate self-selecting treatment with probability = treat.prob
  #   b4=0, ## b4 = seasonal component effect (weight)
  #   b5=0,  ## b5 = Linear Trend:  growth rate
  #   bsts.st.sp=list()
  # ),
  ##----------- ONLY TREND ------------
  ## Correct Spec
  `1a`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=0, b5=.02,
    bsts.state.spec=list(AddLocalLevel)
  ),
  ## Wrong Spec: using student local lienar trend (or other local trend instead of local level)
  `1b`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=0, b5=.02,
    bsts.state.spec=list(AddStudentLocalLinearTrend)
  ),
  ## Wrong Spec: using student local lienar trend (or other local trend instead of local level)
  `1c`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=0, b5=.02,
    bsts.state.spec=list(AddGeneralizedLocalLinearTrend)
  ),
  ##_-------- Seasonality --------------------
  ## True=seasons[12] | Model=none
  `2a`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=1, b5=0,
    nseasons=12, season.frequency=1,
    bsts.state.spec=list() 
  ),
  ## True=seasons[12] | Model=AddTrig[12]
  `2b`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=1, b5=0,
    nseasons=12, season.frequency=1,
    bsts.state.spec=list(
      list(name='AddTrig', nseasons=12, freq=1)
    ) 
  ),
  ## True=seasons[12] | Model=AddTrig[7]**incorrect
  `2c`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=1, b5=0,
    nseasons=20, season.frequency=1,
    bsts.state.spec=list(
      list(name='AddTrig', nseasons=7, freq=1)
    ) 
  ),
  ### True=seasons[20] + growth | Model=AddTrig[7]  + LocalLevel
  `2d`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=1, b5=.03,
    nseasons=20, season.frequency=1,
    bsts.state.spec=list(
      AddLocalLevel,
      list(name='AddTrig', nseasons=7, freq=1)
    ) 
  ),
  ### True=seasons[20] + growth | Model=AddTrig[7]**incorrect** + locallevel
  `2e`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=1, b5=.03,
    nseasons=20, season.frequency=1,
    bsts.state.spec=list(
      AddStudentLocalLinearTrend,
      list(name='AddTrig', nseasons=7, freq=1)
    ) 
  ),
  ##
  `3`=list(
    treat.rule='random', treat.prob=0.5, treat.threshold=NA,
    b4=1, b5=.03,
    nseasons=20, season.frequency=1,
    bsts.state.spec=list(
      AddStudentLocalLinearTrend,
      list(name='AddTrig', nseasons=7, freq=1)
    ),
    b8=.1
  )#,
  ##-------------- AR -------------------------------
  ##
  ##_--------- BOTH AR & TREND ----------------------
  # ## Correct Spec
  # `3a.rand.trendAr.ssTrendAr`=list(
  #   treat.rule='random', treat.prob=0.5, treat.threshold=NA,
  #   b4=1, b5=0,
  #   nseasons=12, season.frequency=1,
  #   bsts.state.spec=list(
  #     list(func=AddTrig, nseasons=12, freq=1, type='seasonal')
  #   ) 
  # ),
  # ## ? Possibly Correct Spec ?
  # `3b.rand.trendAr.ssGenTrend`=list(
  #   treat.rule='random', treat.prob=0.5, treat.threshold=NA,
  #   b4=1, b5=0,
  #   nseasons=12, season.frequency=1,
  #   bsts.state.spec=list(
  #     list(func=AddTrig, nseasons=12, freq=1, type='seasonal')
  #   ) 
  # ),
  # ## Wrong Spec ?  student local lienar trend
  # `3c.rand.trendAr.ssStudentTrend`=list(
  #   treat.rule='random', treat.prob=0.5, treat.threshold=NA,
  #   b4=1, b5=0,
  #   nseasons=12, season.frequency=1,
  #   bsts.state.spec=list(
  #     list(func=AddTrig, nseasons=12, freq=1, type='seasonal')
  #   ) 
  # )
)

## --- TEST ---
# simlist <- list(
#   no.perf.no.gro=list(
#       b4=0, ## b4 = past performance spillover (or persistence)
#       b5=0  ## b5 = growth rate
#     )
# )
##-------------




##================================
##
## MAIN SIMULATION COMPARISON RUN
##
##_-------------------------------
# sim.id <- round(10*as.numeric(Sys.time()))
simlist <- runSimBstsDiDComparison(simlist, n, npds, intpd, effect.types, plot.show = T, plot.save = T)

## Save simulation list as serialized data file
simlist.file <- sprintf('%s_SIMLIST_selection_%s_%s.rds',
                        prefix, length(simlist), simlist$`1a`$sim$id)
saveRDS(simlist, file = file.path(dir_proj, simlist.file))










##=====================================
##
## BUILD META SIMLIST
##  - GRID SEARCH ALL DIMENSIONS OF SIMULATION
##
##=====================================

npds <- 100
# effect.types <- c('constant','quadratic','geometric')

#  BSTS State Space Component Functions 
# bsts::AddAr
# bsts::AddAutoAr
# bsts::AddDynamicRegression
# bsts::AddGeneralizedLocalLinearTrend
# bsts::AddHierarchicalRegressionHoliday
# bsts::AddLocalLevel
# bsts::AddLocalLinearTrend
# bsts::AddMonthlyAnnualCycle
# bsts::AddRandomWalkHoliday
# bsts::AddRegressionHoliday
# bsts::AddSeasonal
# bsts::AddSemilocalLinearTrend
# bsts::AddSharedLocalLevel
# bsts::AddStaticIntercept
# bsts::AddStudentLocalLinearTrend
# bsts::AddTrig
# ## Make BSTS state space configurations
state.configs <- list(
  list(
    AddLocalLevel
  ),
  list(
    AddLocalLinearTrend
  ),
  list(
    AddStudentLocalLinearTrend
  ),
  list(
    AddAr
  ),
  list(
    list(name='AddTrig', nseasons=12, freq=1)
  ),
  list(
    AddLocalLevel, 
    list(name='AddTrig', nseasons=12, freq=1)
  ),
  list(
    AddLocalLevel, 
    AddAr
  ),
  list(
    list(name='AddTrig', nseasons=12, freq=1), 
    AddAr
  ),
  list(
    AddLocalLevel,
    list(name='AddTrig', nseasons=12, freq=1), 
    AddAr
  )
)



actors.sizes <- c(30,300)
noise.levels <- c(1, .01)
self.selections <- c('random', 'below.benchmark')
intpds <- c(2/3, 1/3)   ### round(c(3*npds/4, npds/2, npds/4))

autocors <- c(0, .5)

bsts.nseasons <- c(12, 7) ## correct, incorrect [for DGP_nseasons = 12]
seasonalities <- c(FALSE, TRUE)
linear.trends <- c(0, .05)

##----------------------------------
## MAIN SIMULATION SEARCH LOOP
##----------------------------------
simlist <- list()
# for (l in 1:length(bsts.nseasons)) {
#   bsts.nseason <- bsts.nseasons[l]
for (g in 1:length(actors.sizes)) {
  actor.size <- actor.sizes[g]
  
  for (h in 1:length(noise.levels)) {
    noise.level <- noise.levels[h]
    
    for (i in 1:length(self.selections)) {
      treat.rule <- self.selections[i]
      
      for (ii in 1:length(treat.thresholds)) {
        treat.threshold <- treat.thresholds[ii]
    
        for (j in 1:length(intpds)) {
          intpd <- intpds[j]
          
          for (k in 1:length(autocors)) {
            autocor <- autocors[k]

            
              for (m in 1:length(seasonalities)) {
                seasonality <- seasonalities[m]
                
                for (n in 1:length(linear.trends)) {
                  linear.trend <- linear.trends[n]
                  
                  for (p in 1:length(state.configs)) {
                    state.config <- state.configs[p]
                    
                    key <- sprintf('g%sh%si%sii%sj%sk%sm%sn%sp%s',
                                   g,h,i,ii,j,k,m,n,p)
                    
                    simlist[[key]] <- list(
                      n = actor.size,    ## Number of firms
                      npds = 100,  ## Number of periods
                      intpd = intpd, ## 60% pre-intervention training / 40% post-intervention
                      ##
                      noise.level = noise.level, ## stdev of simulated noise terms
                      ##
                      treat.rule = treat.rule, 
                      treat.prob = ifelse(treat.rule=='random', 0.5, 1), 
                      treat.threshold = ifelse(treat.rule=='random', 1, treat.threshold),
                      b4 = 1,   ## seasonal component weight
                      b5 = linear.trend, ##
                      b9 = autocor  , ## autocorrelation
                      nseasons=12, 
                      season.frequency=1,
                      bsts.state.spec=state.config, 
                      rand.seed = 54321,
                    )
                    
                    
                  }
                  
                }
                
              }
      
            
          }
          
        }
        
      }
      
    }
    
  }
  
}

# }



## PRE DEBUG ###
config1 <- c(g=1,h=1,i=1,j=1,k=1,l=1,m=1,n=1,p=1)
config2 <- c(g=2,h=2,i=2,j=2,k=2,l=2,m=2,n=2,p=2)
###
simlist.debug <- list(
  list(
    
  ),
  list(
    
  )
)
runSimBstsDiDComparison(simlist.debug, n, npds, intpd, effect.types, plot.show = T, plot.save = T)




########################################
##            DEBUG
##
#######################################
# key <- '3a.rand.trendAr.ssTrendAr'
key <- '1b.rand.trend.ssNone'
xsim <- simlist[[key]]$sim
xres <- simlist[[key]]$results


xbstsmod <- xres$quadratic$bsts.CausalImpact$model$bsts.model
plot(xres$quadratic$bsts.CausalImpact)
plot(xbstsmod)







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
# bsts.state.spec <- list(AddAutoAr, AddStudentLocalLinearTrend)
# nss <- length(bsts.state.spec)
# if (nss > 0)
# {
#   for (j in 1:nss) {
#     addStateComponentFunc <- bsts.state.spec[[j]]
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
