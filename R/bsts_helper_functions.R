################################################
##
##  BSTS HELPER FUNCTIONS
##
################################################
require(bsts)
require(Boom)
require(BoomSpikeSlab)










##==============================================
##
##
##
##
##
##
## BSTS STATE SPACE COMPONENTS FUNCTIONS
##
##
##
##
##
##----------------------------------------------
###
##  Add a BSTS state space component to state.space list
##    by applying the arguments in state space component configuration state.conf 
##    in function indicated by state.conf$name
###
updateStateSpaceAddComponentFromConfList <- function(state.space,  y.pre.treat.NAs.post.treat,  state.conf) {
  
  if (class(state.conf) != 'list') {
    cat('\nWARNING: state.conf is not a list. Try wrapping the state component function in a state config list.\n')
    return()
  }
  if (is.null(state.conf$name)) {
    state.conf$name <- ''
  }
  
  if (  state.conf$name == 'AddAr') {
    
    state.space <- AddAr(state.space, y.pre.treat.NAs.post.treat, 
                         lags = state.conf$lags, 
                         sigma.prior = state.conf$sigma.prior, 
                         initial.state.prior = state.conf$initial.state.prior, 
                         sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddAutoAr') {
    
    state.space <- AddAutoAr(state.space, y.pre.treat.NAs.post.treat, 
                             lags = state.conf$lags, 
                             prior = state.conf$prior, 
                             sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddGeneralizedLocalLinearTrend') {
    
    state.space <- AddGeneralizedLocalLinearTrend(state.space, y.pre.treat.NAs.post.treat, 
                                                  level.sigma.prior = state.conf$level.sigma.prior, 
                                                  slope.mean.prior = state.conf$slope.mean.prior, 
                                                  slope.ar1.prior = state.conf$slope.ar1.prior, 
                                                  slope.sigma.prior = state.conf$slope.sigma.prior, 
                                                  initial.level.prior = state.conf$initial.level.prior, 
                                                  initial.slope.prior = state.conf$initial.slope.prior, 
                                                  sdy = state.conf$sdy, 
                                                  initial.y = state.conf$initial.y)
    
  } else if (state.conf$name == 'AddHierarchicalRegressionHoliday') {
    
    state.space <- AddHierarchicalRegressionHoliday(state.space, y.pre.treat.NAs.post.treat, 
                                                    holiday.list = state.conf$holiday.list, 
                                                    coefficient.mean.prior = state.conf$coefficient.mean.prior, 
                                                    coefficient.variance.prior = state.conf$coefficient.variance.prior, 
                                                    time0 = state.conf$time0, 
                                                    sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddLocalLevel') {
    
    state.space <- AddLocalLevel(state.space, y.pre.treat.NAs.post.treat, 
                                 sigma.prior = state.conf$sigma.prior, 
                                 initial.state.prior = state.conf$initial.state.prior, 
                                 sdy = state.conf$sdy, 
                                 initial.y = state.conf$initial.y)
    
  } else if (state.conf$name == 'AddLocalLinearTrend') {
    
    state.space <- AddLocalLinearTrend(state.space, y.pre.treat.NAs.post.treat, 
                                       level.sigma.prior = state.conf$level.sigma.prior, 
                                       slope.sigma.prior = state.conf$slope.sigma.prior, 
                                       initial.level.prior = state.conf$initial.level.prior, 
                                       initial.slope.prior = state.conf$initial.slope.prior, 
                                       sdy = state.conf$sdy, 
                                       initial.y = state.conf$initial.y)
    
  } else if (state.conf$name == 'AddMonthlyAnnualCycle') {
    
    state.space <- AddMonthlyAnnualCycle(state.space, y.pre.treat.NAs.post.treat, 
                                         date.of.first.observation = state.conf$date.of.first.observation, 
                                         sigma.prior = state.conf$sigma.prior, 
                                         initial.state.prior = state.conf$initial.state.prior, 
                                         sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddRandomWalkHoliday') {
    
    state.space <- AddRandomWalkHoliday(state.specification = ,y = ,
                                        holiday = state.conf$holiday,
                                        time0 = state.conf$time0,
                                        sigma.prior = state.conf$sigma.prior,
                                        initial.state.prior = state.conf$initial.state.prior,
                                        sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddRegressionHoliday') {
    
    state.space <- AddRegressionHoliday(state.space, y.pre.treat.NAs.post.treat,
                                        holiday.list = state.conf$holiday.list,
                                        time0 = state.conf$time0,
                                        prior = state.conf$prior,
                                        sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddSeasonal') {
    
    state.space <- AddSeasonal(state.space, y.pre.treat.NAs.post.treat,
                               nseasons = state.conf$nseasons,
                               season.duration = state.conf$season.duration,
                               sigma.prior = state.conf$sigma.prior,
                               initial.state.prior = state.conf$initial.state.prior,
                               sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddSemilocalLinearTrend') {
    
    state.space <- AddSemilocalLinearTrend(state.space, y.pre.treat.NAs.post.treat,
                                           level.sigma.prior = state.conf$level.sigma.prior,
                                           slope.mean.prior = state.conf$slope.mean.prior,
                                           slope.ar1.prior = state.conf$slope.ar1.prior ,
                                           slope.sigma.prior = state.conf$slope.sigma.prior,
                                           initial.level.prior = state.conf$initial.level.prior,
                                           initial.slope.prior = state.conf$initial.slope.prior,
                                           sdy = state.conf$sdy,
                                           initial.y = state.conf$initial.y)
    
  } else if (state.conf$name == 'AddSharedLocalLevel') {
    
    state.space <- AddSharedLocalLevel(state.space, y.pre.treat.NAs.post.treat, 
                                       nfactors = state.conf$nfactors, 
                                       coefficient.prior = state.conf$coefficient.prior, 
                                       initial.state.prior = state.conf$initial.state.prior, 
                                       timestamps = state.conf$timestamps, 
                                       series.id = state.conf$series.id, 
                                       sdy = state.conf$sdy)
    
  } else if (state.conf$name == 'AddStaticIntercept') {
    
    state.space <- AddStaticIntercept(state.space, y.pre.treat.NAs.post.treat, 
                                      initial.state.prior = state.conf$initial.state.prior)
    
  } else if (state.conf$name == 'AddStudentLocalLinearTrend') {
    
    state.space <- AddStudentLocalLinearTrend(state.space, y.pre.treat.NAs.post.treat, 
                                              save.weights = state.conf$save.weights, 
                                              level.sigma.prior = state.conf$level.sigma.prior, 
                                              level.nu.prior = state.conf$level.nu.prior, 
                                              slope.sigma.prior = state.conf$slope.sigma.prior, 
                                              slope.nu.prior = state.conf$slope.nu.prior, 
                                              initial.level.prior = state.conf$initial.level.prior, 
                                              initial.slope.prior = state.conf$initial.slope.prior, 
                                              sdy = state.conf$sdy, 
                                              initial.y = state.conf$initial.y)
    
  } else if (state.conf$name == 'AddTrig') {
    
    state.space <- AddTrig(state.space, y.pre.treat.NAs.post.treat, 
                           period = state.conf$period, 
                           frequencies = state.conf$frequencies,
                           sigma.prior = state.conf$sigma.prior,
                           initial.state.prior = state.conf$initial.state.prior,
                           sdy = state.conf$sdy,
                           method = 'harmonic')
    
  } else {
    ## PASS
    cat('\nNo state components added from this state.conf. Returning original state.space\n')
    print(state.conf)
  }
  
  return(state.space)
  
}



##=======================================
## Get BSTS State Space Component Template
##---------------------------------------
getStateSpaceConf <- function(name, ...) {
  
  args <- list(...)
  ### MAIN CONFIGUATION SELECTION BY COMPONENT NAME
  #
  conf <- if ( name == 'AddAr') {
    
    # list(
    #   name='AddAr',
    #   # func=AddAr,
    #   # state.space = NULL, 
    #   # y.pre.treat.NAs.post.treat = NULL, 
    #   lags = if(!is.null(args$lags)){args$lags} else {1}, 
    #   sigma.prior = if(!is.null(args$sigma.prior)){args$sigma.prior} else {SdPrior(3.0, 1.0)},  ## Boom::SdPrior(sigma.guess, sample.size = .01, initial.value = sigma.guess, fixed = FALSE, upper.limit = Inf)
    #   initial.state.prior = if(!is.null(args$initial.state.prior)){args$initial.state.prior} else {NULL}, ## or use a Boom::MvnPrior(mean, variance)
    #   sdy = args$sdy
    # )
    list(
      name='AddAr',
      # func=AddAr,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL, 
      lags = args$lags, 
      sigma.prior = args$sigma.prior,  ## Boom::SdPrior(sigma.guess, sample.size = .01, initial.value = sigma.guess, fixed = FALSE, upper.limit = Inf)
      initial.state.prior = args$initial.state.prior, ## or use a Boom::MvnPrior(mean, variance)
      sdy = args$sdy
    )
    
  } else if (name == 'AddAutoAr') {
    
    list(
      name='AddAutoAr',
      # func=AddAutoAr,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL, 
      lags = args$lags, 
      prior = args$prior, ##bsts::SpikeSlabArPrior()
      sdy = args$sdy
    )
    
  } else if (name == 'AddHierarchicalRegressionHoliday') {
    
    list(
      name='AddHierarchicalRegressionHoliday',
      # func=AddHierarchicalRegressionHoliday,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL, 
      holiday.list = args$holiday.list, 
      coefficient.mean.prior = args$coefficient.mean.prior, 
      coefficient.variance.prior = args$coefficient.variance.prior, 
      time0 = args$time0, 
      sdy = args$sdy
    )
    
  } else if (name == 'AddLocalLevel') {
    
    list(
      name='AddLocalLevel',
      # func=AddLocalLevel,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL, 
      sigma.prior = args$sigma.prior, 
      initial.state.prior = args$initial.state.prior, 
      sdy = args$sdy, 
      initial.y = args$initial.y
    )
    
  } else if (name == 'AddLocalLinearTrend') {
    
    list(
      name='AddLocalLinearTrend',
      # func=AddLocalLinearTrend,
      # state.space = NULL,
      # y.pre.treat.NAs.post.treat = NULL,
      level.sigma.prior = args$level.sigma.prior, 
      slope.sigma.prior = args$slope.sigma.prior, 
      initial.level.prior = args$initial.level.prior, 
      initial.slope.prior = args$initial.slope.prior, 
      sdy = args$sdy, 
      initial.y = args$initial.y
    )
    
    
  } else if (name == 'AddMonthlyAnnualCycle') {
    
    list(
      name ='AddMonthlyAnnualCycle',
      # func =AddMonthlyAnnualCycle, 
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL, 
      date.of.first.observation = args$date.of.first.observation, 
      sigma.prior = args$sigma.prior, 
      initial.state.prior = args$initial.state.prior, 
      sdy = args$sdy
    )
    
  } else if (name == 'AddRandomWalkHoliday') {
    
    list(
      name = 'AddRandomWalkHoliday',
      # func = AddRandomWalkHoliday,
      # state.specification = NULL,
      # y = NULL,
      holiday = args$holiday,
      time0 = args$time0,
      sigma.prior = args$sigma.prior,
      initial.state.prior = args$initial.state.prior,
      sdy = args$sdy
    )
    
  } else if (name == 'AddRegressionHoliday') {
    
    list(
      name = 'AddRegressionHoliday',
      # func = AddRegressionHoliday,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL,
      holiday.list = args$holiday.list,
      time0 = args$time0,
      prior = args$prior,
      sdy = args$sdy
    )
    
  } else if (name == 'AddSeasonal') {
    
    list(
      name = 'AddSeasonal',
      # func = AddSeasonal,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL,
      nseasons = args$nseasons,
      season.duration = args$season.duration,
      sigma.prior = args$sigma.prior,
      initial.state.prior = args$initial.state.prior,
      sdy = args$sdy
    )
    
  } else if (name %in% c('AddSemilocalLinearTrend','AddGeneralizedLocalLinearTrend')) {
    
    list(
      name = 'AddSemilocalLinearTrend',
      # alias = 'AddGeneralizedLocalLinearTrend',
      # func = AddSemilocalLinearTrend,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL,
      level.sigma.prior = args$level.sigma.prior,
      slope.mean.prior = args$slope.mean.prior,
      slope.ar1.prior = args$slope.ar1.prior ,
      slope.sigma.prior = args$slope.sigma.prior,
      initial.level.prior = args$initial.level.prior,
      initial.slope.prior = args$initial.slope.prior,
      sdy = args$sdy,
      initial.y = args$initial.y
    )
    
  } else if (name == 'AddSharedLocalLevel') {
    
    list(
      name = 'AddSharedLocalLevel',
      # func = AddSharedLocalLevel,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL, 
      nfactors = args$nfactors, 
      coefficient.prior = args$coefficient.prior, 
      initial.state.prior = args$initial.state.prior, 
      timestamps = args$timestamps, 
      series.id = args$series.id, 
      sdy = args$sdy
    )
    
  } else if (name == 'AddStaticIntercept') {
    
    list(
      name = 'AddStaticIntercept',
      # func = AddStaticIntercept,
      # state.space = NULL,
      # y.pre.treat.NAs.post.treat = NULL, 
      initial.state.prior = args$initial.state.prior
    )
    
  } else if (name == 'AddStudentLocalLinearTrend') {
    
    list(
      name = 'AddStudentLocalLinearTrend',
      # func = AddStudentLocalLinearTrend,
      # state.space = NULL, 
      # y.pre.treat.NAs.post.treat = NULL, 
      save.weights = args$save.weights, 
      level.sigma.prior = args$level.sigma.prior, 
      level.nu.prior = args$level.nu.prior, 
      slope.sigma.prior = args$slope.sigma.prior, 
      slope.nu.prior = args$slope.nu.prior, 
      initial.level.prior = args$initial.level.prior, 
      initial.slope.prior = args$initial.slope.prior, 
      sdy = args$sdy, 
      initial.y = args$initial.y
    )
    
  } else if (name == 'AddTrig') {
    
    list(
      name = 'AddTrig',
      # func = AddTrig,
      # state.space = NULL,
      # y.pre.treat.NAs.post.treat = NULL, 
      period = args$period, 
      frequencies = args$frequencies,
      sigma.prior = args$sigma.prior,
      initial.state.prior = args$initial.state.prior,
      sdy = args$sdy,
      method = 'harmonic'
    )
    
  } else {
    list()
  }
  
  return(conf)
  
} ## end func







##=======================================
## Get BSTS State Space Component Template By Simulation Scenario Settings
##---------------------------------------
getStateSpaceConfBySimScenario <- function(name, scenario, ## c('sd.high','sd.low')
                                           x.spikslab.prior=NULL, ## X  values for Boom::SpikeSlabPrior()
                                           y.spikslab.prior=NULL, ## Y  values for Boom::SpikeSlabPrior()
                                           ...) {
  ## Named arguments
  args <- list(...)
  
  ## Default params
  # lags.hi <- 3
  # lags.lo <- 1
  lags.hi <- 3
  lags.lo <- 1
  ## Prior SD
  sig.hi <- .06
  sig.lo <- .015
  ## prior SD initial value
  sig.init.val.hi <- .06
  sig.init.val.lo <- .015
  ## Mean of prior distributions
  mu.hi <- .005
  mu.lo <- .001
  ## Autoregressive component
  ar.mu.hi <- .06
  ar.mu.lo <- .01
  ar.forc.sta <- FALSE
  ar.forc.pos <- FALSE
  ## Weight given to sigma.guess (interpretable as a prior observation count)
  samp.size.prop.hi <- .06
  samp.size.prop.lo <- .015
  ##
  upper.limit <- Inf
  ## Spike and Slab priors on the AddAutoAr (Automatically selected Autoregressive lags)
  spikslab.init.val.hi <- .06
  spikslab.init.val.lo <- .015
  ##
  spikslab.disc.fac.hi <- .5
  spikslab.disc.fac.lo <- .25
  ##
  # spikslab.exp.mod.size.hi <- 3
  # spikslab.exp.mod.size.lo <- 2
  ##
  expected.r2.hi <- .5
  expected.r2.lo <- .1
  ## Latent factors (e.g., AddSharedLocalLevel() )
  n.latent.factors.hi <- 3
  n.latent.factors.lo <- 1
  ## SpikeSlabPrior
  diag.shrinkage.hi <- .5
  diag.shrinkage.lo <- 0 ## set=0 --> gives Zellner's G-prior @see spike.slab.prior {BoomSpikeSlab}
  ##
  bsts.nseasons <- 52  ## 12
  bsts.freq <- 1
  
  ##-------------------------------------------------
  ### MAIN CONFIGUATION SELECTION BY COMPONENT NAME
  ##-------------------------------------------------
  if ( name == 'AddAr') {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {  ## sd.hi, sd.high
      list(name='AddAr', 
           lags = lags.hi, 
           sigma.prior =  SdPrior(sigma.guess=sig.hi, sample.size=samp.size.prop.hi, initial.value=sig.init.val.hi, upper.limit=upper.limit),  ## Boom::SdPrior(sigma.guess, sample.size = .01, initial.value = sigma.guess, fixed = FALSE, upper.limit = Inf)
        initial.state.prior = MvnPrior(mean = rep(sig.init.val.hi, lags.hi), variance = matrix(rep(sig.hi,lags.hi^2),nrow=lags.hi))#, ## or use a Boom::MvnPrior(mean, variance)
        # sdy = args$sdy
      )
    } else if (grepl('sd.lo', scenario, perl = F)) { ## sd.lo, sd.low
      list(name='AddAr', 
           lags = lags.lo, 
           sigma.prior =  SdPrior(sigma.guess=sig.lo, sample.size=samp.size.prop.lo, initial.value=sig.init.val.lo, upper.limit=upper.limit),  ## Boom::SdPrior(sigma.guess, sample.size = .01, initial.value = sigma.guess, fixed = FALSE, upper.limit = Inf)
           initial.state.prior = MvnPrior(mean = rep(sig.init.val.lo, lags.lo), variance = matrix(rep(sig.lo,lags.lo^2),nrow=lags.lo))#, ## or use a Boom::MvnPrior(mean, variance)
           # sdy = args$sdy
      )
    } else {
      list()
    }
    return(conf)
    
  } else if (name == 'AddAutoAr') {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(name='AddAutoAr',
            lags = lags.hi, 
            prior = SpikeSlabArPrior(
              lags = lags.hi,
              prior.inclusion.probabilities = GeometricSequence(length=lags.hi, initial.value=1, discount.factor=spikslab.disc.fac.hi),
              prior.mean = rep(0, lags.hi),
              prior.sd = GeometricSequence(lags.hi, initial.value=1, discount.factor=spikslab.disc.fac.hi),
              sdy=.5,
              prior.df = 1,
              expected.r2 = .5,
              sigma.upper.limit = Inf,
              truncate = TRUE
              )#, ##bsts::SpikeSlabArPrior()
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name='AddAutoAr',
           lags = lags.lo, 
           prior = SpikeSlabArPrior(
             lags = lags.lo,
             prior.inclusion.probabilities = GeometricSequence(length=lags.lo, initial.value=1, discount.factor=spikslab.disc.fac.lo),
             prior.mean = rep(0, lags.lo),
             prior.sd = GeometricSequence(lags.lo, initial.value=1, discount.factor=spikslab.disc.fac.lo),
             sdy=.5,
             prior.df = 1,
             expected.r2 = .5,
             sigma.upper.limit = Inf,
             truncate = TRUE
           )#, ##bsts::SpikeSlabArPrior()
      )
    } else {
      list()
    }
    return(conf)
    
  # } else if (name == 'AddHierarchicalRegressionHoliday') {
  #   
  #   conf <- if (grepl('sd.hi', scenario, perl = F)) {
  #     list(
  #       name='AddHierarchicalRegressionHoliday',
  #       holiday.list = args$holiday.list, 
  #       coefficient.mean.prior = args$coefficient.mean.prior, 
  #       coefficient.variance.prior = args$coefficient.variance.prior, 
  #       time0 = args$time0, 
  #       sdy = args$sdy
  #     )
  #   } else if (grepl('sd.lo', scenario, perl = F)) {
  #     
  #   } else {
  #     list()
  #   }
  #   return(conf)
    
  } else if (name == 'AddLocalLevel') {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(name='AddLocalLevel',
           sigma.prior = SdPrior(sigma.guess=sig.hi, sample.size=samp.size.prop.hi, initial.value=sig.init.val.hi, upper.limit=upper.limit),
           initial.state.prior = NormalPrior(mu=mu.hi, sigma=sig.hi, initial.value=sig.init.val.hi)#,
           # sdy = args$sdy,
           # initial.y = args$initial.y
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name='AddLocalLevel',
           sigma.prior = SdPrior(sigma.guess=sig.lo, sample.size=samp.size.prop.lo, initial.value=sig.init.val.lo, upper.limit=upper.limit),
           initial.state.prior = NormalPrior(mu=mu.lo, sigma=sig.lo, initial.value=sig.init.val.lo)#,
           # sdy = args$sdy,
           # initial.y = args$initial.y
      )
    } else {
      list()
    }
    return(conf)
    
  } else if (name == 'AddLocalLinearTrend') {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(
        name='AddLocalLinearTrend',
        level.sigma.prior= SdPrior(sigma.guess =sig.hi, sample.size = samp.size.prop.hi, initial.value = sig.init.val.hi, upper.limit = upper.limit),
        level.nu.prior= NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi), ### Inherits: DoubleModel(),
        slope.sigma.prior = SdPrior(sigma.guess = sig.hi, sample.size = samp.size.prop.hi, initial.value = sig.init.val.hi, upper.limit = upper.limit),
        slope.nu.prior= NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi), ### Inherits: DoubleModel(),
        initial.level.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi),
        initial.slope.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi)
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(
        name='AddLocalLinearTrend',
        level.sigma.prior= SdPrior(sigma.guess =sig.lo, sample.size = samp.size.prop.lo, initial.value = sig.init.val.lo, upper.limit = upper.limit),
        level.nu.prior= NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo), ### Inherits: DoubleModel(),
        slope.sigma.prior = SdPrior(sigma.guess = sig.lo, sample.size = samp.size.prop.lo, initial.value = sig.init.val.lo, upper.limit = upper.limit),
        slope.nu.prior= NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo), ### Inherits: DoubleModel(),
        initial.level.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo),
        initial.slope.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo)
      )
    } else {
      list()
    }
    return(conf)

    
  # } else if (name == 'AddMonthlyAnnualCycle') {
  #   
  #   conf <- if (grepl('sd.hi', scenario, perl = F)) {
  #     
  #   } else if (grepl('sd.lo', scenario, perl = F)) {
  #     list(
  #       name ='AddMonthlyAnnualCycle',
  #       # func =AddMonthlyAnnualCycle, 
  #       # state.space = NULL, 
  #       # y.pre.treat.NAs.post.treat = NULL, 
  #       date.of.first.observation = args$date.of.first.observation, 
  #       sigma.prior = args$sigma.prior, 
  #       initial.state.prior = args$initial.state.prior, 
  #       sdy = args$sdy
  #     )
  #   } else {
  #     list()
  #   }
  #   return(conf)

    
  # } else if (name == 'AddRandomWalkHoliday') {
  #   
  #   conf <- if (grepl('sd.hi', scenario, perl = F)) {
  #     
  #   } else if (grepl('sd.lo', scenario, perl = F)) {
  #     list(
  #       name = 'AddRandomWalkHoliday',
  #       # func = AddRandomWalkHoliday,
  #       # state.specification = NULL,
  #       # y = NULL,
  #       holiday = args$holiday,
  #       time0 = args$time0,
  #       sigma.prior = args$sigma.prior,
  #       initial.state.prior = args$initial.state.prior,
  #       sdy = args$sdy
  #     )
  #   } else {
  #     list()
  #   }
  #   return(conf)
   
    
  # } else if (name == 'AddRegressionHoliday') {
  #   
  #   conf <- if (grepl('sd.hi', scenario, perl = F)) {
  #     list(
  #       name = 'AddRegressionHoliday',
  #       # func = AddRegressionHoliday,
  #       # state.space = NULL, 
  #       # y.pre.treat.NAs.post.treat = NULL,
  #       holiday.list = args$holiday.list,
  #       time0 = args$time0,
  #       prior = args$prior,
  #       sdy = args$sdy
  #     )
  #   } else if (grepl('sd.lo', scenario, perl = F)) {
  #     
  #   } else {
  #     list()
  #   }
  #   return(conf)

  } else if (name == 'AddSeasonal') {

    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(
        name = 'AddSeasonal',
        nseasons = bsts.nseasons,
        season.duration = bsts.freq,
        sigma.prior = SdPrior(sigma.guess =sig.hi, sample.size =samp.size.prop.hi, initial.value =sig.init.val.hi, upper.limit = upper.limit)#,
        # initial.state.prior = args$initial.state.prior#,
        # sdy = args$sdy
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(
        name = 'AddSeasonal',
        nseasons = bsts.nseasons,
        season.duration = bsts.freq,
        sigma.prior = SdPrior(sigma.guess =sig.lo, sample.size =samp.size.prop.lo, initial.value =sig.init.val.lo, upper.limit = upper.limit)#,
        # initial.state.prior = args$initial.state.prior#,
        # sdy = args$sdy
      )
    } else {
      list()
    }
    return(conf)

    
  } else if (name %in% c('AddSemilocalLinearTrend','AddGeneralizedLocalLinearTrend')) {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(name = 'AddSemilocalLinearTrend',
           level.sigma.prior= SdPrior(sigma.guess=sig.hi, sample.size = samp.size.prop.hi, initial.value = sig.init.val.hi, upper.limit = upper.limit),
           slope.mean.prior = NormalPrior(mu=mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi),
           slope.ar1.prior = Ar1CoefficientPrior(mu = ar.mu.hi, sigma = sig.hi, force.stationary = ar.forc.sta, force.positive = ar.forc.pos, initial.value = ar.mu.hi),
           slope.sigma.prior = SdPrior(sigma.guess =sig.hi, sample.size = samp.size.prop.hi, initial.value = sig.init.val.hi, upper.limit = upper.limit),
           initial.level.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi),
           initial.slope.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi)
        )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name = 'AddSemilocalLinearTrend',
           level.sigma.prior= SdPrior(sigma.guess=sig.lo, sample.size = samp.size.prop.lo, initial.value = sig.init.val.lo, upper.limit = upper.limit),
           slope.mean.prior = NormalPrior(mu=mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo),
           slope.ar1.prior = Ar1CoefficientPrior(mu = ar.mu.lo, sigma = sig.lo, force.stationary = ar.forc.sta, force.positive = ar.forc.pos, initial.value = ar.mu.lo),
           slope.sigma.prior = SdPrior(sigma.guess =sig.lo, sample.size = samp.size.prop.lo, initial.value = sig.init.val.lo, upper.limit = upper.limit),
           initial.level.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo),
           initial.slope.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo)
      )
    } else {
      list()
    }
    return(conf)
    
  } else if (name == 'AddSharedLocalLevel') {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(name = 'AddSharedLocalLevel',
           nfactors = n.latent.factors.hi, 
           # coefficient.prior=SpikeSlabPrior(x = x.spikslab.prior,
           #                                  y = y.spikslab.prior,
           #                                  expected.r2 = expected.r2.hi,
           #                                  prior.df = sig.hi,
           #                                  expected.model.size = spikslab.exp.mod.size.hi,
           #                                  prior.information.weight = samp.size.prop.hi,
           #                                  diagonal.shrinkage = diag.shrinkage.hi,  ## setting=0 --> Zellner's G-prior
           #                                  # optional.coefficient.estimate = NULL,
           #                                  max.flips = -1,  ## <= 0 means all indicators will be sampled
           #                                  # mean.y = mean(y, na.rm = TRUE),
           #                                  # sdy = sd(as.numeric(y), na.rm = TRUE),
           #                                  # prior.inclusion.probabilities = NULL,
           #                                  sigma.upper.limit = upper.limit),
           initial.state.prior = MvnPrior(mean = rep(sig.init.val.hi, lags.hi), variance = matrix(rep(sig.hi, lags.hi^2),nrow=lags.hi))#, 
           # timestamps = args$timestamps, 
           # series.id = args$series.id, 
           # sdy = args$sdy
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name = 'AddSharedLocalLevel',
           nfactors = n.latent.factors.lo, 
           # coefficient.prior=SpikeSlabPrior(x = x.spikslab.prior,
           #                                  y = y.spikslab.prior,
           #                                  expected.r2 = expected.r2.lo,
           #                                  prior.df = sig.lo,
           #                                  expected.model.size = spikslab.exp.mod.size.lo,
           #                                  prior.information.weight = samp.size.prop.lo,
           #                                  diagonal.shrinkage = diag.shrinkage.lo,  ## setting=0 --> Zellner's G-prior
           #                                  # optional.coefficient.estimate = NULL,
           #                                  max.flips = -1,  ## <= 0 means all indicators will be sampled
           #                                  # mean.y = mean(y, na.rm = TRUE),
           #                                  # sdy = sd(as.numeric(y), na.rm = TRUE),
           #                                  # prior.inclusion.probabilities = NULL,
           #                                  sigma.upper.limit = upper.limit),
           initial.state.prior = MvnPrior(mean = rep(sig.init.val.lo, lags.lo), variance = matrix(rep(sig.lo, lags.lo^2),nrow=lags.lo))#, 
           # timestamps = args$timestamps, 
           # series.id = args$series.id, 
           # sdy = args$sdy
      )
    } else {
      list()
    }
    return(conf)
    
  # } else if (name == 'AddStaticIntercept') {
  #   
  #   conf <- if (grepl('sd.hi', scenario, perl = F)) {
  #     list(name = 'AddStaticIntercept',
  #          initial.state.prior = args$initial.state.prior
  #     )
  #   } else if (grepl('sd.lo', scenario, perl = F)) {
  #     
  #   } else {
  #     list()
  #   }
  #   return(conf)

    
  } else if (name == 'AddStudentLocalLinearTrend') {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(name = 'AddStudentLocalLinearTrend',
          level.sigma.prior= SdPrior(sigma.guess=sig.hi, sample.size = samp.size.prop.hi, initial.value = sig.init.val.hi, upper.limit = upper.limit),
          level.nu.prior= LognormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi), ### Inherits: DoubleModel(),
          slope.sigma.prior = SdPrior(sigma.guess =sig.hi, sample.size = samp.size.prop.hi, initial.value = sig.init.val.hi, upper.limit = upper.limit),
          slope.nu.prior= LognormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi), ### Inherits: DoubleModel(),
          initial.level.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi),
          initial.slope.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, initial.value = sig.init.val.hi),
          save.weights = TRUE
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name = 'AddStudentLocalLinearTrend',
           level.sigma.prior= SdPrior(sigma.guess=sig.lo, sample.size = samp.size.prop.lo, initial.value = sig.init.val.lo, upper.limit = upper.limit),
           level.nu.prior= LognormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo), ### Inherits: DoubleModel(),
           slope.sigma.prior = SdPrior(sigma.guess =sig.lo, sample.size = samp.size.prop.lo, initial.value = sig.init.val.lo, upper.limit = upper.limit),
           slope.nu.prior= LognormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo), ### Inherits: DoubleModel(),
           initial.level.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo),
           initial.slope.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, initial.value = sig.init.val.lo),
           save.weights = TRUE
      )
    } else {
      list()
    }
    return(conf)
    
  } else if (name == 'AddTrig') {
    
    conf <- if (grepl('sd.hi', scenario, perl = F)) {
      list(name = 'AddTrig',
           period=bsts.nseasons,
           frequencies=bsts.freq,
           sigma.prior = SdPrior(sigma.guess =sig.hi, sample.size =samp.size.prop.hi, initial.value =sig.init.val.hi, upper.limit = upper.limit),
           # initial.state.prior=NormalPrior(mu=.1, sigma=1, initial.value = .1),
           # initial.state.prior = MvnDiagonalPrior(mean.vector = c(.1, .1), sd.vector = c(1, 1)),
           method='harmonic'
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name = 'AddTrig',
           period=bsts.nseasons,
           frequencies=bsts.freq,
           sigma.prior = SdPrior(sigma.guess =sig.lo, sample.size =samp.size.prop.lo, initial.value =sig.init.val.lo, upper.limit = upper.limit),
           # initial.state.prior=NormalPrior(mu=.1, sigma=1, initial.value = .1),
           # initial.state.prior = MvnDiagonalPrior(mean.vector = c(.1, .1), sd.vector = c(1, 1)),
           method='harmonic'
      )
    } else {
      list()
    }
    return(conf)
    
  } else {
    
    return(list())
  }
  
  
} ## end func






































cat(sprintf('\n\nLoaded BSTS Helper Functions.\n\n'))
