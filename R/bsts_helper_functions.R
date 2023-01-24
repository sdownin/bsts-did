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
  } else if (state.conf$name == 'AddRegression') {
    cat('\nAddRegression called; static regression included in BSTS\n')
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
    
  } else if (name == 'AddRegression') {
    list(name = 'AddRegression') ## static regression betas included in BSTS outside of state space
  } else {
    list()
  }
  
  return(conf)
  
} ## end func







##=======================================
## Get BSTS State Space Component Template By Simulation Scenario Settings
##---------------------------------------
getStateSpaceConfBySimScenario <- function(name, scenario=NA, ## c('sd.high','sd.low')
                                           x.spikslab.prior=NULL, ## X  values for Boom::SpikeSlabPrior()
                                           y.spikslab.prior=NULL, ## Y  values for Boom::SpikeSlabPrior()
                                           ...) {
  scenario <- ifelse(is.na(scenario), 'sd.low', scenario)
  
  ## Named arguments
  args <- list(...)
  
  ## Default params
  # lags.hi <- 3
  # lags.lo <- 1
  lags.hi <- 3
  lags.lo <- 1
  ## Prior SD
  sig.hi <- .1
  sig.lo <- .01
  ## prior SD initial value
  sig.init.val.hi <- .05
  sig.init.val.lo <- .01
  ##
  init.state.prior.sig.hi <- 1
  init.state.prior.sig.lo <- 1
  ## Mean of prior distributions
  mu.hi <- .01
  mu.lo <- 0
  ## Autoregressive component
  ar.mu.hi <- .05
  ar.mu.lo <- .01
  ar.forc.sta <- FALSE
  ar.forc.pos <- FALSE
  ## Weight given to sigma.guess (interpretable as a prior observation count)
  samp.size.prop.hi <- 64
  samp.size.prop.lo <- 32
  ##
  # upper.limit.hi <- Inf
  upper.limit.hi <- Inf
  upper.limit.lo <- 1.2
  ## Spike and Slab priors on the AddAutoAr (Automatically selected Autoregressive lags)
  spikslab.init.val.hi <- .05
  spikslab.init.val.lo <- .01
  ##
  spikslab.disc.fac.hi <- .5
  spikslab.disc.fac.lo <- .25
  ##
  # spikslab.exp.mod.size.hi <- 3
  # spikslab.exp.mod.size.lo <- 2
  ##
  expected.r2.hi <- .8
  expected.r2.lo <- .5
  ## Latent factors (e.g., AddSharedLocalLevel() )
  n.latent.factors.hi <- 2
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
           sigma.prior = SdPrior(sigma.guess=sig.hi, sample.size=samp.size.prop.hi, upper.limit=upper.limit.hi, fixed = FALSE),
           # sigma.prior = NormalInverseGammaPrior(mu.guess = mu.hi, mu.guess.weight = 1, sigma.guess = sig.hi, sigma.guess.weight = 1),
           initial.state.prior = NormalPrior(mu=mu.hi, sigma=init.state.prior.sig.hi, fixed = FALSE)#,
           # sdy = args$sdy,
           # initial.y = args$initial.y
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name='AddLocalLevel',
           sigma.prior = SdPrior(sigma.guess=sig.lo, sample.size=samp.size.prop.lo, upper.limit=upper.limit.lo, fixed = FALSE),
           # sigma.prior = NormalInverseGammaPrior(mu.guess = mu.lo, mu.guess.weight = 1, sigma.guess = sig.lo, sigma.guess.weight = 1),
           initial.state.prior = NormalPrior(mu=mu.lo, sigma=init.state.prior.sig.lo, fixed = FALSE)#,
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
        sigma.prior = SdPrior(sigma.guess =sig.hi, sample.size =samp.size.prop.hi, upper.limit = upper.limit.hi)#,
        # initial.state.prior = args$initial.state.prior#,
        # sdy = args$sdy
      )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(
        name = 'AddSeasonal',
        nseasons = bsts.nseasons,
        season.duration = bsts.freq,
        sigma.prior = SdPrior(sigma.guess =sig.lo, sample.size =samp.size.prop.lo, upper.limit = upper.limit.lo)#,
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
           level.sigma.prior= SdPrior(sigma.guess=sig.hi, sample.size = samp.size.prop.hi, upper.limit = upper.limit.hi, fixed = FALSE),
           slope.mean.prior = NormalPrior(mu=mu.hi, sigma = sig.hi, fixed = FALSE),
           slope.ar1.prior = Ar1CoefficientPrior(mu = ar.mu.hi, sigma = sig.hi, force.stationary = ar.forc.sta, force.positive = ar.forc.pos),
           slope.sigma.prior = SdPrior(sigma.guess =sig.hi, sample.size = samp.size.prop.hi, upper.limit = upper.limit.hi, fixed = FALSE),
           initial.level.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, fixed = FALSE),
           initial.slope.prior = NormalPrior(mu = mu.hi, sigma = sig.hi, fixed = FALSE)
        )
    } else if (grepl('sd.lo', scenario, perl = F)) {
      list(name = 'AddSemilocalLinearTrend',
           level.sigma.prior= SdPrior(sigma.guess=sig.lo, sample.size = samp.size.prop.lo, upper.limit = upper.limit.lo, fixed = FALSE),
           slope.mean.prior = NormalPrior(mu=mu.lo, sigma = sig.lo, fixed = FALSE),
           slope.ar1.prior = Ar1CoefficientPrior(mu = ar.mu.lo, sigma = sig.lo, force.stationary = ar.forc.sta, force.positive = ar.forc.pos),
           slope.sigma.prior = SdPrior(sigma.guess =sig.lo, sample.size = samp.size.prop.lo, upper.limit = upper.limit.lo, fixed = FALSE),
           initial.level.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, fixed = FALSE),
           initial.slope.prior = NormalPrior(mu = mu.lo, sigma = sig.lo, fixed = FALSE)
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
    
  } else if (name == 'AddRegression') {
    
    conf <- list(name = 'AddRegression')
    return(conf)
    
  } else {
    
    return(list(name=NA))
  }
  
  
} ## end func









# #############################################
# ##   BSTS computation on list of simulated time series
# #############################################
# fitBstsUpdateSimlist <- function(simlist,     ## n, npds, intpd moved into simlist elements
#                                  effect.types=c('constant','quadratic','geometric'), 
#                                  sim.id=NA,
#                                  save.items.dir=NA, ## save updated simlist items to seprate RDS files
#                                  bsts.niter=1e3,
#                                  bsts.max.iter=1e5, ## 80000
#                                  plot.show=TRUE, plot.save=TRUE,
#                                  verbose=TRUE, 
#                                  save.sim.rds=TRUE
# ) {
#   ## cache original bsts.niter for dynamic niter updating if MCMC convergence failed
#   bsts.niter.orig <- bsts.niter
#   
#   # print("runSimBstsDiDComparison()::SIMLIST INPUT:")
#   # print(simlist)
#   if (length(simlist) > 0 & length(names(simlist))==0) {
#     names(simlist) <- 1:length(simlist)
#   }
#   
#   ## Simulation ID
#   if (is.na(sim.id)) {
#     sim.id <- simlist[[1]]$sim$id
#     if (is.null(sim.id) | is.na(sim.id)) {
#       sim.id <- round(10*as.numeric(Sys.time()))
#     } 
#   } 
#   
#   ## IF save simlist items is NA, then save images to work_dir
#   ## else save images to save.items.fir
#   save.img.dir <- ifelse(is.na(save.items.dir), getwd(), save.items.dir)
#   
#   ##===============================
#   ##  BSTS State Specification Comparison 
#   ##------------------------------
#   for (i in 1:length(simlist))
#   {
#     key <- names(simlist)[i]
#     key.strip <- gsub('[|]','',key,ignore.case = F, perl = T)
#     if(verbose) cat(sprintf('\n%s, %s\n',i, key))
#     
#     simlist[[key]]$cordf <- data.frame()
#     simlist[[key]]$compare <- list(bsts=list())
#     
#     ## simulation output from simulation scenario = simlist[[key]]
#     npds <- simlist[[key]]$npds
#     intpd <- simlist[[key]]$intpd
#     n <- simlist[[key]]$n
#     noise.level <- simlist[[key]]$noise.level
#     
#     ## Simulation object (containing the simulated timeseries)
#     sim <- simlist[[key]]$sim
#     # sim.id <- simlist[[key]]$sim$id
#     
#     
#     ## list of BSTS State component lists
#     bsts.state.specs <- simlist[[key]]$bsts.state.specs
#     if (length(names(bsts.state.specs))==0) {
#       names(bsts.state.specs) <- 1:length(bsts.state.specs)
#     }
#     
#     # ## BSTS expected model size for spike-and-slab priors
#     expect.mod.size <- ifelse(is.null(simlist[[key]]$expect.mod.size), NA, simlist[[key]]$expect.mod.size)
#     
#     
#     ## Dynamic Treatment Effect Type shapes
#     for (k in 1:length(effect.types)) 
#     {
#       
#       effect.type <- effect.types[k]
#       simdf <- sim$df[sim$df$effect.type == effect.type, ]
#       
#       ## Set group name 'gname' field, where 0 = control, # = period of treatment
#       simdf$match_pd <- as.numeric(simdf$match_pd)
#       simdf$gname <- 0
#       simdf$gname[simdf$group=='treatment'] <- simdf$match_pd[simdf$group=='treatment']
#       ## Remove NAs
#       simdf <- simdf[!is.na(simdf$match_id), ]
#       
#                      
#       ## Init output list wihtin simulations list
#       simlist[[key]]$compare$bsts[[effect.type]] <- list()
#       
#       
#       ##====================
#       ## BSTS Timseries Setup
#       ##--------------------
#       ## Aggregate into timeseries dataframe
#       tsdf <- simdf %>%
#         dplyr::filter( ! is.na(match_id)) %>%
#         group_by(t, t.post.intpd, effect.type, match_pd, gname, group, group.color) %>%
#         dplyr::summarize(
#           n_in_pd = n(),
#           actors = paste(unique(actor), collapse = '|'),
#           y_outcome = mean(y, na.rm=T),
#           y_sum = sum(y, na.rm=T),
#           y_sd = sd(y, na.rm=T),
#           y_min = min(y, na.rm=T),
#           y_max = max(y, na.rm=T),
#           y_skew = skewness(y, na.rm = T, type = 2), ## moment-based distribution
#           y_kurt = kurtosis(y, na.rm = T, type = 2), ## moment-based distribution
#           ##
#           x1_mean = mean(x1, na.rm=T),
#           x2_mean = mean(x2, na.rm=T),
#           x3_mean = mean(x3, na.rm=T),
#           ##
#           c1_mean = mean(c1, na.rm=T),
#           c2_mean = mean(c2, na.rm=T),
#           c3_mean = mean(c3, na.rm=T),
#           #
#           c1_sd = sd(c1, na.rm=T),
#           c2_sd = sd(c2, na.rm=T),
#           c3_sd = sd(c3, na.rm=T),
#           ##
#           c1_skew = skewness(c1, na.rm=T, type = 2),
#           c2_skew = skewness(c2, na.rm=T, type = 2),
#           c3_skew = skewness(c3, na.rm=T, type = 2),
#           #
#           c1_kurt = skewness(c1, na.rm=T, type = 2),
#           c2_kurt = skewness(c2, na.rm=T, type = 2),
#           c3_kurt = skewness(c3, na.rm=T, type = 2),
#           #
#           b1_mean = mean(b1, na.rm=T),
#           b2_mean = mean(b2, na.rm=T),
#           b3_mean = mean(b3, na.rm=T),
#           #
#           u_mean = mean(u, na.rm=T),
#           v_mean = mean(v, na.rm=T)
#         )
#       tsdf$.id <- 1:nrow(tsdf)
#       
#       ## MAKE WIDE TIMESERIES FOR treatment,control groups in n periods
#       val.cols <- c('y_outcome','y_sum','y_min','y_max','y_sd',
#                     'y_skew','y_kurt',
#                     'x1_mean','x2_mean','x3_mean',
#                     'c1_mean','c2_mean','c3_mean',
#                     'c1_sd','c2_sd','c3_sd',
#                     'c1_skew','c2_skew','c3_skew',
#                     'c1_kurt','c2_kurt','c3_kurt',
#                     'b1_mean','b2_mean','b3_mean',
#                     'u_mean','v_mean')
#       ts <- unique(tsdf$t)
#       groups <- unique(tsdf$group)
#       ## init timeseries dataframe - wide
#       tsdfw <- data.frame(t=ts,  stringsAsFactors = F)
#       for (jj in 1:length(groups)) {
#         id.j <- which( tsdf$group == groups[ jj ] ) 
#         for (kk in 1:length(val.cols)) {
#           df.col <- data.frame( tsdf[ id.j , val.cols[ kk ] ] )
#           names(df.col) <- sprintf('%s_%s',groups[ jj ],val.cols[ kk ])
#           tsdfw <- cbind(tsdfw,  df.col)
#         }
#       }
#       
#       
#       # Set up pre- and post-treatment period
#       # pre.period <- as.Date(c("2013-01-01","2016-01-25"))
#       pre.period <- c(1, intpd-1)  
#       # post.period <- as.Date(c("2016-01-26","2018-01-01"))
#       post.period <- c(intpd, npds) 
#       
#       # # BSTS causal effect analysis using CausalImpact package
#       # # CausalImpact option: 
#       # # niter: Number of MCMC samples to draw. More samples lead to more accurate inferences. Defaults to 1000.
#       # # standardize.data: Whether to standardize all columns of the data before fitting the model (i.e., setting Bayes priors), Defaults to TRUE.
#       # # prior.level.sd: Prior standard deviation of the Gaussian random walk of the local level. Defaults to 0.01.
#       # # nseasons: Period of the seasonal components. Default to 1.
#       # # dynamic.regression: Whether to include time-varying regression coefficients. Defaults to FALSE.
#       # impact_amount <- CausalImpact(amount.impact,pre.period,post.period,alpha=0.1, model.args = list(niter = 5000))
#       # summary(impact_amount)
#       # plot(impact_amount)
#       dat <- tsdfw[,c('treatment_y_outcome',
#                       'control_y_outcome',
#                       # 'control_y_sd',#'control_y_sum', 'control_y_min','control_y_max',
#                       # 'control_y_max','control_y_skew', 'control_y_kurt',
#                       'control_c1_mean','control_c2_mean',  'control_c3_mean',
#                       'control_c1_sd','control_c2_sd','control_c3_sd',
#                       'control_c1_skew','control_c2_skew','control_c3_skew'#,
#                       # 'control_c1_kurt','control_c2_kurt',
#                       # 'control_c3_kurt'#,
#                       # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
#                       # 'control_u_mean','control_v_mean'
#       )]
#       ## Train on y pre-treatment but NA's post-treatment
#       y.pre.treat.NAs.post.treat <- c(dat$treatment_y_outcome[1:(intpd-1)], rep(NA,npds-intpd+1))
#       ## Then use the post-treatment response for causal impact estimation
#       post.period.response <- dat$treatment_y_outcome[intpd:npds]
#       ## Covariates (predictors) - Matrix for "formula = y ~ predictors" argument
#       predictors <- dat[, ! names(dat) %in% 'treatment_y_outcome'] ## remove response; convert to matrix
#       # ## Covariates (predictors) - Dataframe for "data" argument
#       # predictors <- as.matrix(predictors) 
#       
#       ## ADD temporal trend to covariates
#       predictors$covar_temp_trend <- (1:npds) + rnorm(npds, 0, noise.level)   ## * simlist$`0.rand.base`$b5
#       
#       ## ***DEBUG***
#       # print("bsts.state.specs:")
#       # print(bsts.state.specs)
#       ##
#       
#       ##----------------------------
#       ## State Space Configuration
#       ##----------------------------
#       ## LOOP OVER BSTS STATE SPECIFICATIONS FOR THE SAME SIMULATION RUN
#       for (h in 1:length(bsts.state.specs)) 
#       {
#         simlist[[key]]$compare$bsts[[effect.type]][[ h ]] <- list()
#         
#         ## h'th BSTS state space configuration (state component list)
#         state.conf <- bsts.state.specs[[ h ]]
#         
#         ## names of 
#         state.comps <- unname(sapply(state.conf, function(x){
#           ifelse(class(x)=='list' & 'name' %in% names(x), x$name, 'ERROR x$name not set in state.conf[x]')
#         }, simplify = T))
#       
#         ## Loop over state space components
#         nss <- length(state.conf)
#         ## State Space Config list
#         st.sp <- list()
#         if (nss > 0) {
#           for (jj in 1:nss) {
#             state.conf.item <- state.conf[[ jj ]]
#             if (class(state.conf.item)=='list') {
#               
#               if(state.conf.item$name=='AddSharedLocalLevel') {
#                 state.conf$coefficient.prior <- SpikeSlabPrior(
#                   x = predictors,
#                   y = dat$treatment_y_mean, ##**NOTE** USING ALL Y VALS (not just y.pre.treat.NAs.post.treat)
#                   expected.r2 = .5, ## [.5]
#                   prior.df = .01, ##[.01]
#                   expected.model.size = 2, ## [1]
#                   prior.information.weight = .01, ## [.01]
#                   diagonal.shrinkage = 0.5,  ## [.5] setting=0 --> Zellner's G-prior
#                   # optional.coefficient.estimate = NULL,
#                   max.flips = -1,  ## <= 0 means all indicators will be sampled
#                   # prior.inclusion.probabilities = NULL,
#                   sigma.upper.limit = Inf
#                 )
# 
#               }
#               ## [[ IF NOT AddSharedLocalLevel(), NOT INCLUDING SPIKESLABPRIOR ]]
#               st.sp <- updateStateSpaceAddComponentFromConfList(st.sp,  
#                                                                 y.pre.treat.NAs.post.treat,  
#                                                                 state.conf.item)
#               if(verbose) cat(sprintf('add to state.space: %s\n',state.conf.item$name))
#             }
#           }
#         } else {
#           ## Default in CausalImpact package
#           st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
#         }
#         # print(st.sp)
#         
#         if(verbose) cat(sprintf('\nRunning BSTS model estimation for state.conf h=%s\n',h))
#         
#         ##-------------------------------
#         ## Regression component of model
#         if ('AddRegression' %in% state.comps) {
#           bsts.input.form <- y.pre.treat.NAs.post.treat ~ . ## with regression
#         } else {
#           bsts.input.form <- y.pre.treat.NAs.post.treat  ## without regression
#         }
#         
#         ##---------------------------------------------------------
#         ## RUN BSTS WITH DYNAMIC niter BASED ON CONVERGENCE 
#         isConverged <- FALSE
#         isMaxIter <- FALSE
#         hasBstsError <- FALSE
#         bsts.niter <- bsts.niter.orig  ## reset to original bsts.niter input value
#         while ( !isConverged  &  !isMaxIter & !hasBstsError  ) {
#           
#           # ##**DEBUG**
#           # print('st.sp')
#           # print(st.sp)
#           # ##
#           
#           ## BSTS model
#           bsts.model <- tryCatch(expr = {
#             bsts(formula = bsts.input.form,
#                  state.specification = st.sp,
#                  data = predictors,
#                  expected.model.size = expect.mod.size,
#                  niter = bsts.niter,
#                  ping = ifelse(verbose, round(bsts.niter/10), 0))
#           },
#           error=function(e) {
#             message(sprintf('bsts() error: %s', as.character(e)))
#           },
#           warning=function(w) {
#             message(sprintf('bsts() warning: %s', as.character(w)))
#           },
#           finally={
#             ##PASS
#           })
#           
#           
#           ## skip if bsts() function threw an error (not return 'bsts' object)
#           if ( class(bsts.model) == 'bsts' ) {
#             ## Use BSTS prediction of counterfactual to estimate CausalImpact
#             impact_amount <- CausalImpact(bsts.model=bsts.model,
#                                           post.period.response = post.period.response,
#                                           alpha=0.05, model.args = list(niter = bsts.niter))
#             ## POSTERIOR PREDICTIVE CHECKS
#             ppcheck.filename <- file.path(save.img.dir,
#                                           sprintf('%s_bsts_post_pred_checks_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
#                                                   prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id))
#             convcheck <- bstsPostPredChecks(bsts.model, filename=ppcheck.filename, return.val = T)
#             ##
#             # print(convcheck)
#             ##
#             if(verbose) cat(sprintf('\nBSTS niter = %s',bsts.niter))
#             if(verbose) cat(convcheck$summary)
#             ## UPDATE CONVERGENCE CHECK FLAG - ELSE RERUN WITH INCREASED bsts.niter
#             # isConverged <- convcheck$converged.all
#             conv.tol <- 0.8
#             conv.min.iter.thresh <- 4e4 ## 40k
#             # isConverged <- convcheck$converged.prop >= conv.tol
#             isConverged <- convcheck$converged.all | (convcheck$converged.prop >= conv.tol & bsts.niter >= conv.min.iter.thresh) ## don't allow incomplete check below minimum threshold of bsts.niter = 10k 
#             if(verbose) print(convcheck$converged)
#             if(verbose) cat(sprintf('Converged proportion = %.3f (tol = %.3f) (min.iter.converg.thresh=%s)\nConverged status = %s\n\n',
#                         convcheck$converged.prop, conv.tol, conv.min.iter.thresh,  isConverged))
#           } else {
#             hasBstsError <- TRUE
#           }
#           
#           
#           ## UPDATE niter --> increase by factor of 2
#           if ( !isConverged ) {
#             # .niter <- round( bsts.niter * (2 - convcheck$converged.prop) ) ## scale niter increment by proportion of failed checks
#             bsts.niter.new <- round( (bsts.niter * 2) / 10) * 10  ## double niter and make divisible by 10 (for bsts ping=10 to divide evenly into niter)
#             if ( bsts.niter.new <= bsts.max.iter ) {
#               bsts.niter <- bsts.niter.new  
#             } else {
#               isMaxIter <- TRUE
#             }
#           }
#           
#           
#         }
#         
#         if (hasBstsError) {
#           next
#         }
#         ##-------------------------------------------------------------
#         
#         if (plot.show) {
#           ##
#           plot(bsts.model, main=sprintf('BSTS Plot: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
#           PlotBstsComponents(bsts.model) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
#           if (bsts.model$has.regression) {
#             PlotBstsSize(bsts.model, main=sprintf('BSTS Size: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
#           }
#           ##
#         }
#         if (plot.save) {
#           #                         key,effect.type,sim.id))
#           p.bsts.impact.all <- plot(impact_amount, c('original','pointwise','cumulative')) # pointwise','cumulative
#           ggsave(filename = file.path(save.img.dir,
#                                       sprintf('%s_bsts_CausalImpact_plot_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
#                                               prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id))
#           )
#         }
# 
#         
#         ##-------------------
#         ## OVERALL CONFIDENCE INTERVALS AND PVALUES
#         ## OVERALL CAUSAL P-VALUE
#         pval.bsts.general <- impact_amount$summary$p[2]
#         # pval.did.general <- ccattgt
#         ci.bsts.general <- c(impact_amount$summary$AbsEffect.lower[1], impact_amount$summary$AbsEffect.upper[1]) 
#         ## MUST USE GROUP AGGREGATION TO GET critical.val.egt
#         ## (we only have one group [i.e., one treatment time])
#         # ci.did.general <- c(agg.group$overall.att - (agg.group$overall.se * agg.group$crit.val.egt),
#         #                     agg.group$overall.att + (agg.group$overall.se * agg.group$crit.val.egt))
#         ##-------------------
#         
#         # ## DID
#         # did.res <- tidy(agg.es)
#         ## BSTS
#         bsts.res <- impact_amount$series
#         # .nidx <- which(names(bsts.res) %in% c('point.effect','point.effect.lower','point.effect.upper'))
#         bsts.res <- bsts.res[ , which(names(bsts.res) %in% c('point.effect','point.effect.lower','point.effect.upper')) ]
#         names(bsts.res) <- c('bsts.point.effect','bsts.point.effect.lower','bsts.point.effect.upper')
#         
#         # plot(did.res)
#         
#         # ## AVERAGE TREATMENT EFFECT USED IN SIMULATION
#         # simdf
#         
#         tr.actors <- unique(simdf$actor[which(simdf$group=='treatment')])  
#         co.actors <- unique(simdf$actor[which(simdf$group=='control')])
#       
#         ##
#         # att.b3 <- mean(b3diff$diff[intpd:npds])
#         # att.did <- agg.es$overall.att
#         att.bsts <- impact_amount$summary$AbsEffect[1]
#         # 
#         
#         ##===============================================================
# 
#         ##
#         simlist[[key]]$compare$bsts[[effect.type]][[ h ]]$CausalImpact <- impact_amount
#         simlist[[key]]$compare$att.bsts <- att.bsts
# 
#       } ## // end h loop over bsts.state components
#       
#       
#     } ## // end k loop over effect types
#     
#     
#     if ( save.sim.rds ) {
#       save.items.dir <- ifelse(is.na(save.items.dir), getwd(), save.items.dir)
#       maxGb <- 6
#       # maxGb <- .001  ## **DEBUG**
#       if ( (object.size(simlist[[key]])/1e9) <= maxGb) {
#         ## IF SMALL ENOUGH, Save simulation list as serialized data file
#         simlist.file <- sprintf('__fitBstsUpdateSimlist__n%s_pd%s_niter%s_%s_%s.rds', n, npds, bsts.niter, sim.id, key.strip)
#         save.file.path <-  file.path(save.items.dir, simlist.file)
#         saveRDS(simlist[[key]], file = save.file.path)
#         ##
#         # ## FREE UP MEMORY
#         # simlist[[key]] <- list(file = save.file.path)
#       } else { 
#         ## IF TOO LARGE, then save BSTS object (with MCMC samples) separately from rest of simulation
#         save.file.paths <- c()
#         ## Save simulation list as serialized data file
#         ## 1.  BSTS
#         simlist.file <- sprintf('__fitBstsUpdateSimlist__n%s_pd%s_niter%s_%s_%s_simlist-compare-bsts.rds', n, npds, bsts.niter, sim.id, key.strip)
#         save.file.path <-  file.path(save.items.dir, simlist.file)
#         saveRDS(simlist[[key]]$compare$bsts, file = save.file.path)
#         simlist[[key]]$compare$bsts <- NULL ## save space
#         save.file.paths[1] <- save.file.path
#         ## 2. rest of simulation object
#         simlist.file <- sprintf('__fitBstsUpdateSimlist__n%s_pd%s_niter%s_%s_%s_simlist-MAIN.rds', n, npds, bsts.niter, sim.id, key.strip)
#         save.file.path <-  file.path(save.items.dir, simlist.file)
#         saveRDS(simlist[[key]], file = save.file.path)
#         save.file.paths[2] <- save.file.path
#         ##
#         # ## FREE UP MEMORY
#         # simlist[[key]] <- list(file = save.file.paths)
#       }
#     } 
#     
#     
#   } ## // end simlist loop i   ##  #; dev.off()
#   
#   return(simlist)
#   
# }





### 
## Plot BSTS State Components - all series compared against observed values
###
plotBstsStateComps <- function(bsts.model, intpd=NA, filename=NA, save.plot=FALSE) {
  # save.plot <- !is.na(filename)
  
  npds <- length(bsts.model$original.series)
  
  if (is.na(intpd)) {
    intpd <- npds
  }
  
  
  niter <- bsts.model$niter
  burn <- round( niter * .2 )
  
  idx.noburn <- (burn+1):niter
  idx.pre    <- 1:(intpd-1)
  
  err.pre.1step <- bsts.model$one.step.prediction.errors[idx.noburn, idx.pre]
  y.pre <- bsts.model$original.series[idx.pre]
  
  ## MEAN ABSOLUTE ERROR (MAE)
  mae <- mean(colMeans(abs(err.pre.1step), na.rm = T), na.rm=T)
  
  
  ## Reverse Extract Prediction (per MCMC iteration) from observed and residual
  ##  @see p.69 in  https://cran.r-project.org/web/packages/bsts/bsts.pdf
  ##              err = y[t] - y_hat[t](y[t-1]) 
  ## y_hat[t](y[t-1]) = y[t] - err
  ## Create replicas of outcome series to subtract each MCMC iteration
  y.pre.rep <- matrix(y.pre, length(y.pre), (niter - burn),  byrow = FALSE)
  ## Get prediction for each MCMC iteration:
  ## t([npds, niter-burn]) -   [niter-burn, npds]
  ## colMeans() returns periodwise means (over MCMC iterations)
  y.pred <-  colMeans( t(y.pre.rep) - err.pre.1step )
  
  # 
  # y.pred.interval <- as.matrix( rbind(q) )
  
  # y.pred <- (cbind(y.pre) - t(err.pre.1step))
  
  
  # ##---------- CHECK REFITTING ON TRUNCATED DATA BEFORE INTERVENTION -----------
  # y.pre <- bsts.model$original.series[idx.pre]
  # bsts.model0 <- bsts(y.pre ~ . , 
  #                    state.specification = bsts.model$state.specification,
  #                    data = as.data.frame(bsts.model$predictors[idx.pre, -1]), ## drop (Intercept) first column 
  #                    timestamps = -1*rev(idx.pre),
  #                    niter = niter
  #                    )
  # bsts.model <- bsts.model0
  # ##
  # response <- bsts.model$original.series
  # ##
  # pred <- if (bsts.model$has.regression) {
  #   newdata <- bsts.model$predictors
  #   newdata <- cbind(response=response, newdata)
  #   predict.bsts(bsts.model, newdata = newdata , burn = burn,
  #                timestamps = idx.pre
  #   ) ## already knows horizon from 
  # } else {
  #   predict.bsts(bsts.model, burn = burn, 
  #                horizon = intpd-1,
  #                olddata = response[1:(intpd-1)] )
  # }
  # # ## DEBUG
  # plot(pred, style='dynamic')
  # # matplot(cbind(t(pred$interval),pred$mean))
  # # ##
  # ##------------------------------------------------------------------------
  
  # # response <- bsts.model$original.series
  # 
  # pred <- if (bsts.model$has.regression) {
  #   newdata <- bsts.model$predictors[idx.pre, ]
  #   newdata <- cbind(response=response[idx.pre], newdata)
  #   predict.bsts(bsts.model, newdata = newdata , burn = burn,
  #                timestamps = idx.pre
  #                ) ## already knows horizon from 
  # } else {
  #   predict.bsts(bsts.model, burn = burn, 
  #                horizon = intpd-1,
  #                olddata = response[1:(intpd-1)] )
  # }
  # # ## DEBUG
  # plot(pred, style='boxplot')
  # # matplot(cbind(t(pred$interval),pred$mean))
  # # ##
  
  sc <- bsts.model$state.contributions
  
  # pred.mean <- colMeans(pred$distribution)
  pred.mean <- y.pred
  
  # bsts.model$final.state
  # dim(sc)
  # par(mfrow=c(2,2))
  components <- dimnames(sc)$component
  ncomps <- length(components)
  y.orig <- as.numeric( bsts.model$original.series[1:(intpd-1)] )
  sc.means.all <-  c(sapply(1:ncomps,function(i)colMeans(sc[,i,])))
  # sc.mins <-  c(sapply(1:ncomps,function(i)min(c(sc[,i,]),na.rm = T)))
  # sc.maxs <-  c(sapply(1:ncomps,function(i)max(c(sc[,i,]),na.rm = T)))
  .vals <- c(y.orig, sc.means.all, pred.mean)  #sc.mins, sc.maxs
  .vals <- .vals[which(!is.null(.vals) & !is.nan(.vals) & !is.na(.vals))]
  .ylim <- range(.vals)
  # .ylim <- c( min(.vals) - .25*diff(range(.vals)),  max(.vals) )
  
  ## Get Mean Absolute Error (MAE) for plot title
  ## y.rep dimensions:  rows [npds before intervention] x cols [draws (niter - burn)]
  # y.rep <- matrix(y.orig, length(y.orig), (niter - burn),  byrow = FALSE)
  # # res <- y.rep - (post.pred.dist - CausalMBSTS$mcmc$eps.samples[, i, ])
  # res <-  y.rep - t(pred$distribution)
  # err.pds <- rowMeans(res)  ## mean over MCMC draws per period
  # mae <- mean( abs(err.pds), na.rm=T) ## mean of absolute per-period errors
  
  if(save.plot) {
    png(filename = filename, width = 10, height = 8, units = 'in', res = 400)
  }
  ##----------- INSIDE PNG PLOT --------------------------------------
  # nrow <- ceiling( sqrt(ncomps) )
  # ncol <- ceiling(ncomps / nrow)
  par(mar=c(2.5,2.5,2.5,1))  ##mfrow=c(nrow,ncol), 
  title <- sprintf('%s (MAE = %.3f)', paste(components,collapse = ' + '), mae)
  plot(x=1:(intpd-1), y.orig, pch=16,
       ylab='Y', xlab='t', ylim=.ylim, main=title)  ## ylim=.ylim
  lines(pred.mean, col='blue', lwd=1.9, ylim=.ylim) ## ylim=.ylim
  for (i in 1:dim(sc)[2]) {
    # plot(colMeans(sc[,i,]), type='l', main=component[i])
    lines(colMeans(sc[,i,]), type='l', col=i, lty=i, lwd=1.5, ylim=.ylim)
  }
  legend('topleft', legend=c('observed', 'predicted', components), 
         lty=c(NA, 1, 1:ncomps), 
         pch=c(16, NA, rep(NA,ncomps)),
         col=c('black', 'blue', 1:ncomps),
         lwd=c(NA, 1.9, rep(1.5,ncomps)))
  ##-------------- END PNG PLOT ----------------------------------
  if (save.plot) {
    dev.off()
  }

}



##
##  Get the  CausalImpact object from the BSTS vs. DiD comparison simlist object
##
getCausalImpactFromSimlist <- function(simlist, key=NA, 
                                       effect.type='quadratic',
                                       state.space.list.id=1
                                       ) {
  if (is.na(key)) {
    if (length(names(simlist))==0) {
      names(simlist) <- as.character( 1:length(simlist) )
    }
    key <- names(simlist)[1]
  }
  bstslist <- simlist[[key]]$compare$bsts[[ effect.type ]]
  return( bstslist[[ state.space.list.id ]]$CausalImpact )
}


##
##  Get the  CausalImpact object from the BSTS vs. DiD comparison simlist object
##
getBstsModelFromSimlist <- function(simlist, key=NA, 
                                     effect.type='quadratic',
                                     state.space.list.id=1
                                    ) {
  
  causimp <- getCausalImpactFromSimlist(simlist, key, 
                                        effect.type, 
                                        state.space.list.id)
  return( causimp$model$bsts.model )
}





##
#
##
heidel.diag.mod <- function (x, eps = 0.1, pvalue=0.05) 
{
  if (is.mcmc.list(x))
    return(lapply(x, heidel.diag.mod, eps))
  x <- as.mcmc(as.matrix(x))
  HW.mat0 <- matrix(0, ncol = 7, nrow = nvar(x))
  dimnames(HW.mat0) <- list(varnames(x),
                            c("stest", "start", "CMV.stat", "pvalue", "htest",
                              "mean", "halfwidth"))
  HW.mat <- HW.mat0
  for (j in 1:nvar(x)) {
    start.vec <- seq(from=start(x), to = end(x)/2, by=niter(x)/10)
    Y <- x[, j, drop = TRUE]    
    n1 <- length(Y)
    ## Schruben's test for convergence, applied sequentially
    ##
    S0 <- spectrum0.ar(window(Y, start=end(Y)/2))$spec
    converged <- FALSE
    for (i in seq(along = start.vec)) {
      Y <- window(Y, start = start.vec[i])
      n <- niter(Y)
      ybar <- mean(Y)
      B <- cumsum(Y) - ybar * (1:n)
      Bsq <- (B * B)/(n * S0)
      I <- sum(Bsq)/n
      if(converged <- !is.na(I) && pcramer(I) < 1 - pvalue)
        break
    }
    ## Recalculate S0 using section of chain that passed convergence test
    S0ci <- spectrum0.ar(Y)$spec
    halfwidth <- 1.96 * sqrt(S0ci/n)
    passed.hw <- !is.na(halfwidth) & (abs(halfwidth/ybar) <= eps)
    if (!converged || is.na(I) || is.na(halfwidth)) {
      nstart <- NA
      passed.hw <- NA
      halfwidth <- NA
      ybar <- NA
    }
    else {
      nstart <- start(Y)
    }
    HW.mat[j, ] <- c(converged, nstart, I, 1 - pcramer(I), 
                     passed.hw, ybar, halfwidth)
  }
  # class(HW.mat) <- "heidel.diag"
  return(HW.mat)
}

###
## POSTERIOR PREDICTIVE CHECKS
###
bstsPostPredChecks <- function(bsts.model, filename=NA, 
                               save.plot=TRUE, return.val=FALSE,
                               burn=NA, conv.alpha=0.05) {
  
  ppcheck.filename <- if (is.na(filename)){
    sprintf('bsts_post_pred_checks_%s.png', round(10*as.numeric(Sys.time())) )
  }  else {
    filename
  }
  
  ## output list of convergence checks (booleans) and residual dataframes
  checklist <- list()
  
  response <-  bsts.model$original.series  
  
  idx.pre <- which( ! is.na(response) )
  idx.post <- which( is.na(response) )
  
  y <- response
  
  # ##**DEBUG**
  # print('as.numeric( response )')
  # print(response)
  # ##**
  
  npds <- length(y)
  
  intpd <- idx.post[1]  ## this is (intpd + 1);    [[should we use (idx.post[1] - 1) to get intpd ?]]
  
  niter <- length(bsts.model$sigma.obs)
  
  if (is.na(burn)) {
    burn <- round( niter * .2 )
  }
  
  idx.noburn <- (burn+1):niter
  
  ## before intervention period bool dummy
  .ind <- (1:npds) %in% idx.pre
  # .ind <- (1:npds) < intpd
  
  ## observed data to predict via BSTS
  hasRegression <- bsts.model$has.regression
  if (hasRegression) {
    newdata <- bsts.model$predictors[idx.pre, ]
    newdata <- cbind(response=response[idx.pre], newdata)
  } else {
    newdata <- NULL
  }
  
  # ##**DEBUG**
  # print('newdata[1:10,]')
  # print(newdata[1:10,])
  # ##**
  # ##**DEBUG**
  # plot(post.pred$mean, col='red',type='l',ylim=c(-2,6)); points(response[1:(intpd-1)], col='black', pch=16, main='DEBUG')
  # print('post.pred:')
  # print(post.pred)
  # ##
  
  ##-- Prediction distribution ----------------------------------------
  ##  in sample forecasts for plotting predicted-observed distributions
  pred.in.samp <- if (hasRegression) {
    predict.bsts(bsts.model, newdata = newdata , burn = burn) ## already knows horizon from
  } else {
    predict.bsts(bsts.model, burn = burn,
                 horizon = intpd-1,
                 olddata = response[idx.pre])
  }
  pred.in.samp.dist <- pred.in.samp$distribution  ## [niter - burn,  1:(intpd-1) ]
  pred.in.samp.mean <- colMeans(pred.in.samp.dist, na.rm = T)
  ##--------------------------------------------------------------------
  
  ##-- Error distributions ------------------------------------
  ##   1-step ahead prediction error 
  # ## matrix [niter - burn,  1:(intpd-1) ]
  # post.pred.dist <- post.pred$distribution
  # ## vector of length(1:(intpd-1))
  # post.pred.mean <- post.pred$mean ##apply(post.pred$distribution, c(1), mean)
  #####
  # standardized errors of 1-step ahead prediction (pre-intervention window)
  err.1step.dist <- bsts.model$one.step.prediction.errors[idx.noburn, idx.pre]
  ##### ^ this is SAME as:
  # # err.1step.std.all <- bsts.prediction.errors(bsts.model, burn = burn, standardize = T)
  # # err.1step.dist <- err.1step.std.all$in.sample[ , idx.pre]
  #####
  ## Periodwise means of 1-step ahead errors from MCMC iterations 
  err.1step.mean <- colMeans(err.1step.dist, na.rm=T)
  ##----------------------------------------------------------
  
  
  ##========================= PLOTTING =========================================
  if (save.plot) {
    png(filename = ppcheck.filename, width = 15, height = 10, units = 'in', res = 400)
  }
  ##----------- INSIDE PNG PLOT --------------------------------------
  par(mfrow=c(2,3), mar=c(2.5,2.5,2.5,1))
  
  ##===================
  ## Posterior Predictive (Y) plots 
  ##
  ##   - use predictions from in-sample window  
  ##      `pred.in.samp.dist`, `pred.in.samp.mean`
  ##
  ##-------------------
  
  ##-----------
  ## Trace plot of Posterior Predictive distribution Markov Chain 
  post.pred.tr <- rowMeans(pred.in.samp.dist)
  ##
  gd <- geweke.diag(post.pred.tr, .1, .5)
  gd.z <- abs(gd$z) ## z statistic
  gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
  gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
  ##
  hd <- heidel.diag.mod(post.pred.tr)
  hd.st.cmv <- hd[1,'CMV.stat'][1]
  hd.st.p <- hd[1,'pvalue'][1]
  hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
  hd.st.start <- hd[1,'start'][1]
  hd.hw.eps <- 0.1 ## default 0.1
  hd.hw <- hd[1,'halfwidth'][1]
  hd.hw.mean <- hd[1,'mean'][1]
  hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
  ##
  rng <- range(post.pred.tr)
  ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
  plot(post.pred.tr, type='l', main='A. Posterior Predicted MCMC Trace, Y' ,
       ylim=ylims
  )
  mtext.postpred <- sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
                            gd.result,gd.z,gd.p,
                            hd.st.result, hd.st.cmv, hd.st.p,
                            hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps)
  mtext(text = mtext.postpred, side = 3, line=-4.5, outer = F)
  ##
  checklist$ck.postpred <- list(
    geweke = list(check=(gd.p >= conv.alpha), z=gd.z, p=gd.p),
    hw.st  = list(check=(hd.st.p >= conv.alpha), cmv=hd.st.cmv, p=hd.st.p),
    hw.hw  = list(check=(abs(hd.hw/hd.hw.mean) < hd.hw.eps), hw.mean=abs(hd.hw/hd.hw.mean), hw.eps=hd.hw.eps)
  )
  ##-----------
  
  ##-----------
  ## DENSITIES
  plot(density(y[.ind]),  main = "B. Density comparison, Y")
  lines(density(pred.in.samp.mean), lwd=2, lty=2, col='red')
  legend('topleft', legend=c('observed','predicted'), lty=c(1,2), col=c('black','red'))
  ##-----------
  
  ##-----------
  ## HISTOGRAMS & BAYESIAN P-VALUES
  # max.distrib <- apply(post.pred, c(2, 3), max)
  max.distrib <- apply(pred.in.samp.dist, 1, max)
  pvalue <- sum(max.distrib >= max(y[.ind]))/length(max.distrib)
  hist(max.distrib, 30, col = "lightblue", border = "grey", 
       main = paste0("C. Bayesian p-val (Max Y) = ", round(pvalue, 2)),
       xlab = "Max. in-sample forecasts")
  abline(v = max(y[.ind]), col = "darkblue", lwd = 3)
  ##-----------
  
  ##===================
  ## Std. Residual plots
  ##
  ##   - use errors from 1-step ahead predictions (pre-intervention window)
  ##      `err.1step.dist`, `err.1step.mean`
  ##
  ##-------------------
  # y.rep <- matrix(y[.ind], length(y[.ind]), (niter - burn),  byrow = FALSE)
  # # res <- y.rep - (post.pred.dist - CausalMBSTS$mcmc$eps.samples[, i, ])
  # res <-  y.rep - t(post.pred.dist)
  # std.res <- res / apply(res,1,sd)   ## [residual i] / [ stdev of residual i]
  # # std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
  err.1step.dist <- err.1step.dist
  
  
  ##-----------
  ## Trace plot of Posterior Predictive distribution Markov Chain 
  # res.tr <- colMeans(std.res)
  res.tr <- rowMeans(err.1step.dist)
  ##
  gd <- geweke.diag(res.tr, .1, .5)
  gd.z <- abs(gd$z) ## z statistic
  gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
  gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
  ##
  hd <- heidel.diag.mod(res.tr)
  hd.st.cmv <- hd[1,'CMV.stat'][1]
  hd.st.p <- hd[1,'pvalue'][1]
  hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
  hd.st.start <- hd[1,'start'][1]
  hd.hw.eps <- 0.1 ## default 0.1
  hd.hw <- hd[1,'halfwidth'][1]
  hd.hw.mean <- hd[1,'mean'][1]
  hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
  ##
  rng <- range(res.tr)
  ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
  plot(res.tr, type='l', main='D. Std.Residual MCMC Trace, Y',
       ylim=ylims)
  mtext.residual <- sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
                            gd.result,gd.z,gd.p,
                            hd.st.result, hd.st.cmv, hd.st.p,
                            hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps)
  mtext(text = mtext.residual, side = 3, line=-4.5, outer = F)
  ##
  checklist$ck.residual <- list(
    geweke = list(check=(gd.p >= conv.alpha), z=gd.z, p=gd.p),
    hw.st  = list(check=(hd.st.p >= conv.alpha), cmv=hd.st.cmv, p=hd.st.p),
    hw.hw  = list(check=(abs(hd.hw/hd.hw.mean) < hd.hw.eps), hw.mean=abs(hd.hw/hd.hw.mean), hw.eps=hd.hw.eps)
  )
  ##-----------
  
  ##-----------
  # Std. Residual Normality plot
  qqnorm(rowMeans(err.1step.dist), main = "E. Std.Residual QQ-plot, Y")
  qqline(rowMeans(err.1step.dist))
  ##-----------
  
  ##-----------
  ## Std Residual ACF
  Acf(rowMeans(err.1step.dist), main = "");title(main='F. Std.Residual ACF, Y')
  ##-----------
  
  ##----------- end PNG PLOT --------------------------------------
  if (save.plot) {
    dev.off()
  }
  
  ##
  if(return.val) {
    checklist$err.1step.dist <- err.1step.dist
    checklist$pred.in.samp.dist <- pred.in.samp.dist
    checklist$summary <- sprintf('\nPosterior Predictive:-----\n%s\nStd. Residuals:-----\n%s\n\n', mtext.postpred, mtext.residual)
    
    ## ALL CONVERGENCE CHECKS
    c1 <- unname( checklist$ck.postpred$geweke$check )
    c2 <- unname( checklist$ck.postpred$hw.st$check )
    c3 <- unname( checklist$ck.postpred$hw.hw$check )
    c4 <- unname( checklist$ck.residual$geweke$check )
    c5 <- unname( checklist$ck.residual$hw.st$check )
    c6 <- unname( checklist$ck.residual$hw.hw$check )
    ##
    ck1 <- ifelse(is.na(c1) | is.nan(c1), FALSE, c1)
    ck2 <- ifelse(is.na(c2) | is.nan(c2), FALSE, c2)
    ck3 <- ifelse(is.na(c3) | is.nan(c3), FALSE, c3)
    ck4 <- ifelse(is.na(c4) | is.nan(c4), FALSE, c4)
    ck5 <- ifelse(is.na(c5) | is.nan(c5), FALSE, c5)
    ck6 <- ifelse(is.na(c6) | is.nan(c6), FALSE, c6)
    ##
    checklist$converged <- c(
      pp.ge = ck1, 
      pp.st = ck2, 
      pp.hw = ck3, 
      r.ge  = ck4, 
      r.st  = ck5, 
      r.hw  = ck6
    )
    checklist$converged.all <- all( checklist$converged )
    checklist$converged.prop <- sum(checklist$converged) / length(checklist$converged)
    
    return(checklist)
  }
  
}
















cat(sprintf('\nLoaded BSTS Helper Functions.\n'))
