library(bsts)

## EX 6
data(iclaims)
ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
model <- bsts(iclaimsNSA ~ ., state.specification = ss, data =
                initial.claims, niter = 1000)
plot(model)
plot(model, "components")
plot(model, "coefficients")
plot(model, "predictors")



## Example 7:  Regressors with multiple time stamps.
number.of.time.points <- 50       ## npds
sample.size.per.time.point <- 10  ## n  (number of actors)
total.sample.size <- number.of.time.points * sample.size.per.time.point
sigma.level <- .1
sigma.obs <- 1

## Simulate some fake data with a local level state component.
trend <- cumsum(rnorm(number.of.time.points, 0, sigma.level))
predictors <- matrix(rnorm(total.sample.size * 2), ncol = 2)
colnames(predictors) <- c("X1", "X2")
coefficients <- c(-10, 10)
regression <- as.numeric(predictors %*% coefficients)
y.hat <- rep(trend, sample.size.per.time.point) + regression
y <- rnorm(length(y.hat), y.hat, sigma.obs)

## Create some time stamps, with multiple observations per time stamp.
first <- as.POSIXct("2013-03-24")
dates <- seq(from = first, length = number.of.time.points, by = "month")
timestamps <- rep(dates, sample.size.per.time.point)

## Run the model with a local level trend, and an unnecessary seasonal component.
ss <- AddLocalLevel(list(), y)
ss <- AddSeasonal(ss, y, nseasons = 7)
model <- bsts(y ~ predictors, ss, niter = 500, timestamps = timestamps,
              seed = 8675309)
plot(model)
plot(model, "components")



#############################################################################








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
prefix <- 'bsts-illus_sensitiv_'
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
##  1. TESTING MULTI-OBSERVATION
##
##
##=======================================

## Static defaults
noise.level <- 1.3
b4 <- 1.0  ## seasonal component weight
b5 <- 0.04 ## yearly growth rate
dgp.nseasons= 52
dgp.freq= 1

## Variables Ranges to grid search (optional)
n <- 100 ## 200 ## list( 100, 200 ) ## ** INTERACT SIZE (N) with GRANULARITY
npds <- 260
intpd <- round( npds * 5/6)
treat.rule <- 'random'  ## 'below.benchmark'
seasonality <- TRUE   ## c(TRUE,  FALSE )
prior.sd.scenario <- 'sd.low' ## list('sd.low','sd.high')  ## sd.low
expect.mod.size <- 1
## FOCAL CONSTRUCT
dgp.ar <- 0  ## list(0, 0.1, 0.5) used in 20230105 run
## STATE SPACE CONFIGURATIONS
# st.sp.lists <- list(
#   `8b`= c('AddLocalLevel','AddSeasonal','AddRegression')#,
#   # `15b`=c('AddSemilocalLinearTrend','AddSeasonal','AddRegression')
# )



simlist <- list()
##--------------------------------------------
## Append simulation configuration to simlist
##--------------------------------------------
key <- 'base'
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
  bsts.state.specs=list(
    level=list(
      getStateSpaceConfBySimScenario('AddLocalLevel', 'sd.low'),
      getStateSpaceConfBySimScenario('AddSeasonal', 'sd.low')
    ),
    reg=list(
      getStateSpaceConfBySimScenario('AddLocalLevel', 'sd.low'),
      getStateSpaceConfBySimScenario('AddSeasonal', 'sd.low'),
      getStateSpaceConfBySimScenario('AddRegression')
    )
  ),
  expect.mod.size=expect.mod.size,
  # bsts.state.specs=list(list(AddSemilocalLinearTrend),list(AddSemilocalLinearTrend,AddStudentLocalLinearTrend)),
  rand.seed = 13579
)


##
effect.types <- c('quadratic') ##  c('constant','geometric','quadratic') ## c('quadratic')  
## Scale bsts iterations in increments (doubling) starting at 10k (more iterations for models that converge more slowly)
bsts.niter.start <- 1000  ## 100k
bsts.niter.max   <- 1000  ## 100k
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
                                      bsts.niter = bsts.niter.start, ## **START at MAX niter for large npds **
                                      bsts.max.iter= bsts.niter.max
)  ## D:\\BSTS_external

# key <- 'reg'
simx <- readRDS(simlist.files$base$file[1])
causimp1 <- getCausalImpactFromSimlist(list(base=simx), 'base', 'quadratic', 1)
causimp2 <- getCausalImpactFromSimlist(list(base=simx), 'base', 'quadratic', 2)
bsts1 <- causimp1$model$bsts.model
bsts2 <- causimp2$model$bsts.model
plot(bsts1)
plot(bsts2)
plot(bsts2, 'coefficients')
plot(bsts2, 'predictors')
summary(bsts1)
summary(bsts2)
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################



simbase <- simlist$base$sim
simdf <- simbase$df
## Set group name 'gname' field, where 0 = control, # = period of treatment
simdf$match_pd <- as.numeric(simdf$match_pd)
simdf$gname <- 0
simdf$gname[simdf$group=='treatment'] <- simdf$match_pd[simdf$group=='treatment']
## Remove NAs
simdf <- simdf[!is.na(simdf$match_id), ]

response <- simdf$y
predictors <- simdf[,c('c1','c2','c3')]


# y.pre.treat.NAs.post.treat <- c(dat$treatment_y_outcome[1:(intpd-1)], rep(NA,npds-intpd+1))
#
st.sp <- list()
st.sp <- AddLocalLevel(st.sp, y = response)
st.sp <- AddSeasonal(st.sp, y = response, nseasons = 52, season.duration = 1)
mobsts <- bsts(formula = response ~ ., 
               state.specification = st.sp, 
               family = 'gaussian', 
               data = predictors, 
               timestamps = predictors$t,
               niter = 500, 
               expected.model.size = 1,
               seed = 12345)

plot(mobsts)
plot(mobsts, 'predictors')
plot(mobsts, 'coefficients')

post.pd.idx <- which(predictors$t >= intpd)

dim(mobsts$one.step.prediction.errors)


err.s <- colMeans(bsts2$one.step.prediction.errors)
err.m <- colMeans(mobsts$one.step.prediction.errors)
plot(density(err.m), 'BSTS 1-step ahead prediction error distributions')
lines(density(err.s), col=2, lty=2)
legend('topleft',legend = c('multi.obs','single.obs'),lty=1:2,col=1:2)
t.test(err.s, err.m, mu=0)


## Use BSTS prediction of counterfactual to estimate CausalImpact
impact_amount <- CausalImpact(bsts.model=mobsts,
                              post.period.response = response[post.pd.idx],
                              alpha=0.05, model.args = list(niter = 500), 
                              )
## POSTERIOR PREDICTIVE CHECKS
# ppcheck.filename <- file.path(getwd(),
#                               sprintf('%s_bsts_post_pred_checks_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
#                                       prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id))
convcheck <- postPredChecks(impact_amount, save.plot = F, return.val = T)
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

plot(bsts.model)
PlotBstsComponents(bsts.model) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
plot(bsts.model, 'predictors')
# PlotBstsState(bsts.model, main=sprintf('BSTS State: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsResiduals(bsts.model, main=sprintf('BSTS Residuals: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsPredictionErrors(bsts.model, main=sprintf('BSTS Pred.Err: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsForecastDistribution(bsts.model, main=sprintf('BSTS Forecast Dist: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
if (bsts.model$has.regression) {
  PlotBstsSize(bsts.model, main=sprintf('BSTS Size: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
}












############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
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
    # ##
    # x1_sum = sum(x1, na.rm=T),
    # x2_sum = sum(x2, na.rm=T),
    # x3_sum = sum(x3, na.rm=T),
    # ##
    # c1_sum = sum(c1, na.rm=T),
    # c2_sum = sum(c2, na.rm=T),
    # c3_sum = sum(c3, na.rm=T),
    # #
    # b1_sum = sum(b1, na.rm=T),
    # b2_sum = sum(b2, na.rm=T),
    # b3_sum = sum(b3, na.rm=T),
    # #
    # u_sum = sum(u, na.rm=T),
    # v_sum = sum(v, na.rm=T),
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
##----------------------------
## State Space Configuration
##----------------------------
h <- 2


## h'th BSTS state space configuration (state component list)
state.conf <- simlist$base$bsts.state.specs[[ h ]]

## names of 
state.comps <- unname(sapply(state.conf, function(x)x$name, simplify = T))

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
        # }  else if (state.conf.item$name=='AddAutoAr') {
        #   .lags <- 3
        #   state.conf$prior <- SpikeSlabArPrior(
        #     lags = .lags,
        #     prior.inclusion.probabilities = GeometricSequence(length= .lags, initial.value= 0.9, discount.factor= 0.3),
        #     prior.mean = rep(0, .lags),
        #     prior.sd = GeometricSequence(length = .lags, initial.value=0.5, discount.factor=0.1),
        #     sdy=.5,
        #     prior.df = 1,
        #     expected.r2 = .5,
        #     sigma.upper.limit = Inf,
        #     truncate = TRUE
        #   )
      }
      ## [[ IF NOT AddSharedLocalLevel(), NOT INCLUDING SPIKESLABPRIOR ]]
      st.sp <- updateStateSpaceAddComponentFromConfList(st.sp,  
                                                        y.pre.treat.NAs.post.treat,  
                                                        state.conf.item)
      if(verbose) cat(sprintf('add to state.space: %s\n',state.conf.item$name))
    }
  }
} else {
  ## Default in CausalImpact package
  st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
}
# print(st.sp)

##-------------------------------
## Regression component of model
if ('AddRegression' %in% state.comps) {
  bsts.input.form <- y.pre.treat.NAs.post.treat ~ . ## with regression formula
} else {
  bsts.input.form <- y.pre.treat.NAs.post.treat  ## without regression vector
}

bsts.niter <- bsts.niter.max

## BSTS model
bsts.model <- bsts(formula = bsts.input.form,
       state.specification = st.sp,
       data = predictors,
       expected.model.size = expect.mod.size,
       niter = bsts.niter,
       ping = ifelse(verbose, round(bsts.niter/10), 0))
 

## Use BSTS prediction of counterfactual to estimate CausalImpact
impact_amount <- CausalImpact(bsts.model=bsts.model,
                              post.period.response = post.period.response,
                              alpha=0.05, model.args = list(niter = bsts.niter))
## POSTERIOR PREDICTIVE CHECKS
# ppcheck.filename <- file.path(getwd(),
#                               sprintf('%s_bsts_post_pred_checks_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
#                                       prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id))
convcheck <- postPredChecks(impact_amount, save.plot = F, return.val = T)
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

plot(bsts.model)
PlotBstsComponents(bsts.model) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
plot(bsts.model, 'predictors')
# PlotBstsState(bsts.model, main=sprintf('BSTS State: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsResiduals(bsts.model, main=sprintf('BSTS Residuals: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsPredictionErrors(bsts.model, main=sprintf('BSTS Pred.Err: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsForecastDistribution(bsts.model, main=sprintf('BSTS Forecast Dist: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
if (bsts.model$has.regression) {
  PlotBstsSize(bsts.model, main=sprintf('BSTS Size: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
}


##





















