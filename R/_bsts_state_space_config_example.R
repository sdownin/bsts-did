##############################
##
##   EXAMPLE:  BSTS STATE SPACE CONFIGS
##
##############################


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
library(forecast)  ## # library(sarima); library(qualV)

set.seed(987654321)  ## reproducibility

## Directories
dir_proj <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'
dir_ext <- 'D:\\BSTS_external'
dir_plot <- file.path(dir_proj, 'plots')
dir_r <- file.path(dir_proj,'R')
##  file prefix for saving images, writing outputs, etc.
prefix <- 'bsts-cconma_'

# ## Load simulation functions - Actor index vectorized simulation
# source(file.path(dir_r,'single_intervention_sim_vec.R')) 
# ## Setting up and adding state space components to state.space list
# source(file.path(dir_r,'bsts_helper_functions.R')) 
# ## BSTS vs. DiD comparison and sensitivity analysis
# source(file.path(dir_r,'bsts_did_comparison_functions.R')) 


########################################
##
##  FUNCTIONS
##
#########################################

##
# Heidel diagnostics modified from source to return the CMV.stat
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
postPredChecks <- function(causimp, filename=NA, 
                           save.plot=TRUE, return.val=FALSE,
                           burn=NA, conv.alpha=0.05) {
  
  ppcheck.filename <- if (is.na(filename)){
    sprintf('bsts_post_pred_checks_%s.png', round(10*as.numeric(Sys.time())) )
  }  else {
    filename
  }
  
  ## output list of convergence checks (booleans) and residual dataframes
  checklist <- list()
  
  response <-  causimp$series$response  
  y <- response
  
  # ##**DEBUG**
  # print('as.numeric( response )')
  # print(response)
  # ##**
  
  npds <- length(y)
  
  intpd <- causimp$model$post.period[1]
  
  niter <- length(causimp$model$bsts.model$sigma.obs)
  
  if (is.na(burn)) {
    burn <- round( niter * .2 )
  }
  
  ## Newdata (post-intervention data) to predict via BSTS
  hasRegression <- causimp$model$bsts.model$has.regression
  if (hasRegression) {
    newdata <- causimp$model$bsts.model$predictors[1:(intpd-1), ]
    newdata <- cbind(response=response[1:(intpd-1)], newdata)
  } else {
    newdata <- NULL
  }
  

  # predict.mbsts()
  post.pred <- if (hasRegression) {
    predict.bsts(causimp$model$bsts.model, newdata = newdata , burn = burn) ## already knows horizon from 
  } else {
    predict.bsts(causimp$model$bsts.model, burn = burn, 
                 horizon = intpd-1,
                 olddata = response[1:(intpd-1)])
  }

  post.pred.dist <- post.pred$distribution
  post.pred.mean <- post.pred$mean ##apply(post.pred$distribution, c(1), mean)
  
  ## before intervention period bool dummy
  .ind <- (1:npds) < intpd
  
  
  if (save.plot) {
    png(filename = ppcheck.filename, width = 15, height = 10, units = 'in', res = 400)
  }
  ##----------- INSIDE PNG PLOT --------------------------------------
  par(mfrow=c(2,3), mar=c(2.5,2.5,2.5,1))
  
  ##============================
  ## Posterior Predictive (Y) plots 
  ##-------------------
  
  ##-----------
  ## Trace plot of Posterior Predictive distribution Markov Chain 
  post.pred.tr <- rowMeans(post.pred.dist)
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
  lines(density(post.pred.mean), lwd=2, lty=2, col='red')
  legend('topleft', legend=c('observed','predicted'), lty=c(1,2), col=c('black','red'))
  ##-----------
  
  ##-----------
  ## HISTOGRAMS & BAYESIAN P-VALUES
  # max.distrib <- apply(post.pred, c(2, 3), max)
  max.distrib <- apply(post.pred$distribution, 1, max)
  pvalue <- sum(max.distrib >= max(y[.ind]))/length(max.distrib)
  hist(max.distrib, 30, col = "lightblue", border = "grey", 
       main = paste0("C. Bayesian p-val (Max Y) = ", round(pvalue, 2)),
       xlab = "Max. in-sample forecasts")
  abline(v = max(y[.ind]), col = "darkblue", lwd = 3)
  ##-----------
  
  
  ##============================
  ## Std. Residual plots
  ##-------------------
  y.rep <- matrix(y[.ind], length(y[.ind]), (niter - burn),  byrow = FALSE)
  # res <- y.rep - (post.pred.dist - CausalMBSTS$mcmc$eps.samples[, i, ])
  res <-  y.rep - t(post.pred.dist)
  std.res <- res / apply(res,1,sd)   ## [residual i] / [ stdev of residual i]
  # std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
  
  ##-----------
  ## Trace plot of Posterior Predictive distribution Markov Chain 
  res.tr <- colMeans(std.res)
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
  qqnorm(rowMeans(std.res), main = "E. Std.Residual QQ-plot, Y")
  qqline(rowMeans(std.res))
  ##-----------
  
  ##-----------
  ## Std Residual ACF
  Acf(rowMeans(std.res), main = "");title(main='F. Std.Residual ACF, Y')
  ##-----------
  
  ##----------- end PNG PLOT --------------------------------------
  if (save.plot) {
    dev.off()
  }
  
  
  ##
  if(return.val) {
    checklist$std.residual <- std.res
    checklist$postpred.dist <- post.pred.dist
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









#########################################
## MAIN SIM SETTINGS
##########################################
n <- 100  ## number of actors (i.e., number of time series)
npds <- 520  #number of periods
intpd <-  npds * .6
effect.types <- c('quadratic') ##  c('constant','geometric','quadratic') ## c('quadratic')  
# bsts.niters <- list(sm.start = 50, sm.max = 100, lg.start = 500, lg.max = 1000)  ## lg=10k for full run
# bsts.ctrl.cats <- 4




##########################################
##
##    LOAD CCONMA DATA HERE
##
###########################################
###
##
##**TODO**
##
##
##




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


# filex <- '____DEBUG__bsts-illus_default-st-sp__GRIDSEARCH_output__n200_pd520_niter1500_covCats1_msize3_16739259770_d1f1g1h1i1j1.rds'
# filex <- '____DEBUG__bsts-illus_default-st-sp__GRIDSEARCH_output__n200_pd520_niter1000_covCatsNA_msize3_16739448493_d1f1g1h1i1j1.rds'
# filex <- '____DEBUG__bsts-illus_default-st-sp__GRIDSEARCH_output__n200_pd520_niter1000_covCatsNA_msize3_16739456005_d1f1g1h1i1j1.rds'
filex <- '____DEBUG__bsts-illus_default-st-sp__GRIDSEARCH_output__n200_pd520_niter1000_covCatsNA_msize3_16739463285_d1f1g1h1i1j1.rds'
#
. <- readRDS(file.path(dir_ext, filex))
# bsts.model <- .$quadratic[[1]]$CausalImpact$model$bsts.model


simdfx <- .$sim$df

### 1 CONTROL GROUP plus covariate series
dfc <- simdfx %>%
  dplyr::filter( 
    ! is.na(match_id), 
    group=='control'
  ) %>% 
  group_by(t) %>%
  dplyr::summarize(
    y = mean(y, na.rm=T),
    # c0_seasonal = mean(season.val, na.rm=T),
    c1_mean = mean(c1, na.rm=T),
    c2_mean = mean(c2, na.rm=T),
    c3_mean = mean(c3, na.rm=T),
    c1_sd = sd(c1, na.rm=T),
    c2_sd = sd(c2, na.rm=T),
    c3_sd = sd(c3, na.rm=T),
    c1_skew = ifelse(length(c1)<=1, NA, skewness(c1, na.rm = T, type = 2)),
    c2_skew = ifelse(length(c2)<=1, NA, skewness(c2, na.rm = T, type = 2)),
    c3_skew = ifelse(length(c3)<=1, NA, skewness(c3, na.rm = T, type = 2))#,
  ) #%>% pivot_wider(names_from, values_from)
dfc$t <- NULL

# dft <- simdfx %>%
#   dplyr::filter( 
#     ! is.na(match_id), 
#     group=='treatment'
#   ) %>% 
#   group_by(t) %>%
#   dplyr::summarize(
#     y_control = mean(y, na.rm=T),
#     # c0_seasonal = mean(season.val, na.rm=T),
#     c1_mean = mean(c1, na.rm=T),
#     c2_mean = mean(c2, na.rm=T),
#     c3_mean = mean(c3, na.rm=T),
#     c1_sd = sd(c1, na.rm=T),
#     c2_sd = sd(c2, na.rm=T),
#     c3_sd = sd(c3, na.rm=T),
#     c1_skew = ifelse(length(c1)<=1, NA, skewness(c1, na.rm = T, type = 2)),
#     c2_skew = ifelse(length(c2)<=1, NA, skewness(c2, na.rm = T, type = 2)),
#     c3_skew = ifelse(length(c3)<=1, NA, skewness(c3, na.rm = T, type = 2))#,
#   ) #%>% pivot_wider(names_from, values_from)
# dft$t <- NULL

##------------- CONTROL GROUP OUTCOME SERIES --------------------------

bsts.niter <- 1000

##**ASSUMES outcome series if the untreated control group series**
##   intpd = intervention period
##   npds = number of periods in time series
##
y.pre <- dfc$y
##
y.pre[intpd:npds] <- NA

## PREDICTORS
predictors <- as.data.frame( dfc[, -1] )  ## exclude outcome series, drop first column


##------------------------------------------------------------
## Initiate empty state space configuration list
st.sp <- list()

## Add local level to trend component of state space
st.sp <- AddLocalLevel(st.sp, y.pre, 
   sigma.prior = SdPrior(
     sigma.guess = 0.01,    ## guess
     sample.size = 0.01,    ## default
     initial.value= 0.01,    ## defaults to sigma.guess
     upper.limit=Inf
   ), 
   initial.state.prior = NormalPrior(
     mu= 0.01,   ##guess
     sigma=0.01,  ## guess
     initial.value=0.01  ## default to mu
   )
)

## Add Seasonality to state space
st.sp <- AddSeasonal(st.sp, y.pre, nseasons = 52, season.duration = 1,
   sigma.prior = SdPrior(
     sigma.guess = 0.01,    ## guess
     sample.size = 0.01,    ## default
     initial.value= 0.01,    ## defaults to sigma.guess
     upper.limit=Inf
   ), 
   initial.state.prior = NormalPrior(
     mu= 0.01,   ##guess
     sigma=0.01,  ## guess
     initial.value=0.01  ## default to mu
   )
)
print(st.sp)


### PREDICTORS FOR REGRESSION
##**TODO**
## predictors data frame is size = [periods, covariates],
##            (and excluding the outcome series)
  

## Fit BSTS model: formula “y.pre ~ .” indicates inclusion of regression  
## component (all covariates in df1[,-1] dataframe)
bsts.fit <- bsts(y.pre ~ . ,  ## This specifies including regression (all columns in pre)
                 state.specification = st.sp,
                 data = predictors,  
                 niter = bsts.niter)

plot(bsts.fit)
plot(bsts.fit, 'components')
plot(bsts.fit, 'coefficients')
plot(bsts.fit, 'predictors')

##---------------------------------------------------------------------
## Initiate empty state space configuration list
st.sp2 <- list()

## Add local level to trend component of state space
# st.sp <- AddLocalLevel(st.sp, y.pre, 
#                        sigma.prior = SdPrior(
#                          sigma.guess = 0.01,    ## guess
#                          sample.size = 0.01,    ## default
#                          initial.value= 0.01,    ## defaults to sigma.guess
#                          upper.limit=Inf
#                        ), 
#                        initial.state.prior = NormalPrior(
#                          mu= 0.01,   ##guess
#                          sigma=0.01,  ## guess
#                          initial.value=0.01  ## default to mu
#                        )
# )
st.sp2 <- AddSemilocalLinearTrend(st.sp2, y.pre)

## Add Seasonality to state space
st.sp2 <- AddSeasonal(st.sp2, y.pre, nseasons = 52, season.duration = 1,
                     sigma.prior = SdPrior(
                       sigma.guess = 0.01,    ## guess
                       sample.size = 0.01,    ## default
                       initial.value= 0.01,    ## defaults to sigma.guess
                       upper.limit=Inf
                     ), 
                     initial.state.prior = NormalPrior(
                       mu= 0.01,   ##guess
                       sigma=0.01,  ## guess
                       initial.value=0.01  ## default to mu
                     )
)
print(st.sp2)


bsts.fit2 <- bsts(y.pre ~ . ,  ## This specifies including regression (all columns in pre)
                 state.specification = st.sp2,
                 data = predictors,  
                 niter = bsts.niter)

plot(bsts.fit2)

CompareBstsModels(list(level=bsts.fit, linearTrend=bsts.fit2))



##------------------------------------------------------------

## Initiate empty state space configuration list
st.sp3 <- list()

## Add local level to trend component of state space
# st.sp <- AddLocalLevel(st.sp, y.pre, 
#                        sigma.prior = SdPrior(
#                          sigma.guess = 0.01,    ## guess
#                          sample.size = 0.01,    ## default
#                          initial.value= 0.01,    ## defaults to sigma.guess
#                          upper.limit=Inf
#                        ), 
#                        initial.state.prior = NormalPrior(
#                          mu= 0.01,   ##guess
#                          sigma=0.01,  ## guess
#                          initial.value=0.01  ## default to mu
#                        )
# )
st.sp3 <- AddLocalLinearTrend(st.sp3, y.pre)

## Add Seasonality to state space
st.sp3 <- AddSeasonal(st.sp3, y.pre, nseasons = 52, season.duration = 1,
                      sigma.prior = SdPrior(
                        sigma.guess = 0.01,    ## guess
                        sample.size = 0.01,    ## default
                        initial.value= 0.01,    ## defaults to sigma.guess
                        upper.limit=Inf
                      ), 
                      initial.state.prior = NormalPrior(
                        mu= 0.01,   ##guess
                        sigma=0.01,  ## guess
                        initial.value=0.01  ## default to mu
                      )
)
print(st.sp3)


bsts.fit3 <- bsts(y.pre ~ . ,  ## This specifies including regression (all columns in pre)
                  state.specification = st.sp3,
                  data = predictors,  
                  niter = bsts.niter)

plot(bsts.fit3)
CompareBstsModels(list(level=bsts.fit, semilocalLinear=bsts.fit2, linear=bsts.fit3))


###-----------------------------------------------------------

## Post period response
y.post <- dfc$y
y.post[1:(intpd-1)] <- NA  ## SETTING values before intervention to missing


## Causal impact estimation: fitted BSTS model  forecasts the counterfactual
causimp <- CausalImpact(bsts.model = bsts.fit,
                        post.period.response = y.post[!is.na(y.post)],
                        alpha=0.05, model.args=list(niter = bsts.niter))

causimp2 <- CausalImpact(bsts.model = bsts.fit2,
                        post.period.response = y.post[!is.na(y.post)],
                        alpha=0.05, model.args=list(niter = bsts.niter))

causimp3 <- CausalImpact(bsts.model = bsts.fit3,
                         post.period.response = y.post[!is.na(y.post)],
                         alpha=0.05, model.args=list(niter = bsts.niter))

causimp
causimp2
causimp3

plot(causimp)
plot(causimp2)
plot(causimp3)

## RUN CONVERGENCE CHECKS AND DIAGNOSTICS
convcheck <- postPredChecks(causimp, save.plot=F, return.val = T)
cat(convcheck$summary)

convcheck <- postPredChecks(causimp2, save.plot=F, return.val = T)










