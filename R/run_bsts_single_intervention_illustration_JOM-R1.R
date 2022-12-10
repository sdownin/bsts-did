#####################################################
##
##   Dynamic Causal Inference Simulations
##
##    Run script for simulations of
##     - internal intervention (self-selection) endogeneity
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
prefix <- 'bsts-illustration_'

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


###
## Aggregate Aggregated Timeseries Simulation Panel Dataframe
###
getAggregatedSimPanelDf <- function(tmpdf, pd.agg, 
                                    na.rm=TRUE) {
  if (is.na(pd.agg) | pd.agg <= 1) {
    cat(sprintf('\npd.agg=%s is either NA or <= 1. Try setting pd.agg >= 2.\n',pd.agg))
  }
  ##
  ts <- sort(unique(tmpdf$t))
  npds.orig <- length(ts)
  npds.new <- npds.orig / pd.agg
  actors <- sort(unique(tmpdf$actor[!is.na(tmpdf$match_id)]))
  
  aggmap <- data.frame(t.old=1:npds.orig,
                       t.new=rep(1:npds.new, each=pd.agg))
  
  ##
  intpd.old <- unique(tmpdf$t[which(tmpdf$t.post.intpd==1)])[1]
  intpd.new <- aggmap$t.new[which(aggmap$t.old == intpd.old)]
  
  ##
  tmpdf$match_pd <- as.numeric(tmpdf$match_pd)
  tmpdf$match_pd[tmpdf$match_pd == intpd.old] <- intpd.new
  
  
  tmpdf$t.agg <- NA
  for (i in 1:nrow(aggmap)) {
    idx <- which( tmpdf$t==aggmap$t.old[i] )
    tmpdf$t.agg[idx] <- aggmap$t.new[i]
  }
  
  if (na.rm) {
    tmpdf <- tmpdf[ !is.na(tmpdf$match_id), ]
  }
  
  aggdf <- tmpdf %>% 
    group_by(effect.type, t.agg, actor) %>% 
    summarize(
      n_actors = n(),
      group = paste(unique(group), collapse = '|'), 
      group.color = paste(unique(group.color), collapse = '|'),
      match_id = paste(unique(match_id), collapse = '|'), 
      match_pd = paste(unique(match_pd), collapse = '|'), 
      ##
      y = mean(y, na.rm=T),
      ##
      x1 = mean(x1, na.rm=T),
      x2 = mean(x2, na.rm=T),
      x3 = mean(x3, na.rm=T),
      #
      b1 = mean(b1, na.rm=T),
      b2 = mean(b2, na.rm=T),
      b3 = mean(b3, na.rm=T),
      ##
      c1 = mean(c1, na.rm=T),
      c2 = mean(c2, na.rm=T),
      c3 = mean(c3, na.rm=T),
      ##
      season.val = mean(season.val, na.rm=T),
      u = mean(u, na.rm=T),
      v = mean(v, na.rm=T)
    )
  
  t.aggs <- sort(unique(aggdf$t.agg))
  aggdf$t.agg.post.intpd <- NA
  for (i in 1:length(t.aggs)) {
    t.new <- t.aggs[ i ]
    aggdf$t.agg.post.intpd[which(aggdf$t.agg==t.new)] <- ( t.new - (intpd.new - 1) )
  }
  
  ## Match names for call in simulation comparison function for DiD vs. BSTS
  aggdf$t <- aggdf$t.agg
  aggdf$t.post.intpd <- aggdf$t.agg.post.intpd
  
  # ##
  # aggdf$gname <- 0
  # aggdf$gname[aggdf$group=='treatment'] <- aggdf$match_pd[aggdf$group=='treatment']
  # aggdf$gname <- as.numeric(aggdf$gname)
  
  ##
  return(aggdf)
}

###
## Update Simlist configurations and simulated panel dataframes for aggregated periods
###
updateSimlistAggregatePd <- function(simlist, pd.agg, na.rm=TRUE) {
  
  for (i in 1:length(simlist)) 
  {
    simx <- simlist[[ i ]]
    
    if (is.null(simx$sim)) {
      cat(sprintf('\nsimlist item i=%s is missing "sim" object. First call runSimUpdateSimlist().\n',i))
    }
    
    aggdf <- getAggregatedSimPanelDf(simx$sim$df, pd.agg=pd.agg, na.rm=na.rm)
    simlist[[i]]$sim$df <- aggdf
    simlist[[i]]$npds <- length(unique(aggdf$t.agg))
    simlist[[i]]$intpd <- unique(aggdf$t.agg[which(aggdf$t.agg.post.intpd==1)])[1]
  }
  
  return(simlist)
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
postPredChecks <- function(causimp, filename=NA, 
                           save.plot=TRUE, return.val=FALSE,
                           burn=NA, conv.alpha=0.05) {
  
  ppcheck.filename <- if (is.na(filename)){
    sprintf('bsts_post_pred_checks_%s.png', round(10*as.numeric(Sys.time())) )
  }  else {
    filename
  }
  
  response <- causimp$series$response
  y <- response
  
  npds <- length(y)
  
  intpd <-causimp$model$post.period[1]
  
  niter <- nrow(causimp$model$bsts.model$coefficients)
  
  if (is.na(burn)) {
    burn <- round( niter * .2 )
  }
  
  ## Newdata (post-intervention data) to predict via BSTS
  newdata <- causimp$model$bsts.model$predictors[1:(intpd-1), ]
  newdata <- cbind(response=response[1:(intpd-1)], newdata)
  
  # predict.mbsts()
  post.pred <- predict.bsts(causimp$model$bsts.model, newdata = newdata , burn = burn)
  post.pred.dist <- post.pred$distribution
  post.pred.mean <- post.pred$mean ##apply(post.pred$distribution, c(1), mean)
  
  ## before intervention period bool dummy
  .ind <- (1:npds) < intpd
  
  if (save.plot) {
    png(filename = ppcheck.filename, width = 15, height = 10, units = 'in', res = 400)
  }
  ##----------- INSIDE PNG PLOT --------------------------------------
    par(mfrow=c(2,3), mar=c(2.5,2.5,2.5,1))
  
    ## DENSITIES
    plot(density(y[.ind]),  main = "A. Density comparison, Y")
    lines(density(post.pred.mean), lwd=2, lty=2, col='red')
    legend('topleft', legend=c('observed','predicted'), lty=c(1,2), col=c('black','red'))
    
    ## HISTOGRAMS & BAYESIAN P-VALUES
    # max.distrib <- apply(post.pred, c(2, 3), max)
    max.distrib <- apply(post.pred$distribution, 1, max)
    pvalue <- sum(max.distrib >= max(y[.ind]))/length(max.distrib)
    hist(max.distrib, 30, col = "lightblue", border = "grey", 
         main = paste0("B. Bayesian p-val (Max Y) = ", round(pvalue, 2)),
         xlab = "Max. in-sample forecasts")
    abline(v = max(y[.ind]), col = "darkblue", lwd = 3)
    
    ## Trace plot of Posterior Predictive distribution Markov Chain 
    post.pred.tr <- rowMeans(post.pred.dist)
    ##----------
    gd <- geweke.diag(post.pred.tr, .1, .5)
    gd.z <- abs(gd$z) ## z statistic
    gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
    gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
    ##
    hd <- heidel.diag.mod(post.pred.tr)
    hd.st.cmv <- hd[1,'CMV.stat']
    hd.st.p <- hd[1,'pvalue']
    hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
    hd.st.start <- hd[1,'start']
    hd.hw.eps <- 0.1 ## default 0.1
    hd.hw <- hd[1,'halfwidth']
    hd.hw.mean <- hd[1,'mean']
    hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
    ##----------
    rng <- range(post.pred.tr)
    ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
    plot(post.pred.tr, type='l', main='C. Posterior Predicted MCMC Trace, Y' ,
         ylim=ylims
         )
    mtext(text = sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
                         gd.result,gd.z,gd.p,
                         hd.st.result, hd.st.cmv, hd.st.p,
                         hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps), 
          side = 3, line=-4.5, outer = F)
    
    # Residual plots
    y.rep <- matrix(y[.ind], length(y[.ind]), (niter - burn),  byrow = FALSE)
    # res <- y.rep - (post.pred.dist - CausalMBSTS$mcmc$eps.samples[, i, ])
    res <-  y.rep - t(post.pred.dist)
    std.res <- res / apply(res,1,sd)   ## [residual i] / [ stdev of residual i]
    # std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
    qqnorm(rowMeans(std.res), main = "D. Std.Residual QQ-plot, Y")
    qqline(rowMeans(std.res))
    Acf(rowMeans(std.res), main = "");title(main='E. Std.Residual ACF, Y')
    
    ## Trace plot of Posterior Predictive distribution Markov Chain 
    res.tr <- colMeans(std.res)
    ##----------
    gd <- geweke.diag(res.tr, .1, .5)
    gd.z <- abs(gd$z) ## z statistic
    gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
    gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
    ##
    hd <- heidel.diag.mod(res.tr)
    hd.st.cmv <- hd[1,'CMV.stat']
    hd.st.p <- hd[1,'pvalue']
    hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
    hd.st.start <- hd[1,'start']
    hd.hw.eps <- 0.1 ## default 0.1
    hd.hw <- hd[1,'halfwidth']
    hd.hw.mean <- hd[1,'mean']
    hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
    ##----------
    rng <- range(res.tr)
    ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
    plot(res.tr, type='l', main='F. Std.Residual MCMC Trace, Y',
         ylim=ylims)
    mtext(text = sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
                         gd.result,gd.z,gd.p,
                         hd.st.result, hd.st.cmv, hd.st.p,
                         hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps), 
          side = 3, line=-4.5, outer = F)
  ##----------- end PNG PLOT --------------------------------------
  if (save.plot) {
    dev.off()
  }
    
  ##
  if(return.val) {
    return(list(
      std.res=std.res, 
      post.pred.dist=post.pred.dist
    ))
  }

}


###
## GGPLOT OF DYNMAIC DID from ATTGT object
##  - MODIFIED FROM Callaway & Sant'Anna 2021
###
ggdid.agg.es <- function(attgt,
                         ylim=NULL,
                         xlab=NULL,
                         ylab=NULL,
                         title="",
                         # xgap=NA,
                         legend=TRUE,
                         ref_line = 0,
                         theming = FALSE,
                         alpha=0.05,
                         ...) {
  
  ## DiD pre-test parallel trends Wald Chi-Sq  test stat and p-val
  W <- attgt$W[1]
  Wpval <- attgt$Wpval[1]
  W.result <- ifelse(Wpval < alpha, 'FAIL', 'PASS')
  
  ## DYNAMIC EFFECTS AND EVENT STUDIES
  object <- aggte(attgt, type = "dynamic", bstrap = TRUE)
  
  if ( !(object$type %in% c("dynamic","group","calendar")) ) {
    stop(paste0("Plot method not available for this type of aggregation"))
  }

  ##
  post.treat <- 1*(object$egt >= 0)
  results <- cbind.data.frame(year=object$egt,
                              att=object$att.egt,
                              att.se=object$se.egt,
                              post=as.factor(post.treat))
  results <- rbind(results[1,], results)
  results[1, c('att','att.se')] <- NA
  results$year[1] <- min(results$year, na.rm = T) - 1
  results$c <- ifelse(is.null(object$crit.val.egt), abs(qnorm(.025)), object$crit.val.egt)
  
  npds <- nrow(results)
  # xgap <- ifelse(is.na(xgap), round(npds/4), xgap)
  
  if (title == "") {
    # get title right depending on which aggregation
    title <- ifelse(object$type=="group", "Average Effect by Group", ifelse(object$type=="dynamic", "Average Effect by Length of Exposure", "Average Effect by Time Period"))
  }
  
  # p <- gplot(results, ylim, xlab, ylab, title, xgap, legend, ref_line, theming)
  # 
  # p
  rng.y <- range(results$att, na.rm = T)
  rng.x <- range(results$year, na.rm=T)
  an.y <- max( (results$att + results$c * results$att.se), na.rm = T) + 0.3*diff(rng.y)
  an.x <- min(as.numeric(results$year), na.rm=T) + ifelse(npds >= 7, 0.3*diff(rng.x), 0.4*diff(rng.x))
  
  p <- ggplot(results,
              aes(x=as.numeric(year), y=att, ymin=(att-c*att.se),
                  ymax=(att+c*att.se))) +
    geom_point(colour='black', size=1.8) +
    #geom_ribbon(aes(x=as.numeric(year)), alpha=0.2) +
    geom_errorbar(colour='darkgray', width=0.3) +
    scale_y_continuous(limits=ylim) +
    #scale_x_discrete(breaks=dabreaks, labels=as.character(dabreaks)) +
    # scale_x_continuous(breaks=as.numeric(dabreaks), labels=as.character(dabreaks)) +
    # scale_color_manual(drop=FALSE, values=c("#e87d72","#56bcc2"), breaks = c(0, 1), labels = c('Pre','Post')) +
    labs(x = xlab, y = ylab, title = title, color = NULL) +
    # scale_color_manual(values = c('black','darkgray')) +
    geom_vline(xintercept= -0.5, lty=2) +  ## -1 to align intervention time
    geom_hline(yintercept = 0) +
    ylab('ATT') + xlab('Event Time') +
    theme_bw() +
    theme(legend.position='none') +
    annotate('text', x=an.x, y=an.y, label=sprintf('Parallel Trends Pretest: %s (Wald Chisq=%.1f, p=%.3f)',W.result,W,Wpval)) +
    ggtitle(' DiD: Average Effect by Length of Exposure')
  
  p
}





####################################
##  RUN SIMULATION IN LOOP OVER EFFECT TYPES
##  - 1. runSimSingleInterventionEffectComparison() on simlist
##  - 2. DiD & BSTS results
##  - 3. comparison of DiD & BSTS performance
##   returns full simlist
######################################
runSimUpdateSimlist <- function(simlist,     ## n, npds, intpd moved into simlist elements
                                effect.types=c('constant','quadratic','geometric'), 
                                sim.id=round(10*as.numeric(Sys.time())),
                                plot.show=F, plot.save=F) {
  
  # print("runSimBstsDiDComparison()::SIMLIST INPUT:")
  # print(simlist)
  if (length(simlist) > 0 & length(names(simlist))==0) {
    names(simlist) <- 1:length(simlist)
  }
  ##----------------------------
  ## Run Simulation List
  ##----------------------------
  # sim.id <- ifelse( is.na(sim.id), round(as.numeric(Sys.time())), sim.id)
  for (i in 1:length(simlist)) {
    key <- names(simlist)[i]
    sim <- simlist[[key]]
    cat(sprintf('\nScenario label: %s\n\n', key))
    ##  
    set.seed( ifelse(is.null(sim$rand.seed), 54321, sim$rand.seed) )
    noise.level <- ifelse(is.null(sim$noise.level), 0, sim$noise.level)
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
      noise.level = noise.level,
      ##
      b4 = ifelse(is.null(sim$b4), 0, sim$b4), ## past performance
      b5 = ifelse(is.null(sim$b5), 0, sim$b5), ## growth (linear function of time t)
      b9 = ifelse(is.null(sim$b9), 0, sim$b9), ## Autocorrelation
      ## Dynamic treatment effect function parameters
      w0 = 2.0,  ## constant
      w1 = 0.18,  ## linear
      w2 = -1.25 / (.1*sim$npds)^2  ,  ##-0.009,  ## quadratic  ## -0.005, ## ***** made steeper curve= -0.008 *****
      ## TODO: CHECK w2.shift SENSITIVITY - this shifts quadratic curve several periods to the right so that treatment effect increases slowly  
      w2.shift = -round( sqrt(sim$npds)*.8 ),  ## optimal value here is likely a function of the combination of treatment effect function parameters
      ##
      nseasons  = ifelse(is.null(sim$dgp.nseasons), NA, sim$dgp.nseasons), 
      season.frequency = ifelse(is.null(sim$dgp.freq), NA, sim$dgp.freq),
      ## Plotting
      plot.show = plot.show, ## TRUE
      plot.save = plot.save  ## TRUE
    )
  }
  
  return(simlist)
}

  
# ## Save simulation list as serialized data file
# simlist.file <- sprintf('single_intervention_SIMLIST_selection_endog_%s_%s.rds',
#                         length(simlist), sim.id)
# saveRDS(simlist, file = file.path(dir_plot, simlist.file))





####################################
##  MAIN COMPARISON FUNCTION
##  - 1. runSimSingleInterventionEffectComparison() on simlist
##  - 2. DiD & BSTS results
##  - 3. comparison of DiD & BSTS performance
##   returns full simlist
######################################
runSimUpdateCompareBstsDiD <- function(simlist,     ## n, npds, intpd moved into simlist elements
                                    effect.types=c('constant','quadratic','geometric'), 
                                    sim.id=NA,
                                    save.items.dir=NA, ## save updated simlist items to seprate RDS files
                                    bsts.niter=5e3
                                    ) {
  
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
  save.img.dir <- ifelse(is.na(save.items.dir), getwd(), save.items.dir)
  
  ##===============================
  ##  BSTS State Specification Comparison 
  ##------------------------------
  for (i in 1:length(simlist))
  {
    key <- names(simlist)[i]
    key.strip <- gsub('[|]','',key,ignore.case = F, perl = T)
    cat(sprintf('\n%s, %s\n',i, key))
    
    simlist[[key]]$cordf <- data.frame()
    simlist[[key]]$compare <- list(did=list(), bsts=list(), res.tbl=list(), 
                                   att.err.tbl=list(), 
                                   att.err.mean.bsts=list(),
                                   att.err.mean.did=list(),
                                   att.err.sd.bsts=list(),
                                   att.err.sd.did=list())
    
    ## simulation output from simulation scenario = simlist[[key]]
    npds <- simlist[[key]]$npds
    intpd <- simlist[[key]]$intpd
    n <- simlist[[key]]$n
    sim <- simlist[[key]]$sim
    # sim.id <- simlist[[key]]$sim$id

    ## list of BSTS State component lists
    bsts.state.specs <- simlist[[key]]$bsts.state.specs
    if (length(names(bsts.state.specs))==0) {
      names(bsts.state.specs) <- 1:length(bsts.state.specs)
    }
  
    
    ## Dynamic Treatment Effect Type shapes
    for (k in 1:length(effect.types)) 
    {
    
      effect.type <- effect.types[k]
      simdf <- sim$df[sim$df$effect.type == effect.type, ]
  
      simlist[[key]]$compare$bsts[[effect.type]] <- list()
      simlist[[key]]$compare$res.tbl[[effect.type]] <- list()
      simlist[[key]]$compare$att.err.tbl[[effect.type]] <- list()
      simlist[[key]]$compare$att.err.mean.bsts[[effect.type]] <- list()
      simlist[[key]]$compare$att.err.mean.did[[effect.type]] <- list()
      simlist[[key]]$compare$att.err.sd.bsts[[effect.type]] <- list()
      simlist[[key]]$compare$att.err.sd.did[[effect.type]] <- list()
      
      
      ##------------------------------
      ## DiD
      ##------------------------------
      # simdf <- simdf[simdf$effect.type==effect.type, ]
      ## Set group name 'gname' field, where 0 = control, # = period of treatment
      simdf$match_pd <- as.numeric(simdf$match_pd)
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
      
      ## PLOT DID DYNAMIC EFFECT from ATTGT object
      p.agg.es <- ggdid.agg.es(ccattgt) 
      ggsave(filename = file.path(save.img.dir,
                                  sprintf('%s_did_dynamic_effect_n%s_pd%s_ss%s_%s_%s_%s.png',
                                          prefix,n,npds,h,key.strip,effect.type,sim.id))
      )
      
      
      ## Get first treatment group actor
      tr.actor.1 <- simdf$actor[which(simdf$group=='treatment')[1]]
      
      ## SIMPLE AGGREGATION (OVERALL EFFECT) ATT
      agg.simple <- aggte(ccattgt, type='simple', bstrap = TRUE)
      # summary(agg.simple)
      
      ## GROUP EFFECT AGGREGATION (for reporting overall confidence intervals)
      agg.group <- aggte(ccattgt, type = "group", bstrap = TRUE)
      
      ## DYNAMIC EFFECTS AND EVENT STUDIES
      agg.es <- aggte(ccattgt, type = "dynamic", bstrap = TRUE)
      # summary(agg.es)
      # tidy(agg.es)
      
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
      self.select.cor <- sigmamat[2,1]  ## either off-diagonal element
      
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
                    'x1_mean','x2_mean','x3_mean',
                    'c1_mean','c2_mean','c3_mean',
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
                      'control_c1_mean','control_c2_mean','control_c3_mean'#,
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
      ## LOOP OVER BSTS STATE SPECIFICATIONS FOR THE SAME SIMULATION RUN
      for (h in 1:length(bsts.state.specs)) 
      {
        simlist[[key]]$compare$bsts[[effect.type]][[ h ]] <- list()
        
        ## h'th BSTS state space configuration (state component list)
        state.conf <- bsts.state.specs[[ h ]]
        
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
              cat(sprintf('add to state.space: %s\n',state.conf.item$name))
            }
          }
        } else {
          ## Default in CausalImpact package
          st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
        }
        # print(st.sp)
        
        cat(sprintf('\nRunning BSTS model estimation for state.conf h=%s\n',h))
        
        ## BSTS model
        bsts.model <- tryCatch(expr = {
            bsts(y.pre.treat.NAs.post.treat ~ . ,
                 state.specification = st.sp,
                 data = predictors,
                 niter = bsts.niter)
          },
          error=function(e) {
            message(sprintf('bsts() error: %s', as.character(e)))
            # message(cond)
            # # Choose a return value in case of error
            # return(NA)
          },
          warning=function(w) {
            message(sprintf('bsts() warning: %s', as.character(w)))
            # message(paste("URL caused a warning:", url))
            # return(NA)
          },
          finally={
            ##PASS
          })
        # bsts.model <- bsts(y.pre.treat.NAs.post.treat ~ . ,
        #                    state.specification = st.sp,
        #                    data = predictors,
        #                    niter = bsts.niter)
        if ( class(bsts.model) != 'bsts' ) {
          next
        }
        ##
        plot(bsts.model, main=sprintf('BSTS Plot: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
        PlotBstsComponents(bsts.model) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
        # PlotBstsState(bsts.model, main=sprintf('BSTS State: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
        # PlotBstsResiduals(bsts.model, main=sprintf('BSTS Residuals: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
        # PlotBstsPredictionErrors(bsts.model, main=sprintf('BSTS Pred.Err: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
        # PlotBstsForecastDistribution(bsts.model, main=sprintf('BSTS Forecast Dist: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
        PlotBstsSize(bsts.model, main=sprintf('BSTS Size: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
        ##
        
        # ## BSTS model for Dynamic Regression
        # bsts.model <- bsts(y.pre.treat.NAs.post.treat,
        #                    state.specification = st.sp,
        #                    niter = 5000)
        
        ## Use BSTS prediction of counterfactual to estimate CausalImpact
        impact_amount <- CausalImpact(bsts.model=bsts.model,
                                      post.period.response = post.period.response,
                                      alpha=0.05, model.args = list(niter = bsts.niter))
        ## POSTERIOR PREDICTIVE CHECKS
        ppcheck.filename <- file.path(save.img.dir,
                                      sprintf('%s_bsts_post_pred_checks_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
                                              prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id))
        postPredChecks(impact_amount, filename=ppcheck.filename)
        
        # ##
        # summary(impact_amount)
        # summary(impact_amount$model$bsts.model)
        # plot(impact_amount)
        
        # summary(impact_amount)
        # png(filename=sprintf('single_intervention_BSTS_CausalImpact_plot_%s_%s_%s.png',
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
        ci.did.general <- c(agg.group$overall.att - (agg.group$overall.se * agg.group$crit.val.egt),
                            agg.group$overall.att + (agg.group$overall.se * agg.group$crit.val.egt))
        ##-------------------
        
        ## DID
        did.res <- tidy(agg.es)
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
        
        b3diff <- data.frame(
          treat=simdf %>% dplyr::filter(group=='treatment' & actor==tr.actors[1]) %>% mutate(treat=b3) %>% dplyr::select(treat),
          ctrl=simdf %>% dplyr::filter(group=='control' & actor==co.actors[1]) %>% mutate(ctrl=b3) %>% dplyr::select(ctrl),
          diff=NA
        )
        b3diff$diff <- b3diff$treat - b3diff$ctrl
        
        ##
        att.b3 <- mean(b3diff$diff[intpd:npds])
        att.did <- agg.es$overall.att
        att.bsts <- impact_amount$summary$AbsEffect[1]
        
        # simdf %>% 
        #   filter(group=='control' & actor==tr.actors[1]) %>% 
        #   dplyr::mutate(b3.diff=b3.treat-b3.ctrl) %>% 
        #   dplyr::select(b3.diff) 
        
        ## time df add empty row at top of event time columns
        time.df <- did.res[ ,c('term','event.time')]
        # time.df <- rbind(time.df,time.df[nrow(time.df),])
        # time.df[nrow(time.df),] <- NA
        ## add placeholder row at top and set to NAs
        time.df <- rbind(time.df[1,], time.df)
        time.df[1,] <- NA
        time.df$event.time[1] <-  min(time.df$event.time, na.rm = T) - 1
        
        # ## *** FIX MISALIGNED INTERVENTION TIME INDEX BETWEEN DID & BSTS
        # ##     - add empty row on top of DiD (which indexed treatment at t=0) so now treatment is t=1
        # ##       which matches BSTS (treatment at t=1)
        did.res.adj <- did.res[ ,c('estimate','point.conf.low','point.conf.high')]
        names(did.res.adj) <- c('did.estimate','did.point.conf.low','did.point.conf.high')
        .na.row <- did.res.adj[1, ]
        .na.row[1:nrow(.na.row), ] <- NA
        did.res.adj <- rbind(.na.row, did.res.adj)
        # did.res <- rbind(did.res[-1, ], .na.row, .na.row)
        # ## ##    - remove first row to pull index t=0 up to t=1
        # did.res <- did.res[-1, ]
        
        ## Results comparison table
        res.tbl <- cbind(
          time.df,
          bsts.res,
          did.res.adj,
          b3.treat=b3diff$treat,
          b3.ctrl=b3diff$ctrl,
          b3.att=b3diff$diff 
        )
        ## ROUND RESULTS TABLE (round numeric columns)
        # res.tbl <- res.tbl[ ! is.na(res.tbl$estimate), ]
        num.cols <- names(res.tbl)[ ! names(res.tbl) %in% c('term','event.time') ]
        # res.tbl[ , num.cols] <- round( as.numeric(res.tbl[ , num.cols]), 4)  ## ** TODO: Find out why this as.numeric stopped working on 'list'
        res.tbl[ , num.cols] <- round( res.tbl[ , num.cols], 4)
        # for (i in 1:length(num.cols)) {
        #   res.tbl[ , num.cols[i] ] <- round( as.numeric( res.tbl[ , num.cols[i] ] ), 4)
        # }
        # ## MOVE ID COLUMNS TO FRONT
        # .col.idx <- which(names(res.tbl) %in% c('term','event.time'))
        # res.tbl4 <- cbind(res.tbl[, .col.idx],  res.tbl[, -.col.idx] )
        # # View(res.tbl4)
        
        ##PLOT INCLUSION PROBABILITIES
        png(filename = file.path(save.img.dir,
                                 sprintf('%s_BSTS_inclusion_probs_n%s_pd%s_ss%s_%s_%s_%s.png',
                                         prefix,n,npds,h,key.strip,effect.type,sim.id)))
        plot(impact_amount$model$bsts.model,'coefficients', main=sprintf('%s %s', key,effect.type))
        dev.off()
        
        ## PLOT DYNAMIC EFFECTS COMPARISON - DID vs. BSTS vs. DGP
        dyndf <- rbind(
          data.frame(event.time=res.tbl$event.time, series='1.BSTS', ATT=res.tbl$bsts.point.effect),
          data.frame(event.time=res.tbl$event.time, series='2.DiD', ATT=res.tbl$did.estimate),
          data.frame(event.time=res.tbl$event.time, series='3.DGP', ATT=res.tbl$b3.att)
        )
        hue2 <- hue_pal()(2)
        p.err1 <- ggplot(dyndf, aes(x=event.time, y=ATT, color=series,fill=series,linetype=series,shape=series)) + 
          geom_line(na.rm=T, size=.9) + geom_point(na.rm=T, size=2) +
          theme_bw() + xlab('Event Time') + 
          ggtitle(sprintf('Comparison of Mean ATT Estimates (DGP = %.3f):  BSTS = %.3f; DiD = %.3f',att.b3,att.bsts,att.did)) + 
          geom_vline(xintercept=-0.5, linetype='dotted')+
          geom_hline(yintercept = 0) +
          scale_color_manual(values=c(hue2[1],hue2[2],'black')) + 
          scale_shape_manual(values=c(17,19,NA)) +
          scale_linetype_manual(values=c(2,3,1))  + 
          theme(legend.position='top')
        # png(filename = file.path(save.img.dir,
        #                          sprintf('%s_BSTS_dynamic_treatment_effect_comparison_n%s_pd%s_ss%s_%s_%s_%s.png',
        #                                  prefix,n,npds,h,key.strip,effect.type,sim.id)))
        # matplot.main <- sprintf('%s: DGP = %.3f; DiD = %.3f;  BSTS = %.3f',
        #                             key.strip, mean(b3diff$diff[intpd:nrow(b3diff)]), agg.es$overall.att, impact_amount$summary$AbsEffect[1])
        # matplot(x = res.tbl$event.time, y=res.tbl[,c('point.effect','estimate','b3.att')],,
        #         type='o',lty=c(2,3,1),pch=c(1,20,NA),lwd=c(1,1,1),
        #         col=c('red','blue','black'),
        #         main=matplot.main, ylab='ATT',xlab='t')
        # legend('topright',legend=c('BSTS','DiD','DGP'),col=c('red','blue','black'),lty=c(2,3,1),pch=c(1,20,NA),lwd=c(1,1,1)) 
        # dev.off()
        
        ## PLOT ATT ESTIMATE ERROR DISTRIBUTIONS COMPARISON 
        errdf <- rbind(
          data.frame(method='BSTS', error=(res.tbl$bsts.point.effect[1:(intpd-1)] - res.tbl$b3.att[1:(intpd-1)])),
          data.frame(method='DiD', error=(res.tbl$did.estimate[1:(intpd-1)] - res.tbl$b3.att[1:(intpd-1)]))
        )
        vline.dat <- errdf %>% dplyr::group_by(method) %>% dplyr::summarize(grp.mean=mean(error,na.rm=T))
        p.err2 <- ggplot(errdf, aes(x=error, colour=method,fill=method)) + 
          geom_vline(xintercept=0)+
          geom_density(alpha=0.3, na.rm=T) +
          geom_vline(data=vline.dat, aes(xintercept=grp.mean, color=method), linetype="dashed",size=1.2) +
          # geom_histogram(alpha=0.2, position = 'identity', na.rm = T) +
          ggtitle(sprintf('Pre-Intervention Pointwise Error:\n Mean: BSTS = %.3f; DiD = %.3f\n SD:   BSTS = %.3f; DiD = %.3f',
                          mean(errdf$error[errdf$method=='BSTS'],na.rm=T), 
                          mean(errdf$error[errdf$method=='DiD'],na.rm=T),
                          sd(errdf$error[errdf$method=='BSTS'],na.rm=T), 
                          sd(errdf$error[errdf$method=='DiD'],na.rm=T)  )
                  ) + theme_bw() + xlab('Residuals')
        # ggsave(filename = file.path(save.img.dir,
        #                             sprintf('%s_ATT_est_pointwise_error_distributions_n%s_pd%s_ss%s_%s_%s_%s.png',
        #                                     prefix,n,npds,h,key.strip,effect.type,sim.id)))
        
        
        ## PLOT ATT ESTIMATE ERROR DISTRIBUTIONS COMPARISON 
        errdf <- rbind(
          data.frame(method='BSTS', error=(res.tbl$bsts.point.effect[intpd:npds] - res.tbl$b3.att[intpd:npds])),
          data.frame(method='DiD', error=(res.tbl$did.estimate[intpd:npds] - res.tbl$b3.att[intpd:npds]))
        )
        vline.dat <- errdf %>% dplyr::group_by(method) %>% dplyr::summarize(grp.mean=mean(error,na.rm=T))
        p.err3 <- ggplot(errdf, aes(x=error, colour=method,fill=method)) + 
          geom_vline(xintercept=0)+
          geom_density(alpha=0.3, na.rm=T) +
          geom_vline(data=vline.dat, aes(xintercept=grp.mean, color=method), linetype="dashed",size=1.2) +
          # geom_histogram(alpha=0.2, position = 'identity', na.rm = T) +
          ggtitle(sprintf('Post-Intervention ATT Bias:\n Mean: BSTS = %.3f; DiD = %.3f\n SD:   BSTS = %.3f; DiD = %.3f',
                          mean(errdf$error[errdf$method=='BSTS'],na.rm=T), 
                          mean(errdf$error[errdf$method=='DiD'],na.rm=T),
                          sd(errdf$error[errdf$method=='BSTS'],na.rm=T), 
                          sd(errdf$error[errdf$method=='DiD'],na.rm=T))
                  ) + theme_bw() + xlab('Residuals (Distance of Pointwise ATT Estimate from DGP)')
        # ggsave(filename = file.path(save.img.dir,
        #                             sprintf('%s_ATT_est_pointwise_error_distributions_n%s_pd%s_ss%s_%s_%s_%s.png',
        #                                     prefix,n,npds,h,key.strip,effect.type,sim.id)))
        
        # ## SIMULATED DATA GROUPS FOR EFFECT.TYPE=k
        # df.group.series <-  ddply(sim$df %>% dplyr::filter(effect.type==effect.type), .(t,effect.type,group), summarize,
        #                           min=min(y, na.rm=T),
        #                           cl=quantile(y, probs=0.025, na.rm=T),
        #                           med=median(y, na.rm=T),
        #                           cu=quantile(y, probs=0.975, na.rm=T),
        #                           max=max(y, na.rm=T))
        # ## ## sim$df.summary %>% dplyr::filter(effect.type==effect.type)
        # p.group.series <- ggplot( df.group.series, aes(x=as.numeric(t), y=med, color=group)) +
        #   geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
        #   geom_line(size=1.2) +
        #   geom_point(aes(x=as.numeric(t), y=min),pch=1,alpha=.3) + geom_point(aes(x=as.numeric(t),y=max),pch=1,alpha=.3) +
        #   geom_hline(yintercept=0) + geom_vline(xintercept=intpd, lty=2) +
        #   ylab('Y') +
        #   # facet_grid( effect.type ~ . ) +
        #   theme_bw() + theme(legend.position='top') # + ggtitle(plot.main)
        ## ## ALL SIMULATED ACTOR TIMESERIRES
        # df.ind.series <- simlist[[key]]$sim$df.plot %>% dplyr::filter(effect.type == effect.type)
        df.plot <- simlist[[key]]$sim$df
        df.ind.series <- df.plot[which(df.plot$effect.type == effect.type),]
        # df.ind.series$t0 <-df.ind.series$t - df.ind.series$t.post.intpd
        p.ind.series <- ggplot(data=df.ind.series, mapping=aes(x=t,y=y, color=group, group=actor)) +
          geom_line(size=1.05, alpha=0.2) +
          # geom_point(, color=rgb(0.2, 0.2, 0.8, 0.1))+
          geom_hline(yintercept=0)  + # facet_grid( effect.type ~ . ) +
          geom_vline(xintercept= (intpd - 0.5), linetype='dotted')+
          scale_color_manual(values = c(rgb(.2,.2,.8,.3), rgb(.8,.2,.2,.3))) +
          theme_bw() + theme(legend.position='bottom') + ggtitle('Simulated Time Series') #+
        # guides(color=guide_legend(nrow=3,byrow=TRUE)) #+
        
        ## BSTS
        # ## GROUP SUMMARY TIME SERIES
        bsts.wide <- as.data.frame(impact_amount$series)
        bsts.wide$t <- 1:nrow(bsts.wide)
        bsts.wide$t0 <- bsts.wide$t - intpd 
        # bsts.long <- gather(bsts.wide, point, val, point.effect:point.effect.upper, factor_key = T)
        # ##
        # print(bsts.long)
        # ia <- as.data.frame(impact_amount$series)
        # bsts.long <- rbind(
        #   data.frame(variable='point.effect', value=ia$point.effect),
        #   data.frame(variable='point.effect.lower', value=ia$point.effect.lower),
        #   data.frame(variable='point.effect.upper', value=ia$point.effect.upper)
        # )
        p.bsts.impact <- ggplot(bsts.wide, aes(x=t0, y=point.effect)) +
          geom_ribbon(aes(ymin=point.effect.lower,ymax=point.effect.upper), alpha=.25, size=.01, lty=1) +
          geom_line(size=1.1, lty=1) +
          # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
          geom_hline(yintercept=0) + geom_vline(xintercept=-0.5, lty=2) +
          ylab('ATT') +  xlab('Event Time') + # facet_grid( effect.type ~ . ) +
          theme_bw() + #theme(legend.position='bottom') + 
          ggtitle('BSTS CausalImpact: Pointwise Effect By Length of Exposure')
          # ggtitle(sprintf('BSTS CausalImpact (p = %.3f)',
          #                 ,pval.bsts.general))
        
        ## COMBINE TIMESERIES COMPARISON AND ERROR DISTRIBUTION PLOTS OUTPUT
        # ggarrange(p.err2, p.err3, p.err1, ncol=2, nrow=2, widths = , common.legend = F)
        ggdraw() +
          draw_plot(p.ind.series, x= 0 , y= 4/5, width=1, height=1/5) +
          draw_plot(p.agg.es, x=0, y=3/5, width=1, height=1/5) +  ## *** DiD ***
          draw_plot(p.bsts.impact, x=0, y=2/5, width=1, height=1/5) +  ## *** BSTS ***
          draw_plot(p.err1, x = 0,  y = 1/5, width = 1, height = 1/5) +
          draw_plot(p.err2, x = 0,  y = 0, width = .5, height = 1/5) +
          draw_plot(p.err3, x = .5, y = 0, width = .5, height = 1/5) +
          draw_plot_label(label = c("A", "B", "C",'D','E','F'), size = 15,
                          x = c(0, 0, 0, 0, 0, .5), y = c(5/5, 4/5, 3/5, 2/5, 1/5, 1/5))
        ggsave(filename = file.path(save.img.dir,
                                    sprintf('%s_ATT_pointwise_error_distribution_compare_n%s_pd%s_ss%s_niter%s_%s_%s_%s.png',
                                            prefix,n,npds,h,bsts.niter,key.strip,effect.type,sim.id)),
               height=15, width=9, units = 'in', dpi = 300)
        
        ##===============================================================
        ## 1-step ahead prediction error
        bsts.pred.er <- bsts.prediction.errors(impact_amount$model$bsts.model)$in.sample[,(1:(intpd-1))]
        
        ## Append results to output list
        simlist[[key]]$compare$did[[effect.type]]$attgt <- ccattgt
        simlist[[key]]$compare$did[[effect.type]]$agg.simple <- agg.simple
        simlist[[key]]$compare$did[[effect.type]]$agg.es <- agg.es
        simlist[[key]]$compare$did[[effect.type]]$self.select.cor <- self.select.cor
        ##
        simlist[[key]]$compare$bsts[[effect.type]][[ h ]]$CausalImpact <- impact_amount
        simlist[[key]]$compare$bsts[[effect.type]][[ h ]]$cumu.pred.error <-  cumsum(colSums(abs(bsts.pred.er)))
        ##
        simlist[[key]]$compare$res.tbl[[effect.type]][[ h ]] <- res.tbl
        simlist[[key]]$compare$att.err.tbl[[effect.type]][[ h ]] <- errdf
        simlist[[key]]$compare$att.err.mean.bsts[[effect.type]][[ h ]] <- mean(errdf$error[errdf$method=='BSTS'],na.rm = T)
        simlist[[key]]$compare$att.err.mean.did[[effect.type]][[ h ]]  <- mean(errdf$error[errdf$method=='DiD'],na.rm = T)
        simlist[[key]]$compare$att.err.sd.bsts[[effect.type]][[ h ]] <- sd(errdf$error[errdf$method=='BSTS'],na.rm = T)
        simlist[[key]]$compare$att.err.sd.did[[effect.type]][[ h ]]  <- sd(errdf$error[errdf$method=='DiD'],na.rm = T)
      
      } ## // end h loop over bsts.state components
      

    } ## // end k loop over effect types
    
    
    if ( ! is.na(save.items.dir) ) {
      ## Save simulation list as serialized data file
      simlist.file <- sprintf('__GRIDSEARCH_output__%s_%s.rds', sim.id, key.strip)
      save.file.path <-  file.path(save.items.dir, simlist.file)
      saveRDS(simlist[[key]], file = save.file.path)
      ## FREE UP MEMORY
      simlist[[key]] <- list(file = save.file.path)
    } 
    
    
  } ## // end simlist loop i   ##  #; dev.off()
  
  return(simlist)

  }


################################################################################################################
################################################################################################################





##=======================================
##
##  STATESPACE CONFIG GRIDSEARCH
##   1. search over priors/args in state space component function
##
##
##=======================================


## Static Defaults
# n <- 100  ## NUmber of actors (i.e., number of timeseries)
# npds <- 60 #100
# intpd <- round( npds * 2/3 )
noise.level <- 1.4
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
ns <- list(100) ## 400
sim.lengths <- list(120)
treat.rules <- list('random')  ## 'below.benchmark'
seasonalities <- list(TRUE)   ## c(TRUE,  FALSE )
prior.sd.scenarios <- list('sd.low')  ## sd.low
## FOCAL CONSTRUCT
dgp.ars <- list(0)  ## 0.6  ## .1,.2,.4
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
  # `6`=c('AddAutoAr'),
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
            
            ##--------------------------------------------
            ## Append simulation configuration to simlist
            ##--------------------------------------------
            key <- sprintf('d%s|f%s|g%s|h%s|i%s|j%s', d,f,g,h,i,j)
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
              # bsts.state.specs=list(list(AddSemilocalLinearTrend),list(AddSemilocalLinearTrend,AddStudentLocalLinearTrend)),
              rand.seed = 7531
            )
            
          } ## // end prior.sd.scenarios loop
          
        } ## // end treat.rules loop
        
      } ## // end seasonalities loop
      
    } ## // end ARs loop
    
  } ## // end nps loop  (sim lengths)
  
} ## // end ns loop (number of actors)


effect.types = c('quadratic')

simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
                               plot.show = FALSE, plot.save = FALSE )


simlist.files <- runSimUpdateCompareBstsDiD(simlist, 
                                            effect.types = effect.types,
                                            save.items.dir= dir_ext,
                                            bsts.niter=1000)  ## D:\\BSTS_external


## LOAD INDIVIDUAL SIMULATION COMPARISON LIST
key <- sprintf('d%sf%sg%sh%si%sj%s', 1,1,1,1,1,1)
simx <- readRDS(file.path(dir_ext,sprintf('__GRIDSEARCH_output__%s_%s.rds',simlist[[1]]$sim$id, key )))

## ORIGINAL DATA 120 pds = Monthly data over 10 years
##----------------------
## Aggregate at 4 periods (quarterly)
##---------------------
# simlist <- runSimUpdateSimlist(simlist, effect.types = effect.types,
#                                plot.show = F, plot.save = F )
simlist <- updateSimlistAggregatePd(simlist, pd.agg = 4 )
simlist.files <- runSimUpdateCompareBstsDiD(simlist, 
                                            effect.types = effect.types,
                                            save.items.dir= dir_ext,
                                            bsts.niter=1000) 
##----------------------
## Aggregate at 12 periods (yearly)
##---------------------
simlist <- updateSimlistAggregatePd(simlist, pd.agg = 12 )
simlist.files <- runSimUpdateCompareBstsDiD(simlist, 
                                            effect.types = effect.types,
                                            save.items.dir= dir_ext,
                                            bsts.niter=30000) 






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
# updateSimlistAggregatePd <- function(simlist, pd.agg, na.rm=TRUE) {
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
