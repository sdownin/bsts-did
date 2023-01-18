##
##   Internal intervention simulation - Vectorized 
##
##     - vectorized actor loop i=1:n  to improve performance
##
##

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(scales)
library(foreach)



##==============================
## Helper Functions
##------------------------------

# logistic <- function(x, shift=1, scale=1) {
#   ## 1  / ( 1    + exp(-x) )
#   scale / (shift + exp(-x))
# }
# 
# getLogisticProportion <- function(x, shift=1, scale=1) {
#   hi <- logistic( 9e8, shift, scale) ## upper bound
#   lo <- logistic(-9e8, shift, scale) ## lower bound
#   rng <- max(hi,lo) - min(hi,lo)     ## range lower->upper
#   (logistic(x) - lo) / rng           ## x proportion of range
# }





##             Treat      Post-int   Treat-post-int   Past Perf.  Growth  noise   covariates (age/type/level)
## y[t] = b0 + b1*x1[t] + b2*x2[t] + b3*x1[t]*x2[t] + b4*y[t-1] + b5*t  + u[t]  + b6*c1[t] + b7*c2[t] + b8*c3[t]
## x1 = rule[[? f(y[t-1])]]


## Two cases to simulate: 
## one time for intervention?
## staggered over time ?


##
## Produce seasonality component for simulation 
##   based on nseasons argument
##
getSinBySeasons <- function(period.values, nseasons, freq=1,
                            noise.mean=0, noise.sd=0, vert.scale=1) {
  freq.scale <- freq * 2*pi / nseasons
  season.effect <- vert.scale * sin( period.values * freq.scale )
  noise <- rnorm(length(period.values), noise.mean, noise.sd)
  return(season.effect + noise)
}


##================================
## MAIN VARIABLE FUNCTIONS
##--------------------------------
## outcome
yFunc <- function(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9, x1, x2, season.val, t, u, c1, c2, c3, y.tm1, localLevel) {
  # b0 + b1*x1 + b2*x2 + b3*x1*x2 + b4*season.val + b5*t + u + b6*c1 + b7*c2 + b8*c3 + b9*y.tm1
  b0 + b1*x1 + b2*x2 + b3*x1*x2 + b4*season.val + b5*localLevel + u + b6*c1 + b7*c2 + b8*c3 + b9*y.tm1
}



## Self-select into treatment as function of past performance
## [char] treat.rule Rule for actor to select into treatment group
##                    'past'=below past performance treat.threshold (focal or all actors)  --> intervention
##                    'random'=Bernoulli draw with p=treat.threshold
##                    'worst'=worst performing actor goes for treatment with p=treat.threshold
##                    'below.benchmark'=actors with past performance below benchmark (eg, median) select into treatment with p=treat.threshold
## 
x1Func <- function(x3.treat.post.int, y.tm1, v, ## lengths must match
                   g0=0.01, g1=0.01, 
                   # c1=1, c2=0, c3=1, ## NOT IN SELF_TREATMENT FUNCTION 
                   treat.rule='below.benchmark', 
                   treat.prob=0.5, ## probability of self-selecting into treatment (from set of actors below treat.threshold of treatment criterion)
                   treat.threshold=0.8, ## treatment threshold criterion - below this quantile, self-select into treatment with probability = treat.prob
                   records=c(), ## dataframe with past records used for 'self' or 'all' treat.rule
                   shift=1, scale=1) {
  
  n <- length(x3.treat.post.int) ### x3.treat.post.int, y.tm1, v, <-- lengths must match
  p <- treat.prob ## treatment self-selection threshold
  ## output vector
  out <- rep(0, n)
  
  ## treatment is monotonically increasing (cannot untreat after being treated)
  idx.tr <- which( x3.treat.post.int == 1 )
  out[idx.tr] <- rep(1, length(idx.tr))
  
  ## TODO: see if transform carries any implications (?)
  ##       comparison of z vs. records should be invariant to linear transform (?)
  z <- g0 + g1*y.tm1 + v
  records <- g0 + g1*records + v

  ## negative z flips the s-curve
  ## so higher performance makes actor
  ## LESS likely (approaching zero)
  ## to take the intervention
  ##  (ie, endogenously self-select into treatment group)
  # .x1.t <- getLogisticProportion( -z, shift, scale )
  # out <- sample(1:0, size = 1, prob = c(.x1.t, 1-.x1.t))
  idx <- c()
  if (treat.rule == 'below.benchmark') {
    benchmark <- quantile(records, probs=treat.threshold, na.rm=TRUE)
    if (length(benchmark) == 0 | is.na(benchmark)){
      cat("\nmissing benchmark in treat.type=bottom.half. Records return NULL benchmark\n")
      cat(records)
    }
    idx <- which( z < benchmark & out == 0 )
    
    
  } else if (treat.rule == 'past') {
    ## **NOTE** MUST BE common records across all firms
    benchmark <- quantile(records, probs=treat.threshold, na.rm=TRUE)
    
    if (length(benchmark) == 0 | is.na(benchmark)){
      print("missing benchmark. records return NULL threshold quantile")
      print(records)
    }
    idx <- which( z < benchmark & out == 0 )
    
  } else if (treat.rule == 'worst') {
  
    if (length(records) > 0) {
      idx <- which( z == min(records) & out == 0 )
    }
    
  } else if (treat.rule == 'random') {
    
    idx <- which( out == 0 )

  } else {
    print(sprintf('treatment rule "%s" not known',treat.rule))
    return()
  }
  
  ## Update output x1 vector
  if (length(idx) > 0)   {
    ##  for those actors who were not treated (out[t-1]==0)
    out[idx] <- sample(0:1, size = length(idx), prob = c(1-p, p), replace = TRUE)
  }
  
  # cat(sprintf('prop=%.3f \n\n',out))
  return(out)
} 


## post-intervention dummy
x2Func <- function(t, intpd, n) {
  x2 <- ifelse( t >= intpd, 1, 0)
  rep(x2, n)
} 

## Treatment effect
b3Func <- function(t.post.intpd, x1, x2, ## must be same length vectors
                   type='quadratic', 
                   w0=1.0, w1=.3, w2=-.012, #w2.shift=0,
                   non.negative=TRUE) {
  
  n <- length(t.post.intpd) ## t.post.intpd, x2 <-- must be same length vectors
  # print('b3Func() t.post.intpd')
  # print(t.post.intpd)
  # print('b3Func() n')
  # print(n)
  t <- t.post.intpd ## time period 't' is number of periods after intervention
  
  if (length(x2) != n){
    cat(sprintf('\nb3Func() length of x2 (%s) != length of t.post.intpd (%s)\n',length(x2), n))
    return()
  }

  ## init default treatment effect values as zeros
  b3 <- rep(0, n)
  
  ## if pre-treatment period (x==0) then treatment effect is defined as 0
  idx <- which( x1 * x2  == 1 )
  n.idx <- length(idx)
  
  if (length(n.idx) > 0)
  {
    ## Compute treatment effect by period and dynamic shape
    ## new replacement values vector length == n.idx
    if(type == 'constant') {
      new <- rep(w0, n.idx)  ## length == n.idx
    }
    else if (type == 'linear') {
      new <- w1*t[idx] + w0    ## length == n.idx
    }
    else if (type == 'quadratic') {
      # cat('b3Func() t[idx]', t[idx])
      # cat('\n')
      tx <- t[idx] # + w2.shift ## ? TODO: check sensitivity for gradual treatment effect increase (negative number shifts curve rightward)
      new <-  w2*(tx^2) + w1*tx + 0 ## replace w0 with 0    ## length == n.idx
    }
    else if (type == 'geometric') {
      denom <- t[idx]  ## length == n.idx
      denom[denom==0] <- 1  ## prevent division by zero
      geom.effect.multiplier <- 7  ## how much larger is geom shape initial effect than constant effect level
      new <- geom.effect.multiplier * w0 / denom
    }
    else {
      cat(sprintf('\neffect type %s not known\n',type))
      return()
    }
    
    ## udpate new dynamic treatment effects for period t
    b3[idx] <- new
  }
  
  # cat(sprintf('\nb3Func: b3 before max() = %.3f\n',b3))
  if (non.negative) {
    b3[b3<0] <- 0
  }
  
  # print('b3Func() FINISHED.')
  
  return(b3)
}

##=================================
##
## Internal Intervention Simulation Function
## 
##=================================
runSimSingleIntervention <- function(
    n, ##200, ## NUMBER OF actors  
    npds, ##520, ## NUMBER OF PERIODS
    intpd, ##round( 520 * (5/6) ),
    ystart, 
    effect.type,##'constant',  ## c('constant','quadratic','geometric'), ## treatment effect shapes
    treat.rule, ##'below.benchmark',
    treat.prob, ##0.5, ## probability of self-selecting into treatment, if below treat.threshold 
    treat.threshold, ##0.8, ## threshold of treatment criterion, below which actor self-selects into treatment with probability = treat.prob
    benchmark.type, ##'self', ## 'all', 'self'
    rand.seed, ##54321,  ## pseudo-random number generator seed for replication
    ##
    sim.id, ## Index or timestamp of simulation for  saving, etc.
    ##
    noise.level, #1,
    ## # PEFORMANCE [Y] FUNCTION PARAMETERS
    b0, ##.001, ## intercept
    b1, ##.001, ## treatment dummy
    b2, ##.001, ## post intervention dummy
    # b3 = .001, ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
    b4, ##0, ## Weight of seasonality effect
    # b5, ##.04, ## LEVEL: Random Walk 
    ## Covariates
    # b6, ##.2, ##  c1  ## small variance (.01) random walk
    # b7, ##.2, ##  c2  ## large variance (.1) random walk
    # b8, ##.2, ##  c3  ## gamma var
    ## Autocorrelation
    b9, ##0,
    ## # TREATMENT EFFECT FUNCTION WEIGHTS 
    w0, ##1.2, ## constant
    w1, ##0.3, ## linear
    w2, ##-0.12, ## quadratic
    # w2.shift, ##-12, ## shift quadratic curve to allow gradual treatment effect increase (negative shifts curve rightward)
    ## SEASONALITY 
    nseasons, ##  NA omits seasonality; integer values incluce seasonality with nseasons
    season.frequency, ##1, ## number of completed sin waves within 1 cycle of nseasons
    ## # ENDOGENOUS TREATMENT SELECTION [x1[t](y[t-1])] FUNCTION PARAMETERS
    g0 = .1,
    g1 = .1,
    ## # EXPONENTIAL FUNCTION PARAMETERS -- CURRENTLY NOT USED
    logit.shift = 1,
    logit.scale = 1,
    # ## PLOTTING
    plot.save = TRUE,
    plot.show = TRUE,
    plot.wide.w=10,
    plot.wide.h=7,
    plot.tall.w=7,
    plot.tall.h=10,
    plot.sq.w=9,
    plot.sq.h=9,
    plot.dpi=300,
    verbose=TRUE,
    dgp.prior.sd.weight=.01,
    ...
  ) {
  
  # print('DEBUG runSimSingleIntervention()::')
  # # print(as.list(match.call()))
  # print(c(as.list(environment()), list(...)))
  
  ## --- BEGIN SIM ---
  
  SIMID <- ifelse(is.na(sim.id), round(10*as.numeric(Sys.time())), sim.id)
  
  ntypes <- 1 ## length(effect.types)
  
  sim.count <- ntypes * npds  ## n
  
  if(verbose) cat(sprintf('\nRunning Single Intervention Simulation for %s iterations:\n  %s periods (vec.len.=%s),  effect type = %s\n',
              sim.count, npds, n, effect.type))
  
  ## Initialize progress bar
  counter <- 0
  if(verbose) {
      pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = sim.count, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  }

  
  set.seed(rand.seed)
  
  ##================================
  ## MAIN RUN
  ##--------------------------------
  ## VARIABLES IN SIMULATION PER actor PER PERIOD
  # y <- NA
  # x1 <- NA
  # x2 <- NA
  # b3 <- NA
  # u <- NA
  # v <- NA
  ##------
  ## OUTPUT OBJECTS (to be output as list)
  df <- data.frame(stringsAsFactors = F)
  # df.t0 <- NA
  # df.summary=NA
  # df.t0.summary=NA
  # df.plot.treat.map=NA
  # plot.allseries=NA
  # plot.group=NA
  # plot.group.staggered=NA
  
  
  for (t in 1:npds) ## PERIOD LOOP
  # for (t in 2:23)
  {
    # cat(sprintf('\n\n\n----- t = %s -----\n\n\n',t))
    
    ## TODO: Vectorize the actor loop to improve performance (slow loops in R)
    # for (i in 1:n) ## actors LOOP
    # {
      # cat(sprintf(' i = %s \n',i))
    
    ## noise terms
    v <- rnorm(n=n,mean=0, sd=noise.level )  ##- .5 ##
    u <- rnorm(n=n,mean=0, sd=noise.level ) ## 
    
    ##----------  Past Period Variables---------------------
    ## indices of last period observations
    idx.tm1 <- if (t==1) { c() } else { which( df$t == t-1 & df$effect.type==effect.type ) }
    ## treated dummy (post-intervention-treatment-effect)
    x3.tm1 <- if (t==1) { rep(0, n) } else { df$x3[idx.tm1] }
    ## focal actor past performance
    y.tm1 <- if (t==1) { rep(ystart, n) } else { df$y[idx.tm1] }
    ###
    
    ## State Space local level 
    # ## LEVEL TERM RANDOM WALK
    localLevel.tm1 <- if(t==1) {rep(0, n)} else { mean(df$localLevel[idx.tm1], na.rm=T) } 
    
    ## DYNAMIC REGRESSION BETAS --------------
    ##   (same for all n; using mean( cov[t-1] ) should be same as selection first item [1])
    ## FULL SEASONALITY COMPONENT IS INCLUDED IN OUTCOME ( b5 := 1 )
    ## cov_season (seasonality = sinusoid with 52 pds, 1 frequency)
    b4.tm1 <- if (t==1) { 1 } else { mean(df$b4[idx.tm1], na.rm=T) }
    ## OTHER COVARIATES ----------------------
    ## c1  
    b6.tm1 <- if (t==1) { 1 } else { mean(df$b6[idx.tm1], na.rm=T) }
    ## c2  
    b7.tm1 <- if (t==1) { 1 } else { mean(df$b7[idx.tm1], na.rm=T) }
    ## c3  Temporal drift  (noise on a function of time)
    b8.tm1 <- if (t==1) { 1 } else { mean(df$b8[idx.tm1], na.rm=T) }
    
    ##------------------------------------------------------
    ## Get past performance records for treat.rule == 'past'
    ##------------------------------------------------------
    # perf.records <- if (i==1) { ## if first actor (df not updated yet)
    #   if (t==1){ ystart } else { y } ## first actor must use default start value of performance or past performance in 'y' vector
    # } else {
    #   if (t==1) {
    #     switch(benchmark.type,
    #            "self"=ystart,  ## focal actor i not in data.frame yet when t==1, so just use ystart
    #            "all"=df$y)
    #   } else {
    #     switch(benchmark.type,
    #            "self" = df$y[which(df$actor == i)],
    #            "all" = df$y)
    #   }
    # }
    ##------------------------------------------------
    ## TODO: Fix past performance criterion for treatment self-selection 
    # len.memory <- 3
    # .ts <- sapply(1:len.memory, function(a)max(t-a,0))
    # t.in.memory <- .ts[.ts > 0]
    # perf.records <- if (t==1) { ## if first actor record
    #   switch(benchmark.type,
    #          "self" = ystart, ## first time for actor i must use default start value of performance or past performance in 'y' vector
    #          "all" = df$y)
    # } else {
    #   switch(benchmark.type,
    #          "self" = df$y[which(df$actor == i & df$t %in% t.in.memory)],
    #          "all" = df$y[which(df$t %in% t.in.memory)])
    # }
    perf.records <- y.tm1
    # print('perf.recrods')
    # print(perf.records)
    # Treatment dummy (actor self-select into treatment based on treat.rule)

    ## Intervention period
    if ( t == 1 ) {
      
        x1 <- rep(0,n)
        
    } else if (t == intpd) {
      ## ACTORS SELECT INTO TREATMENT GROUP BY CRITERIA
      x1 <- x1Func(x3.tm1, y.tm1, v,
                   g0, g1, 
                   # c1, c2, c3,  ## NOT IN SELF_TREATMENT FUNCTION 
                   treat.rule = treat.rule, 
                   treat.prob = treat.prob,
                   treat.threshold = treat.threshold, #0.5 for 'worst'; 0.1 for 'past' (?)
                   records=perf.records
      )
      
    } else {
      ## Keep previous period
      x1 <- df$x1[ which( df$t == (t-1) ) ]
    }
    

    
    ## post-intervention dummy var
    x2 <- x2Func(t, intpd, n)
    
    # ## Treatment Effect
    # x3 <- x1 * x2

    
    
    ## DYNAMIC TREATMENT EFFECT
    ## periods after intervention
    t.post.intpd <- rep(0, n)
    if( length(idx.tm1) > 0 ) { ## after first period, can take previous period values
      df.tm1 <- df[idx.tm1, ]
      tm1.post.intpd <- df.tm1$t.post.intpd
      # idx.tr <- which( df.tm1$x3 == 1 ) ## has been treated already
      idx.tr <- which( x1*x2 == 1 ) ## has been treated by this period
      tm1.post.intpd[idx.tr] <- tm1.post.intpd[idx.tr] + 1 ## increment treated actors' post-intervention period by 1
      # cat('tm1.post.intpd: ', tm1.post.intpd)
      t.post.intpd <- tm1.post.intpd  # set this year as updated values from last year
    }
    # cat('t.post.intpd: ', t.post.intpd)
    ##
    # shift.t.post.intpd <- t.post.intpd - intpd
    b3 <- b3Func(t.post.intpd, x1, x2, type=effect.type,
                 w0=w0, w1=w1, w2=w2) ## w2.shift=w2.shift
    
    
    ##-------------------------------------------
    ## ALL VALUES ARE THE SAME (indexed only at time; not by individual )
    # localLevel.mu <- mean(localLevel.tm1, na.rm=T) ## should be the same as taking item [1]
    localLevel <- rnorm(n, localLevel.tm1, sd= noise.level * dgp.prior.sd.weight )  ## save level value u[t] for all actors (not u[i,t])
    
    
    ##-------------------------------------------
    ## Seasonal Component and scaling yearly growth to period-growth
    if ( is.null(nseasons) | is.null(season.frequency) | is.na(nseasons) | is.na(season.frequency) ) {
      season.val <- 0
      ## *** DEFAULT SETTING assumes yearly growth rate = (b5 / (pds/yr)), where pds/yr = 52 weeks
      ##     That is, the simulation assumes weekly periods for growth rate, even if seasonality component is missing from DGP.
      b4 <-  rnorm(1, ( b4.tm1 / (npds / 52)  ), sd = noise.level * .5 )
    } else {
      season.vals <- getSinBySeasons(1:npds, nseasons, freq=season.frequency,
                                     noise.mean=0, noise.sd = 0, # add noise below
                                     vert.scale =  1.5  ) 
      season.val <- season.vals[t]
      ## Scale linear growth to equal b5 value for every completed seasonal cycle (e.g., 1 year)
      # season.frequency <- ifelse(season.frequency == 0, 1, season.frequency)
      # nseasons <- ifelse( nseasons == 0 , 1, nseasons )
      # growth.scale <-  npds /  ( nseasons * season.frequency )
      # b5 <- rnorm(1, (b5.tm1 / growth.scale), sd = noise.level * 0.5 )
      b4 <- rnorm(1, season.val, sd = noise.level * .5 )
    }
    
    ## Weight of local level
    b5 <- 1 ## default to unit weight of the local level (see Brodersen et al 2015)
    # b5 <- rnorm(1, b6.tm1, sd = noise.level * dgp.prior.sd.weight)
    ## Weight of covariates: Dynamic covariate betas (see Brodersen et al 2015)
    b6 <- rnorm(1, b6.tm1, sd = noise.level * dgp.prior.sd.weight)
    b7 <- rnorm(1, b7.tm1, sd = noise.level * dgp.prior.sd.weight)
    b8 <- rnorm(1, b8.tm1, sd = noise.level * dgp.prior.sd.weight)
    
    
    ## --------------- Covariates -------------------------
    ## INITIAL VALUES OF COVARIATES
    c1.tm1 <- if (t==1) { rep(.1, n) } else { df$c1[idx.tm1] }
    ## START HIGHER = 10
    c2.tm1 <- if (t==1) { rep(.1, n) } else { df$c2[idx.tm1] }
    ## 
    c3.tm1 <- if (t==1) { rep(.1, n) } else { df$c3[idx.tm1] }
    
    # ## RANDOM WALKW WITH DRIFT (noise in the local level)
    # # c1 <- rpois(n, lambda = noise.level*0.8) + 1
    # # c1.tm1.drifted.mean <- rnorm(n, c1.tm1, sd = noise.level * 0.1 )
    # # c1.tm1.drifted.mean <- rnorm(n, c1.tm1, sd = 0 )  ##***CHANGED***
    # # c1 <- rnorm(n, c1.tm1, sd = noise.level * 0.1 )
    # c1 <- rnorm(n, .0025*t, sd = noise.level * .2)
    # 
    # ### COVARIATE SERIES 2
    # ### OPTION A
    # # c2 <- rnorm(n, -c2.tm1, sd=noise.level * 0.5 )
    # ## OPTION B - direct comparison with different correlations with y (outcome affects by b7, y(b7) )
    # # c2.tm1.drifted.mean <- rnorm(n, c2.tm1, sd = noise.level * 0.1 )
    # # c2 <- rnorm(n, c2.tm1, sd = noise.level * 0.01 )
    # # c2 <- rnorm(n,  c2.tm1, sd = noise.level * 0.1 )
    # ######
    # c2.sin.vals <- getSinBySeasons(1:npds,  nseasons = 52,
    #                                freq= ifelse(is.na(season.frequency),1, season.frequency * 8),
    #                                noise.mean=0, noise.sd = 0, # add noise below
    #                                vert.scale =  runif(1, .01, .1)  )
    # .id <- ifelse(t<2, 1, t-1)
    # c2 <- rnorm(n, c2.sin.vals[.id], sd = noise.level * 0.1)
    # 
    # # # ## COVARIATE SERIES 3
    # # ## OPTION A
    # # c3 <- rgamma(n, shape = log(t + 1) + .1, scale = noise.level * 0.1 )
    # # ## OPTION B - direct comparison with different correlations with y (outcome affects by b8, y(b8) )
    # # c3.tm1.drifted.mean <- rnorm(n, c3.tm1, sd = noise.level * 0.1 )
    # # c3 <- rnorm(n, -.001*(t/2), sd = noise.level * 0.1 )
    # c3 <- rnorm(n, c3.tm1, sd = noise.level )
    # # c3 <- rpois(n, lambda = max(.1, c3.tm1)  )
    
    ##**MULTIVARIATE NORMAL CORRELATED RANDOM NOISE COVARIATES**
    c1 <- rnorm(n, .5, noise.level * 1)
    c2 <- rnorm(n, .02 * t, noise.level * 1.5)
    c3 <- rnorm(n, -.01 * t, noise.level * .5)
    ## 
    sig.mat <- matrix( c( 1,.3,.1,
                         .3, 1,.2,
                         .1,.2, 1), ncol=3, byrow = T)
    mu.vec <- c(.1, .2, .15)
    rmv.mat <- mvtnorm::rmvnorm(n, 
                                mean = mu.vec, 
                                sigma = sig.mat) 
    # correlated random variables
    Cmvt <- cbind(c1, c2, c3) %*% sig.mat
    
    ##-------------------------------------
    ## PERFORMANCE
    y <- yFunc(b0=b0, 
               b1=b1, 
               b2=b2, 
               b3=b3, 
               b4=b4, 
               b5=b5, 
               b6=b6, 
               b7=b7, 
               b8=b8, 
               b9=b9, 
               x1=x1, 
               x2=x2, 
               season.val=season.val, 
               t=t, 
               u=u, 
               c1=Cmvt[,1],
               c2=Cmvt[,2], 
               c3=Cmvt[,3],
               y.tm1=y.tm1,
               localLevel=localLevel
              )
    ##--------------------------------------
  
    
    ##--------------------------------------
    ## UPDATE PERIOD DATAFRAME
    ##--------------------------------------
    df.t <- data.frame(actor=1:n, t=rep(t, n), 
                       t.post.intpd=t.post.intpd,
                       effect.type=rep(effect.type, n),
                       y=y, x1=x1, x2=x2, x3=x1*x2,
                       season.val=season.val,
                       b1=rep(b1, n), 
                       b2=rep(b2, n), 
                       b3=b3, ## DYNAMIC TREATMENT EFFECT by GROUP (differs between treatment, control; needs time-varying covariates storage)
                       b4=rep(b4, n), 
                       b5=rep(b5, n),
                       b6=rep(b6, n),
                       b7=rep(b7, n),
                       b8=rep(b8, n),
                       localLevel=localLevel,
                       c1=c1, c2=c2, c3=c3,
                       u=u, v=v,
                       stringsAsFactors = F)
    
    if (nrow(df.t) != n) {
      cat(sprintf('nrow(df.t) != n: %s != %s\n', nrow(df.t), n))
      break
    }
    
    ## APPEND TO OUTPUT DATAFRAME
    df <- rbind(df, df.t)
    
    ## Update Progress Bar
    ## starting at 0 --> add 1 after completing first item / before progress bar update
    counter <- counter + 1 
    setTxtProgressBar(pb, counter)
      
  } ## /end period loop
    
    
  # }
  
  ## DEBUG 
  df <- unique(df)
  df$actor <- as.numeric(as.character(df$actor)) ## avoid possible factor issue (?)
  
  ## Which actors are treated 
  idx.x3.eq1 <- which(df$x3==1)
  if (length(idx.x3.eq1) == 0) {
    cat(sprintf('\n\nWARNING: No treated actors. All actors have x3[t]==0 for t=%s\n\n',t))
  }
  
  
  ## COHORT GROUP ASSIGNMENT 
  ## ASSIGN DEFAULT GROUP AND COLOR (before updating treated cohort)
  df$group <- 'control'
  ggcolors <- hue_pal()(2)
  df$group.color <- ggcolors[1] ## control red color
  ## Assign cohort group (treatment, control)
  ## based on rule: treatment=treated by end; control=not treated by end
  actor.eff.tr <- sort(unique(df$actor[idx.x3.eq1]))  ## & df$effect.type==effect.type
  
  if (length(actor.eff.tr)==0) {
    cat(sprintf('\n\nWARNING: zero actors have x3[t]==1 for t=%s\n\n',t))
  }
  
  for (f in actor.eff.tr) {
    # cat(' ',f)
    ## all observation indices for actor f with effect type=effect.type
    idx.act.eff.tr <- which( df$actor==f)  ## & df$effect.type==effect.type
    ## if actor is treated at any time (any x3=0) then apply treatment group to all actor observations (periods) for this effect type
    df$group[idx.act.eff.tr] <- 'treatment'
    df$group.color[idx.act.eff.tr] <- ggcolors[2]  ## treatment green color
  }
   
  
  # ## Matching Treatment subject to Control subject
  # ##  match control subject according to:
  # ##  - unmatched control subject at same time period as treated subject in treatment period
  # ##  - if multiple unmatched control subjects, match on covariates (c1,c2,c3)
  # # a.trea <- unique(df$actor[df$group=='treatment' & df$effect.type==effect.type])
  # ## Simple random selection of unmatched control group members by period
  # a.ctrl <- sort(unique(df$actor[df$group=='control']))  ## & df$effect.type==effect.type
  # df.ctrl <- df[which(df$actor %in% a.ctrl), ]
  # min.idx <- if(length(a.ctrl)<length(actor.eff.tr)) {a.ctrl} else {actor.eff.tr}##min(c(length(actor.eff.tr),length(a.ctrl)))
  # matches <- data.frame(treat=actor.eff.tr, control=NA, stringsAsFactors = F) 
  df$match_id <- NA
  df$match_pd <- NA
  if (length(actor.eff.tr)>0)
  {
    for (i in 1:length(actor.eff.tr)) {
      # cat(sprintf('%s ',i))
      actor <- actor.eff.tr[i]
      
      ## Unmatched controls
      um.ctrl <- unique(df$actor[which(df$group=='control' & is.na(df$match_id))]) ## & df$effect.type==effect.type
      if (length(um.ctrl)==0) {
        if(verbose) cat(sprintf(' skipping i=%s: actor=%s\n',i,actor.eff.tr[i]))
        next
      }
      
      pd.tr <- min(df$t[which(df$t.post.intpd>0 & df$actor==actor)])  ## & df$effect.type==effect.type
      
      a.ctrl <- sample(um.ctrl, size = 1, replace = T)
      ## MATCH ID:  effect_pd__treated_control
      match_id <- sprintf('%s_%s__%s_%s',effect.type,pd.tr,actor,a.ctrl) ##paste(c(effect.type,actor,pd.tr), collapse = '_')
      
      idx.match <- which(df$actor %in% c(actor,a.ctrl))
      df$match_id[ idx.match ] <- match_id   # &  df$effect.type==effect.type
      df$match_pd[ idx.match ] <- pd.tr
      
      # batch <- df[which(df$actor==idx & df$effect.type==effect.type), ]
      # grp <- 'treatment' #unique(batch$group) ## should only be one element as all rows should be from same subject in 1 group
      
      # matches <- rbind(matches, data.frame(treat=actor, control=.control))
    }
  }

  ## close progress bar
  close(pb)
  
  
  
  # ##--------------------------------
  # ## STAGGERED DID ARRANEMENT:
  # ##--------------------------------
  # ## Reset each series period as t=0 when intervention selected
  # df.t0 <- df
  # df.t0$t0 <- NA
  # if (length(actor.eff.tr)>0)
  # {
  #   for (i in 1:length(actor.eff.tr)) {
  #     actor <- actor.eff.tr[i]
  #     actor.match.id <- unique( df$match_id[which(df$actor==actor)] )
  #     actor.ctrl <- unique( df$actor[which(df$match_id==actor.match.id & df$actor!=actor)] )
  #     ##
  #     x3.idx <- which(df.t0$actor==actor  & df.t0$x3 > 0)   ## & df.t0$effect.type==effect.type
  #     actor.x3.t <- ifelse(length(x3.idx)>0, min(df.t0$t[x3.idx]), NA)
  #     isTreated <- length(actor.x3.t) > 0 & any( ! is.na(actor.x3.t))
  #     ## apply staggered timing to both treated actor and corresponding matched control (actor with same match_id)
  #     df.t0$t0[which(df.t0$actor==actor)]      <- if (isTreated) { (1:npds) - actor.x3.t + 1 } else { (1:npds) - intpd + 1 }
  #     df.t0$t0[which(df.t0$actor==actor.ctrl)] <- if (isTreated) { (1:npds) - actor.x3.t + 1 } else { (1:npds) - intpd + 1 }
  #   } 
  # }
  # if (length(df.t0$t0)==1 & is.na(df.t0$t0[1])){
  #   cat(sprintf('\n\nWARNING: df.t0 has zero rows (only default t0=NA)\n\n'))
  # }
  # # for (f in 1:n) {
  # #   x3.idx <- which(df.t0$actor==f  & df.t0$x3 > 0)   ## & df.t0$effect.type==effect.type
  # #   f.x3.t <- ifelse(length(x3.idx)>0, min(df.t0$t[x3.idx]), NA)
  # #   df.t0$t0[which(df.t0$actor==f)] <- if (length(f.x3.t)==0 | all(is.na(f.x3.t))) {   ## & df.t0$effect.type==effect.type
  # #     (1:npds) - intpd + 1
  # #   } else {
  # #     (1:npds) - f.x3.t + 1 ## t0 = f.x3.t period
  # #   }
  # # }
  

  ##============================
  ## OUTPUT LIST
  ##----------------------------
  return(list(
    df=df,
    # df.t0=df.t0,
    # df.summary=df.summary,
    # df.t0.summary=df.t0.summary,
    # df.plot=df.plot,
    # plot.allseries=p1,
    # plot.group=p2,
    # plot.group.staggered=p3,
    id=SIMID
  ))

}




# 
# ##### DEBUG ###########
# n = 100 ## NUMBER OF actors
# npds = 80 ## NUMBER OF PERIODS
# intpd = round( npds / 4 ) ## intervention after first fifth
# ystart = 0
# effect.types = c('constant','quadratic','geometric') ## treatment effect shapes
# treat.rule = 'below.benchmark'
# treat.threshold = 0.025 ## treatment threshold (either probability or quantile of performance below which low performer(s) self-select into treatment)
# benchmark.type = 'self' ## 'all', 'self'
# seed = 54321  ## pseudo-random number generator seed for replication
# ##
# sim.id = "DEBUG_SIMID" ## Index or timestamp of simulation for  saving etc.
# ## # PEFORMANCE [Y] FUNCTION PARAMETERS
# b0 = .001 ## intercept
# b1 = .001 ## treatment dummy
# b2 = .001 ## post intervention dummy
# # b3 = .001 ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
# b4 = 1/3 ## spillover of past performance on current performance (how much of treatment effect persists across periods)
# b5 = .01 ## growth rate (linear effect of time: proportion of time t added to linear combination in yFunc() performance )
# ## Covariates
# b6 = 1 ## Age [123...]
# b7 = 0 ## type [01]
# b8 = 1 ## level [012]
# ## # TREATMENT EFFECT FUNCTION WEIGHTS
# w0 = 1.3 ## constant
# w1 = 0.3 ## linear
# w2 = -0.12 ## quadratic
# w2.shift = -10 ## shift quadratic curve to allow gradual treatment effect increase (negative shifts curve rightward)
# ## SEASONALITY
# nseasons = NA ##  NA omits seasonality; integer values incluce seasonality with nseasons
# season.frequency = 1 ## number of completed sin waves within 1 cycle of nseasons
# ## # ENDOGENOUS TREATMENT SELECTION [x1[t](y[t-1])] FUNCTION PARAMETERS
# g0 = .1
# g1 = .1
# ## # EXPONENTIAL FUNCTION PARAMETERS -- CURRENTLY NOT USED
# logit.shift = 1
# logit.scale = 1
# # ## PLOTTING
# plot.save = TRUE
# plot.show = TRUE
# plot.wide.w=10
# plot.wide.h=5
# plot.tall.w=5
# plot.tall.h=10
# plot.sq.w=9
# plot.sq.h=9
# plot.dpi=300
# ############



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



###
## Internal Intervention Simulation Treatment Effect Comparison
## Run Internal Intervention Simulation over Multiple Treatment Effects
##  to output results and plots 
###
runSimSingleInterventionEffectComparison <- function(
    ## Effect types vector
    effect.types = c('constant','quadratic','geometric'),
    ## Simulation base settings
    n = 100, ## NUMBER OF actors  
    npds = 520, ## NUMBER OF PERIODS
    intpd = 520 * 0.6,
    ystart = 0,
    treat.rule = 'below.benchmark',
    treat.prob = 1.0, ## probability of self-selecting into treatment, if below treat.threshold 
    treat.threshold = 0.5, ## threshold of treatment criterion, below which actor self-selects into treatment with probability = treat.prob
    # cov.scenario= list(c1=.1,c2=.2,c3=.3), ## list of vectors of covariances of c1,c2,c3 with outcome
    benchmark.type = 'self', ## 'all', 'self'
    rand.seed = 54321,  ## pseudo-random number generator seed for replication
    ##
    sim.id = NA, ## Index or timestamp of simulation for  saving, etc.
    ## # PEFORMANCE [Y] FUNCTION PARAMETERS
    b0 = .01, ## intercept
    b1 = .01, ## treatment dummy
    b2 = .01, ## post intervention dummy
    # b3 = .001, ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
    b4 = 1, ## 
    b5 = 1, ## 
    ## Covariates 
    # b6 = .1, ## c1  
    # b7 = .2, ## c2  
    # b8 = .3, ## c3  
    ## Autocorrelation
    b9 = 0,
    ## # TREATMENT EFFECT FUNCTION WEIGHTS 
    w0 = 1.5, ## constant
    w1 = 0.13, ## linear
    w2 =  -0.007/sqrt(520), ## = -0.01/sqrt(npds) , ## quadratic
    # w2.shift = -16, ## =  -round( sqrt(npds)*.7 ) shift quadratic curve to allow gradual treatment effect increase (negative shifts curve rightward)
    ## SEASONALITY 
    nseasons = NA, ##  NA omits seasonality; integer values incluce seasonality with nseasons
    season.frequency = 1, ## number of completed sin waves within 1 cycle of nseasons
    ## # ENDOGENOUS TREATMENT SELECTION [x1[t](y[t-1])] FUNCTION PARAMETERS
    g0 = .1,   ## TODO CHECK SENSITIVITY TO THESE PARAMS
    g1 = .1,
    ## # EXPONENTIAL FUNCTION PARAMETERS -- CURRENTLY NOT USED
    logit.shift = 1,
    logit.scale = 1,
    # ## PLOTTING
    plot.save = TRUE,
    plot.show = TRUE,
    plot.wide.w=10,
    plot.wide.h=7,
    plot.tall.w=7,
    plot.tall.h=10,
    plot.sq.w=9,
    plot.sq.h=9,
    plot.dpi=300,
    ## 
    verbose=TRUE,
    dgp.prior.sd.weight=.01,
    ...
  ) {   #effect.types, sim.id, plot.show, plot.save, 
  
  SIMID <- ifelse(is.na(sim.id), round(10*as.numeric(Sys.time())), sim.id)
  
  ##========================================
  ## Loop Simulation over effect types list
  ##----------------------------------------
  df <- data.frame()
  # df.t0 <- data.frame()
  for (k in 1:length(effect.types)) {
    cat(sprintf(' k=%s ',k))
    effect.type <- effect.types[k]
    ## Main Simulation for one effect type
    out <- runSimSingleIntervention(
      n = n, 
      npds = npds, 
      intpd = intpd, 
      ystart = ystart, 
      effect.type = effect.type, 
      treat.rule = treat.rule, 
      treat.prob = treat.prob,
      treat.threshold = treat.threshold, 
      benchmark.type = benchmark.type, 
      rand.seed = rand.seed, 
      sim.id=SIMID, 
      # cov.scenario=cov.scenario,
      b0=b0, 
      b1=b1, 
      b2=b2,     
      b4=b4, 
      b5=b5, 
      # b6=cov.scenario[['c1']], ## covariate c1 weight in y (correlation with outcome)
      # b7=cov.scenario[['c2']], ## covariate c2 weight in y (correlation with outcome)
      # b8=cov.scenario[['c3']], ## covariate c3 weight in y (correlation with outcome)
      b9=b9,
      w0=w0, 
      w1=w1, 
      w2=w2, 
      # w2.shift=w2.shift, 
      nseasons=nseasons, 
      season.frequency=season.frequency,
      g0=g0, 
      g1=g1, 
      logit.shift=logit.shift, 
      logit.scale=logit.scale,
      plot.save=plot.save, 
      plot.show=plot.show, 
      plot.wide.w=plot.wide.w, 
      plot.wide.h=plot.wide.h, 
      plot.tall.w=plot.tall.w, 
      plot.tall.h=plot.tall.h, 
      plot.sq.w=plot.sq.w, 
      plot.sq.h=plot.sq.h, 
      plot.dpi=plot.dpi,
      verbose=verbose,
      dgp.prior.sd.weight=dgp.prior.sd.weight,
       ...
    )
    ##
    df <- rbind(df, out$df)
    # df.t0 <- rbind(df.t0, out$df.t0)
  }
  
  # ## ** DEBUG **
  # return(list(
  #   df=df,
  #   df.t0=df.t0,
  #   # df.plot=df.plot,
  #   # plot.allseries=p1,
  #   # plot.group=p2,
  #   # plot.group.staggered=p3,
  #   id=SIMID
  # ))
  
  df.summary <- NA
  # df.t0.summary <- NA
  df.att <- NA
  df.plot <- NA
  plot.allseries <- NA
  plot.group <- NA
  # plot.group.staggered <- NA
  p1 <- NA
  p2 <- NA
  # p3 <- NA
  
  ##=============================================
  ## PLOTS
  ##---------------------------------------------
  if (plot.show | plot.save)
  {
    PLOTID <- round(10*as.numeric(Sys.time()))
    
    if(verbose) cat('\n rendering plots...')
    ## firm-period observations count to find proportion of treated firms by effect.type
    cnt <- ddply(df, .(effect.type), summarize,
                 count.tr=sum(group=='treatment'),
                 count.ct=sum(group=='control'),
                 prop.tr=round(sum(group=='treatment')/length(actor), 2))
    plot.main <- sprintf('past perf. = %.2f; growth = %.2f,  treat.threshold = %.2f;  prop.treat = %s',
                         b4, b5, treat.threshold, paste(cnt$prop.tr,collapse = '|'))
    
    
    ## default ggplot2 colors
    # ggcolors <- c("#F8766D", "#00BA38")
    ggcolors <- hue_pal()(2)
    
    
    ## ----- 1. ALL INDIVIDUAL TIME SERIES ------------------
    # color <- rgb(.3,.3,.8,.7)
    df.plot <- df[which(df$actor %in% 1:n),] ## filter actors
    df.plot$actor <- as.factor(df.plot$actor)
    p1 <- ggplot(data=df.plot, mapping=aes(x=t,y=y, color=group, group=actor)) +
      geom_line(size=1.05, alpha=0.3) +
      # geom_point(, color=rgb(0.2, 0.2, 0.8, 0.1))+
      geom_hline(yintercept=0) + facet_grid( effect.type ~ . ) +
      theme_bw() + theme(legend.position='none') + ggtitle(plot.main) #+
    # guides(color=guide_legend(nrow=3,byrow=TRUE)) #+
    #
    # facet_wrap(.~factor(df$x1) + factor(df$x2))
    if(plot.show) print(p1)
    if(plot.save) {
      p1.file <- sprintf('internal_intervention_all_actor_timeseries_%s_%s.png',SIMID,PLOTID)
      ggsave(filename=file.path(dir_plot, p1.file), plot=p1,
             width=plot.tall.w,heigh=plot.tall.h,dpi=plot.dpi,units='in')
    }
    
    
    ## ---- 2. GROUP SUMMARY TIME SERIES -----------------------
    df.summary <- ddply(df, .(t,effect.type,group), summarize,
                        min=min(y, na.rm=T),
                        cl=quantile(y, probs=0.025, na.rm=T),
                        med=median(y, na.rm=T),
                        cu=quantile(y, probs=0.975, na.rm=T),
                        max=max(y, na.rm=T))
    # df.summary.long <- reshape2::melt(df.summary, id.vars=c('t','effect.type'),
    #                                   variable.name='series', value.name = 'y')
    ## GROUP SUMMARY TIME SERIES
    p2 <- ggplot(df.summary, aes(x=t, y=med, color=group)) +
      geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
      geom_line(size=1.2) +
      geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
      geom_hline(yintercept=0) + geom_vline(xintercept=intpd, lty=2) +
      ylab('Y') +
      facet_grid( effect.type ~ . ) +
      theme_bw() + theme(legend.position='top') + ggtitle(plot.main)
    if(plot.show) print(p2)
    if(plot.save) {
      p2.file <- sprintf('internal_intervention_group_summary_series_%s_%s.png',SIMID,PLOTID)
      ggsave(filename=file.path(dir_plot, p2.file), plot=p2,
             width=plot.tall.w,heigh=plot.tall.h,dpi=plot.dpi,units='in')
    }
    
    
    # ## ---- 3. STAGGERED DID DESIGN BY GROUPS---------------------
    # df.t0.summary <- ddply(df.t0, .(t0,effect.type,group), summarize,
    #                        min=min(y, na.rm=T),
    #                        cl=quantile(y, probs=0.025, na.rm=T),
    #                        med=median(y, na.rm=T),
    #                        cu=quantile(y, probs=0.975, na.rm=T),
    #                        max=max(y, na.rm=T))
    # # df.summary.long <- reshape2::melt(df.summary, id.vars=c('t','effect.type'),
    # #                                   variable.name='series', value.name = 'y')
    # p3 <- ggplot(df.t0.summary, aes(x=t0, y=med, color=group)) +
    #   geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
    #   geom_line(size=1.2) +
    #   # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
    #   geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
    #   ylab('Y') +
    #   xlim(c(-intpd+2, npds-intpd-2)) +
    #   facet_grid( . ~ effect.type ) +
    #   theme_bw() + theme(legend.position='top') + ggtitle(plot.main)
    # if(plot.show) print(p3)
    # if(plot.save) {
    #   p3.file <- sprintf('internal_intervention_staggered_DiD_%s_%s.png',SIMID,PLOTID)
    #   ggsave(filename=file.path(dir_plot, p3.file), plot=p3,
    #          width=plot.wide.w,heigh=plot.wide.h,dpi=plot.dpi,units='in')
    # }
    
    ## ---- 4. Median Treatment Effect on Treated -  STAGGERED DID DESIGN ---------------------
    
    # df.att <- ddply(df.t0, .(t0,effect.type,group), summarize,
    #                 # min=min(y, na.rm=T),
    #                 # cl=quantile(y, probs=0.025, na.rm=T),
    #                 med=median(y, na.rm=T)#,
    #                 # cu=quantile(y, probs=0.975, na.rm=T),
    #                 # max=max(y, na.rm=T)
    #                 )
    # df.att[is.nan(df.att)] <- NA
    # .treat <- ddply(df.t0[df.t0$group=='treatment', ], .(t0,effect.type), summarize,
    #                        min=min(y, na.rm=T),
    #                        cl=quantile(y, probs=0.025, na.rm=T),
    #                        med=median(y, na.rm=T),
    #                        cu=quantile(y, probs=0.975, na.rm=T),
    #                        max=max(y, na.rm=T))
    # .control <- ddply(df.t0[df.t0$group=='control', ], .(t0,effect.type), summarize,
    #                 min=min(y, na.rm=T),
    #                 cl=quantile(y, probs=0.025, na.rm=T),
    #                 med=median(y, na.rm=T),
    #                 cu=quantile(y, probs=0.975, na.rm=T),
    #                 max=max(y, na.rm=T))
    # if (nrow(.treat) != nrow(.control))
    #   cat(sprintf('WARNING: treatment (%s) != control (%s) rows\n',nrow(.treat), nrow(.control)))
    # df.att <- .treat
    # df.att$min <- 
    # df.att$cl <- 
    # df.att$med <- 
    # df.att$cu <- 
    # df.att$max <- 
    # df.summary.long <- reshape2::melt(df.summary, id.vars=c('t','effect.type'),
    #                                   variable.name='series', value.name = 'y')
    ##########
    # df.att <- data.frame()
    # t0s <- sort(unique(df.t0.summary$t0)) 
    # for (t in 1:length(t0s)) {
    #   t0 <- t0s[t]
    #   for (k in 1:length(effect.types)) {
    #     batch.t0 <- df[which(df.t0$t0 == t0 & df$effect.type==effect.types[k]), ]
    #     batch.t0.treat <- batch.t0[which(batch.t0$group=='treatment'),]
    #     batch.t0.ctrl  <- batch.t0[which(batch.t0$group=='control'),]
    #     rowdf <- data.frame(t0=t0,effect.type=effect.types[k], 
    #                         y.treat=mean(batch.t0.treat$y), 
    #                         y.control=mean(batch.t0.ctrl$y), 
    #                         y.att=mean(batch.t0.treat$y)-mean(batch.t0.ctrl$y)#,
    #                         # y.att=mean(batch.t0.treat$y - batch.t0.ctrl$y)
    #                         )
    #     df.att <- rbind(df.att, rowdf)
    #   }
    # }
    # ##
    # p4 <- ggplot(df.att, aes(x=t0, y=y.att)) +
    #   # geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
    #   geom_line(size=1.2) +
    #   # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
    #   geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
    #   ylab('ATT (Mean[Treated - Control] by staggered period)') +
    #   xlim(c(-intpd+2, npds-intpd-2)) +
    #   facet_grid( . ~ effect.type ) +
    #   theme_bw() + theme(legend.position='top') + ggtitle(plot.main)
    # if(plot.show) print(p4)
    # if(plot.save) {
    #   p4.file <- sprintf('internal_intervention_staggered_DiD_ATT_%s_%s.png',SIMID,PLOTID)
    #   ggsave(filename=file.path(dir_plot, p4.file), plot=p4,
    #          width=plot.wide.w,heigh=plot.wide.h,dpi=plot.dpi,units='in')
    # }
    
    if(verbose) cat('done.\n\n')
  }
  
  
  ##=================
  ## Return
  ##-----------------
  
  return(list(
    df=df,
    # df.t0=df.t0,
    df.summary=df.summary,
    # df.t0.summary=df.t0.summary,
    df.att=df.att,
    df.plot=df.plot,
    plot.allseries=p1,
    plot.group=p2,
    # plot.group.staggered=p3,
    id=SIMID
  ))
  
  
}



####################################
##  RUN INTERVENTION SIMULATION IN LOOP OVER EFFECT TYPES
######################################
runSimUpdateSimlist <- function(simlist,     ## n, npds, intpd moved into simlist elements
                                effect.types=c('constant','quadratic','geometric'), 
                                sim.id=round(10*as.numeric(Sys.time())),
                                dgp.prior.sd.weight=.01,
                                plot.show=F, plot.save=F, verbose=TRUE) {
  
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
    if(verbose) cat(sprintf('\nScenario label: %s\n\n', key))
    ##  
    set.seed( ifelse(is.null(sim$rand.seed), 54321, sim$rand.seed) )
    noise.level <- ifelse(is.null(sim$noise.level), 0, sim$noise.level)
    ##
    # cov.scenario <- if(is.null(sim$cov.scenario)){list(c1=.1, c2=.2, c3=.3)} else{ sim$cov.scenario }
    # cov.scenario.key <- ifelse(is.null(sim$cov.scenario.key), 'mid', sim$cov.scenario.key)
    ##
    simlist[[key]]$sim <- runSimSingleInterventionEffectComparison(
      effect.types = effect.types,
      n = sim$n, ## NUMBER OF FIRMS
      npds = sim$npds, ## NUMBER OF PERIODS
      intpd = sim$intpd, ## intervention after first section
      ystart = 0.1,
      treat.rule = ifelse(is.null(sim$treat.rule), NA, sim$treat.rule),
      treat.prob =  ifelse(is.null(sim$treat.prob), NA, sim$treat.prob), #0.95,  ## 0.6
      treat.threshold = ifelse(is.null(sim$treat.threshold), NA, sim$treat.threshold),  # 0.015
      sim.id = sim.id, ## defaults to timestamp
      # cov.scenario = cov.scenario,
      # cov.scenario.key = cov.scenario.key,
      ##
      noise.level = noise.level,
      ##
      b4 = ifelse(is.null(sim$b4), 0, sim$b4), ## past performance
      b5 = ifelse(is.null(sim$b5), 1, sim$b5), ## weight of level component in y outcome
      b9 = ifelse(is.null(sim$b9), 0, sim$b9), ## Autocorrelation
      ## Dynamic treatment effect polynomial parameters
      w0 = ifelse(is.null(sim$w0), 1.5 , sim$w0), ## constant
      w1 = ifelse(is.null(sim$w1), 0.13 , sim$w1), ## linear
      w2 = ifelse(is.null(sim$w2), -.01 / sqrt(sim$npds) , sim$w2), ## ,  ##-0.009,  ## quadratic  ## -0.005, ## ***** made steeper curve= -0.008 *****
      # ## TODO: CHECK w2.shift SENSITIVITY - this shifts quadratic curve several periods to the right so that treatment effect increases slowly  
      # w2.shift = ifelse(is.null(sim$w2.shift), -round( sqrt(sim$npds)*.8 ) , sim$w2.shift),
      # w2.shift = -round( sqrt(sim$npds)*.8 ),  ## optimal value here is likely a function of the combination of treatment effect function parameters
      ##
      nseasons  = ifelse(is.null(sim$dgp.nseasons), NA, sim$dgp.nseasons), 
      season.frequency = ifelse(is.null(sim$dgp.freq), NA, sim$dgp.freq),
      ##
      dgp.prior.sd.weight=dgp.prior.sd.weight,
      ## Plotting
      plot.show = plot.show, ## TRUE
      plot.save = plot.save, ## TRUE
      verbose = verbose ## echo messages to console
    )
  }
  
  return(simlist)
}





cat('\nLoaded Single Intervention Simulation - vectorized.\n')



# runSimInternalInterventionEffectComparison(effect.types=c('constant','geometric','quadratic'))



################################################################
################################################################
######################## DEBUG / TESTING #######################
################################################################
################################################################

# ## DEFAULT NO GROWTH NO PERFORMANCE SPILLOVER
# sim1 <- runInternalInterventionSim(
#   n = n, ## NUMBER OF actors  
#   npds = npds,  ## NUMBER OF PERIODS
#   intpd = round( npds / 5 ), ## intervention after first section
#   ystart = 0,
#   effect.types = c('constant','quadratic','geometric'),
#   benchmark.type = 'self', ## 'all', 'self'
#   treat.threshold = 0.02,
#   ##
#   b4 = 0, ## past performance
#   b5 = 0, ## growth (linear function of time t)
#   ##
#   w0 = 1.3,
#   w1 = 0.3,
#   w2 = 0.12,
# )
# 


# sim1 <- simlist[[1]]
# 
# treat.map <- sim1$df.plot.treat.map
# df <- sim1$df
# df$actor <- as.numeric(as.character(df$actor))
# 
# ( ctr.map <- treat.map$actor[which(treat.map$group == 'control')] )
# ( ctr.df <- unique(df$actor[which(df$group=='control')]) )





# ##---------------------------------------
# ## LOW GROWTH LOW PERFORMANCE SPILLOVER
# sim2 <- runInternalInterventionSim(
#   n = 100, ## NUMBER OF actors  
#   npds = 100, ## NUMBER OF PERIODS
#   intpd = round( npds / 5 ), ## intervention after first section
#   ystart = 0,
#   effect.types = c('constant','quadratic','geometric'),
#   benchmark.type = 'self', ## 'all', 'self'
#   treat.threshold = 0.02,
#   ##
#   b4 = 0.1, ## past performance
#   b5 = 0.01, ## growth (linear function of time t)
#   ##
#   w0 = 1.3,
#   w1 = 0.3,
#   w2 = 0.12,
# )
# 
# 
# ##---------------------------------------
# ## LOW GROWTH MID PERFORMANCE SPILLOVER
# sim3 <- runInternalInterventionSim(
#   n = 100, ## NUMBER OF actors  
#   npds = 100, ## NUMBER OF PERIODS
#   intpd = round( npds / 5 ), ## intervention after first section
#   ystart = 0,
#   effect.types = c('constant','quadratic','geometric'),
#   benchmark.type = 'self', ## 'all', 'self'
#   treat.threshold = 0.02,
#   ##
#   b4 = 0.3, ## past performance
#   b5 = 0.01, ## growth (linear function of time t)
#   ##
#   w0 = 1.3,
#   w1 = 0.3,
#   w2 = 0.12,
# )
# 
# ##---------------------------------------
# ## LOW GROWTH HI PERFORMANCE SPILLOVER
# sim4 <- runInternalInterventionSim(
#   n = 100, ## NUMBER OF actors  
#   npds = 100, ## NUMBER OF PERIODS
#   intpd = round( npds / 5 ), ## intervention after first section
#   ystart = 0,
#   effect.types = c('constant','quadratic','geometric'),
#   benchmark.type = 'self', ## 'all', 'self'
#   treat.threshold = 0.02,
#   ##
#   b4 = 0.6, ## past performance
#   b5 = 0.01, ## growth (linear function of time t)
#   ##
#   w0 = 1.3,
#   w1 = 0.3,
#   w2 = 0.12,
# )

# #
# yvec <- c()
# x1vec <- c()
# x2vec <- c()
# b3vec <- c()
# for (t in 1:npds) {
#
#   vt <- runif(1)
#   ut <- runif(1)
#    
#   y.tm1 <- ifelse(t==1, ystart, yvec[t-1])
#  
#   x2vec[t] <- x2(t, intpd)
#  
#   .x1.t <- x1(g0, g1, vt,  y.tm1, logit.shift, logit.scale)
#   .x1.mean <- ifelse(t==1, 0, mean(x1vec) )
#   x1vec[t] <- ifelse(.x1.t > .x1.mean, 1, 0)
#  
#   b3vec[t] <- b3(t)
#   yvec[t] <- y(b0, b1, b2, b3vec[t], ut, x1vec[t], x2vec[t])
#
# }
#
# par(mfrow=c(2,2), mar=c(3,4,2,1))
# plot(yvec, main='y', type='o');abline(h=0)
# plot(x1vec, main='x1',type='o');abline(h=0)
# plot(b3vec, main='b3',type='o');abline(h=0)
#






# ## Self-select into treatment as function of past performance
# ## [char] treat.rule Rule for actor to select into treatment group
# ##                    'self'=below own actor's past performance treat.threshold --> intervention
# ##                    'all'=below all actors' past performance records treat.threshold --> intervention
# ##                    'random'=Bernoulli draw with p=treat.threshold
# ## 
# x1Func <- function(g0, g1, y.tm1, v, 
#                    treat.rule='self', 
#                    treat.threshold=0.5, ## above x proportion --> select intervention
#                    df=NA, ## dataframe with past records used for 'self' or 'all' treat.rule
#                    self=NA, ## name of focal actor (used with 'self' treat.rule)
#                    shift=1, scale=1) {
#   # z <- g0 + g1*y.tm1 + v
#   cat(sprintf('z=%.3f | ',z))
#   ## Default to zero unless catching treatment rule
#   out <- 0
#   ## negative z flips the s-curve
#   ## so higher performance makes actor
#   ## LESS likely (approaching zero)
#   ## to take the intervention
#   ##  (ie, endogenously self-select into treatment group)
#   # .x1.t <- getLogisticProportion( -z, shift, scale )
#   # out <- sample(1:0, size = 1, prob = c(.x1.t, 1-.x1.t))
#   if (treat.rule == 'self') {
#     x <- mean( df$y[which(df$actor==self)], na.rm=TRUE)
#     if (y.tm1 <= x) {
#       out <- 1
#     }
#   } else if (treat.rule == 'all') {
#     x <- mean( df$y, na.rm=TRUE)
#     
#   } else if (treat.rule == 'random') {
#     p <- treat.threshold  ## P(treatment==1)
#     sample(0:1, size = 1, prob = c(1-p, p))
#   } else {
#     print(sprintf('treatment rule "%s" not known',treat.rule))
#     return()
#   }
#   cat(sprintf('prop=%.3f \n\n',out))
#   return(out)
# } 
