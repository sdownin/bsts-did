##
##   Internal intervention simulation - Vectorized 
##
##     - vectorized actor loop i=1:n  to improve performance
##
##

library(ggplot2)
library(plyr)
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





##             Treat      Post-int   Treat-post-int   Past Perf.  Growth  Age/sex/marr.   noise
## y[t] = b0 + b1*x1[t] + b2*x2[t] + b3*x1[t]*x2[t] + b4*y[t-1] + b5*t  + [covariates]  + u[t] ##TODO: covariates
## x1 = rule[[? f(y[t-1])]]


## Two cases to simulate: 
## one time for intervention?
## staggered over time ?


##================================
## MAIN VARIABLE FUNCTIONS
##--------------------------------
## outcome
yFunc <- function(b0, b1, b2, b3, b4, b5, x1, x2, y.tm1, t, u) {
  b0 + b1*x1 + b2*x2 + b3*x1*x2 + b4*y.tm1 + b5*t + u
}


## Self-select into treatment as function of past performance
## [char] treat.rule Rule for actor to select into treatment group
##                    'past'=below past performance treat.threshold (focal or all actors)  --> intervention
##                    'random'=Bernoulli draw with p=treat.threshold
##                    'worst'=worst performing actor goes for treatment with p=treat.threshold
##                    'below.median'=actors with past performance below median select into treatment with p=treat.threshold
## 
x1Func <- function(x3.treat.post.int, y.tm1, v, ## lengths must match
                   g0=0.001, g1=0.001, 
                   treat.rule='below.median', 
                   treat.threshold=0.1, ## above x proportion --> select intervention
                   records=c(), ## dataframe with past records used for 'self' or 'all' treat.rule
                   shift=1, scale=1) {
  
  n <- length(x3.treat.post.int) ### x3.treat.post.int, y.tm1, v, <-- lengths must match
  p <- treat.threshold ## treatment self-selection threshold
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
  if (treat.rule == 'below.median') {
    
    benchmark <- quantile(records, probs=0.5, na.rm=TRUE)
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
b3Func <- function(t.post.intpd, x2, ## must be same length vectors
                   type='quadratic', 
                   w0=1.3, w1=.3, w2=-.012, w2.shift=0,
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
  idx <- which( x2 == 1 )
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
      tx <- t[idx] + w2.shift ## ? TODO: check sensitivity for gradual treatment effect increase (negative number shifts curve rightward)
      new <-  w2*(tx^2) + w1*tx + w0    ## length == n.idx
    }
    else if (type == 'geometric') {
      denom <- t[idx]  ## length == n.idx
      denom[denom==0] <- 1  ## prevent division by zero
      new <- w0 / denom
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


# # #### DEBUG ###########
# n = 50 ## NUMBER OF actors
# npds = 60 ## NUMBER OF PERIODS
# intpd = round( npds / 5 ) ## intervention after first fifth
# ystart = 0
# effect.types = c('constant','quadratic','geometric') ## treatment effect shapes
# treat.rule = 'below.median'
# treat.threshold = 0.05 ## treatment threshold (either probability or quantile of performance below which low performer(s) self-select into treatment)
# benchmark.type = 'self' ## 'all', 'self'
# seed = 54321  ## pseudo-random number generator seed for replication
# ##
# sim.id = NA ## Index or timestamp of simulation for  saving, etc.
# ## # PEFORMANCE [Y] FUNCTION PARAMETERS
# b0 = .001 ## intercept
# b1 = .001 ## treatment dummy
# b2 = .001 ## post intervention dummy
# b3 = .001 ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
# b4 = 1/3 ## spillover of past performance on current performance (how much of treatment effect persists across periods)
# b5 = .01 ## growth rate (linear effect of time: proportion of time t added to linear combination in yFunc() performance )
# ## # TREATMENT EFFECT FUNCTION WEIGHTS
# w0 = 1.3
# w1 = 0.3
# w2 = 0.12
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


##=================================
##
## Internal Intervention Simulation Function
## 
##=================================
runInternalInterventionSim <- function(
    n = 50, ## NUMBER OF actors  
    npds = 60, ## NUMBER OF PERIODS
    intpd = round( npds / 5 ), ## intervention after first fifth
    ystart = 0,
    effect.types = c('constant','quadratic','geometric'), ## treatment effect shapes
    treat.rule = 'below.median',
    treat.threshold = 0.05, ## treatment threshold (either probability or quantile of performance below which low performer(s) self-select into treatment)
    benchmark.type = 'self', ## 'all', 'self'
    seed = 54321,  ## pseudo-random number generator seed for replication
    ##
    sim.id = NA, ## Index or timestamp of simulation for  saving, etc.
    ## # PEFORMANCE [Y] FUNCTION PARAMETERS
    b0 = .001, ## intercept
    b1 = .001, ## treatment dummy
    b2 = .001, ## post intervention dummy
    # b3 = .001, ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
    b4 = 1/3, ## spillover of past performance on current performance (how much of treatment effect persists across periods)
    b5 = .01, ## growth rate (linear effect of time: proportion of time t added to linear combination in yFunc() performance )
    ## # TREATMENT EFFECT FUNCTION WEIGHTS 
    w0 = 1.3, ## constant
    w1 = 0.3, ## linear
    w2 = -0.12, ## quadratic
    w2.shift = -10, ## shift quadratic curve to allow gradual treatment effect increase (negative shifts curve rightward)
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
    plot.wide.h=5,
    plot.tall.w=5,
    plot.tall.h=10,
    plot.sq.w=9,
    plot.sq.h=9,
    plot.dpi=300,
    ...
  ) {
  
  ## --- BEGIN SIM ---

  SIMID <- ifelse(is.na(sim.id), round(10*as.numeric(Sys.time())), sim.id)
  
  ntypes <- length(effect.types)
  
  sim.count <- ntypes * npds  ## n
  
  cat(sprintf('\nRunning Internal Intervention Simulation for %s iterations:\n  %s periods, %s effect types\n',
              sim.count, npds, ntypes))
  
  ## Initialize progress bar
  counter <- 0
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = sim.count, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 80,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  set.seed(seed)
  
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
  
  for (k in 1:length(effect.types))
  {
    effect.type <- effect.types[k]
      
    for (t in 1:npds) ## PERIOD LOOP
    {
      # cat(sprintf('\n\n\n----- t = %s -----\n\n\n',t))
      
      ## TODO: Vectorize the actor loop to improve performance (slow loops in R)
      # for (i in 1:n) ## actors LOOP
      # {
        # cat(sprintf(' i = %s \n',i))
      
      ## noise terms
      v <- rnorm(n=n,mean=0,sd=.5)  ##- .5 ##
      u <- rnorm(n=n,mean=0,sd=.5) ## 
      
      ##----------  Past Period Variables---------------------
      ## indices of last period observations
      idx.tm1 <- if (t==1) { c() } else { which( df$t == t-1 & df$effect.type==effect.type ) }
      ## treated dummy (post-intervention-treatment-effect)
      x3.tm1 <- if (t==1) { rep(0, n) } else { df$x3[idx.tm1] }
      ## focal actor past performance
      y.tm1 <- if (t==1) { rep(ystart, n) } else { df$y[idx.tm1] }
      
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

      ## TREATMENT SELF SELECTION DUMMY
      x1 <- x1Func(x3.tm1, y.tm1, v,
                   g0, g1, 
                   treat.rule = 'below.median', 
                   treat.threshold = treat.threshold, #0.5 for 'worst'; 0.1 for 'past' (?)
                   records=perf.records
                   )
      
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
      b3 <- b3Func(t.post.intpd, x2, type=effect.type,
                   w0=ifelse(effect.type=='geometric', w0*4, w0),
                   w1=w1, w2=w2, w2.shift=w2.shift)
      
      ## PERFORMANCE
      y <- yFunc(b0, b1, b2, b3, b4, b5, x1, x2, y.tm1, t, u)

      
      ##--------------------------------------
      ## UPDATE PERIOD DATAFRAME
      ##--------------------------------------
      df.t <- data.frame(actor=1:n, t=rep(t, n), 
                         t.post.intpd=t.post.intpd,
                         effect.type=rep(effect.type, n),
                         y=y, x1=x1, x2=x2, x3=x1*x2,
                         b1=rep(b1, n), b2=rep(b2, n), 
                         b3=b3, u=u, v=v,
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
    
  } ## /end effect loop
    
  # }
  
  ## DEBUG 
  df <- unique(df)
  df$actor <- as.numeric(as.character(df$actor)) ## avoid possible factor issue (?)
  
  
  ## COHORT GROUP ASSIGNMENT 
  ## ASSIGN DEFAULT GROUP AND COLOR (before updating treated cohort)
  df$group <- 'control'
  ##
  ggcolors <- hue_pal()(2)
  df$group.color <- ggcolors[1] ## control red color
  ##
  ## Assign cohort group (treatment, control)
  ## based on rule: treatment=treated by end; control=not treated by end
  for (k in 1:ntypes) {
    effect.type <- effect.types[k]
    actor.eff.tr <- sort(unique(df$actor[which(df$x3==1 & df$effect.type==effect.type)]))
    for (f in actor.eff.tr) {
      # cat(' ',f)
      ## all observation indices for actor f with effect type=effect.type
      idx.act.eff.tr <- which( df$actor==f & df$effect.type==effect.type )
      ## if actor is treated at any time (any x3=0) then apply treatment group to all actor observations (periods) for this effect type
      df$group[idx.act.eff.tr] <- 'treatment'
      df$group.color[idx.act.eff.tr] <- ggcolors[2]  ## treatment green color
    }
  }

  ## close progress bar
  close(pb)
  
  
  ## STAGGERED DID ARRANEMENT:
  ## Reset each series period as t=0 when intervention selected
  df.t0 <- df
  df.t0$t0 <- NA
  for (effect.type in effect.types) {
    for (f in 1:n) {
      x3.idx <- which(df.t0$actor==f & df.t0$effect.type==effect.type  & df.t0$x3 > 0)
      f.x3.t <- ifelse(length(x3.idx)>0, min(df.t0$t[x3.idx]), NA)
      df.t0$t0[which(df.t0$actor==f & df.t0$effect.type==effect.type)] <- if (length(f.x3.t)==0 | all(is.na(f.x3.t))) {
        (1:npds) - intpd + 1
      } else {
        (1:npds) - f.x3.t + 1 ## t0 = f.x3.t period
      }
    }
  }



  ##=============================================
  ## PLOTS
  ##---------------------------------------------
  if (plot.show | plot.save)
  {
    PLOTID <- round(10*as.numeric(Sys.time()))

    cat('\n rendering plots...')
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
    df.plot$actor <-  as.factor(df.plot$actor)
    p1 <- ggplot(data=df.plot, mapping=aes(x=t,y=y, color=group, group=actor)) +
      geom_line(size=1.05, alpha=0.3) +
      # geom_point(, color=rgb(0.2, 0.2, 0.8, 0.1))+
      geom_hline(yintercept=0) + facet_grid( . ~ effect.type) +
      theme_bw() + theme(legend.position='none') + ggtitle(plot.main) #+
      # guides(color=guide_legend(nrow=3,byrow=TRUE)) #+
    #
    # facet_wrap(.~factor(df$x1) + factor(df$x2))
    if(plot.show) print(p1)
    if(plot.save) {
      p1.file <- sprintf('internal_intervention_all_actor_timeseries_%s_%s.png',SIMID,PLOTID)
      ggsave(filename=file.path(dir_plot, p1.file), plot=p1,
             width=plot.wide.w,heigh=plot.wide.h,dpi=plot.dpi,units='in')
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
      facet_grid( . ~ effect.type ) +
      theme_bw() + theme(legend.position='top') + ggtitle(plot.main)
    if(plot.show) print(p2)
    if(plot.save) {
      p2.file <- sprintf('internal_intervention_group_summary_series_%s_%s.png',SIMID,PLOTID)
      ggsave(filename=file.path(dir_plot, p2.file), plot=p2,
             width=plot.wide.w,heigh=plot.wide.h,dpi=plot.dpi,units='in')
    }


    ## ---- 3. STAGGERED DID DESIGN BY GROUPS---------------------
    df.t0.summary <- ddply(df.t0, .(t0,effect.type,group), summarize,
                           min=min(y, na.rm=T),
                           cl=quantile(y, probs=0.025, na.rm=T),
                           med=median(y, na.rm=T),
                           cu=quantile(y, probs=0.975, na.rm=T),
                           max=max(y, na.rm=T))
    # df.summary.long <- reshape2::melt(df.summary, id.vars=c('t','effect.type'),
    #                                   variable.name='series', value.name = 'y')
    p3 <- ggplot(df.t0.summary, aes(x=t0, y=med, color=group)) +
      geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
      geom_line(size=1.2) +
      # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
      geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
      ylab('Y') +
      xlim(c(-intpd+2, npds-intpd-2)) +
      facet_grid( . ~ effect.type ) +
      theme_bw() + theme(legend.position='top') + ggtitle(plot.main)
    if(plot.show) print(p3)
    if(plot.save) {
      p3.file <- sprintf('internal_intervention_staggered_DiD_%s_%s.png',SIMID,PLOTID)
      ggsave(filename=file.path(dir_plot, p3.file), plot=p3,
             width=plot.wide.w,heigh=plot.wide.h,dpi=plot.dpi,units='in')
    }

    cat('done.\n\n')
  }
 
  
  ##============================
  ## OUTPUT LIST
  ##----------------------------
  return(list(
    df=df,
    df.t0=df.t0,
    df.summary=df.summary,
    df.t0.summary=df.t0.summary,
    df.plot=df.plot,
    plot.allseries=p1,
    plot.group=p2,
    plot.group.staggered=p3,
    id=SIMID
  ))

}


cat('\nLoaded Internal Intervention Simulation - vectorized.\n')






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
