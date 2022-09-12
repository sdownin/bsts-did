##
##   Internal intervention simulation
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
yFunc <- function(b0, b1, b2, b3, b4, b5, x1, x2, y.t.minus1, t, u) {
  b0 + b1*x1 + b2*x2 + b3*x1*x2 + b4*y.t.minus1 + b5*t + u
}


## Self-select into treatment as function of past performance
## [char] treat.rule Rule for firm to select into treatment group
##                    'past'=below past performance treat.threshold (focal or all firms)  --> intervention
##                    'random'=Bernoulli draw with p=treat.threshold
##                     'worst'=worst performing firm goes for treatment with p=treat.threshold
## 
x1Func <- function(x3.treat.post.int,
                   g0, g1, y.t.minus1, v, 
                   treat.rule='worst', 
                   treat.threshold=0.5, ## above x proportion --> select intervention
                   records=c(), ## dataframe with past records used for 'self' or 'all' treat.rule
                   shift=1, scale=1) {
  
  if (x3.treat.post.int==1)
    return(1) ## treatment is monotonically increasing (cannot untreat after being treated)
  
  ## TODO: see if transform carries any implications (?)
  ##       comparison of z vs. records should be invariant to linear transform (?)
  z <- g0 + g1*y.t.minus1 + v
  records <- g0 + g1*records + v
  # cat(sprintf('z=%.3f | rule=%s | ',z, treat.rule))
  ##
  # cat(sprintf('rule=%s \n\n',treat.rule))
  ## Default to zero unless catching treatment rule
  out <- 0
  ## negative z flips the s-curve
  ## so higher performance makes firm
  ## LESS likely (approaching zero)
  ## to take the intervention
  ##  (ie, endogenously self-select into treatment group)
  # .x1.t <- getLogisticProportion( -z, shift, scale )
  # out <- sample(1:0, size = 1, prob = c(.x1.t, 1-.x1.t))
  if (treat.rule == 'past') {
    
    benchmark <- quantile(records, probs=treat.threshold, na.rm=TRUE)
    
    if (length(benchmark) == 0 | is.na(benchmark)){
      print("missing benchmark. records return NULL threshold quantile")
      print(records)
    }
    if (z < benchmark) {
      # cat(sprintf('\nbenchmark = %.3f\n',benchmark))
      out <- 1 
    }
    
  } else if (treat.rule == 'worst') {
    
    if (z == min(records)) { #focal firm's past performance is worst in last period
      p <- treat.threshold  ## P(treatment==1)
      out <- sample(0:1, size = 1, prob = c(1-p, p))
      # cat(sprintf('\nz == min(records) | out=%s\n',out))
    }
    
  } else if (treat.rule == 'random') {
    
    p <- treat.threshold  ## P(treatment==1)
    out <- sample(0:1, size = 1, prob = c(1-p, p))
    
  } else {
    print(sprintf('treatment rule "%s" not known',treat.rule))
    return()
  }
  
  # cat(sprintf('prop=%.3f \n\n',out))
  return(out)
} 


## post-intervention dummy
x2Func <- function(t, intpd) {
  ifelse( t >= intpd, 1, 0)
} 

## Treatment effect
b3Func <- function(t.post.intpd, x2, type='quadratic', 
                   w0=1.3, w1=.3, w2=.12, 
                   non.negative=TRUE) {
  t <- t.post.intpd ## time period 't' is number of periods after intervention
  if (x2==0) ## zero in pre-intervention periods
    return(0)
  ## constant  effect
  if(type == 'constant') {
    b3 <- w0
  }
  else if (type == 'linear') {
    b3 <- w1*t + w0
  }
  else if (type == 'quadratic') {
    b3 <-  -w2*t^2 + w1*t + w0
  }
  else if (type == 'geometric') {
    b3 <- w0 / ifelse(t!=0, t, 1)
  }
  else {
    print(sprintf('effect type %s not known',type))
    return()
  }
  # cat(sprintf('\nb3Func: b3 before max() = %.3f\n',b3))
  if (non.negative) {
    b3 <- max(b3, 0)
  }
  return(b3)
}






##=================================
##
## Internal Intervention Simulation Function
## 
##=================================
runInternalInterventionSim <- function(
    n = 50, ## NUMBER OF FIRMS  
    npds = 60, ## NUMBER OF PERIODS
    intpd = round( npds / 5 ), ## intervention after first fifth
    ystart = 0,
    effect.types = c('constant','quadratic','geometric'), ## treatment effect shapes
    benchmark.type = 'self', ## 'all', 'self'
    treat.threshold = 0.02, ## treatment threshold (either probability or quantile of performance below which low performer(s) self-select into treatment)
    seed = 54321,  ## pseudo-random number generator seed for replication
    ##
    sim.id = NA, ## Index or timestamp of simulation for  saving, etc.
    ## # PEFORMANCE [Y] FUNCTION PARAMETERS
    b0 = .001, ## intercept
    b1 = .001, ## treatment dummy
    b2 = .001, ## post intervention dummy
    b3 = .001, ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
    b4 = 1/3, ## spillover of past performance on current performance (how much of treatment effect persists across periods)
    b5 = .01, ## growth rate (linear effect of time: proportion of time t added to linear combination in yFunc() performance )
    ## # TREATMENT EFFECT FUNCTION WEIGHTS 
    w0 = 1.3,
    w1 = 0.3,
    w2 = 0.12,
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
  
  sim.count <- ntypes * npds * n
  
  cat(sprintf('\nRunning Internal Intervention Simulation for %s iterations:\n  %s firms, %s periods, %s effect types\n',
              sim.count, n, npds, ntypes))
  
  ## Initialize progress bar
  counter <- 0
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = sim.count, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 60,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  set.seed(seed)
  
  ##================================
  ## MAIN RUN
  ##--------------------------------
  ## VARIABLES IN SIMULATION PER FIRM PER PERIOD
  # y <- NA
  # x1 <- NA
  # x2 <- NA
  # b3 <- NA
  # u <- NA
  # v <- NA
  ##------
  ## OUTPUT OBJECTS (to be output as list)
  df <- data.frame(stringsAsFactors = F)
  df.t0 <- NA
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
      
      ## TODO: Vectorize the firm loop to improve performance (slow loops in R)
      for (i in 1:n) ## FIRMS LOOP
      {
        # cat(sprintf(' i = %s \n',i))
        
        ## noise terms
        v <- rnorm(n=1,mean=0,sd=.5)  ##- .5 ##
        u <- rnorm(n=1,mean=0,sd=.5) ## 
        
        ## focal firm past performance
        y.t.minus1 <- if (t==1) { ystart } else { df$y[which(df$firm==i & df$t==t-1 & df$effect.type==effect.type)] }
        
        ## post-intervention dummy var
        x2 <- x2Func(t, intpd)
        
        ## periods after intervention
        t.post.intpd <- if(t <= intpd) {
          0
        } else {
          ifelse(max(df$x3[which(df$firm==i & df$effect.type==effect.type)]) > 0, ## has had post-intervention-pd treatment
                 df$firm.t.post.intpd[which(df$firm==i & df$t==t-1 & df$effect.type==effect.type)] + 1,
                 0)
        }
        
        ##------------------------------------------------------
        ## Get past performance records for treat.rule == 'past'
        ##------------------------------------------------------
        # perf.records <- if (i==1) { ## if first firm (df not updated yet)
        #   if (t==1){ ystart } else { y } ## first firm must use default start value of performance or past performance in 'y' vector
        # } else {
        #   if (t==1) {
        #     switch(benchmark.type,
        #            "self"=ystart,  ## focal firm i not in data.frame yet when t==1, so just use ystart
        #            "all"=df$y)
        #   } else {
        #     switch(benchmark.type,
        #            "self" = df$y[which(df$firm == i)],
        #            "all" = df$y)
        #   }
        # }
        len.memory <- 3
        .ts <- sapply(1:len.memory, function(a)max(t-a,0))
        t.in.memory <- .ts[.ts > 0] 
        perf.records <- if (t==1) { ## if first firm record 
          switch(benchmark.type,
                 "self" = ystart, ## first time for firm i must use default start value of performance or past performance in 'y' vector
                 "all" = df$y)
        } else {
          switch(benchmark.type,
                 "self" = df$y[which(df$firm == i & df$t %in% t.in.memory  & df$effect.type==effect.type)],
                 "all" = df$y[which(df$t %in% t.in.memory  & df$effect.type==effect.type)])
        }
        # print('perf.recrods')
        # print(perf.records)
        ## Treatment dummy (firm self-select into treatment based on treat.rule)
        x3.t.minus1 <- if (t==1) { 0 } else { df$x3[which(df$firm==i & df$t==t-1 & df$effect.type==effect.type)] }
        x1 <- x1Func(x3.t.minus1,
                     g0, g1, y.t.minus1, v,
                     treat.rule = 'worst', 
                     treat.threshold = treat.threshold, #0.5 for 'worst'; 0.1 for 'past' (?)
                     records=perf.records,
                     logit.shift, logit.scale)
        # .x1.t <- x1Func(g0, g1, v[t],  y.t.minus1,
        #                 logit.shift, logit.scale)
        # x1[t] <- .x1.t # round(.x1.t)
        # x1[t] <- sample(0:1, size = 1, prob = c(.x1.t, 1-.x1.t))
        # .x1.mean <- ifelse(t==1, 0, mean(x1vec) )
        # x1vec[t] <- ifelse(.x1.t > .x1.mean, 1, 0)
        
        b3 <- b3Func(t.post.intpd, x2, type=effect.type,
                     w0=ifelse(effect.type=='geometric', w0*6, w0))
        
        y <- yFunc(b0, b1, b2, b3, b4, b5, x1, x2, y.t.minus1, t, u)
        
        
        df.t <- data.frame(firm=as.numeric(i), t=t, firm.t.post.intpd=t.post.intpd,
                           effect.type=effect.type,
                           y=y, x1=x1, x2=x2, x3=x1 * x2,
                           b1=b1, b2=b2, b3=b3, 
                           u=u, v=v,
                           group=NA, group.color=NA,
                           stringsAsFactors = F)
        df <- rbind(df, df.t)
        
        ## Update Progress Bar
        ## starting at 0 --> add 1 after completing first item / before progress bar update
        counter <- counter + 1 
        setTxtProgressBar(pb, counter)
        
      }
      
    }
    
    ## Assign cohort group label (treatment, control)
    ## based on rule: treatment=treated by end; control=not treated by end
    for (f in 1:n) {
      idx <- which(df$firm==f & df$effect.type==effect.type)
      df$group[idx] <- ifelse(any(df$x3[idx] == 1 | df$firm.t.post.intpd[idx] > 1), 
                              'treatment', 'control')
    }
    ggcolors <- hue_pal()(2)
    df$group.color <- ggcolors[1]
    df$group.color[which(df$group=='treatment')] <- ggcolors[1]
    
  }
 
  
  
  close(pb)
  
  df <- unique(df)
  df$firm <- as.numeric(as.character(df$firm))
  
  
  ## STAGGERED DID ARRANEMENT:
  ## Reset each series period as t=0 when intervention selected
  df.t0 <- df
  df.t0$t0 <- NA
  for (f in 1:n) {
    x3.idx <- which(df.t0$firm==f & df.t0$x3 > 0)
    f.x3.t <- ifelse(length(x3.idx)>0, min(df.t0$t[x3.idx]), NA)
    df.t0$t0[which(df.t0$firm==f)] <- if (length(f.x3.t)==0 | is.na(f.x3.t)) {
      (1:npds) - intpd + 1
    } else {
      (1:npds) - f.x3.t + 1 ## t0 = f.x3.t period
    }
  }


  ##=============================================
  ## PLOTS
  ##---------------------------------------------
  if (plot.show | plot.save)
  {
    PLOTID <- round(10*as.numeric(Sys.time()))

    cat('\n rendering plots...')
    cnt <- plyr::count(df$group)
    prop.treat.obs <- cnt$freq[which(cnt$x == 'treatment')] / sum(cnt$freq)
    prop.treat.firms <- length(unique(df$firm[which(df$group == 'treatment')])) / n
    plot.main <- sprintf('past perf. = %.2f; growth = %.2f,  treat.threshold = %.2f;  prop.treat.obs = %.2f;  prop.treat.firms = %.2f',
                         b4, b5, treat.threshold, prop.treat.obs, prop.treat.firms)


    ## default ggplot2 colors
    # ggcolors <- c("#F8766D", "#00BA38")
    ggcolors <- hue_pal()(2)


    ## ----- 1. ALL INDIVIDUAL TIME SERIES ------------------
    # color <- rgb(.3,.3,.8,.7)
    df$firm <- as.factor(df$firm)
    df.plot <- df[which(df$firm %in% 1:n),] ## filter firms
    treat.map <- unique(df.plot[,c('firm','group')] ) ##***
    treat.map$firm <- as.numeric(as.character(treat.map$firm))
    p1 <- ggplot(data=df.plot, mapping=aes(x=t,y=y, color=firm)) +
      geom_line(aes(color=firm), size=1.05, alpha=0.3) +
      # geom_point(, color=rgb(0.2, 0.2, 0.8, 0.1))+
      geom_hline(yintercept=0) + facet_grid( . ~ effect.type) +
      theme_bw() + theme(legend.position='none') + ggtitle(plot.main) +
      # guides(color=guide_legend(nrow=3,byrow=TRUE)) +
      scale_color_manual(values = sapply(1:nrow(treat.map),function(a){
        ifelse(treat.map$group[a]=='treatment', ggcolors[2], ggcolors[1])
        }))  #rep(rgb(.2,.2,.8,.2),n)
    #
    # facet_wrap(.~factor(df$x1) + factor(df$x2))
    if(plot.show) print(p1)
    if(plot.save) {
      p1.file <- sprintf('internal_intervention_all_firm_timeseries_%s_%s.png',SIMID,PLOTID)
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
      geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
      geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
      ylab('Y') +
      xlim(c(-intpd + 1, round(npds*0.7))) +
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
    df.plot.treat.map=treat.map,
    plot.allseries=p1,
    plot.group=p2,
    plot.group.staggered=p3,
    id=SIMID
  ))

}


cat('\nLoaded Internal Intervention Simulation.\n')






################################################################
################################################################
######################## DEBUG / TESTING #######################
################################################################
################################################################

# ## DEFAULT NO GROWTH NO PERFORMANCE SPILLOVER
# sim1 <- runInternalInterventionSim(
#   n = n, ## NUMBER OF FIRMS  
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
# df$firm <- as.numeric(as.character(df$firm))
# 
# ( ctr.map <- treat.map$firm[which(treat.map$group == 'control')] )
# ( ctr.df <- unique(df$firm[which(df$group=='control')]) )





# ##---------------------------------------
# ## LOW GROWTH LOW PERFORMANCE SPILLOVER
# sim2 <- runInternalInterventionSim(
#   n = 100, ## NUMBER OF FIRMS  
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
#   n = 100, ## NUMBER OF FIRMS  
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
#   n = 100, ## NUMBER OF FIRMS  
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
#   y.t.minus1 <- ifelse(t==1, ystart, yvec[t-1])
#  
#   x2vec[t] <- x2(t, intpd)
#  
#   .x1.t <- x1(g0, g1, vt,  y.t.minus1, logit.shift, logit.scale)
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
# ## [char] treat.rule Rule for firm to select into treatment group
# ##                    'self'=below own firm's past performance treat.threshold --> intervention
# ##                    'all'=below all firms' past performance records treat.threshold --> intervention
# ##                    'random'=Bernoulli draw with p=treat.threshold
# ## 
# x1Func <- function(g0, g1, y.t.minus1, v, 
#                    treat.rule='self', 
#                    treat.threshold=0.5, ## above x proportion --> select intervention
#                    df=NA, ## dataframe with past records used for 'self' or 'all' treat.rule
#                    self=NA, ## name of focal firm (used with 'self' treat.rule)
#                    shift=1, scale=1) {
#   # z <- g0 + g1*y.t.minus1 + v
#   cat(sprintf('z=%.3f | ',z))
#   ## Default to zero unless catching treatment rule
#   out <- 0
#   ## negative z flips the s-curve
#   ## so higher performance makes firm
#   ## LESS likely (approaching zero)
#   ## to take the intervention
#   ##  (ie, endogenously self-select into treatment group)
#   # .x1.t <- getLogisticProportion( -z, shift, scale )
#   # out <- sample(1:0, size = 1, prob = c(.x1.t, 1-.x1.t))
#   if (treat.rule == 'self') {
#     x <- mean( df$y[which(df$firm==self)], na.rm=TRUE)
#     if (y.t.minus1 <= x) {
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
