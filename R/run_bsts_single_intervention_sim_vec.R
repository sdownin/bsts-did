#####################################################
##
##   Dynamic Causal Inference Simulations
##
##    Run script for simulations of
##     - internal intervention (self-selection) endogeneity
##     - ...
##
#####################################################
library(tidyr)

## Directories
dir_proj <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'   
dir_plot <- file.path(dir_proj, 'plots')
dir_r <- file.path(dir_proj,'R')

setwd(dir_proj)



##==============================
## Load simulation functions
##------------------------------
# source(file.path(dir_r,'internal_intervention_sim.R'))
source(file.path(dir_r,'single_intervention_sim_vec.R'))  ## vectorized -- in development


##=================================
## Simulation Function Arguments
##=================================
# runInternalInterventionSim <- function(
#   n = 50, ## NUMBER OF actors  
#   npds = 60, ## NUMBER OF PERIODS
#   intpd = round( npds / 5 ), ## intervention after first fifth
#   ystart = 0,
#   effect.types = c('constant','quadratic','geometric'), ## treatment effect shapes
#   treat.rule = 'below.benchmark',
#   treat.threshold = 0.05, ## treatment threshold (either probability or quantile of performance below which low performer(s) self-select into treatment)
#   benchmark.type = 'self', ## 'all', 'self'
#   seed = 54321,  ## pseudo-random number generator seed for replication
#   ##
#   sim.id = NA, ## Index or timestamp of simulation for  saving, etc.
#   ## # PEFORMANCE [Y] FUNCTION PARAMETERS
#   b0 = .001, ## intercept
#   b1 = .001, ## treatment dummy
#   b2 = .001, ## post intervention dummy
#   # b3 = .001, ## treatment effect (replaced by function b3Func() for dynamic treatment effect)
#   b4 = 1/3, ## spillover of past performance on current performance (how much of treatment effect persists across periods)
#   b5 = .01, ## growth rate (linear effect of time: proportion of time t added to linear combination in yFunc() performance )
#   ## # TREATMENT EFFECT FUNCTION WEIGHTS 
#   w0 = 1.3, ## constant
#   w1 = 0.3, ## linear
#   w2 = -0.12, ## quadratic
#   w2.shift = -10, ## shift quadratic curve to allow gradual treatment effect increase (negative shifts curve rightward)
#   ## # ENDOGENOUS TREATMENT SELECTION [x1[t](y[t-1])] FUNCTION PARAMETERS
#   g0 = .1,
#   g1 = .1,
#   ## # EXPONENTIAL FUNCTION PARAMETERS -- CURRENTLY NOT USED
#   logit.shift = 1,
#   logit.scale = 1,
#   # ## PLOTTING
#   plot.save = TRUE,
#   plot.show = TRUE,
#   plot.wide.w=10,
#   plot.wide.h=5,
#   plot.tall.w=5,
#   plot.tall.h=10,
#   plot.sq.w=9,
#   plot.sq.h=9,
#   plot.dpi=300,
#   ...
# ) {
# ##=================================





#######################################################
## SIMULATION RUNS
#######################################################
n <- 200    ## Number of firms
npds <- 100 ## 100  ## Number of periods
intpd <- round( npds / 4 )

## dynamic treatment effect types (shapes)
effect.types <- c('constant','quadratic','geometric')

# Specify simulation parameters
simlist <- list(
  `1.no.perf.no.gro`=list(
    b4=0, ## b4 = past performance spillover (or persistence)
    b5=0  ## b5 = growth rate
  ),
  `2.no.perf.lo.gro`=list(
    b4=0, b5=.02
  ),
  `3.mi.perf.lo.gro`=list(
    b4=3/4, b5=.02
  )
)
## --- TEST ---
# simlist <- list(
#   no.perf.no.gro=list(
#       b4=0, ## b4 = past performance spillover (or persistence)
#       b5=0  ## b5 = growth rate
#     )
# )
##-------------

## Run Simulation List
sim.id <- round(as.numeric(Sys.time()))
for (i in 1:length(simlist)) {
  
  key <- names(simlist)[i]
  sim <- simlist[[key]]
  
  cat(sprintf('\nScenario label: %s\n\n', key))
  # cat(str(simlist[[key]]))
    
  simlist[[key]]$sim <- runSimSingleInterventionEffectComparison(
    n = n, ## NUMBER OF FIRMS
    npds = npds, ## NUMBER OF PERIODS
    intpd = intpd, ## intervention after first section
    ystart = 0,
    effect.types = effect.types,
    # benchmark.type = 'self', ## 'all', 'self'
    # treat.rule = 'random',
    treat.rule = 'below.benchmark',
    treat.prob = .95, #0.95,  ## 0.6
    treat.threshold = 1/3,  # 0.015
    sim.id = sim.id, ## defaults to timestamp
    ##
    b4 = sim$b4, ## past performance
    b5 = sim$b5, ## growth (linear function of time t)
    ## Dynamic treatment effect function parameters
    w0 = 1.7,  ## constant
    w1 = 0.18,  ## linear
    w2 = -0.005,  ## quadratic
    ## TODO: CHECK w2.shift SENSITIVITY - this shifts quadratic curve several periods to the right so that treatment effect increases slowly  
    w2.shift = -round(.4*intpd),  ## optimal value here is likely a function of the combination of treatment effect function parameters
    ## Plotting
    plot.show = F, ## TRUE
    plot.save = F  ## TRUE
  )
  
}

## Save simulation list as serialized data file
simlist.file <- sprintf('single_intervention_SIMLIST_%s.rds', simlist[[1]]$sim$id)
saveRDS(simlist, file = file.path(dir_plot, simlist.file))



##====================
##  DEBUG
##--------------
View(simlist$`3.mi.perf.lo.gro`$sim$df)

library(tibble)
library(did)

# library(dplyr)
simdf <- simlist$`1.no.perf.no.gro`$sim$df
simdf <- simdf[which(simdf$effect.type=='quadratic'), ]
## Set group name 'gname' field, where 0 = control, # = period of treatment
simdf$gname <- 0
simdf$gname[simdf$group=='treatment'] <- simdf$match_pd[simdf$group=='treatment']

ccattgt <- att_gt(yname = "y", ## "Y",
                  tname = "t",
                  idname = "actor",
                  gname = "gname",
                  xformla = ~c1 + c2 + c3,
                  data = simdf #,
                  # panel = F
                  
)

ccattgt

agg.simple <- aggte(ccattgt, type='simple')
summary(agg.simple)

## DYNAMIC EFFECTS AND EVENT STUDIES
agg.es <- aggte(ccattgt, type = "dynamic")
summary(agg.es)
tidy(agg.es)
ggdid(agg.es)

## CALENDAR TIME EFFECT
agg.ca <- aggte(ccattgt, type = "calendar")
summary(agg.ca)
tidy(agg.es)
ggdid(agg.ca)


## Correlation of simulated to inferred
cor(cbind(simdf$b3[simdf$actor==1][-1], agg.es$att.egt))
## MATPLOT of 
matplot(cbind(simdf$b3[simdf$actor==60][-1], agg.es$att.egt), type='o',pch=1:2)


cor(cbind(simdf$x1[simdf$t==intpd], simdf$y[simdf$t==(intpd-1)]))

# endog <- data.frame(threshold=c(1/2, 1/3,1/4,1/5,1/6),
#                     cor=c(-0.2523,-.2030,-.2086,-.2107,-.2214))

###=========================================
##  Aggregate actor series into 1 total
##------------------------------------------
# actors <- sort(unique(simdf$actor))
# for (i in 1:length(actors)) {
#   actor <- actors[i]
#   
# }

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
  
##
##
##  ***** START HERE *****
##
##

## Timeseries Dataframe - Wide
tsdfw <- tsdf %>% pivot_wider(# id_cols = .(t,group),
                     # id_expand = FALSE,
                     names_from = c(group),
                     # names_prefix = "",
                     # names_sep = "_",
                     # names_glue = NULL,
                     # names_sort = FALSE,
                     # names_vary = "fastest",
                     # names_expand = FALSE,
                     # names_repair = "check_unique",
                     values_from = c('y_mean','y_max','y_min'),   #,
                     # values_fill = NULL,
                     values_fn = mean,
                     # unused_fn = NULL
                     )

# id.cols <- c('t','group')
val.cols <- c('y_mean','y_sum','y_min','y_max','y_sd',
              'x1_mean','x2_mean','x3_mean',
              'c1_mean','c2_mean','c3_mean',
              'b1_mean','b2_mean','b3_mean',
              'u_mean','v_mean'
              )
ts <- unique(tsdf$t)
groups <- unique(tsdf$group)
tsdfw <- data.frame(t=ts,  stringsAsFactors = F)
for (j in 1:length(groups)) {
  id.j <- which( tsdf$group == groups[j] ) 
  # sprintf('%s_%s_%s',ts[i],groups[j],)
  # z <- tsdf[id.ij, ]
  # tsdfw <- cbind(tsdfw, group=rep(groups[j],length(ts)))
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
dat <- tsdfw[,c('treatment_y_mean','control_y_mean','control_y_sd','control_y_sum', 'control_y_min','control_y_max',
                'control_c1_mean','control_c2_mean','control_c3_mean',
                # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
                'control_u_mean','control_v_mean'
                )]
impact_amount <- CausalImpact(dat, pre.period,post.period,
                              alpha=0.05, model.args = list(niter = 5000))
summary(impact_amount)
plot(impact_amount)
bstsfit <- impact_amount$series

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


res.tbl <- cbind(bstsfit[ ,c('point.effect','point.effect.lower','point.effect.upper')],
                 did.res[ ,c('term','event.time','estimate','point.conf.low','point.conf.high')],
                 b3.treat=b3diff$treat,
                 b3.ctrl=b3diff$ctrl,
                 b3.att=b3diff$diff 
                 )
# res.tbl <- res.tbl[ ! is.na(res.tbl$estimate), ]
num.cols <- names(res.tbl)[ ! names(res.tbl) %in% c('term','event.time') ]
for (i in 1:length(num.cols)) {
  res.tbl[ , num.cols[i] ] <- round( as.numeric( res.tbl[ , num.cols[i] ] ), 4)
}
res.tbl4 <- cbind(res.tbl[,4:5],  res.tbl[,-(4:5)] )


matplot(res.tbl4[,c('point.effect','estimate')],type='l',col=c('black','red'),
        main=sprintf('ATT[DGP] = %.3f;  ATT[DiD] = %.3f;  ATT[BSTS] = %.3f',
                     mean(b3diff$diff[intpd:nrow(b3diff)]), agg.es$overall.att, impact_amount$summary$AbsEffect[1]))
legend('topright',legend=c('BSTS','DiD'),col=c('black','red'),lty=c(1,2)) 

plot(impact_amount$model$bsts.model,'coefficients')

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

  



##===========================================
## CHECK BSTS FOR HOW TO HANDLE 
##------------------------------------------

library(bsts)
data(iclaims)

ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
model <- bsts(initial.claims$iclaimsNSA,
               state.specification = ss,
               niter = 1000)

ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
model1 <- bsts(initial.claims$iclaimsNSA,
               state.specification = ss,
               niter = 1000)
plot(model)
plot(model1, "components")  # plot(model1, "comp") works too!
plot(model1, "help")

pred1 <- predict(model1, horizon = 12)
plot(pred1, plot.original = 156)

# Fit a bsts model with expected model size 1, the default.
model2 <- bsts(iclaimsNSA ~ .,
               state.specification = ss,
               niter = 1000,
               data = initial.claims)


# Fit a bsts model with expected model size 5, to include more coefficients.
model3 <- bsts(iclaimsNSA ~ .,
               state.specification = ss,
               niter = 1000,
               data = initial.claims,
               expected.model.size = 5)  # Passed to SpikeSlabPrior.

plot(model2, 'comp')
plot(model3, 'comp')

plot(model2, 'coef')
plot(model3, 'coef')


# ## 
# agg.ca <- aggte(ccattgt, type = "calendar")
# summary(agg.ca)
# ggdid(agg.ca)
# 
# ## 
# agg.gr <- aggte(ccattgt, type = "group")
# summary(agg.gr)
# ggdid(agg.gr)

##---------------
## // end debug
##======================




##--------------------------------------
##  COMBINED PLOT FACET GRID
##   - Simulation Scenario by dynamic treatment effect shape
##--------------------------------------
dfx.t0.summary <- data.frame(stringsAsFactors = F)
dfx.att <- data.frame(stringsAsFactors = F)
for (i in 1:length(simlist)) {
  dfx.sim.i <- simlist[[i]]$sim$df.t0.summary
  dfx.sim.i$scenario <- names(simlist)[i]
  dfx.t0.summary <- rbind(dfx.t0.summary,  dfx.sim.i)
  ##
  dfx.att.i <- simlist[[i]]$sim$df.att
  dfx.att.i$scenario <- names(simlist)[i]
  dfx.att <- rbind(dfx.att,  dfx.att.i)
}

pall <- ggplot(dfx.t0.summary, aes(x=t0, y=med, color=group)) +
  geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
  geom_line(size=1.2) +
  # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
  ylab('Y') +
  xlim(c(-intpd+2, npds-intpd-2)) +
  facet_grid( scenario ~ effect.type ) +
  theme_bw() + theme(legend.position='top') 
print(pall)
pall.file <- sprintf('single_intervention_staggered_DiD_COMBINED_%s.png',sim.id)   ### ******* ??? *************
ggsave(filename=file.path(dir_plot, pall.file), plot=pall,
       width=10,heigh=10,dpi=300,units='in')

## Average Treatment Effect on the Treated (empirical difference between treated and controlled by day from simulation)
patt <- ggplot(dfx.att, aes(x=t0, y=y.att)) +
  # geom_ribbon(aes(ymin=cl,ymax=cu,fill=group), alpha=.15, size=.01, lty=1) +
  geom_line(size=1.2) +
  # geom_point(aes(x=t, y=min),pch=1,alpha=.3) + geom_point(aes(x=t,y=max),pch=1,alpha=.3) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
  ylab('ATT') +
  xlim(c(-intpd+2, npds-intpd-2)) +
  facet_grid( scenario ~ effect.type ) +
  theme_bw() + theme(legend.position='top') 
print(patt)

########################## END ##########################################
