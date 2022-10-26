#####################################################
##
##   Dynamic Causal Inference Simulations
##
##    Run script for simulations of
##     - internal intervention (self-selection) endogeneity
##     - ...
##
#####################################################


## Directories
dir_proj <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'   
dir_plot <- file.path(dir_proj, 'plots')
dir_r <- file.path(dir_proj,'R')



##==============================
## Load simulation functions
##------------------------------
# source(file.path(dir_r,'internal_intervention_sim.R'))
source(file.path(dir_r,'internal_intervention_sim_vec.R'))  ## vectorized -- in development


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
n <- 150    ## Number of firms
npds <- 80  ## Number of periods
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
    
  simlist[[key]]$sim <- runSimInternalInterventionEffectComparison(
    n = n, ## NUMBER OF FIRMS
    npds = npds, ## NUMBER OF PERIODS
    intpd = intpd, ## intervention after first section
    ystart = 0,
    effect.types = effect.types,
    # benchmark.type = 'self', ## 'all', 'self'
    treat.rule = 'below.benchmark',
    treat.threshold = 0.02,  # 0.015
    sim.id = sim.id, ## defaults to timestamp
    ##
    b4 = sim$b4, ## past performance
    b5 = sim$b5, ## growth (linear function of time t)
    ## Dynamic treatment effect function parameters
    w0 = 1.7,  ## constant
    w1 = 0.18,  ## linear
    w2 = -0.005,  ## quadratic
    ## TODO: CHECK w2.shift SENSITIVITY - this shifts quadratic curve several periods to the right so that treatment effect increases slowly  
    w2.shift = -round(.2*intpd),  ## optimal value here is likely a function of the combination of treatment effect function parameters
    ## Plotting
    plot.show=TRUE,
    plot.save=TRUE
  )
  
}

## Save simulation list as serialized data file
simlist.file <- sprintf('internal_intervention_SIMLIST_%s.rds', simlist[[1]]$sim$id)
saveRDS(simlist, file = file.path(dir_plot, simlist.file))



##--------------------------------------
##  COMBINED PLOT FACET GRID
##   - Simulation Scenario by dynamic treatment effect shape
##--------------------------------------
dfx.t0.summary <- data.frame(stringsAsFactors = F)
for (i in 1:length(simlist)) {
  dfx.sim.i <- simlist[[i]]$sim$df.t0.summary
  dfx.sim.i$scenario <- names(simlist)[i]
  dfx.t0.summary <- rbind(dfx.t0.summary,  dfx.sim.i)
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
pall.file <- sprintf('internal_intervention_staggered_DiD_COMBINED_%s.png',sim.id)
ggsave(filename=file.path(dir_plot, pall.file), plot=pall,
       width=10,heigh=10,dpi=300,units='in')



########################## END ##########################################
