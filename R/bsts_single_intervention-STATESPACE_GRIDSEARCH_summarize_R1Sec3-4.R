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

## Directories
dir_proj <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'
dir_ext <- 'D:\\BSTS_external'
dir_plot <- file.path(dir_proj, 'plots')
dir_r <- file.path(dir_proj,'R')

setwd(dir_proj)


##==============================
##  file prefix for saving images, writing outputs, etc.
##-----------------------------
prefix <- 'single-int-st-sp-summ_'

##==============================
## Load simulation functions
##------------------------------
# source(file.path(dir_r,'internal_intervention_sim.R'))
source(file.path(dir_r,'single_intervention_sim_vec.R'))  ## vectorized -- in development


# ## DEBUG ###
# effect.types=c('constant','quadratic','geometric') 
# sim.id=round(10*as.numeric(Sys.time()))
# plot.show=F
# plot.save=F
# ##----------

# '__GRIDSEARCH_output__16685974280_j1g1h2k1i2ii2m2n1.rds'
# 'single_intervention-statespace-config__BSTS_inclusion_probs_ss8_j1g1h2k1i2ii2m2n1_geometric_16685974280'

##==================================
## GET SIMLIST FROM EXTERNAL STORAGE
##----------------------------------
# '__GRIDSEARCH_output__16690657090_g1h1i1j1.rds'
sim.id <- 16690657090
##
g <- 1
h <- 1
i <- 1
j <- 1
# i <- 1
# ii<- 1
# m <- 1
# n <- 1
key.vec <- sprintf('g%s|h%s|i%s|j%s',
                    g,  h,  i, j)
key <- paste(key.vec, collapse = '|')
key.strip <- gsub('[|]','',key.vec)
# key <- 'j1g1h1k1i1ii1m1n1'
##
file.name <- sprintf('__GRIDSEARCH_output__%s_%s.rds', sim.id, key.strip)

file.exists(file.path(dir_ext, file.name))

# ( sim.files <- dir(dir_ext, pattern = as.character(sim.id)) )
# 
# k2files <- sim.files[ grep(key, x = sim.files, ignore.case = T, perl = T) ]
simi <- readRDS(file.path(dir_ext, file.name))

#### DEBUG #####
(i.bsts <- simi$compare$bsts$quadratic[[1]]$CausalImpact$model$bsts.model)

par(mfrow=c(3,3));for(i in 1:ncol(i.bsts$predictors))hist(i.bsts$predictors[,i],main=colnames(i.bsts$predictors)[i])




##########################################
##
##  Summarize Simulation Scenarios
##
##########################################
##
compdf <- data.frame(stringsAsFactors = F)
##
sim.files <- dir(dir_ext, pattern =sprintf('__GRIDSEARCH_output__%s',sim.id))
##
sim.files.todo <- sim.files  #[ ! dd %in% ]


# graphics.off()

effect.types <- c('constant','quadratic','geometric')
# simlist <- list()
for (i in 1:length(sim.files.todo)) 
{
  sim.file.i <- sim.files.todo[i]
  cat(sprintf('\n%s_%s\n',i,sim.file.i))
  simi <- readRDS(file.path(dir_ext, sim.file.i))
  
  # simlist[[i]] <- simi
  # simi <- simlist[[i]]
  intpd <- simi$intpd
  

  ## LOOP OVER DYNAMIC EFFECT SHAPES
  for (k in 1:length(effect.types)) 
  {   
      effect.type <- effect.types[k]

      ## LOOP OVER BSTS STATE SPACE CONFIGURATIONS
      bsts.res.list <- simi$compare$bsts[[ effect.type ]]
      for (j in 1:length(bsts.res.list)) {
        
        res.tbl <- simi$compare$res.tbl[[ effect.type ]][[ j ]]
        
        if (is.null(res.tbl) | length(res.tbl)==0) {
          cat(sprintf('\nNo res.tbl: Skipping j=%s\n',j))
          next
        }
        
        bsts.res.j <- bsts.res.list[[ j ]]
      
        agg.es <- simi$compare$did[[ effect.type ]]$agg.es
        impact_amount <- simi$compare$bsts[[ effect.type ]][[ j ]]$CausalImpact
        
        ## Treatment Effect
        gatt.b3 <- mean( as.numeric( res.tbl$b3.att[intpd:nrow(res.tbl)] ) )
        gatt.did <- agg.es$overall.att
        gatt.bsts <- impact_amount$summary$AbsEffect[1]
        
        ## TODO: Change to loop j over BSTS configs list
        ## BSTS State Specifications
        # bsts.state.spec.j <- simi$compare$bsts[[ j ]][[ effect.type ]]$CausalImpact$model$bsts.model$state.specification  ## with all bsts state specs
        bsts.state.spec.j <- simi$compare$bsts[[ effect.type ]][[j]]$CausalImpact$model$bsts.model$state.specification
        bsts.state.comps <- sapply(bsts.state.spec.j,function(z)class(z)[1])
        
        
        # # ############################################
        # # ########################### TODO ### ### ###
        # # ############################################
        # ## Append results to output list
        # simlist[[key]]$compare$did[[effect.type]]$attgt <- ccattgt
        # simlist[[key]]$compare$did[[effect.type]]$agg.simple <- agg.simple
        # simlist[[key]]$compare$did[[effect.type]]$agg.es <- agg.es
        # simlist[[key]]$compare$did[[effect.type]]$self.select.cor <- self.select.cor
        # ##
        # simlist[[key]]$compare$bsts[[effect.type]][[ h ]]$CausalImpact <- impact_amount
        # simlist[[key]]$compare$bsts[[effect.type]][[ h ]]$cumu.pred.error <-  cumsum(colSums(abs(bsts.pred.er)))
        # ##
        # simlist[[key]]$compare$res.tbl[[effect.type]][[ h ]] <- res.tbl
        # simlist[[key]]$compare$att.err.tbl[[effect.type]][[ h ]] <- errdf
        # simlist[[key]]$compare$att.err.mean.bsts[[effect.type]][[ h ]] <- mean(errdf$error[errdf$method=='BSTS'],na.rm = T)
        # simlist[[key]]$compare$att.err.mean.did[[effect.type]][[ h ]]  <- mean(errdf$error[errdf$method=='DiD'],na.rm = T)
        # simlist[[key]]$compare$att.err.sd.bsts[[effect.type]][[ h ]] <- sd(errdf$error[errdf$method=='BSTS'],na.rm = T)
        # simlist[[key]]$compare$att.err.sd.did[[effect.type]][[ h ]]  <- sd(errdf$error[errdf$method=='DiD'],na.rm = T)
        # # ###############################
        # # ############################
        # # #########################
        
        ## Simulation scenario (row) dataframe
        df.ijk <- data.frame(
          n = simi$n,
          npds = simi$npds,
          intpd = simi$intpd,
          effect.type = effect.type,
          noise.level = simi$noise.level,
          prior.sd.scenario = simi$prior.sd.scenario,
          treat.rule = simi$treat.rule,
          treat.prob = simi$treat.prob,
          treat.threshold = simi$treat.threshold,
          b4 = simi$b4,
          b5 = simi$b5,
          b9 = simi$b9,
          dgp.nseasons = simi$dgp.nseasons,
          dgp.freq = simi$dgp.freq,
          seasonality = simi$seasonality,
          ##
          bsts.state.comps = paste(bsts.state.comps, collapse = '|'),
          bsts.state.comps.len = length(bsts.state.comps),
          ## Performance Measures (OVERALL GENERAL ATT)
          gatt.b3 = gatt.b3,
          gatt.did = gatt.did,
          gatt.bsts = gatt.bsts,
          ##
          bsts.cumu.pred.err.mean =  mean( simi$compare$bsts[[effect.type]][[ j ]]$cumu.pred.error, na.rm = T),
          bsts.cumu.pred.err.sum =  sum( simi$compare$bsts[[effect.type]][[ j ]]$cumu.pred.error, na.rm = T),
          # bsts.cumu.pred.err.mean =  mean( simi$compare$bsts[[ j ]][[effect.type]]$cumu.pred.error, na.rm = T), ## with all bsts state specs
          # bsts.cumu.pred.err.sum =  sum( simi$compare$bsts[[ j ]][[effect.type]]$cumu.pred.error, na.rm = T), ## with all bsts state specs
          ###
          did.b3.diff = (gatt.did - gatt.b3),
          bsts.b3.diff = (gatt.bsts - gatt.b3),
          ##
          att.err.mean.bsts = simi$compare$att.err.mean.bsts[[effect.type]][[ j ]],
          att.err.mean.did = simi$compare$att.err.mean.did[[effect.type]][[ j ]] ,
          att.err.sd.bsts = simi$compare$att.err.sd.bsts[[effect.type]][[ j ]],
          att.err.sd.did = simi$compare$att.err.sd.did[[effect.type]][[ j ]],
          ##
          self.select.cor.did = simi$compare$did[[effect.type]]$self.select.cor,
          ##
          rand.seed = simi$rand.seed,
          file = sim.file.i
        )
        
        # ## AUTOCORRELATION VALUES
        # for (g in 1:length(dgp.ars)) {
        #   dgp.ar <- dgp.ars[[ g ]]
        #   
        #   ## SEASONALITY
        #   for (h in 1:length(seasonalities)) {
        #     seasonality <- seasonalities[[ h ]]
        #     
        #     ## ENDOGENEITY (SELF-SELECTION)
        #     for (i in 1:length(treat.rules)) {
        #       treat.rule <- treat.rules[[ i ]]
        #       
        #       ##--------------------------------------------
        #       ## Setup state space configurations
        #       ##--------------------------------------------
        #       bsts.state.specs <- list()
        #       ## PRIOR NOISE / UNCERTAINTY (PRIOR STDEV)
        #       for (j in 1:length(prior.sd.scenarios)) {
        #         prior.sd.scenario <- prior.sd.scenarios[[ j ]]
        #         
        #         ## STATE SPACE COMPONENTS CONFIGURATION
        #         for (k in 1:length(st.sp.lists)) {
        #           st.sp.vec <- st.sp.lists[[ k ]]
        
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
    
  }
  
}

compdf$dgp.freq[is.na(compdf$dgp.freq)] <- 0
compdf$dgp.nseasons[is.na(compdf$dgp.nseasons)] <- 0

compdf$bsts.state.comps.f <- as.factor(compdf$bsts.state.comps)
compdf$bsts.state.comps.i <- as.numeric(compdf$bsts.state.comps.f)
compdf$bsts.state.comps.c <- LETTERS[ compdf$bsts.state.comps.i ]

write.csv(compdf, file = 'STATESPACE_GRIDSEARCH_ATT_simulation_summary_df.csv',row.names = F)
# write.csv(compdf, file = 'STATESPACE_GRIDSEARCH_ATT_simulation_summary_df_v2.csv',row.names = F)
saveRDS(compdf, file = 'STATESPACE_GRIDSEARCH_ATT_simulation_summary_df_v2.rds')


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
# heat.x.cols <- c('bsts.state.comps','effect.type')
heat.x.cols <- c('bsts.state.comps.c','bsts.state.comps.i') ## prior.sd.scenario
heat.y.cols <- c('b9', 'effect.type','dgp.nseasons','prior.sd.scenario','treat.threshold','treat.rule','b5')
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
##
col.scale.vals <- c(
  min(compdf.stack$b3.diff[compdf.stack$b3.diff < 0], na.rm=T),
  median(compdf.stack$b3.diff[compdf.stack$b3.diff < 0], na.rm=T),
  0,
  median(compdf.stack$b3.diff[compdf.stack$b3.diff > 0], na.rm=T),
  max(compdf.stack$b3.diff[compdf.stack$b3.diff > 0], na.rm=T)
)
## PLOT FACET HEATMAP
compdf.stack$b3.diff.alog <- sapply(compdf.stack$b3.diff, function(x) {
  log(abs(x) + .001) * ifelse(x>0, 1, -1) 
} )
ggplot(compdf.stack, aes(factor(heat.x), factor(heat.y), fill= b3.diff.alog )) +
  geom_tile() + facet_wrap( . ~ factor(stats.type)) +
  xlab(paste(heat.x.cols,collapse = ' | ')) + 
  ylab(paste(heat.y.cols,collapse=' | ')) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6)) +
  scale_fill_gradient2() + theme_bw()
  # scale_fill_gradientn(colours=c('red','yellow','white','cyan','blue'), 
  #                      # values=rescale(c(-1,0-.Machine$double.eps,0,0+.Machine$double.eps,1)),
  #                      values=rescale(col.scale.vals)
  #                      )
ggsave(filename = sprintf('grid_search_heatmap_did_bsts_%s.png',sim.id), width = 20, height = 20, units = 'in', dpi = 300)
  # scale_fill_gradient(low="white", high="blue") 






#################################
## SUMMARY OF SCENARIO GRIDSEARCH
##  COMPARISON OF BSTS VS DID
#################################
heat.ys <- sort(unique(compdf$heat.y))
heat.xs <- sort(unique(compdf$heat.x))
z <- data.frame()
for (i in 1:length(heat.ys)) {
  cat(sprintf('\n i=%s ',i))
  scenario <- heat.ys[i]
  
  best.bsts.bias <- NA
  best.bsts.state.space <- NA
  best.heat.x <- NA
  best.j <- NA
  for (j in 1:length(heat.xs)) {
    
    stconfig <- heat.xs[j]
    idx <- which( compdf$heat.y == scenario & compdf$heat.x == stconfig )
    
    if (length(idx)>0) {
      dfij <- compdf[idx[1], ]
      
      if ( is.na(best.bsts.bias) | abs(dfij$bsts.b3.diff) < abs(best.bsts.bias) ) {
        best.bsts.bias <- dfij$bsts.b3.diff
        best.bsts.state.space <- dfij$bsts.state.comps
        best.j <- j
        best.heat.x <- stconfig
      }

    }
    
  }
  
  
  tmpdf <- data.frame(i=i, j=best.j, heat.y=scenario, heat.x=best.heat.x,
                      bsts.state.space=best.bsts.state.space,
                      bsts.bias=best.bsts.bias, did.bias=dfij$did.b3.diff )
  z <- rbind(z, tmpdf)
  
}

View(z)
z$bsts.better <- apply(z[,c('bsts.bias','did.bias')],1,function(x)abs(x[1]) < abs(x[2]))


df.gs.summary <- z %>% 
  group_by(bsts.better, bsts.state.space) %>% 
  summarize(
    n=n(),
    bsts.mean= round( mean(bsts.bias, na.rm=T), 3),
    bsts.sd= round( sd(bsts.bias, na.rm=T), 3),
    did.mean= round( mean(did.bias, na.rm=T), 3),
    did.sd= round( sd(did.bias, na.rm=T), 3),
    abs.bsts.did = round( mean(abs(bsts.bias)) / mean(abs(did.bias)), 3)#,
    # bsts.did = mean(bsts.bias) / mean(did.bias)
  ) %>%
  arrange(abs.bsts.did)

write.csv(df.gs.summary, file = 'bsts_did_simulation_gridsearch_ATT_bias_summary.csv',row.names = F)







##==================================
## GET SIMLIST FROM EXTERNAL STORAGE
##----------------------------------
# '__GRIDSEARCH_output__16690657090_g1h1i1j1.rds'
sim.id <- 16690657090
##
g <- 1  ## autocorrelation  0, .1, .2, .4
# h <- 1  ## Seasonality: TRUE (12 pds), FALSE
h <- 2
i <- 1  ## self-selection:  No (random); Yes (below median)
j <- 1  ## priod SD: 0.01, 0.1
# i <- 1
# ii<- 1
# m <- 1
# n <- 1
key.vec <- sprintf('g%s|h%s|i%s|j%s',
                   g,  h,  i, j)
key <- paste(key.vec, collapse = '|')
key.strip <- gsub('[|]','',key.vec)
# key <- 'j1g1h1k1i1ii1m1n1'
##
file.name <- sprintf('__GRIDSEARCH_output__%s_%s.rds', sim.id, key.strip)

file.exists(file.path(dir_ext, file.name))

# ( sim.files <- dir(dir_ext, pattern = as.character(sim.id)) )
# 
# k2files <- sim.files[ grep(key, x = sim.files, ignore.case = T, perl = T) ]
simi <- readRDS(file.path(dir_ext, file.name))

effect.types <- c('constant','geometric','quadratic')
st.sp.id <- 5
for (effect.type in effect.types) {
  res.tbl.filename <- sprintf('bsts_did_gridsearch_res-tbl_%s_%s.csv',effect.type,key.strip)
  xdf <- simi$compare$res.tbl[[ effect.type ]][[ st.sp.id ]]
  print(xdf)
  write.csv(xdf, res.tbl.filename)
}


#### DEBUG #####
(i.bsts <- simi$compare$bsts$quadratic[[1]]$CausalImpact$model$bsts.model)

par(mfrow=c(3,3));for(i in 1:ncol(i.bsts$predictors))hist(i.bsts$predictors[,i],main=colnames(i.bsts$predictors)[i])
#################



simi$compare$did$quadratic$attgt


intpd <- 67
res5 <- residuals( simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model)[,1:(intpd-1)]
res8 <- residuals( simi$compare$bsts$quadratic[[8]]$CausalImpact$model$bsts.model)[,1:(intpd-1)]
qqdist(res5); qqdist(res8)
par(mfrow=c(1,2))
AcfDist(res5)
AcfDist(res8)

burnin  <- 200
bsts.iter <- 1000
err5 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model$one.step.prediction.errors[burnin:bsts.iter, 1:(intpd-1)]
err8 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model$one.step.prediction.errors[burnin:bsts.iter, 1:(intpd-1)]
err15 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model$one.step.prediction.errors[burnin:bsts.iter, 1:(intpd-1)]
par(mfrow=c(1,2))
plot(rowMeans(err5 ^ 2), type = "l"); plot(rowMeans(err8 ^ 2), type = "l")
plot(rowMeans(err15 ^ 2), type = "l")


coefmc5 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model$coefficients[burnin:bsts.iter,]
par(mfrow=c(3,3))
for (i in 1:ncol(coefmc5)) {
  attr <- dimnames(coefmc5)[[2]][i]
  hist(coefmc5[,i], main=attr, xlab='')
}

## CODA DIAGNOSTICS
library(coda)
geweke.diag(rowMeans(err5 ^ 2), .1, .5)
geweke.diag(coefmc5, .1, .5)

heidel.diag(rowMeans(err5 ^ 2))
heidel.diag(coefmc5)


CompareBstsModels(list(
  # simi$compare$bsts$quadratic[[1]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[2]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[3]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[4]]$CausalImpact$model$bsts.model,
  `Trig`=simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[6]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[7]]$CausalImpact$model$bsts.model,
  `Trig|LocalLevel`=simi$compare$bsts$quadratic[[8]]$CausalImpact$model$bsts.model#,
  # simi$compare$bsts$quadratic[[9]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[10]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[11]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[12]]$CausalImpact$model$bsts.model,
  # simi$compare$bsts$quadratic[[13]]$CausalImpact$model$,,
  # simi$compare$bsts$quadratic[[14]]$CausalImpact$model$bsts.model,
  # `Trig|SemilocalLinearTrend`=simi$compare$bsts$quadratic[[15]]$CausalImpact$model$bsts.model
), burn = 100)



bsts5 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model
bsts8 <- simi$compare$bsts$quadratic[[8]]$CausalImpact$model$bsts.model
bsts15 <- simi$compare$bsts$quadratic[[15]]$CausalImpact$model$bsts.model

plot(bsts5)

PlotBstsComponents(bsts5, burn = 100)
PlotBstsComponents(bsts8, burn = 100)
PlotBstsComponents(bsts15, burn = 100) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
# PlotBstsState(bsts.model, main=sprintf('BSTS State: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsResiduals(bsts.model, main=sprintf('BSTS Residuals: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsPredictionErrors(bsts.model, main=sprintf('BSTS Pred.Err: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsForecastDistribution(bsts.model, main=sprintf('BSTS Forecast Dist: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsSize(bsts5)
# PlotBstsSize(bsts8)
# PlotBstsSize(bsts15)
# 
# PlotBstsState(bsts5)
# PlotBstsState(bsts8)
# PlotBstsState(bsts15)
# 
# PlotBstsResiduals(bsts5)
# PlotBstsPredictionErrors(bsts5)
# PlotBstsForecastDistribution(bsts5)
##

# ## BSTS model for Dynamic Regression
# bsts.model <- bsts(y.pre.treat.NAs.post.treat,
#                    state.specification = st.sp,
#                    niter = 5000)

## Use BSTS prediction of counterfactual to estimate CausalImpact
impact_amount <- CausalImpact(bsts.model=bsts.model,
                              post.period.response = post.period.response,
                              alpha=0.05, model.args = list(niter = bsts.niter))
# ##
# summary(impact_amount)
# summary(impact_amount$model$bsts.model)
# plot(impact_amount)

# summary(impact_amount)
# png(filename=sprintf('single_intervention_BSTS_CausalImpact_plot_%s_%s_%s.png',
#                         key,effect.type,sim.id))
p.bsts.impact.all <- plot(impact_amount, c('original','pointwise','cumulative')) # pointwise','cumulative







































######################################################
######################################################
# '__GRIDSEARCH_output__16690657090_g1h1i1j1.rds'
sim.id <- 16690657090
##
g <- 1  ## autocorrelation  0, .1, .2, .4
# h <- 1  ## Seasonality: TRUE (12 pds), FALSE
h <- 2
i <- 1  ## self-selection:  No (random); Yes (below median)
j <- 1  ## priod SD: 0.01, 0.1
# i <- 1
# ii<- 1
# m <- 1
# n <- 1
key.vec <- sprintf('g%s|h%s|i%s|j%s',
                   g,  h,  i, j)
key <- paste(key.vec, collapse = '|')
key.strip <- gsub('[|]','',key.vec)
# key <- 'j1g1h1k1i1ii1m1n1'
##
file.name <- sprintf('__GRIDSEARCH_output__%s_%s.rds', sim.id, key.strip)

file.exists(file.path(dir_ext, file.name))

# ( sim.files <- dir(dir_ext, pattern = as.character(sim.id)) )
# 
# k2files <- sim.files[ grep(key, x = sim.files, ignore.case = T, perl = T) ]
simi <- readRDS(file.path(dir_ext, file.name))

npds <- 100
noise.level <- 1.2

simdf <- simi$sim$df

simdf$gname <- 0
simdf$gname[simdf$group=='treatment'] <- simdf$match_pd[simdf$group=='treatment']


tsdf <- simdf %>%
  dplyr::filter( ! is.na(match_id)) %>%
  dplyr::group_by(t, t.post.intpd, match_pd, gname, group) %>%
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

## MAKE WIDE TIMESERIES FOR treatment,control groups in n periods
val.cols <- c('y_mean','y_sum','y_min','y_max','y_sd',
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
dat <- tsdfw[,c('treatment_y_sum','control_y_sum','control_y_sd',#'control_y_sum', 'control_y_min','control_y_max',
                'control_c1_mean','control_c2_mean','control_c3_mean'#,
                # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
                # 'control_u_mean','control_v_mean'
)]
## Train on y pre-treatment but NA's post-treatment
y.pre.treat.NAs.post.treat <- c(dat$treatment_y_sum[1:(intpd-1)], rep(NA,npds-intpd+1))
## Then use the post-treatment response for causal impact estimation
post.period.response <- dat$treatment_y_sum[intpd:npds]
## Covariates (predictors) - Matrix for "formula = y ~ predictors" argument
predictors <- dat[, ! names(dat) %in% 'treatment_y_sum'] ## remove response; convert to matrix
# ## Covariates (predictors) - Dataframe for "data" argument
# predictors <- as.matrix(predictors) 

## ADD temporal trend to covariates
predictors$covar_temp_trend <- (1:npds) + rnorm(npds, 0, noise.level)   ## * simlist$`0.rand.base`$b5


st.sp <- list()
st.sp <- AddTrig(st.sp, y.pre.treat.NAs.post.treat, period = 12, frequencies = 1)
# st.sp <- AddAr(st.sp, y.pre.treat.NAs.post.treat, lags=2)
st.sp <- AddAutoAr(st.sp, y.pre.treat.NAs.post.treat)

bsts.fit <- bsts(y.pre.treat.NAs.post.treat ~ . ,
                 state.specification = st.sp,
                 data = predictors,
                 niter = 10000)



  
##
plot(bsts.fit)
PlotBstsComponents(bsts.fit) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
# PlotBstsState(bsts.model, main=sprintf('BSTS State: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsResiduals(bsts.model, main=sprintf('BSTS Residuals: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsPredictionErrors(bsts.model, main=sprintf('BSTS Pred.Err: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
# PlotBstsForecastDistribution(bsts.model, main=sprintf('BSTS Forecast Dist: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
PlotBstsSize(bsts.fit)
##


intpd <- 67
res <- residuals( bsts.fit )[,1:(intpd-1)]
# res8 <- residuals( simi$compare$bsts$quadratic[[8]]$CausalImpact$model$bsts.model)[,1:(intpd-1)]
qqdist(res); qqdist(res)
par(mfrow=c(1,1))
AcfDist(res)
# AcfDist(res8)

burnin  <- 2000
bsts.iter <- 10000
err <- bsts.fit$one.step.prediction.errors[burnin:bsts.iter, 1:(intpd-1)]
# err8 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model$one.step.prediction.errors[burnin:bsts.iter, 1:(intpd-1)]
# err15 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model$one.step.prediction.errors[burnin:bsts.iter, 1:(intpd-1)]
par(mfrow=c(1,1))
plot(rowMeans(err ^ 2), type = "l") #; plot(rowMeans(err8 ^ 2), type = "l")
# plot(rowMeans(err15 ^ 2), type = "l")


coefmc <- bsts.fit$coefficients[burnin:bsts.iter,]
par(mfrow=c(3,3))
for (i in 1:ncol(coefmc5)) {
  attr <- dimnames(coefmc)[[2]][i]
  hist(coefmc[,i], main=attr, xlab='')
}

## CODA DIAGNOSTICS
library(coda)
geweke.diag(rowMeans(err ^ 2), .1, .5)
geweke.diag(coefmc, .1, .5)

heidel.diag(rowMeans(err ^ 2))
heidel.diag(coefmc)


CompareBstsModels(list(
  bsts.fit,
  bsts.fit
), burn = 2000)



bsts5 <- simi$compare$bsts$quadratic[[5]]$CausalImpact$model$bsts.model
bsts8 <- simi$compare$bsts$quadratic[[8]]$CausalImpact$model$bsts.model
bsts15 <- simi$compare$bsts$quadratic[[15]]$CausalImpact$model$bsts.model

plot(bsts5)

PlotBstsComponents(bsts.fit, burn = 2000)


######################################################
######################################################





# compdf.scaled <- compdf.stack
# col.idx <- which(names(compdf) %in% c('b3.diff')) #which(unname(sapply(compdf.scaled,function(x)is.numeric(x))))
# compdf.scaled[,col.idx] <- scale(compdf.scaled[,col.idx])
# compdf.scaled.w <- spread(compdf.scaled, 'effect.type', 'b3.diff', fill = NA, convert = FALSE, drop = TRUE, sep = NULL)
# cmat <- as.matrix(compdf.scaled.w[,c('bsts.state.comps.i','b3.diff')])

# heat.ys <- unique(compdf$heat.y)
# compdf.wide <- data.frame()
# ## Loop over all Scenarios (simulation data-structure combinations )
# for (i in 1:length(heat.ys)) {
#   cat(sprintf('\ni=%s, %s\n',i,heat.ys[i]))
#   heat.y <- heat.ys[ i ]
#   compdf.i <- compdf[compdf$heat.y == heat.y,  ]
#   heat.xs <- unique(compdf.i$heat.x)
#   
#   z <- rbind(
#     data.frame(heat.x=compdf.i$heat.x, heat.y=heat.y, stat.type='bsts', b3.diff=compdf.i$bsts.b3.diff),
#     data.frame(heat.x=compdf.i$heat.x, heat.y=heat.y, stat.type='did', b3.diff=compdf.i$did.b3.diff)
#   )
#   names(z)[ncol(z)] <- heat.y
#   
#   if (i == 1) {
#     compdf.wide <- z
#   } else {
#     .dfi <- z[,c('heat.x','heat.y')] ## as data.frame
#     # names(.dfi)[1] <- heat.y
#     # compdf.wide <- cbind(compdf.wide, .dfi )
#     compdf.wide <- merge(x=compdf.wide, y=.dfi, by='heat.x', all=T)
#   }
#   # ## Loop over BSTS State Space Config for each DGP Scenario
#   # for (j in 1:length(heat.xs)) {
#   #   heat.x.j <- heat.xs[ j ]
#   #   compdf.y.i[compdf.y.i$heat.x==heat.x.j, ]
#   # }
#   # compdf.wide <- cbind(compdf.wide, data.frame(compdf.y.i$bsts.b3.diff) ) ## ***TODO***
#   # names(compdf.wide)[ncol(compdf.wide)] <- st.sps[i]
# }



# 
# ##===========================================================
# ##
# ## PLOT HEATMAP OF STATE SPACE CONFIGURATIONS
# ##
# ##------------------------------------------------------------
# 
# heat.xs <- unique(compdf$heat.x)
# compdf.wide <- data.frame()
# ## Loop over all Scenarios (simulation data-structure combinations )
# for (i in 1:length(heat.xs)) {
#   cat(sprintf('\ni=%s, %s\n',i,heat.xs[i]))
#   
#   heat.x <- heat.xs[ i ]
#   compdf.i <- compdf[compdf$heat.x == heat.x,  ]
#   heat.ys <- unique(compdf.i$heat.y)
# 
#   z <- rbind(
#     data.frame(heat.y=compdf.i$heat.y, stat.type='bsts', b3.diff=compdf.i$bsts.b3.diff),
#     data.frame(heat.y=compdf.i$heat.y, stat.type='did', b3.diff=compdf.i$did.b3.diff)
#   )
#   names(z)[ncol(z)] <- heat.x
# 
#   if (i == 1) {
#     compdf.wide <- z
#   } else {
#     .dfi <- z ## as data.frame
#     # names(.dfi)[which(names(.dfi)=='b3.diff')] <- 
#     # names(.dfi)[1] <- heat.y
#     # compdf.wide <- cbind(compdf.wide, .dfi )
#     compdf.wide <- merge(x=compdf.wide, y=.dfi, by=c('heat.y','stat.type'), all=T)
#   }
#   # ## Loop over BSTS State Space Config for each DGP Scenario
#   # for (j in 1:length(heat.xs)) {
#   #   heat.x.j <- heat.xs[ j ]
#   #   compdf.y.i[compdf.y.i$heat.x==heat.x.j, ]
#   # }
#   # compdf.wide <- cbind(compdf.wide, data.frame(compdf.y.i$bsts.b3.diff) ) ## ***TODO***
#   # names(compdf.wide)[ncol(compdf.wide)] <- st.sps[i]
# }
# 
# # View(compdf.wide)
# stat.type <- 'bsts'
# ##
# compdf.wide.stat <- compdf.wide[compdf.wide$stat.type==stat.type, ] 
# mat <- as.matrix( compdf.wide.stat[,-c(1:2)] )
# mat.scale <- scale(mat)
# rownames(mat.scale) <- compdf.wide.stat$heat.y
# xdist <- dist(mat.scale)
# xcl <- hclust(xdist)
# xdend <-as.dendrogram(xcl)
# xdendord <- order.dendrogram(xdend)
# compdf$heat.y.ord <- factor(x=compdf$heat.y,
#                               levels=compdf.wide.stat$heat.y[ xdendord ], #xdendord
#                               ordered = TRUE)
# ##
# txdist <- dist(t(mat.scale))
# txcl <- hclust(txdist)
# txdend <-as.dendrogram(txcl)
# txdendord <- order.dendrogram(txdend)
# compdf$heat.x.ord <- factor(x=compdf$heat.x,
#                             levels=levels(factor(compdf$heat.x))[ txdendord ], #xdendord
#                             ordered = TRUE)
# ##
# compdf$bsts.b3.diff.alog <- sapply(compdf$bsts.b3.diff, function(x) {
#   log(abs(x) + .01) * ifelse(x>0, 1, -1) 
# } )
# ##
# 
# # hmout <- heatmap(mat.scale, keep.dendro = T)
# # hmout
# 
# ##
# p.bsts.heatmap <- ggplot(compdf, aes(y=heat.y.ord, x=heat.x.ord, fill= bsts.b3.diff.alog )) + 
#   geom_tile() + 
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
#   ylab('Data Generating Process (DGP) Scenarios')  + xlab('BSTS State Space Configurations') + 
#   scale_fill_gradient2() + theme_bw()
# print(p.bsts.heatmap)
# ggsave(filename='heatmap_statespace_gridsearch_hclust__ggplot.png', width = 10, height = 10, units = 'in', dpi = 300)
# 
# library(RColorBrewer)
# 
# # par(mar=c(3,3,3,6))
# png(filename = 'heatmap_hclust_statespace_gridsarch__heatmap-func.png', width = 12, height = 10, units = 'in', res=300)
# hdend <- heatmap(mat.scale,col= RColorBrewer::brewer.pal(11,'RdBu'), keep.dendro = T, margins=c(5,12))
# dev.off()
# 
# 
# ## COMPARE STATE SPACE CONFIGURATIONS
# cSum <- colSums(mat.scale, na.rm = T)
# cSd <- apply(mat.scale,1,function(x)sd(x,na.rm=T),simplify = T)
# ggplot(compdf, aes(bsts.b3.diff.alog, colour=effect.type, fill=effect.type, line=effect.type)) + 
#   geom_density(alpha=0.2,size=.6) + facet_wrap( . ~ bsts.state.comps) +
#   geom_vline(xintercept=0, lty=2, size=.2) +
#   theme_bw() + ggtitle('Robustness of Pointwise ATT Error to Dynamic Treatment Effect Shape')
# ggsave(filename='statespace_gridsearch_compare_ATT_error_by_dynamic_treatment_effect_shape_bsts.png')
# ### ALL separate densities
# ggplot(compdf, aes(bsts.b3.diff.alog, colour=bsts.state.comps, fill=bsts.state.comps, line=bsts.state.comps)) + 
#   geom_density(alpha=0.2,size=.6) + facet_grid( effect.type ~ bsts.state.comps) +
#   geom_vline(xintercept=0, lty=2, size=.2) +
#   theme_bw() + ggtitle('Robustness of Pointwise ATT Error to Dynamic Treatment Effect Shape')
# ggsave(filename='statespace_gridsearch_compare_ATT_error_by_dynamic_treatment_effect_shape_ALL_DENSITIES_FACETGRID.png',
#        width=20,height = 8,units = 'in',dpi = 500)
# 
# 
# compdf$has.seasonality <- sapply(compdf$bsts.state.comps,function(x)ifelse(grepl('trig',x,ignore.case = T,perl = T),'Season','No Season'), simplify = T)
# compdf$has.level <- sapply(compdf$bsts.state.comps,function(x)ifelse(grepl('locallevel',x,ignore.case = T,perl = T),'Level','No Level'), simplify = T)
# ggplot(compdf, aes(bsts.b3.diff.alog, colour=effect.type, fill=effect.type, line=effect.type)) + 
#   geom_density(alpha=0.2,size=.6) + facet_grid(has.level  ~ has.seasonality) +
#   geom_vline(xintercept=0, lty=2, size=.2) +
#   theme_bw() + ggtitle('Robustness of Pointwise ATT Error to Dynamic Treatment Effect Shape')
# ggsave(filename='statespace_gridsearch_compare_ATT_error_by_dynamic_treatment_effect_shape_2x2_season_locallevel.png')
# 
# 
# ##===========================================
# ##
# ##  COMPARE  BSTS vs DID in HEATMAP AND DENSITIY FACETS
# ##
# ##-------------------------------------------
# 
# ## LONG DF
# compdf$bsts.did.b3.absdiff <- abs(compdf$bsts.b3.diff) -  abs(compdf$did.b3.diff)
# compdf$bsts.did.b3.absdiff.alog <- sapply(compdf$bsts.did.b3.absdiff, function(x) {
#   log(abs(x) + .001) * ifelse(x>0, 1, -1)
# } )
# 
# ## WIDE DF
# # View(compdf.wide)
# # stat.type <- 'bsts'
# ##
# compdf.wide.stat.bsts <- compdf.wide[compdf.wide$stat.type=='bsts', ] 
# compdf.wide.stat.did <- compdf.wide[compdf.wide$stat.type=='did', ] 
# compdf.wide.bsts.did.absdiff <- compdf.wide.stat.bsts
# for (i in 1:nrow(compdf.wide.bsts.did.absdiff)) {
#   compdf.wide.bsts.did.absdiff[i,-c(1:2)] <- abs(compdf.wide.stat.bsts[i, -c(1:2)]) - abs(compdf.wide.stat.did[i, -c(1:2)])
# }
# compdf.wide.bsts.did.absdiff.alog <- compdf.wide.bsts.did.absdiff
# for (i in 1:nrow(compdf.wide.bsts.did.absdiff)) {
#   for (j in 3:ncol(compdf.wide.bsts.did.absdiff)) {
#     x <- compdf.wide.bsts.did.absdiff[i,j]
#     compdf.wide.bsts.did.absdiff.alog[i,j] <-  log(abs(x) + .001) * ifelse(x>0, 1, -1) 
#   }
# }
# 
# mat <- as.matrix( compdf.wide.bsts.did.absdiff.alog[,-c(1:2)] )
# mat.scale <- scale(mat)
# rownames(mat.scale) <- compdf.wide.bsts.did.absdiff.alog$heat.y
# xdist <- dist(mat.scale)
# xcl <- hclust(xdist)
# ##
# xdend <-as.dendrogram(xcl)
# xdendord <- order.dendrogram(xdend)
# compdf$heat.y.ord <- factor(x=compdf$heat.y,
#                             levels=compdf.wide.stat$heat.y[ xdendord ], #xdendord
#                             ordered = TRUE)
# ##
# txdist <- dist(t(mat.scale))
# txcl <- hclust(txdist)
# txdend <-as.dendrogram(txcl)
# txdendord <- order.dendrogram(txdend)
# compdf$heat.x.ord <- factor(x=compdf$heat.x,
#                             levels=levels(factor(compdf$heat.x))[ txdendord ], #xdendord
#                             ordered = TRUE)
# # ##
# # compdf$bsts.b3.diff.alog <- sapply(compdf$bsts.b3.diff, function(x) {
# #   log(abs(x) + .01) * ifelse(x>0, 1, -1) 
# # } )
# ##
# # hmout <- heatmap(mat.scale, keep.dendro = T)
# # hmout
# 
# ##
# p.bsts.heatmap <- ggplot(compdf, aes(y=heat.y.ord, x=heat.x.ord, fill= bsts.did.b3.absdiff.alog )) + 
#   geom_tile() + 
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
#   ylab('Data Generating Process (DGP) Scenarios')  + xlab('BSTS State Space Configurations') + 
#   scale_fill_gradient2() + theme_bw()
# print(p.bsts.heatmap)
# ggsave(filename='heatmap_statespace_gridsearch_hclust_ABS-DIFF__ggplot.png', width = 10, height = 10, units = 'in', dpi = 300)
# 
# library(RColorBrewer)
# 
# # par(mar=c(3,3,3,6))
# png(filename = 'heatmap_hclust_statespace_gridsarch_ABS-DIFF__heatmap-func.png', width = 12, height = 10, units = 'in', res=300)
# hdend <- heatmap(mat.scale,col= RColorBrewer::brewer.pal(11,'RdBu'), keep.dendro = T, margins=c(5,12))
# dev.off()
# 
# 
# ## COMPARE STATE SPACE CONFIGURATIONS
# cSum <- colSums(mat.scale, na.rm = T)
# cSd <- apply(mat.scale,1,function(x)sd(x,na.rm=T),simplify = T)
# ggplot(compdf, aes(bsts.b3.diff.alog, colour=effect.type, fill=effect.type, line=effect.type)) + 
#   geom_density(alpha=0.2,size=.6) + facet_wrap( . ~ bsts.state.comps) +
#   geom_vline(xintercept=0, lty=2, size=.2) +
#   theme_bw() + ggtitle('Robustness of Pointwise ATT Error to Dynamic Treatment Effect Shape')
# ggsave(filename='statespace_gridsearch_compare_ABS-DIFF_ATT_error_by_dynamic_treatment_effect_shape_bsts.png')
# ### ALL separate densities
# ggplot(compdf, aes(bsts.b3.diff.alog, colour=bsts.state.comps, fill=bsts.state.comps, line=bsts.state.comps)) + 
#   geom_density(alpha=0.2,size=.6) + facet_grid( effect.type ~ bsts.state.comps) +
#   geom_vline(xintercept=0, lty=2, size=.2) +
#   theme_bw() + ggtitle('Robustness of Pointwise ATT Error to Dynamic Treatment Effect Shape')
# ggsave(filename='statespace_gridsearch_compare_ABS-DIFF_ATT_error_by_dynamic_treatment_effect_shape_ALL_DENSITIES_FACETGRID.png',
#        width=20,height = 8,units = 'in',dpi = 500)
# 
# 
# compdf$has.seasonality <- sapply(compdf$bsts.state.comps,function(x)ifelse(grepl('trig',x,ignore.case = T,perl = T),'Season','No Season'), simplify = T)
# compdf$has.level <- sapply(compdf$bsts.state.comps,function(x)ifelse(grepl('locallevel',x,ignore.case = T,perl = T),'Level','No Level'), simplify = T)
# ggplot(compdf, aes(bsts.b3.diff.alog, colour=effect.type, fill=effect.type, line=effect.type)) + 
#   geom_density(alpha=0.2,size=.6) + facet_grid(has.level  ~ has.seasonality) +
#   geom_vline(xintercept=0, lty=2, size=.2) +
#   theme_bw() + ggtitle('Robustness of Pointwise ATT Error to Dynamic Treatment Effect Shape')
# ggsave(filename='statespace_gridsearch_compare_ABS-DIFF_ATT_error_by_dynamic_treatment_effect_shape_2x2_season_locallevel.png')
# 
# 
# 
# 
# st.sp.perf <- compdf %>% 
#   dplyr::group_by( bsts.state.comps,effect.type) %>%   ## effect.type
#   dplyr::summarize(perf.sd=sd(bsts.did.b3.absdiff.alog,na.rm=T),
#                    perf.mean=mean(bsts.did.b3.absdiff.alog,na.rm=T)) %>% 
#   dplyr::arrange(perf.mean)
# print(st.sp.perf, n=50)
# plot(st.sp.perf$perf.sd, st.sp.perf$perf.mean)
# # for (i in 3:ncol(compdf.wide.stat)) {
# #   x <- compdf.wide.stat[,i] 
# #   x[x<0] <- -1 * log( abs(x[x<0]) + .001)
# #   x[x>0] <- log( x[x>0] + .001 )
# #   compdf.wide.stat[,i] <- x
# # }

## HEATMAP() FUNCTION
# hdend <- heatmap(as.matrix(compdf.wide.stat[,-c(1:2)]), col=RColorBrewer::brewer.pal(11,'RdBu'), 
#                  keep.dendro = T#, 
#                  # RowColColors = rainbow(nrow(compdf.wide.stat))
#                  )

# compdf.wide.stat <- compdf.wide[compdf.wide$stat.type==stat.type, ] 
# mat <- as.matrix( compdf.wide.stat[,-c(1:2)] )
# for (i in 1:ncol(mat)) {
#   x <- mat[,i]
#   x[is.na(x)] <- 99
#   x[x<0] <- -1 * log( abs(x[x<0]) + .001)
#   x[x>0] <- log( x[x>0] + .001 )
#   mat[,i] <- x
# }
# mat.scale <- scale(mat)
# # compdf$bsts.b3.diff.alog <- sapply(compdf$bsts.b3.diff, function(x) {
# #   log(abs(x) + .001) * ifelse(x>0, 1, -1) 
# # } )
# hdend <- heatmap(mat.scale,col= RColorBrewer::brewer.pal(11,'RdBu'), keep.dendro = T)
# 
# compdf$heat.y.ord2 <- factor(x=compdf$heat.y,
#                               levels=compdf.wide.stat$heat.y[ hdend$rowInd ],
#                               ordered = TRUE)
# # compdf$heat.x.ord2 <- factor(x=compdf$heat.x,
# #                              levels=levels(factor(compdf$heat.x))[ hdend$colInd ],
# #                              ordered = TRUE)
# ##
# # compdf$bsts.b3.diff.alog <- sapply(compdf$bsts.b3.diff, function(x) {
# #   log(abs(x) + .001) * ifelse(x>0, 1, -1) 
# # } )
# p.bsts.heatmap2 <- ggplot(compdf, aes(y=heat.y.ord2, x=heat.x.ord2, fill= bsts.b3.diff.alog )) + 
#   geom_tile() + 
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
#   ylab('Data Generating Process (DGP) Scenarios')  + xlab('BSTS State Space Configurations') + 
#   scale_fill_gradient2() + theme_bw()
# print(p.bsts.heatmap2)





# mat2 <- as.matrix( compdf.wide.stat[,-c(1:2)] )
# mat2 <- 1/mat2
# mat2[mat2==Inf] <- 1/1e-3
# mat.scale2 <- scale(mat2)
# hdend2 <- heatmap(mat.scale2, col= RColorBrewer::brewer.pal(11,'RdBu'), keep.dendro = T)
# 
# compdf$heat.y.ord2 <- factor(x=compdf$heat.y,
#                              levels=compdf.wide.stat$heat.y[ hdend2$rowInd ],
#                              ordered = TRUE)
# compdf$heat.x.ord2 <- factor(x=compdf$heat.x,
#                              levels=levels(factor(compdf$heat.x))[ hdend2$colInd ],
#                              ordered = TRUE)
# ##
# compdf$bsts.b3.diff.alog <- sapply(compdf$bsts.b3.diff, function(x) {
#   log(abs(x) + .001) * ifelse(x>0, 1, -1) 
# } )
# p.bsts.heatmap2 <- ggplot(compdf, aes(y=heat.y.ord2, x=heat.x.ord2, fill= bsts.b3.diff.alog )) + 
#   geom_tile() + 
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
#   ylab('Data Generating Process (DGP) Scenarios')  + xlab('BSTS State Space Configurations') + 
#   scale_fill_gradient2() + theme_bw()
# print(p.bsts.heatmap2)





  # ggplot(compdf.stack, aes(x=heat.x, y=hclust.ord, fill= b3.diff.alog )) +
  # geom_tile() + facet_wrap( . ~ factor(stats.type)) +
  # xlab(paste(heat.x.cols,collapse = ' | ')) + 
  # ylab(paste(heat.y.cols,collapse=' | ')) + 
  # scale_x_discrete(labels = function(x) str_wrap(x, width = 6)) +
  # scale_fill_gradient2() + theme_bw()
# compdf.stack$hclust.ord <- factor(x=compdf.stack$heat.x,
#                                   levels=compdf.stack$heat.x[ xdendord ],
#                                   ordered = TRUE)


# xdist <- dist(compdf.stack$b3.diff.alog)
# xcl <- hclust(xdist)
# xdend <-as.dendrogram(xcl)
# xdendord <- order.dendrogram(xdend)
# compdf.stack$hclust.ord <- factor(x=compdf.stack$heat.x,
#                                   levels=compdf.stack$heat.x[ xdendord ],
#                                   ordered = TRUE)
# 
# compdf.stack.bsts <- compdf.stack[compdf.stack$stats.type=='bsts',]
# compdf.stack.did  <- compdf.stack[compdf.stack$stats.type=='did',]

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





# 
# heatmap(
#   as.matrix(dat), Rowv=NA,
#   Colv=as.dendrogram(hclust(dist(t(as.matrix(dat)))))
# )

# 
# ##=====================================
# ##  Structural Changes
# ##-------------------------------------
# # ## devtools::install_github("KevinKotze/tsm")
# # ## devtools::install_github("cran/fArma")
# # library(tsm)
# # library(fArma)
# library(strucchange)
# i <- 3
# simi <- simlist[[ i ]]
# simidfs <- simi$sim$df.summary
# par(mfrow=c(3,2))
# for (.effect.type in effect.types) {
#   for (.group in c('treatment','control')) {
#     y <- simidfs$med[simidfs$effect.type==.effect.type & simidfs$group==.group]
#     dat <- data.frame(cbind(y[-1], y[-(length(y))]))
#     colnames(dat) <- c("ylag0", "ylag1")
#     fs <- Fstats(ylag0 ~ ylag1, data = dat)
#     print( breakpoints(fs) ) # where the breakpoint is
#     print( sctest(fs, type = "supF") )  # the significance of this breakpoint
#     plot(fs, main=sprintf('%s %s',.group,.effect.type), ylim=c(0,100)); abline(v=(simi$intpd - 1)/simi$npds)
#   }
# }
# 
# 
# 
# 
# 
# 
# # ## Results comparison table
# # res.tbl <- cbind(bsts.res[ ,c('point.effect','point.effect.lower','point.effect.upper')],
# #                  did.res[ ,c('term','event.time','estimate','point.conf.low','point.conf.high')],
# #                  b3.treat=b3diff$treat,
# #                  b3.ctrl=b3diff$ctrl,
# #                  b3.att=b3diff$diff 
# # )
# # ## ROUND RESULTS TABLE (round numeric columns)
# # # res.tbl <- res.tbl[ ! is.na(res.tbl$estimate), ]
# # num.cols <- names(res.tbl)[ ! names(res.tbl) %in% c('term','event.time') ]
# # res.tbl[ , num.cols] <- round( as.numeric(res.tbl[ , num.cols]), 4)
# # # for (i in 1:length(num.cols)) {
# # #   res.tbl[ , num.cols[i] ] <- round( as.numeric( res.tbl[ , num.cols[i] ] ), 4)
# # # }
# # ## MOVE ID COLUMNS TO FRONT
# # .col.idx <- which(names(res.tbl) %in% c('term','event.time'))
# # res.tbl4 <- cbind(res.tbl[, .col.idx],  res.tbl[, -.col.idx] )
# # # View(res.tbl4)
# 
# ##PLOT INCLUSION PROBABILITIES
# # png(filename = sprintf('%s_BSTS_inclusion_probs_ss%s_%s_%s_%s.png',
# #                        prefix,h,key.strip,effect.type,sim.id))
# # plot(impact_amount$model$bsts.model,'coefficients', main=sprintf('%s %s', key,effect.type))
# # dev.off()
# 
# ## PLOT DYNAMIC EFFECTS COMPARISON - DID vs. BSTS vs. DGP
# png(filename = sprintf('%s_BSTS_dynamic_treatment_effect_comparison_ss%s_%s_%s_%s.png',
#                        prefix,h,key.strip,effect.type,sim.id))
# res.tbl.filename <- sprintf('%s: DGP = %.3f; DiD = %.3f;  BSTS = %.3f',
#                             key.strip, mean(res.tbl$b3.att[intpd:nrow(b3.att)]), agg.es$overall.att, impact_amount$summary$AbsEffect[1])
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
