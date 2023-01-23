##=======================================================================
##
##
##
##
##
##
##        BSTS vs. DiD 
##        COMPARISON and SENSITIVITY ANALYSIS FUNCTIONS 
##
##
##
##
##
##=======================================================================
library(e1071)


###
## Aggregate Timeseries Simulation Panel Dataframe
###
getAggregatedSimPanelDf <- function(tmpdf, pd.agg, 
                                    na.rm=TRUE) {
  if (is.na(pd.agg) | pd.agg <= 1) {
    cat(sprintf('\npd.agg=%s is either NA or <= 1. Try setting pd.agg >= 2.\n',pd.agg))
  }
  ##
  ts <- sort(unique(tmpdf$t))
  npds.orig <- length(ts)
  npds.new <- round( npds.orig / pd.agg )
  actors <- sort(unique(tmpdf$actor[!is.na(tmpdf$match_id)]))
  
  aggmap <- data.frame(t.old=1:npds.orig,
                       t.new=rep(1:npds.new, each=pd.agg))
  
  ##
  intpd.old <- unique(tmpdf$t[which(tmpdf$t.post.intpd==1)])[1]
  intpd.new <- aggmap$t.new[which(aggmap$t.old == intpd.old)]
  
  ##
  tmpdf$match_pd <- as.numeric(tmpdf$match_pd)
  tmpdf$match_pd[tmpdf$match_pd == intpd.old] <- intpd.new
  
  
  tmpdf$t.agg <- rep(NA, nrow(tmpdf))
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
      # ##
      # c1_sd = sd(c1, na.rm=T),
      # c2_sd = sd(c2, na.rm=T),
      # c3_sd = sd(c3, na.rm=T),
      # ##
      # c1_skew = skewness(c1, na.rm=T, type = 2),
      # c2_skew = skewness(c2, na.rm=T, type = 2),
      # c3_skew = skewness(c3, na.rm=T, type = 2),
      # #
      # c1_kurt = skewness(c1, na.rm=T, type = 2),
      # c2_kurt = skewness(c2, na.rm=T, type = 2),
      # c3_kurt = skewness(c3, na.rm=T, type = 2),
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
  ## remove unnecessary columns
  aggdf$t.agg <- NULL
  aggdf$t.agg.post.intpd <- NULL
  
  # ##
  # aggdf$gname <- 0
  # aggdf$gname[aggdf$group=='treatment'] <- aggdf$match_pd[aggdf$group=='treatment']
  # aggdf$gname <- as.numeric(aggdf$gname)
  
  ##
  return(aggdf)
}

##=======================================================================###
## Update Simlist configurations and simulated panel dataframes for aggregated periods
###
updateSimlistAggregateSimDfPd <- function(simlist, pd.agg, na.rm=TRUE) {
  
  for (i in 1:length(simlist)) 
  {
    simx <- simlist[[ i ]]
    
    if (is.null(simx$sim)) {
      cat(sprintf('\nsimlist item i=%s is missing "sim" object. First call runSimUpdateSimlist().\n',i))
    }
    
    aggdf <- getAggregatedSimPanelDf(simx$sim$df, pd.agg=pd.agg, na.rm=na.rm)
    simlist[[i]]$sim$df <- aggdf
    simlist[[i]]$npds <- length(unique(aggdf$t))
    simlist[[i]]$intpd <- unique(aggdf$t[which(aggdf$t.post.intpd==1)])[1]
  }
  
  return(simlist)
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
  object <- aggte(attgt, type = "dynamic", bstrap = TRUE, na.rm = F)
  
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
  an.x <- mean(as.numeric(results$year), na.rm=T) - 0.1*diff(rng.x)
  
  p <- ggplot(results,
              aes(x=as.numeric(year), y=att, ymin=(att-c*att.se),
                  ymax=(att+c*att.se))) +
    geom_errorbar(colour='darkgray', width=0.3) +
    geom_point(colour='black', size=1.8) +
    #geom_ribbon(aes(x=as.numeric(year)), alpha=0.2) +
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
##  MAIN COMPARISON FUNCTION
##  - 1. runSimSingleInterventionEffectComparison() on simlist
##  - 2. DiD & BSTS results
##  - 3. comparison of DiD & BSTS performance
##   returns full simlist
######################################
runSimCompareBstsDiD <- function(simlist,     ## n, npds, intpd moved into simlist elements
                                 effect.types=c('constant','quadratic','geometric'), 
                                 sim.id=NA,
                                 save.items.dir=NA, ## save updated simlist items to separate RDS files
                                 bsts.niter=5000, 
                                 bsts.max.iter=8e4, ## 80000 
                                 bsts.ctrl.cats=2,  ## bins for each covariate (c1,c2,c3)
                                 bsts.expect.mod.size=1,
                                 plot.show=TRUE, plot.save=TRUE,
                                 verbose = TRUE
) {
  ## cache original bsts.niter for dynamic niter updating if MCMC convergence failed
  bsts.niter.orig <- bsts.niter
  
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
    if(verbose) cat(sprintf('\n%s, %s\n',i, key))
    
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
    noise.level <- simlist[[key]]$noise.level
    covariates.type <- simlist[[key]]$covariates.type
    rand.seed <- simlist[[key]]$rand.seed
    
    ## Simulation object (containing the simulated timeseries)
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
      ## Remove NAs
      simdf <- simdf[!is.na(simdf$match_id), ]
      
      ## Rows
      nrow.sdf <- nrow(simdf)
      
      ## Compute Multiperiod DiD for Avg Treatment Effect on Treated
      set.seed( rand.seed )
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
      .filename <- sprintf('%s_did_dynamic_effect_n%s_pd%s_cov-%s_%s_%s_%s.png',
              prefix,n,npds,covariates.type,key.strip,effect.type,sim.id)
      ggsave(filename = file.path(save.img.dir,  .filename),  
             width = 9, height = 6, units = 'in',dpi = 300)
      
      
      ## Get first treatment group actor
      tr.actor.1 <- simdf$actor[which(simdf$group=='treatment')[1]]
      
      ## SIMPLE AGGREGATION (OVERALL EFFECT) ATT
      agg.simple <- aggte(ccattgt, type='simple', bstrap = TRUE, na.rm = F)
      # summary(agg.simple)
      
      ## GROUP EFFECT AGGREGATION (for reporting overall confidence intervals)
      agg.group <- aggte(ccattgt, type = "group", bstrap = TRUE, na.rm = F)
      
      ## DYNAMIC EFFECTS AND EVENT STUDIES
      agg.es <- aggte(ccattgt, type = "dynamic", bstrap = TRUE, na.rm = F)
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
      
      
      ##===========================================
      ##
      ##
      ## BSTS Timseries Setup
      ##
      ##
      ##----------------------------------------

      ## AGGREGATE LONG DATAFRAME OF COVARIATES 
      ## TREATMENT
      bsts.df <- simdf %>%
        dplyr::filter( 
          ! is.na(match_id), 
          group=='treatment'
        ) %>%
        group_by(t) %>%
        summarize(
          y_treatment = mean(y, na.rm=T)
        ) 
      # %>%
      #   mutate(t=t,y_outcome=mean, .keep='used')
      
      if ( is.na(bsts.ctrl.cats) ) { 
        
        ## ONLY 1 WHOLE-CONTROL-GROUP AGGREGATED AS ONE SERIES
        cov.df.wide <- simdf %>%
          dplyr::filter( 
            ! is.na(match_id), 
            group=='control'
          ) %>%
          group_by(t) %>%
          dplyr::summarize(
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
        
      } else if (bsts.ctrl.cats < 2) {
        
        ### 1 CONTROL GROUP plus covariate series
        cov.df.wide <- simdf %>%
          dplyr::filter( 
            ! is.na(match_id), 
            group=='control'
          ) %>% 
          group_by(t) %>%
          dplyr::summarize(
            y_control = mean(y, na.rm=T),
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
        
      } else {
        
        ### SYNTHETIC CONTROL GROUPS (intersection of all covariate factors)
        cat('\ncreating synthetic control groups as covariate series...')
        ## create quantile categories (bins) 
        # ## quantile upper limits; Get equal probability sequence from number of categories to create bins
        # quant.up.lims <- (1:ncats)*(1/ncats)
        c1cats <- cut(simdf$c1, bsts.ctrl.cats)
        simdf$c1.f <- LETTERS[as.integer(as.factor(c1cats))]
        # c1lvls <- sort(unique(c1cats))
        # simdf$c1.f <- sapply(1:nrow.sdf, function(i)sprintf('%s%s',LETTERS[which(c1lvls==c1cats[i])],c1cats[i]))
        ##
        c2cats <- cut(simdf$c2, bsts.ctrl.cats)
        simdf$c2.f <- LETTERS[as.integer(as.factor(c2cats))]
        # c2lvls <- sort(unique(c2cats))
        # simdf$c2.f <- sapply(1:nrow.sdf, function(i)sprintf('%s%s',LETTERS[which(c2lvls==c2cats[i])],c2cats[i]))
        ##
        c3cats <- cut(simdf$c3, bsts.ctrl.cats)
        simdf$c3.f <- LETTERS[as.integer(as.factor(c3cats))]
        # c3lvls <- sort(unique(c3cats))
        # simdf$c3.f <- sapply(1:nrow.sdf, function(i)sprintf('%s%s',LETTERS[which(c3lvls==c3cats[i])],c3cats[i]))
        ### SYNTHETIC CONTROL GROUPS
        .cov.df <- simdf %>%
          dplyr::filter( 
            ! is.na(match_id), 
            group=='control'
          ) %>%
          group_by(t, c3.f, c2.f, c1.f) %>%
          dplyr::summarize(
            cov_mean = mean(y, na.rm=T)#,
            # cov_sum = sum(y, na.rm=T),
            # cov_sd = sd(y, na.rm=T)
          )  
        .cov.df$cov_cat.f <- apply(.cov.df[,c('c1.f','c2.f','c3.f')], 1, function(x)paste(x,collapse = ''))
        .cov.df$c1.f <- NULL
        .cov.df$c2.f <- NULL
        .cov.df$c3.f <- NULL
        ###
        cov.df.wide <- .cov.df %>% 
          pivot_wider(names_from = cov_cat.f, values_from=c(cov_mean))
        max.cov.missing <- 0.7
        cov.cols.keep <- apply(cov.df.wide[,-1], 2, function(x){ 
           ( ( sum(!is.na(x)) / length(x)  ) > max.cov.missing ) & ## have enough non-missing values
           ( !is.na(x[1]) | !is.na(x[length(x)]) )    ## has either first or last row non-missing
        })
        ## KEEP IF First or last row is not NA (for fill() below) 
        cov.cols.keep.names <- names(cov.df.wide[,-1])[cov.cols.keep]  ##exclude the 1st column 't'
        cov.df.wide <- cov.df.wide %>% dplyr::select(all_of(c('t', cov.cols.keep.names)))
        ##
        cov.cols.need.fill.bool <- apply(cov.df.wide[,-1], 2, function(x){ 
          ( sum(!is.na(x)) / length(x)  ) < 1
        })
        # cov.cols.nonNApct.keep <- (cov.cols.nonNApct2)
        cov.cols.fill <- names(cov.df.wide)[ cov.cols.need.fill.bool ] 
        cov.df.wide <- cov.df.wide %>% 
          ungroup() %>% 
          tidyr::fill(all_of(cov.cols.fill), .direction = 'downup') #%>%
          # tidyr::fill(all_of(cov.cols.fill), .direction = 'up')
        # for (key in cov.cols.keep.names) {
        # cov.df.wide <- cov.df.wide %>% dplyr::fill('B|A|B', .direction='down')
          # cov.df.wide <- rbind(
          #   .cov.df %>% fill(.data[[key]], .direction = 'downup')
          # )
        # }
        ### FILL REMAINING NAs (after 'downup' fill) with zeros
        cov.df.wide[is.na(cov.df.wide)] <- 0
        ##
        cat('done.\n')
        
      }
 
      
      ## JOIN TREATMENT AND CONTROL (WIDE) DATAFRAME into BSTS INPUT DATA
      bsts.df <- bsts.df %>% full_join(cov.df.wide, by='t')
      bsts.df$t <- NULL
      ##############################
      # ## Aggregate into timeseries dataframe
      # tsdf <- simdf %>%
      #   dplyr::filter( ! is.na(match_id)) %>%
      #   group_by(t, t.post.intpd, effect.type, match_pd, gname, group, group.color) %>%
      #   dplyr::summarize(
      #     n_in_pd = n(),
      #     actors = paste(unique(actor), collapse = '|'),
      #     y_outcome = mean(y, na.rm=T),
      #     y_sum = sum(y, na.rm=T),
      #     y_sd = sd(y, na.rm=T),
      #     y_min = min(y, na.rm=T),
      #     y_max = max(y, na.rm=T),
      #     y_skew = skewness(y, na.rm = T, type = 2), ## moment-based distribution
      #     y_kurt = kurtosis(y, na.rm = T, type = 2), ## moment-based distribution
      #     # ##
      #     # x1_sum = sum(x1, na.rm=T),
      #     # x2_sum = sum(x2, na.rm=T),
      #     # x3_sum = sum(x3, na.rm=T),
      #     # ##
      #     # c1_sum = sum(c1, na.rm=T),
      #     # c2_sum = sum(c2, na.rm=T),
      #     # c3_sum = sum(c3, na.rm=T),
      #     # #
      #     # b1_sum = sum(b1, na.rm=T),
      #     # b2_sum = sum(b2, na.rm=T),
      #     # b3_sum = sum(b3, na.rm=T),
      #     # #
      #     # u_sum = sum(u, na.rm=T),
      #     # v_sum = sum(v, na.rm=T),
      #     ##
      #     x1_mean = mean(x1, na.rm=T),
      #     x2_mean = mean(x2, na.rm=T),
      #     x3_mean = mean(x3, na.rm=T),
      #     ## 
      #     c1_mean = mean(c1, na.rm=T),
      #     c2_mean = mean(c2, na.rm=T),
      #     c3_mean = mean(c3, na.rm=T),
      #     #
      #     c1_sd = sd(c1, na.rm=T),
      #     c2_sd = sd(c2, na.rm=T),
      #     c3_sd = sd(c3, na.rm=T),
      #     ##
      #     c1_skew = skewness(c1, na.rm=T, type = 2),
      #     c2_skew = skewness(c2, na.rm=T, type = 2),
      #     c3_skew = skewness(c3, na.rm=T, type = 2),
      #     #
      #     b1_mean = mean(b1, na.rm=T),
      #     b2_mean = mean(b2, na.rm=T),
      #     b3_mean = mean(b3, na.rm=T),
      #     #
      #     u_mean = mean(u, na.rm=T),
      #     v_mean = mean(v, na.rm=T)
      #   )
      # tsdf$.id <- 1:nrow(tsdf)
      
      # ## MAKE WIDE TIMESERIES FOR treatment,control groups in n periods
      # val.cols <- c('y_outcome','y_sum','y_min','y_max','y_sd',
      #               'y_skew','y_kurt',
      #               'x1_mean','x2_mean','x3_mean',
      #               'c1_mean','c2_mean','c3_mean',
      #               'c1_sd','c2_sd','c3_sd',
      #               'c1_skew','c2_skew','c3_skew',
      #               'c1_kurt','c2_kurt','c3_kurt',
      #               'b1_mean','b2_mean','b3_mean',
      #               'u_mean','v_mean')
      # ts <- unique(tsdf$t)
      # groups <- unique(tsdf$group)
      # ## init timeseries dataframe - wide
      # tsdfw <- data.frame(t=ts,  stringsAsFactors = F)
      # for (jj in 1:length(groups)) {
      #   id.j <- which( tsdf$group == groups[ jj ] ) 
      #   for (kk in 1:length(val.cols)) {
      #     df.col <- data.frame( tsdf[ id.j , val.cols[ kk ] ] )
      #     names(df.col) <- sprintf('%s_%s',groups[ jj ],val.cols[ kk ])
      #     tsdfw <- cbind(tsdfw,  df.col)
      #   }
      # }
      
      
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
      # dat <- tsdfw[,c('y_treatment',
      #                 'y_control',
      #                 'c1_mean','c2_mean','c3_mean',
      #                 'c1_sd','c2_sd','c3_sd'
      # )]
      dat <- bsts.df
      ## Train on y pre-treatment but NA's post-treatment
      y.pre.treat.NAs.post.treat <- c(dat$y_treatment[1:(intpd-1)], rep(NA,npds-intpd+1))
      ## Then use the post-treatment response for causal impact estimation
      post.period.response <- dat$y_treatment[intpd:npds]
      ## Covariates (predictors) - Matrix for "formula = y ~ predictors" argument
      predictors <- dat[, ! names(dat) %in% 'y_treatment'] ## remove response; convert to matrix
      # ## Covariates (predictors) - Dataframe for "data" argument
      # predictors <- as.matrix(predictors) 
      
      # ## ADD temporal trend to covariates
      # predictors$covar_temp_trend <- (1:npds) + rnorm(npds, 0, noise.level)   ## * simlist$`0.rand.base`$b5
      
      
      ##----------------------------
      ## State Space Configuration
      ##----------------------------
      ## LOOP OVER BSTS STATE SPECIFICATIONS FOR THE SAME SIMULATION RUN
      if (length(bsts.state.specs)<1){
        stop('\nbsts.state.specs list length = 0. Add BSTS state space specifications.\n')
      }
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
              if(verbose) cat(sprintf('add to state.space: %s\n',state.conf.item$name))
            }
          }
        } else {
          ## Default in CausalImpact package
          st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
        }
        # print(st.sp)
        
        if(verbose) cat(sprintf('\nRunning BSTS model estimation for state.conf h=%s\n',h))
        
        ##-------------------------------
        ## Regression component of model
        if ('AddRegression' %in% state.comps) {
          bsts.input.form <- y.pre.treat.NAs.post.treat ~ . ## with regression formula
        } else {
          bsts.input.form <- y.pre.treat.NAs.post.treat  ## without regression vector
        }
        
        
        ##---------------------------------------------------------
        ## RUN BSTS WITH DYNAMIC niter BASED ON CONVERGENCE 
        isConverged <- FALSE
        isMaxIter <- FALSE
        hasBstsError <- FALSE
        bsts.niter <- bsts.niter.orig  ## reset to original bsts.niter input value
        while ( !isConverged  &  !isMaxIter & !hasBstsError  ) {
          set.seed( rand.seed )
          ## BSTS model
          bsts.model <- tryCatch(expr = {
            bsts(formula = bsts.input.form,
                 state.specification = st.sp,
                 data = predictors,
                 expected.model.size = bsts.expect.mod.size,
                 niter = bsts.niter,
                 ping = ifelse(verbose, round(bsts.niter/10), 0))
          },
          error=function(e) {
            message(sprintf('bsts() error: %s', as.character(e)))
          },
          warning=function(w) {
            message(sprintf('bsts() warning: %s', as.character(w)))
          },
          finally={
            ##PASS
          })
          
          
          ## skip if bsts() function threw an error (not return 'bsts' object)
          if ( class(bsts.model) == 'bsts' ) {
            ## Use BSTS prediction of counterfactual to estimate CausalImpact
            impact_amount <- CausalImpact(bsts.model=bsts.model,
                                          post.period.response = post.period.response,
                                          alpha=0.05, model.args = list(niter = bsts.niter))
            ## POSTERIOR PREDICTIVE CHECKS
            ppcheck.filename <- file.path(save.img.dir,
                                          sprintf('%s_bsts_post_pred_checks_n%s_pd%s_ss%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_%s.png',
                                                  prefix,n,npds,h,bsts.niter,bsts.ctrl.cats,bsts.expect.mod.size,covariates.type,key.strip,effect.type,sim.id))
            convcheck <- bstsPostPredChecks(bsts.model, filename=ppcheck.filename, return.val = T)
            ##
            # print(convcheck)
            ##
            if(verbose) cat(sprintf('\nBSTS niter = %s',bsts.niter))
            if(verbose) cat(convcheck$summary)
            ## UPDATE CONVERGENCE CHECK FLAG - ELSE RERUN WITH INCREASED bsts.niter
            # isConverged <- convcheck$converged.all
            conv.tol <- 0.8
            conv.min.iter.thresh <- 4e4 ## 40k
            # isConverged <- convcheck$converged.prop >= conv.tol
            isConverged <- convcheck$converged.all | (convcheck$converged.prop >= conv.tol & bsts.niter >= conv.min.iter.thresh) ## don't allow incomplete check below minimum threshold of bsts.niter = 10k 
            if(verbose) print(convcheck$converged)
            if(verbose) cat(sprintf('Converged proportion = %.3f (tol = %.3f) (min.iter.converg.thresh=%s)\nConverged status = %s\n\n',
                        convcheck$converged.prop, conv.tol, conv.min.iter.thresh,  isConverged))
          } else {
            hasBstsError <- TRUE
          }
          

          ## UPDATE niter --> increase by factor of 2
          if ( !isConverged ) {
            # .niter <- round( bsts.niter * (2 - convcheck$converged.prop) ) ## scale niter increment by proportion of failed checks
            bsts.niter.new <- round( (bsts.niter * 2) / 10) * 10  ## double niter and make divisible by 10 (for bsts ping=10 to divide evenly into niter)
            if ( bsts.niter.new <= bsts.max.iter ) {
              bsts.niter <- bsts.niter.new  
            } else {
              isMaxIter <- TRUE
            }
          }

       
        }
        
        if (hasBstsError) {
          next
        }
        ##-------------------------------------------------------------
        
        if (plot.show) {
          ##
          # plot(bsts.model, main=sprintf('BSTS Plot: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
          # PlotBstsComponents(bsts.model) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
          # PlotBstsState(bsts.model, main=sprintf('BSTS State: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
          # PlotBstsResiduals(bsts.model, main=sprintf('BSTS Residuals: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
          # PlotBstsPredictionErrors(bsts.model, main=sprintf('BSTS Pred.Err: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
          # PlotBstsForecastDistribution(bsts.model, main=sprintf('BSTS Forecast Dist: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
        }

        ##
        # ## BSTS model for Dynamic Regression
        # bsts.model <- bsts(y.pre.treat.NAs.post.treat,
        #                    state.specification = st.sp,
        #                    niter = 5000)
        
        # ##
        # summary(impact_amount)
        # summary(impact_amount$model$bsts.model)
        # plot(impact_amount)
        
        if (plot.save) 
        {
          ## ----------- BSTS COMPONENTS PLOT --------------------------------------
          .filename <- sprintf('%s_bsts_state_components_n%s_pd%s_ss%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_%s.png',
                               prefix,n,npds,h,bsts.niter,bsts.ctrl.cats,bsts.expect.mod.size, covariates.type,key.strip,effect.type,sim.id)
          png(file.path(save.img.dir, .filename), width = 10, height = 6, units = 'in', res = 300)
          PlotBstsComponents(bsts.model, burn = bsts.niter*.2) 
          dev.off()
          
          ## ----------- BSTS MODEL OUTPUT --------------------------------------
          .filename <- sprintf('%s_bsts_model_pred-obs-circles_n%s_pd%s_ss%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_%s.png',
                               prefix,n,npds,h,bsts.niter,bsts.ctrl.cats,bsts.expect.mod.size, covariates.type,key.strip,effect.type,sim.id)
          png(file.path(save.img.dir, .filename), width = 7, height = 6, units = 'in', res = 300)
          plot(bsts.model, main=sprintf('BSTS Model: %s',paste(state.comps,collapse = ' + ')))
          dev.off()
          
          ## ----------- BSTS CAUSAL IMPACT --------------------------------------
          # summary(impact_amount)
          # png(filename=sprintf('single_intervention_BSTS_CausalImpact_plot_%s_%s_%s.png',
          #                         key,effect.type,sim.id))
          .filename <- sprintf('%s_bsts_CausalImpact_plot_n%s_pd%s_ss%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_%s.png',
                              prefix,n,npds,h,bsts.niter,bsts.ctrl.cats,bsts.expect.mod.size, covariates.type,key.strip,effect.type,sim.id)
          p.bsts.impact.all <- plot(impact_amount, c('original','pointwise','cumulative')) # pointwise','cumulative
          ggsave(filename = file.path(save.img.dir, .filename), width = 8, height = 12, units = 'in',dpi = 300)
          # dev.off()
          
          
          ## ----------- BSTS REGRESSION  ----------------
          if (bsts.model$has.regression) {
            
            ## ----------- BSTS MODEL SIZE DISTRIBUTION  ----------------
            .filename <- sprintf('%s_bsts_size_dist_n%s_pd%s_ss%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_%s.png',
                                 prefix,n,npds,h,bsts.niter,bsts.ctrl.cats,bsts.expect.mod.size, covariates.type,key.strip,effect.type,sim.id)
            png(file.path(save.img.dir, .filename), width = 7, height = 6, units = 'in', res = 300)
            PlotBstsSize(bsts.model, main=sprintf('BSTS Model Size Distribution (expected = %s)',
                                                  bsts.expect.mod.size))
            dev.off()
            
            ## ----------- BSTS MODEL SIZE DISTRIBUTION  ----------------
            .filename <- sprintf('%s_BSTS_inclusion_probs_n%s_pd%s_ss%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_%s.png',
                              prefix,n,npds,h,bsts.niter,bsts.ctrl.cats,bsts.expect.mod.size, covariates.type,key.strip,effect.type,sim.id)
            png(filename = file.path(save.img.dir, .filename),  width = 7, height = 6, units = 'in', res = 300)
            plot(bsts.model,'coefficients', main=sprintf('BSTS Inclusion Probabilities (expected size = %s)',
                                                         bsts.expect.mod.size))
            dev.off()
            
          }
          
        }

        
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
          treat=simdf %>% 
            dplyr::filter(group=='treatment' & actor==tr.actors[1]) %>% 
            ungroup() %>%
            dplyr::mutate(treat=b3) %>% 
            dplyr::select(treat),  ##ungroup() %>% 
          ctrl=simdf %>% 
            dplyr::filter(group=='control' & actor==co.actors[1]) %>% 
            ungroup() %>%
            dplyr::mutate(ctrl=b3) %>% 
            dplyr::select(ctrl),
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
        
        if (plot.save) 
        {
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
          ##TODO
          ##***CHANGE TO 1-STEP AHEAD PREDICTION ERROR BEF
          errdf <- rbind(
            data.frame(method='BSTS', 
                       error=(res.tbl$bsts.point.effect[1:(intpd-1)] - res.tbl$b3.att[1:(intpd-1)])
                       ),  ## MAPE - percentage
            data.frame(method='DiD',  
                       error=(res.tbl$did.estimate[1:(intpd-1)] - res.tbl$b3.att[1:(intpd-1)])
                       )
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
            data.frame(method='BSTS', 
                       error=(res.tbl$bsts.point.effect[intpd:npds] - res.tbl$b3.att[intpd:npds])
                       ),
            data.frame(method='DiD', 
                       error=(res.tbl$did.estimate[intpd:npds] - res.tbl$b3.att[intpd:npds])
                       )
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
          
          ##----------------------------------
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
                                        sprintf('%s_ATT_pointwise_error_distribution_compare_n%s_pd%s_ss%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_%s.png',
                                                prefix,n,npds,h,bsts.niter,bsts.ctrl.cats,bsts.expect.mod.size, covariates.type,key.strip,effect.type,sim.id)),
                   height=15, width=9, units = 'in', dpi = 300)
        }

        
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
        simlist[[key]]$compare$bsts[[effect.type]][[ h ]]$convcheck <- convcheck
        ##
        simlist[[key]]$compare$att.b3 <- att.b3
        simlist[[key]]$compare$att.did <- att.did
        simlist[[key]]$compare$att.bsts <- att.bsts
        simlist[[key]]$compare$bias.bsts<- att.bsts - att.b3
        simlist[[key]]$compare$bias.did <- att.did - att.b3
        simlist[[key]]$compare$perf.adv.bsts <- abs(att.did - att.b3) / abs(att.bsts - att.b3)
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
      maxGb <- 6
      # maxGb <- .001  ## **DEBUG**
      if ( (object.size(simlist[[key]])/1e9) <= maxGb) {
        ## IF SMALL ENOUGH, Save simulation list as serialized data file
        simlist.file <- sprintf('__%s_GRIDSEARCH_output__n%s_pd%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s.rds', 
                                prefix, n, npds, bsts.niter, bsts.ctrl.cats, bsts.expect.mod.size, covariates.type, sim.id, key.strip)
        save.file.path <-  file.path(save.items.dir, simlist.file)
        saveRDS(simlist[[key]], file = save.file.path)
        ## FREE UP MEMORY
        simlist[[key]] <- list(file = save.file.path)
      } else { 
        ## IF TOO LARGE, then save BSTS object (with MCMC samples) separately from rest of simulation
        save.file.paths <- c()
        ## Save simulation list as serialized data file
        ## 1.  BSTS
        simlist.file <- sprintf('__%s_GRIDSEARCH_output__n%s_pd%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_simlist-compare-bsts.rds', 
                                prefix, n, npds, bsts.niter, bsts.ctrl.cats, bsts.expect.mod.size, covariates.type, sim.id, key.strip)
        save.file.path <-  file.path(save.items.dir, simlist.file)
        saveRDS(simlist[[key]]$compare$bsts, file = save.file.path)
        simlist[[key]]$compare$bsts <- NULL ## save space
        save.file.paths[1] <- save.file.path
        ## 2. rest of simulation object
        simlist.file <- sprintf('__%s_GRIDSEARCH_output__n%s_pd%s_niter%s_covCats%s_msize%s_cov-%s_%s_%s_simlist-MAIN.rds', 
                                prefix, n, npds, bsts.niter, bsts.ctrl.cats, bsts.expect.mod.size, covariates.type, sim.id, key.strip)
        save.file.path <-  file.path(save.items.dir, simlist.file)
        saveRDS(simlist[[key]], file = save.file.path)
        save.file.paths[2] <- save.file.path
        ## FREE UP MEMORY
        simlist[[key]] <- list(file = save.file.paths)
      }
    } 
    
    
  } ## // end simlist loop i   ##  #; dev.off()
  
  return(simlist)
  
}
















cat(sprintf('\nLoaded BSTS vs. DiD Comparison and Sensitivity Analysis Functions.\n'))








#####################################################################
##*****MOVED*** TO  single_intervention_sim_vec.R file
# ####################################
# ##  RUN INTERVENTION SIMULATION IN LOOP OVER EFFECT TYPES
# ######################################
# runSimUpdateSimlist <- function(simlist,     ## n, npds, intpd moved into simlist elements
#                                 effect.types=c('constant','quadratic','geometric'), 
#                                 sim.id=round(10*as.numeric(Sys.time())),
#                                 plot.show=F, plot.save=F, verbose=TRUE) {
#   
#   # print("runSimBstsDiDComparison()::SIMLIST INPUT:")
#   # print(simlist)
#   if (length(simlist) > 0 & length(names(simlist))==0) {
#     names(simlist) <- 1:length(simlist)
#   }
#   ##----------------------------
#   ## Run Simulation List
#   ##----------------------------
#   # sim.id <- ifelse( is.na(sim.id), round(as.numeric(Sys.time())), sim.id)
#   for (i in 1:length(simlist)) {
#     key <- names(simlist)[i]
#     sim <- simlist[[key]]
#     if(verbose) cat(sprintf('\nScenario label: %s\n\n', key))
#     ##  
#     set.seed( ifelse(is.null(sim$rand.seed), 54321, sim$rand.seed) )
#     noise.level <- ifelse(is.null(sim$noise.level), 0, sim$noise.level)
#     ##
#     cov.scenario <- if(is.null(sim$cov.scenario)){list(c1=.1, c2=.2, c3=.3)} else{ sim$cov.scenario }
#     ##
#     simlist[[key]]$sim <- runSimSingleInterventionEffectComparison(
#       effect.types = effect.types,
#       n = sim$n, ## NUMBER OF FIRMS
#       npds = sim$npds, ## NUMBER OF PERIODS
#       intpd = sim$intpd, ## intervention after first section
#       ystart = 0.1,
#       treat.rule = ifelse(is.null(sim$treat.rule), NA, sim$treat.rule),
#       treat.prob =  ifelse(is.null(sim$treat.prob), NA, sim$treat.prob), #0.95,  ## 0.6
#       treat.threshold = ifelse(is.null(sim$treat.threshold), NA, sim$treat.threshold),  # 0.015
#       sim.id = sim.id, ## defaults to timestamp
#       cov.scenario = cov.scenario,
#       ##
#       noise.level = noise.level,
#       ##
#       b4 = ifelse(is.null(sim$b4), 0, sim$b4), ## past performance
#       b5 = ifelse(is.null(sim$b5), 0, sim$b5), ## growth (linear function of time t)
#       b9 = ifelse(is.null(sim$b9), 0, sim$b9), ## Autocorrelation
#       ## Dynamic treatment effect polynomial parameters
#       w0 = ifelse(is.null(sim$w0), 2.0 , sim$w0), ## constant
#       w1 = ifelse(is.null(sim$w1), 0.18 , sim$w1), ## linear
#       w2 = ifelse(is.null(sim$w2), -.1 / (sim$npds^.6) , sim$w2), ## ,  ##-0.009,  ## quadratic  ## -0.005, ## ***** made steeper curve= -0.008 *****
#       ## TODO: CHECK w2.shift SENSITIVITY - this shifts quadratic curve several periods to the right so that treatment effect increases slowly  
#       w2.shift = ifelse(is.null(sim$w2.shift), -round( sqrt(sim$npds)*.8 ) , sim$w2.shift),
#       # w2.shift = -round( sqrt(sim$npds)*.8 ),  ## optimal value here is likely a function of the combination of treatment effect function parameters
#       ##
#       nseasons  = ifelse(is.null(sim$dgp.nseasons), NA, sim$dgp.nseasons), 
#       season.frequency = ifelse(is.null(sim$dgp.freq), NA, sim$dgp.freq),
#       ## Plotting
#       plot.show = plot.show, ## TRUE
#       plot.save = plot.save, ## TRUE
#       verbose = verbose ## echo messages to console
#     )
#   }
#   
#   return(simlist)
# }


# ## Save simulation list as serialized data file
# simlist.file <- sprintf('single_intervention_SIMLIST_selection_endog_%s_%s.rds',
#                         length(simlist), sim.id)
# saveRDS(simlist, file = file.path(dir_plot, simlist.file))






########### MOVED TO bsts_helper_functions ############################
# 
# ##
# #
# ##
# heidel.diag.mod <- function (x, eps = 0.1, pvalue=0.05) 
# {
#   if (is.mcmc.list(x))
#     return(lapply(x, heidel.diag.mod, eps))
#   x <- as.mcmc(as.matrix(x))
#   HW.mat0 <- matrix(0, ncol = 7, nrow = nvar(x))
#   dimnames(HW.mat0) <- list(varnames(x),
#                             c("stest", "start", "CMV.stat", "pvalue", "htest",
#                               "mean", "halfwidth"))
#   HW.mat <- HW.mat0
#   for (j in 1:nvar(x)) {
#     start.vec <- seq(from=start(x), to = end(x)/2, by=niter(x)/10)
#     Y <- x[, j, drop = TRUE]    
#     n1 <- length(Y)
#     ## Schruben's test for convergence, applied sequentially
#     ##
#     S0 <- spectrum0.ar(window(Y, start=end(Y)/2))$spec
#     converged <- FALSE
#     for (i in seq(along = start.vec)) {
#       Y <- window(Y, start = start.vec[i])
#       n <- niter(Y)
#       ybar <- mean(Y)
#       B <- cumsum(Y) - ybar * (1:n)
#       Bsq <- (B * B)/(n * S0)
#       I <- sum(Bsq)/n
#       if(converged <- !is.na(I) && pcramer(I) < 1 - pvalue)
#         break
#     }
#     ## Recalculate S0 using section of chain that passed convergence test
#     S0ci <- spectrum0.ar(Y)$spec
#     halfwidth <- 1.96 * sqrt(S0ci/n)
#     passed.hw <- !is.na(halfwidth) & (abs(halfwidth/ybar) <= eps)
#     if (!converged || is.na(I) || is.na(halfwidth)) {
#       nstart <- NA
#       passed.hw <- NA
#       halfwidth <- NA
#       ybar <- NA
#     }
#     else {
#       nstart <- start(Y)
#     }
#     HW.mat[j, ] <- c(converged, nstart, I, 1 - pcramer(I), 
#                      passed.hw, ybar, halfwidth)
#   }
#   # class(HW.mat) <- "heidel.diag"
#   return(HW.mat)
# }
# 
########### MOVED TO bsts_helper_functions ############################
# ###
# ## POSTERIOR PREDICTIVE CHECKS
# ###
# postPredChecks <- function(causimp, filename=NA, 
#                            save.plot=TRUE, return.val=FALSE,
#                            burn=NA, conv.alpha=0.05) {
#   
#   ppcheck.filename <- if (is.na(filename)){
#     sprintf('bsts_post_pred_checks_%s.png', round(10*as.numeric(Sys.time())) )
#   }  else {
#     filename
#   }
#   
#   ## output list of convergence checks (booleans) and residual dataframes
#   checklist <- list()
#   
#   response <-  causimp$series$response  
#   y <- response
# 
#   # ##**DEBUG**
#   # print('as.numeric( response )')
#   # print(response)
#   # ##**
#     
#   npds <- length(y)
#   
#   intpd <- causimp$model$post.period[1]
#   
#   niter <- length(causimp$model$bsts.model$sigma.obs)
#   
#   if (is.na(burn)) {
#     burn <- round( niter * .2 )
#   }
#   
#   ## Newdata (post-intervention data) to predict via BSTS
#   hasRegression <- causimp$model$bsts.model$has.regression
#   if (hasRegression) {
#     newdata <- causimp$model$bsts.model$predictors[1:(intpd-1), ]
#     newdata <- cbind(response=response[1:(intpd-1)], newdata)
#   } else {
#     newdata <- NULL
#   }
#   
#   # ##**DEBUG**
#   # print('newdata[1:10,]')
#   # print(newdata[1:10,])
#   # ##**
#   
#   # predict.mbsts()
#   post.pred <- if (hasRegression) {
#     predict.bsts(causimp$model$bsts.model, newdata = newdata , burn = burn) ## already knows horizon from 
#   } else {
#     predict.bsts(causimp$model$bsts.model, burn = burn, 
#                  horizon = intpd-1,
#                  olddata = response[1:(intpd-1)])
#   }
# 
#   # ##**DEBUG**
#   # plot(post.pred$mean, col='red',type='l',ylim=c(-2,6)); points(response[1:(intpd-1)], col='black', pch=16, main='DEBUG')
#   # print('post.pred:')
#   # print(post.pred)
#   # ##
#   
#   post.pred.dist <- post.pred$distribution
#   post.pred.mean <- post.pred$mean ##apply(post.pred$distribution, c(1), mean)
#   
#   ## before intervention period bool dummy
#   .ind <- (1:npds) < intpd
#   
#   if (save.plot) {
#     png(filename = ppcheck.filename, width = 15, height = 10, units = 'in', res = 400)
#   }
#   ##----------- INSIDE PNG PLOT --------------------------------------
#   par(mfrow=c(2,3), mar=c(2.5,2.5,2.5,1))
#   
#   ##===================
#   ## Posterior Predictive (Y) plots 
#   ##-------------------
#   
#   ##-----------
#   ## Trace plot of Posterior Predictive distribution Markov Chain 
#   post.pred.tr <- rowMeans(post.pred.dist)
#   ##
#   gd <- geweke.diag(post.pred.tr, .1, .5)
#   gd.z <- abs(gd$z) ## z statistic
#   gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
#   gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
#   ##
#   hd <- heidel.diag.mod(post.pred.tr)
#   hd.st.cmv <- hd[1,'CMV.stat'][1]
#   hd.st.p <- hd[1,'pvalue'][1]
#   hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
#   hd.st.start <- hd[1,'start'][1]
#   hd.hw.eps <- 0.1 ## default 0.1
#   hd.hw <- hd[1,'halfwidth'][1]
#   hd.hw.mean <- hd[1,'mean'][1]
#   hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
#   ##
#   rng <- range(post.pred.tr)
#   ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
#   plot(post.pred.tr, type='l', main='A. Posterior Predicted MCMC Trace, Y' ,
#        ylim=ylims
#   )
#   mtext.postpred <- sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
#                             gd.result,gd.z,gd.p,
#                             hd.st.result, hd.st.cmv, hd.st.p,
#                             hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps)
#   mtext(text = mtext.postpred, side = 3, line=-4.5, outer = F)
#   ##
#   checklist$ck.postpred <- list(
#     geweke = list(check=(gd.p >= conv.alpha), z=gd.z, p=gd.p),
#     hw.st  = list(check=(hd.st.p >= conv.alpha), cmv=hd.st.cmv, p=hd.st.p),
#     hw.hw  = list(check=(abs(hd.hw/hd.hw.mean) < hd.hw.eps), hw.mean=abs(hd.hw/hd.hw.mean), hw.eps=hd.hw.eps)
#   )
#   ##-----------
#   
#   ##-----------
#   ## DENSITIES
#   plot(density(y[.ind]),  main = "B. Density comparison, Y")
#   lines(density(post.pred.mean), lwd=2, lty=2, col='red')
#   legend('topleft', legend=c('observed','predicted'), lty=c(1,2), col=c('black','red'))
#   ##-----------
#   
#   ##-----------
#   ## HISTOGRAMS & BAYESIAN P-VALUES
#   # max.distrib <- apply(post.pred, c(2, 3), max)
#   max.distrib <- apply(post.pred$distribution, 1, max)
#   pvalue <- sum(max.distrib >= max(y[.ind]))/length(max.distrib)
#   hist(max.distrib, 30, col = "lightblue", border = "grey", 
#        main = paste0("C. Bayesian p-val (Max Y) = ", round(pvalue, 2)),
#        xlab = "Max. in-sample forecasts")
#   abline(v = max(y[.ind]), col = "darkblue", lwd = 3)
#   ##-----------
#   
#   ##===================
#   ## Std. Residual plots
#   ##-------------------
#   y.rep <- matrix(y[.ind], length(y[.ind]), (niter - burn),  byrow = FALSE)
#   # res <- y.rep - (post.pred.dist - CausalMBSTS$mcmc$eps.samples[, i, ])
#   res <-  y.rep - t(post.pred.dist)
#   std.res <- res / apply(res,1,sd)   ## [residual i] / [ stdev of residual i]
#   # std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
#   
#   ##-----------
#   ## Trace plot of Posterior Predictive distribution Markov Chain 
#   res.tr <- colMeans(std.res)
#   ##
#   gd <- geweke.diag(res.tr, .1, .5)
#   gd.z <- abs(gd$z) ## z statistic
#   gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
#   gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
#   ##
#   hd <- heidel.diag.mod(res.tr)
#   hd.st.cmv <- hd[1,'CMV.stat'][1]
#   hd.st.p <- hd[1,'pvalue'][1]
#   hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
#   hd.st.start <- hd[1,'start'][1]
#   hd.hw.eps <- 0.1 ## default 0.1
#   hd.hw <- hd[1,'halfwidth'][1]
#   hd.hw.mean <- hd[1,'mean'][1]
#   hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
#   ##
#   rng <- range(res.tr)
#   ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
#   plot(res.tr, type='l', main='D. Std.Residual MCMC Trace, Y',
#        ylim=ylims)
#   mtext.residual <- sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
#                             gd.result,gd.z,gd.p,
#                             hd.st.result, hd.st.cmv, hd.st.p,
#                             hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps)
#   mtext(text = mtext.residual, side = 3, line=-4.5, outer = F)
#   ##
#   checklist$ck.residual <- list(
#     geweke = list(check=(gd.p >= conv.alpha), z=gd.z, p=gd.p),
#     hw.st  = list(check=(hd.st.p >= conv.alpha), cmv=hd.st.cmv, p=hd.st.p),
#     hw.hw  = list(check=(abs(hd.hw/hd.hw.mean) < hd.hw.eps), hw.mean=abs(hd.hw/hd.hw.mean), hw.eps=hd.hw.eps)
#   )
#   ##-----------
#   
#   ##-----------
#   # Std. Residual Normality plot
#   qqnorm(rowMeans(std.res), main = "E. Std.Residual QQ-plot, Y")
#   qqline(rowMeans(std.res))
#   ##-----------
#   
#   ##-----------
#   ## Std Residual ACF
#   Acf(rowMeans(std.res), main = "");title(main='F. Std.Residual ACF, Y')
#   ##-----------
#   
#   ##----------- end PNG PLOT --------------------------------------
#   if (save.plot) {
#     dev.off()
#   }
#   
#   ##
#   if(return.val) {
#     checklist$std.residual <- std.res
#     checklist$postpred.dist <- post.pred.dist
#     checklist$summary <- sprintf('\nPosterior Predictive:-----\n%s\nStd. Residuals:-----\n%s\n\n', mtext.postpred, mtext.residual)
#     
#     ## ALL CONVERGENCE CHECKS
#     c1 <- unname( checklist$ck.postpred$geweke$check )
#     c2 <- unname( checklist$ck.postpred$hw.st$check )
#     c3 <- unname( checklist$ck.postpred$hw.hw$check )
#     c4 <- unname( checklist$ck.residual$geweke$check )
#     c5 <- unname( checklist$ck.residual$hw.st$check )
#     c6 <- unname( checklist$ck.residual$hw.hw$check )
#     ##
#     ck1 <- ifelse(is.na(c1) | is.nan(c1), FALSE, c1)
#     ck2 <- ifelse(is.na(c2) | is.nan(c2), FALSE, c2)
#     ck3 <- ifelse(is.na(c3) | is.nan(c3), FALSE, c3)
#     ck4 <- ifelse(is.na(c4) | is.nan(c4), FALSE, c4)
#     ck5 <- ifelse(is.na(c5) | is.nan(c5), FALSE, c5)
#     ck6 <- ifelse(is.na(c6) | is.nan(c6), FALSE, c6)
#     ##
#     checklist$converged <- c(
#       pp.ge = ck1, 
#       pp.st = ck2, 
#       pp.hw = ck3, 
#       r.ge  = ck4, 
#       r.st  = ck5, 
#       r.hw  = ck6
#     )
#     checklist$converged.all <- all( checklist$converged )
#     checklist$converged.prop <- sum(checklist$converged) / length(checklist$converged)
#     
#     return(checklist)
#   }
#   
# }



# 
# ###
# ## POSTERIOR PREDICTIVE CHECKS
# ###
# postPredChecksPreIntervention <- function(causimp, filename=NA, 
#                                            save.plot=TRUE, return.val=FALSE,
#                                            burn=NA, conv.alpha=0.05) {
#   
#   ppcheck.filename <- if (is.na(filename)){
#     sprintf('bsts_post_pred_checks_%s.png', round(10*as.numeric(Sys.time())) )
#   }  else {
#     filename
#   }
#   
#   ## output list of convergence checks (booleans) and residual dataframes
#   checklist <- list()
#   
#   response <-  causimp$series$response  
#   y <- response
#   
#   # ##**DEBUG**
#   # print('as.numeric( response )')
#   # print(response)
#   # ##**
#   
#   npds <- length(y)
#   
#   intpd <- causimp$model$post.period[1]
#   
#   niter <- length(causimp$model$bsts.model$sigma.obs)
#   
#   if (is.na(burn)) {
#     burn <- round( niter * .2 )
#   }
#   
#   ## Newdata (post-intervention data) to predict via BSTS
#   if (causimp$model$bsts.model$has.regression) {
#     newdata <- causimp$model$bsts.model$predictors[1:(intpd-1), ]
#     newdata <- cbind(response=response[1:(intpd-1)], newdata)
#   } else {
#     # newdata <- cbind(response=response[1:(intpd-1)], 
#     #                  intercept=rep(1, intpd-1))
#     newdata <- NULL
#   }
#   
#   # ##**DEBUG**
#   # print('newdata[1:10,]')
#   # print(newdata[1:10,])
#   # ##**
#   
#   # predict.mbsts()
#   post.pred <- predict.bsts(causimp$model$bsts.model, newdata = newdata , burn = burn, horizon=intpd-1)
#   ##**DEBUG**
#   print('post.pred:')
#   print(post.pred)
#   ##
#   post.pred.dist <- post.pred$distribution
#   post.pred.mean <- post.pred$mean ##apply(post.pred$distribution, c(1), mean)
#   
#   ## before intervention period bool dummy
#   .ind <- (1:npds) < intpd
#   
#   if (save.plot) {
#     png(filename = ppcheck.filename, width = 15, height = 10, units = 'in', res = 400)
#   }
#   ##----------- INSIDE PNG PLOT --------------------------------------
#   par(mfrow=c(2,2), mar=c(2.5,2.5,2.5,1))
#   
#   ##===================
#   ## Posterior Predictive (Y) plots 
#   ##-------------------
#   
#   ##-----------
#   ## Trace plot of Posterior Predictive distribution Markov Chain 
#   post.pred.tr <- rowMeans(post.pred.dist)
#   ##
#   gd <- geweke.diag(post.pred.tr, .1, .5)
#   gd.z <- abs(gd$z) ## z statistic
#   gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
#   gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
#   ##
#   hd <- heidel.diag.mod(post.pred.tr)
#   hd.st.cmv <- hd[1,'CMV.stat'][1]
#   hd.st.p <- hd[1,'pvalue'][1]
#   hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
#   hd.st.start <- hd[1,'start'][1]
#   hd.hw.eps <- 0.1 ## default 0.1
#   hd.hw <- hd[1,'halfwidth'][1]
#   hd.hw.mean <- hd[1,'mean'][1]
#   hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
#   ##
#   rng <- range(post.pred.tr)
#   ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
#   plot(post.pred.tr, type='l', main='A. Posterior Predicted MCMC Trace, Y' ,
#        ylim=ylims
#   )
#   mtext.postpred <- sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
#                             gd.result,gd.z,gd.p,
#                             hd.st.result, hd.st.cmv, hd.st.p,
#                             hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps)
#   mtext(text = mtext.postpred, side = 3, line=-4.5, outer = F)
#   ##
#   checklist$ck.postpred <- list(
#     geweke = list(check=(gd.p >= conv.alpha), z=gd.z, p=gd.p),
#     hw.st  = list(check=(hd.st.p >= conv.alpha), cmv=hd.st.cmv, p=hd.st.p),
#     hw.hw  = list(check=(abs(hd.hw/hd.hw.mean) < hd.hw.eps), hw.mean=abs(hd.hw/hd.hw.mean), hw.eps=hd.hw.eps)
#   )
#   ##-----------
#   
#   ##-----------
#   ## DENSITIES
#   plot(density(y[.ind]),  main = "B. Density comparison, Y")
#   lines(density(post.pred.mean), lwd=2, lty=2, col='red')
#   legend('topleft', legend=c('observed','predicted'), lty=c(1,2), col=c('black','red'))
#   ##-----------
#   
#   ##-----------
#   ## HISTOGRAMS & BAYESIAN P-VALUES
#   # max.distrib <- apply(post.pred, c(2, 3), max)
#   max.distrib <- apply(post.pred$distribution, 1, max)
#   pvalue <- sum(max.distrib >= max(y[.ind]))/length(max.distrib)
#   hist(max.distrib, 30, col = "lightblue", border = "grey", 
#        main = paste0("C. Bayesian p-val (Max Y) = ", round(pvalue, 2)),
#        xlab = "Max. in-sample forecasts")
#   abline(v = max(y[.ind]), col = "darkblue", lwd = 3)
#   ##-----------
#   
#   ##===================
#   ## Std. Residual plots
#   ##-------------------
#   y.rep <- matrix(y[.ind], length(y[.ind]), (niter - burn),  byrow = FALSE)
#   # res <- y.rep - (post.pred.dist - CausalMBSTS$mcmc$eps.samples[, i, ])
#   res <-  y.rep - t(post.pred.dist)
#   std.res <- res / apply(res,1,sd)   ## [residual i] / [ stdev of residual i]
#   # std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
#   
#   ##-----------
#   ## Trace plot of Posterior Predictive distribution Markov Chain 
#   res.tr <- colMeans(std.res)
#   ##
#   gd <- geweke.diag(res.tr, .1, .5)
#   gd.z <- abs(gd$z) ## z statistic
#   gd.p <- 2*pnorm(gd.z, mean = 0, sd = 1, lower.tail = F)  ## P-value for 2-sided test
#   gd.result <- ifelse(gd.p < conv.alpha, 'FAIL', 'PASS')
#   ##
#   hd <- heidel.diag.mod(res.tr)
#   hd.st.cmv <- hd[1,'CMV.stat'][1]
#   hd.st.p <- hd[1,'pvalue'][1]
#   hd.st.result <- ifelse( hd.st.p < conv.alpha, 'FAIL', 'PASS')
#   hd.st.start <- hd[1,'start'][1]
#   hd.hw.eps <- 0.1 ## default 0.1
#   hd.hw <- hd[1,'halfwidth'][1]
#   hd.hw.mean <- hd[1,'mean'][1]
#   hd.hw.result <- ifelse( abs(hd.hw/hd.hw.mean) < hd.hw.eps, 'PASS', 'FAIL' )
#   ##
#   rng <- range(res.tr)
#   ylims <- rng + c( -.05*diff(rng), .3*diff(rng) )
#   plot(res.tr, type='l', main='D. Std.Residual MCMC Trace, Y',
#        ylim=ylims)
#   mtext.residual <- sprintf('Geweke: %s (z=%.2f, p=%.2f)\nH&W Stationarity: %s (CMV=%.2f, p=%.2f)\nH&W Halfwidth: %s (hw/mean=%.2f < eps=%.2f)',
#                             gd.result,gd.z,gd.p,
#                             hd.st.result, hd.st.cmv, hd.st.p,
#                             hd.hw.result, abs(hd.hw/hd.hw.mean), hd.hw.eps)
#   mtext(text = mtext.residual, side = 3, line=-4.5, outer = F)
#   ##
#   checklist$ck.residual <- list(
#     geweke = list(check=(gd.p >= conv.alpha), z=gd.z, p=gd.p),
#     hw.st  = list(check=(hd.st.p >= conv.alpha), cmv=hd.st.cmv, p=hd.st.p),
#     hw.hw  = list(check=(abs(hd.hw/hd.hw.mean) < hd.hw.eps), hw.mean=abs(hd.hw/hd.hw.mean), hw.eps=hd.hw.eps)
#   )
#   ##-----------
#   
#   ##-----------
#   # Std. Residual Normality plot
#   qqnorm(rowMeans(std.res), main = "E. Std.Residual QQ-plot, Y")
#   qqline(rowMeans(std.res))
#   ##-----------
#   
#   ##-----------
#   ## Std Residual ACF
#   Acf(rowMeans(std.res), main = "");title(main='F. Std.Residual ACF, Y')
#   ##-----------
#   
#   ##----------- end PNG PLOT --------------------------------------
#   if (save.plot) {
#     dev.off()
#   }
#   
#   ##
#   if(return.val) {
#     checklist$std.residual <- std.res
#     checklist$postpred.dist <- post.pred.dist
#     checklist$summary <- sprintf('\nPosterior Predictive:-----\n%s\nStd. Residuals:-----\n%s\n\n', mtext.postpred, mtext.residual)
#     
#     ## ALL CONVERGENCE CHECKS
#     c1 <- unname( checklist$ck.postpred$geweke$check )
#     c2 <- unname( checklist$ck.postpred$hw.st$check )
#     c3 <- unname( checklist$ck.postpred$hw.hw$check )
#     c4 <- unname( checklist$ck.residual$geweke$check )
#     c5 <- unname( checklist$ck.residual$hw.st$check )
#     c6 <- unname( checklist$ck.residual$hw.hw$check )
#     ##
#     ck1 <- ifelse(is.na(c1) | is.nan(c1), FALSE, c1)
#     ck2 <- ifelse(is.na(c2) | is.nan(c2), FALSE, c2)
#     ck3 <- ifelse(is.na(c3) | is.nan(c3), FALSE, c3)
#     ck4 <- ifelse(is.na(c4) | is.nan(c4), FALSE, c4)
#     ck5 <- ifelse(is.na(c5) | is.nan(c5), FALSE, c5)
#     ck6 <- ifelse(is.na(c6) | is.nan(c6), FALSE, c6)
#     ##
#     checklist$converged <- c(
#       pp.ge = ck1, 
#       pp.st = ck2, 
#       pp.hw = ck3, 
#       r.ge  = ck4, 
#       r.st  = ck5, 
#       r.hw  = ck6
#     )
#     checklist$converged.all <- all( checklist$converged )
#     checklist$converged.prop <- sum(checklist$converged) / length(checklist$converged)
#     
#     return(checklist)
#   }
#   
# }

