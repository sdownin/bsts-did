##
## Dynamic Treatment Effects via Difference-in-Differences
##  - multiple period DiD
##  - staggered DiD
##  - - vs. static DiD
##
## Problem that BSTS can solve: 
##  How to bin event study design
##


library(tibble)
library(did)
library(bsts)
library(CausalImpact)
library(lubridate)
library(plyr);library(dplyr)
library(ggplot2)

work_dir <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS'
setwd(work_dir)

# ########## DEBUG ##########
# ref <- readRDS(file.path('D:\\BSTS_external','__GRIDSEARCH_output__16690657090_g1h1i1j1.rds'))
# ###########################


# dat <- read.csv(file.path(work_dir,'cconma_purchase_data.csv'))
dat <- read.csv(file.path(work_dir,'cconma_data_sellers.csv'))

dat <- as_tibble(dat)

names(dat)
View(dat)

dat$follower_reg_date <- ymd(dat$follower_reg_date)

## ADD GROUP LABELS for TREATMENT and CONTROL
dat$group <- NA
dat$group[ ! is.na(dat$follower_reg_date) ] <- 'treatment'
dat$group[   is.na(dat$follower_reg_date) ] <- 'control'


# ## weekly batches
# nbatch <- ceiling(as.numeric( 
#   (max(dat$follower_reg_date, na.rm=T) - min(dat$follower_reg_date, na.rm=T)) / 7
# )) ## weeks

# dat$follower_reg_day  <- day(dat$follower_reg_date)
dat$follower_reg_week  <- week(dat$follower_reg_date)
dat$follower_reg_month <- month(dat$follower_reg_date)
dat$follower_reg_quarter  <- quarter(dat$follower_reg_date)
dat$follower_reg_year  <- year(dat$follower_reg_date)

## Remove unreasonable age
dat <- dat[which(dat$age > 10 & dat$age < 110),]

# ## Make category of ages
# idx.age <- which(dat$age < 120 & dat$age > 0)
# summary(dat$age[idx.age])
# sd(dat$age[idx.age])
# hist(dat$age[idx.age]); abline(v=10*(1:10), col='red')
# ## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# ## 6.0    42.0    48.0    47.9    53.0   101.0 


## Use decade ("50s" as age category)
dat$age_cat <- sprintf('%ss', 10 * floor(dat$age / 10) )


# # dat$cohort_pd_week_num <- NA
# # dat$cohort_pd_week_t0 <- NA
# dat$order_pd_week_num <- NA
# dat$order_pd_week_t0 <- NA
# dat$cohort_pd_week <- NA
# dat$cohort_pd_week_num <- NA
# dat$match_id_pd_week <- NA

dat$order_pd_num <- NA
dat$order_pd_t0 <- NA
dat$cohort_pd <- NA
dat$cohort_pd_num <- NA
dat$match_id <- NA
dat$match_pd <- NA



# dat$order_day <- day(dat$order_date)
dat$order_week <- week(dat$order_date)
dat$order_month <- month(dat$order_date)
dat$order_quarter <- quarter(dat$order_date)
dat$order_year <- year(dat$order_date)


# ##======================================================
# ## Cohort and order period columns by choice of period length (month, week, ...)
# pd_type <- 'month'
# ##------------------------------------------------------

##
## RUN DID FOR SINGLE INTERVENTION AGGREGATION AT pd_type
##
runSingleInterventionDiDPdAgg <- function(dat, pd_type)
{

  
  dat$pd_type <- pd_type
  
  if (pd_type == 'month') {

    dat$order_pd <- apply(data.frame(yr=dat$order_year,w=dat$order_month),1,function(x)sprintf('%s_%02d',x[1],x[2]))
    dat$cohort_pd <- apply(data.frame(yr=dat$follower_reg_year,w=dat$follower_reg_month),1,function(x){
      if(any(is.na(x))) {
        return(NA)
      } else {
        return(sprintf('%s_%02d',x[1],x[2]))
      }
    })

  } else if (pd_type == 'week') {

    dat$order_pd <- apply(data.frame(yr=dat$order_year,w=dat$order_week),1,function(x)sprintf('%s_%02d',x[1],x[2]))
    dat$cohort_pd <- apply(data.frame(yr=dat$follower_reg_year,w=dat$follower_reg_week),1,function(x){
      if(any(is.na(x))) {
        return(NA)
      } else {
        return(sprintf('%s_%02d',x[1],x[2]))
      }
    })

  # } else if (pd_type == 'day') {
  #
  #   dat$order_pd <- apply(data.frame(yr=dat$order_year,w=dat$order_day),1,function(x)sprintf('%s_%02d',x[1],x[2]))
  #   dat$cohort_pd <- apply(data.frame(yr=dat$follower_reg_year,w=dat$follower_reg_day),1,function(x){
  #     if(any(is.na(x))) {
  #       return(NA)
  #     } else {
  #       return(sprintf('%s_%02d',x[1],x[2]))
  #     }
  #   })

  } else if (pd_type == 'quarter') {

    dat$order_pd <- apply(data.frame(yr=dat$order_year,w=dat$order_quarter),1,function(x)sprintf('%s_%02d',x[1],x[2]))
    dat$cohort_pd <- apply(data.frame(yr=dat$follower_reg_year,w=dat$follower_reg_quarter),1,function(x){
      if(any(is.na(x))) {
        return(NA)
      } else {
        return(sprintf('%s_%02d',x[1],x[2]))
      }
    })

  } else if (pd_type == 'year') {
    
    dat$order_pd <- apply(data.frame(yr=dat$order_year),1,function(x)sprintf('%s',x[1]))
    dat$cohort_pd <- apply(data.frame(yr=dat$follower_reg_year),1,function(x){
      if(any(is.na(x))) {
        return(NA)
      } else {
        return(sprintf('%s',x[1]))
      }
    })
    
  }
  
  # ## MONTH
  # dat$order_pd_month <- apply(data.frame(yr=dat$order_year,w=dat$order_month),1,function(x)sprintf('%s_%02d',x[1],x[2]))
  # dat$cohort_pd_month <- apply(data.frame(yr=dat$follower_reg_year,w=dat$follower_reg_month),1,function(x){
  #   if(any(is.na(x))) {
  #     return(NA)
  #   } else {
  #     return(sprintf('%s_%02d',x[1],x[2]))
  #   }
  # })
  # ## WEEK 
  # dat$order_pd_week <- apply(data.frame(yr=dat$order_year,w=dat$order_week),1,function(x)sprintf('%s_%02d',x[1],x[2]))
  # dat$cohort_pd_week <- apply(data.frame(yr=dat$follower_reg_year,w=dat$follower_reg_week),1,function(x){
  #   if(any(is.na(x))) {
  #     return(NA)
  #   } else {
  #     return(sprintf('%s_%02d',x[1],x[2]))
  #   }
  # })
  # ## YEAR  
  # dat$order_pd_year <- apply(data.frame(yr=dat$order_year,w=dat$order_year),1,function(x)sprintf('%s_%02d',x[1],x[2]))
  # dat$cohort_pd_year <- apply(data.frame(yr=dat$follower_reg_year,w=dat$follower_reg_year),1,function(x){
  #   if(any(is.na(x))) {
  #     return(NA)
  #   } else {
  #     return(sprintf('%s_%02d',x[1],x[2]))
  #   }
  # })
  
  
  
  
  # ## Cohort and order period columns by choice of period length (month, week, ...)
  # pd_type <- 'month'
  # order_pd_num_col <- sprintf('order_pd_%s_num', pd_type)
  # order_pd_col <- sprintf('order_pd_%s', pd_type)
  # order_pd_t0_col <- sprintf('order_pd_%s_t0', pd_type)
  # cohort_pd_col <- sprintf('cohort_pd_%s', pd_type)
  # cohort_pd_num_col <- sprintf('cohort_pd_%s_num', pd_type)
  
  ##=============================
  ## BUILD MATCHED COHORTS DATAFRAME
  ##  - MATCH TREATED-CONTROL SUBJECTS TO MAKE OUT DATAFRAME
  ##----------------------------

  ##===================================
  ## control population
  ##-----------------------------------
  co <- dat[ which(dat$group=='control'), ]
  
  ##====================================
  ## Treatment Group
  ##------------------------------------
  ## DROP missing week rows (NAs)
  df <- dat[ which(dat$group=='treatment'), ]
  
  
  ## PERIOD COHORTS
  pds.cohort <- sort(unique(df$cohort_pd), decreasing = F)
  npds.cohort <- length(pds.cohort)
  

  # for (t in 1:npds.cohort) {
  #   df[which(df[[cohort_pd_col]]==pds.cohort[t]), cohort_pd_num_col] <- t
  # }
  
  # 
  # df$order_week <- week(df$order_date)
  # df$order_year <- year(df$order_date)
  # 
  # df$order_pd_week <- apply(data.frame(yr=df$order_year,w=df$order_week),1,function(x)sprintf('%s_%02d',x[1],x[2]))
  # 
  
  ## Loop to create periods
  cohorts <- data.frame()
  
  ## Loop over treated cohort to get orders in period t
  ##  and match a control subject to add control orders in period t 
  for (t in 1:npds.cohort) 
  {
    cat(sprintf('\nt=%s\n',t))
    pd <- pds.cohort[t]
    
    ## indices of orders in this pd (week t)
    idx <- which(df$cohort_pd == pd)
    # ## members who made orders in this pd (week t)
    # mem.t <- sort(unique(df$mem_no[idx]))
    ## orders by members who made orders in this pd (ie, orders in this period)
    x <- df[idx,]
    
    # ## new cohort members in this pd (week t)
    # new <- unique(x$mem_no[which( ! x$mem_no %in% c(cohorts$treatment,cohorts$control))])
    ## new cohort members in this pd (week t)
    new <- unique(x$seller_mem_no[which( ! x$seller_mem_no %in% c(cohorts$treatment,cohorts$control))])
    
    # ## Update cohort time period
    # # df$cohort_pd_week_num[idx] <- t
    # # df$cohort_pd_week_t0[idx] <- t - 
    # ## Clock time
    # df$order_pd_week_num[idx] <- t 
    # ## Event time (for dynamic effect, event study study)
    # df$order_pd_week_t0[idx] <- t - df$cohort_pd_week_num[idx] + 1
    
  
    ## IF NEW MEMBERS 
    if (length(new)>0)
    {
      # 
      ##-----------------------
      ## Matched pair control group 
      ## - add matched pairs from co into df
      ## - else remove unmatched rows from df (if no matches remaining in co)
      ## - - 
      ##-----------------------
      ## loop though the members who made orders in this pd (week t)
      # for (j in 1:length(mem.t))
      for (j in 1:length(new))
      {
        cat(sprintf('j=%s ',j))
        MISSING_CONTROL_J <- TRUE
        
        ## j'th new member
        new.mem.j <- new[j]
        
        # ## j'th member who made order in this pd
        # xj <- x[ which(x$mem_no == new.mem.j),  ]
        ## j'th member who made order in this pd
        xj <- x[ which(x$seller_mem_no == new.mem.j),  ]
        
        ## Natch on observables
        married <- unique(xj$married)
        # age <- unique(xj$age)
        age_cat <- unique(xj$age_cat)
        sex <- unique(xj$sex)
        ## row indices of control orders data for members with matched observed attributes (married,age,sex) in this pd t
        id.co <- which(co$married == married & co$age_cat == age_cat & co$sex == sex & co$order_pd == pd)
        
        # ## Control members who are potential matches for treated subject j
        # mem.co <- sort(unique(co$mem_no[id.co]))
        ## Control members who are potential matches for treated subject j
        mem.co <- sort(unique(co$seller_mem_no[id.co]))
        
        if (length(mem.co)>0)
        {
          for (k in 1:length(mem.co)) {
            ## first matched member in control orders data
            mem.co.k <- mem.co[k]
            
            if (MISSING_CONTROL_J &  ! mem.co.k %in% c(cohorts$treatment,cohorts$control) ) {
  
                ## ADD COHORT PAIRS
                ## add new cohort for this pd to the cohort dataframe
              batch <- data.frame(treatment=new.mem.j, control=mem.co.k, pd=pd, 
                                  pd_match_id=sprintf('%s_%s__%s_%s',pd_type,pd,new.mem.j,mem.co.k))
              cohorts <- rbind(cohorts, batch)
              
              ## SWITCH FLAG TO EXIT CONTROL MATCH SEARCH LOOP
              MISSING_CONTROL_J <- FALSE 
            }
  
          } ## // end k loop (over control members for j'th treated member who ordered in this period)
  
        } ## // end if (has matched control mem_no)
        
      } ## // end j loop (over members who ordered in this period)
      
      
    }  else { ## // end if (has orders in this period)
      cat(sprintf('\nzero new for j=%s\n',j))
    }
  
  } ## // end t loop (over weekly periods)
  
  
  ## Combined output dataframe with order data for matched treatment and control subjects
  # out <- data.frame()
  
  ##----------------------
  ## MATCHED SAMPLE 
  ##----------------------
  # ms <- dat[which(dat$mem_no %in% c(cohorts$treatment,cohorts$control)), ]
  ms <- dat[which(dat$seller_mem_no %in% c(cohorts$treatment,cohorts$control)), ]
  
  
  
  
  
  # pds <- sort(unique(dat$order_pd_week),decreasing = F)
  pds <- sort(unique(ms$order_pd),decreasing = F)
  npds <- length(pds)
  
  match_pd_num_single <- which( pds == pds.cohort[1] )
  ms$gname <- 0
  ms$gname[ms$group=='treatment'] <- match_pd_num_single
  # ms$gname <- as.factor(ms$gname)
  
  ## Add to cohort period label to control group 
  ##  (which don't have cohort groups initially bc no follower_reg_date)
  for (i in 1:nrow(cohorts)) {
    cat(sprintf('i=%s ',i))
    ctrl_mem <- cohorts$control[i]
    trt_mem <- cohorts$treatment[i]
    pd <- cohorts$pd[i]
    # # ms$cohort_pd[which(ms$mem_no==ctrl_mem)] <- cohorts$pd
    # ms$cohort_pd[which(ms$mem_no==ctrl_mem)] <- pd
    # ms$match_id[which(ms$mem_no %in% c(trt_mem,ctrl_mem))] <- cohorts$pd_match_id[i]
    # ms$match_pd[which(ms$mem_no %in% c(trt_mem))] <- pd
    ms$cohort_pd[which(ms$seller_mem_no==ctrl_mem)] <- pd
    ms$match_id[which(ms$seller_mem_no %in% c(trt_mem,ctrl_mem))] <- cohorts$pd_match_id[i]
    ms$match_pd[which(ms$seller_mem_no %in% c(trt_mem))] <- pd
  }
  
  ## Add cohort_pd_num for all treatment and control subjects 
  ## (after having adding control cohort period in cohort loop above)
  for (t in 1:npds.cohort) {
    pd <- pds.cohort[t]
    ms$cohort_pd_num[which(ms$cohort_pd == pd)] <- t
  }
  
  # 
  # ## Add cohort time period numbering to control matches for their treated subject
  # for (t in 1:npds.cohort) {
  #   pd <- pds.cohort[t]
  # 
  #   id.c <- which( matchsamp[[cohort_pd_col]] == pd & matchsamp$group == 'control' )
  #   id.t <- which( matchsamp[[cohort_pd_col]] == pd & matchsamp$group == 'treatment')
  # 
  #   matchsamp[idx,]
  # 
  # 
  #   idx <- which(dat[[cohort_pd_col]] == pd)
  #   ## *** UPDATE cohort numbering in original dataframe 'dat' ***
  #   dat[ which(dat[[cohort_pd_col]] == pd)  , cohort_pd_num_col] <- t
  # 
  # 
  # }
  # 
  # 
  # 
  # # pds <- sort(unique(dat$order_pd_week),decreasing = F)
  # pds <- sort(unique(ms$order_pd),decreasing = F)
  # npds <- length(pds)
  # 
  
  
  ## Add order period week number and staggered pd number - event time
  for (t in 1:npds) {
    cat(sprintf('\nt=%s\n',t))
    pd <- pds[t]
    
    ## indices of orders in this pd (t)
    idx <- which(ms$order_pd == pd)
    # ## members who made orders in this pd (t)
    # mem.t <- sort(unique(dat$mem_no[idx]))
    # ## orders by members who made orders in this pd (ie, orders in this period)
    # x <- dat[idx,]
    # 
    # ## new cohort members in this pd (t)
    # new <- unique(x$mem_no[which( ! x$mem_no %in% cohorts$mem_no)])
    
    ## Update cohort time period
    # df$cohort_pd_week_num[idx] <- t
    # df$cohort_pd_week_t0[idx] <- t - 
    ## Clock time
    ms$order_pd_num[idx] <- t 
    ## Event time (for dynamic effect, event study study)
    ms$order_pd_t0[idx] <- t - ms$cohort_pd_num[idx] #+ 1 
    
  }
  
  # plyr::count( ms$order_pd_t0)
  
  
  
  # ## Save output to serialized list (long runtime for out dataframe combining matched control and )
  # saveRDS(list(matchsamp=ms, cohorts=cohorts), 
  #         file = sprintf('cconma_purchase_data_intervention_matched_sample_%spairs.rds',nrow(cohorts)))
  
  
  
  # ## SET GROUP = 0 if untreated (control)
  # ms$cohort_pd_num[which(ms$group == 'control')] <- 0
  # ##
  # msatt <- att_gt(yname = 'order_sum', ## "Y",
  #                   tname = 'order_pd_num',
  #                   idname = 'mem_no',
  #                   gname = 'cohort_pd_num',
  #                   xformla = ~age, #+ factor(married) + factor(sex), ##age_mean + married_y_prop + sex_f_prop,
  #                   data = ms,
  #                 control_group = 'notyettreated'
  # )
  
  ms$order_pd_t_single <- 0
  ms$order_pd_t_single[ms$order_pd_t0 >= 1] <- 1
  
  # plyr::count(ms$order_pd_num)

    
  # ms$match_pd_num_single <- 0
  # ms$match_pd_num_single[ms$group=='treatment'] <- 
    
  # #########
  # ##------------------------------
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
  # ##################
  
  
  
  
  ##===============================================
  ## MONTHLY PERIOD
  ##-----------------------------------------------
  ## Summarize period panel dataframe for period: MONTH
  matchdata <- ms %>%   
    # group_by(order_pd_num, cohort_pd_num, group, age_cat, sex, married) %>%    ## order_pd_t0
    group_by(order_pd_num, gname, group, age_cat, sex, married) %>%    ## order_pd_t_single
    summarise(
      n_in_pd = n(),
      id_covar = paste(c(unique(age_cat),unique(sex),unique(married)),collapse = '_'),
      order_cnt_percap_mean = sum(order_cnt,na.rm=T)/n(),
      order_sum_percap_mean = sum(order_sum,na.rm=T)/n(),
      ##
      order_cnt_mean = mean(order_cnt, na.rm=T),
      order_cnt_max = max(order_cnt, na.rm=T),
      order_cnt_min = min(order_cnt, na.rm=T),
      order_cnt_sd = sd(order_cnt, na.rm=T),
      order_cnt_tot = sum(order_cnt, na.rm=T),
      order_sum_mean = mean(order_sum, na.rm=T),
      order_sum_max = max(order_sum, na.rm=T),
      order_sum_min = min(order_sum, na.rm=T),
      order_sum_sd = sd(order_sum, na.rm=T),
      order_sum_tot = sum(order_sum, na.rm=T),
      age_mean = mean(age, na.rm=T),
      age_cat_mode = mode(age_cat),
      married_y_prop = sum(married=='Y',na.rm=T)/n(),
      sex_f_prop = sum(sex=='F',na.rm=T)/n()
    )# %>%
    #filter(
    #  order_sum_mean >  quantile(order_sum_mean, .005),  ##Windsorize at 99%
    #  order_sum_mean <= quantile(order_sum_mean, .995)
    #)
  
  ## PLOT
  p0 <- ggplot(data=matchdata, aes(x=order_pd_num, y=order_sum_mean, colour=group)) + geom_line() + 
    facet_wrap( . ~ gname) +
    theme_bw() ## +  geom_vline(xintercept=0, lty=2)
    # facet_wrap( . ~ cohort_pd_num) + theme_bw() ## +  geom_vline(xintercept=0, lty=2)
  p0
  
  # matchdata <- out %>%   
  #   group_by(order_pd_week_num, group) %>%
  #   summarise(
  #     n_in_pd = n(),
  #     order_pd_week_cnt_mean = mean(order_cnt, na.rm=T),
  #     order_pd_week_cnt_max = max(order_cnt, na.rm=T),
  #     order_pd_week_cnt_min = min(order_cnt, na.rm=T),
  #     order_pd_week_cnt_sd = sd(order_cnt, na.rm=T),
  #     order_pd_week_sum_mean = mean(order_sum, na.rm=T),
  #     order_pd_week_sum_max = max(order_sum, na.rm=T),
  #     order_pd_week_sum_min = min(order_sum, na.rm=T),
  #     order_pd_week_sum_sd = sd(order_sum, na.rm=T)
  #   ) %>%
  #   filter(
  #     # n > 1,
  #     # mass > 50
  #   )
  
  
  
  # ## Loop over treated cohort to get orders in period t
  # ##  and match a control subject to add control orders in period t 
  # for (i in 1:nrow(cohorts)) 
  # {
  #   
  #   treat <- 
  #   ctrl <- co[which(co$mem_no %in% cohorts$control[i]), ]
  #   
  # } ## // end t loop (over weekly periods)
  
  ## SET GROUP = 0 if untreated (control)
  # matchdata$cohort_pd_num[which(matchdata$group == 'control')] <- 0
  matchdata$id_covar_factor <- as.factor(matchdata$id_covar)
  matchdata$id_covar_num <- as.numeric(matchdata$id_covar_factor)
  
  ##-----------------------------------------------------
  ### FILL IN zeros FOR ALL 'actor' series missing by period
  # for () {
  #   for () {
  #     
  #   }
  # }
  # ------------------------------------------------------
  
  # matchdata$gname <- as.numeric(as.character(matchdata$gname))
  # 
  # attgt <- att_gt(yname = "order_sum_mean", ## "order_sum_percap_mean", ## "Y",  
  #                   tname = "order_pd_num",
  #                   idname = "id_covar_num",
  #                   gname = "gname",
  #                   xformla = ~age_mean + married_y_prop + sex_f_prop,
  #                   data = matchdata , 
  #                   panel = F  ## NOT BALANCED
  #                   
  # )


  return(list(#attgt=attgt, 
              cohorts=cohorts,
              matchdata=matchdata,
              matchsample=ms#, 
              # p0=p0
              ))
}  ## END FUNCTION

###################################################################################

l.y <- runSingleInterventionDiDPdAgg(dat, pd_type = 'year')
l.q <- runSingleInterventionDiDPdAgg(dat, pd_type = 'quarter')
l.m <- runSingleInterventionDiDPdAgg(dat, pd_type = 'month')
# l.w <- runSingleInterventionDiDPdAgg(dat, pd_type = 'week')


.yname <- 'order_sum_percap_mean'
attgt.y <- att_gt(yname = .yname, ## "order_sum_percap_mean", ## "Y",  
                tname = "order_pd_num",
                idname = "id_covar_num",
                gname = "gname",
                xformla = ~age_mean + married_y_prop + sex_f_prop,
                data = l.y$matchdata , 
                panel = F)  ## NOT BALANCED
attgt.q <- att_gt(yname = .yname, ## "order_sum_percap_mean", ## "Y",  
                  tname = "order_pd_num",
                  idname = "id_covar_num",
                  gname = "gname",
                  xformla = ~age_mean + married_y_prop + sex_f_prop,
                  data = l.q$matchdata , 
                  panel = F)  ## NOT BALANCED
attgt.m <- att_gt(yname = .yname, ## "order_sum_percap_mean", ## "Y",  
                  tname = "order_pd_num",
                  idname = "id_covar_num",
                  gname = "gname",
                  xformla = ~age_mean + married_y_prop + sex_f_prop,
                  data = l.m$matchdata , 
                  panel = F)  ## NOT BALANCED



# attgt.q <- l.q$attgt
# attgt.m <- l.m$attgt

# agg.q.s <- aggte(attgt.q, type='simple')
# summary(agg.q.s)

## DYNAMIC EFFECTS AND EVENT STUDIES
agg.y.e <- aggte(attgt.y, type = "dynamic")
summary(agg.y.e)
ggdid(agg.y.e)

agg.q.e <- aggte(attgt.q, type = "dynamic")
summary(agg.q.e)
ggdid(agg.q.e)

agg.m.e <- aggte(attgt.m, type = "dynamic")
summary(agg.m.e)
ggdid(agg.m.e)


# ## Group-Specific Effects
# agg.gs <- aggte(example_attgt, type = "group")
# summary(agg.gs)
# ggdid(agg.gs)


# ## Calendar Time Effects
# agg.ct <- aggte(example_attgt, type = "calendar")
# summary(agg.ct)
# ggdid(agg.ct)


.yname <- 'order_cnt'
df2 <- l.x$matchsample
df2$order_pd_num_2 <- 0
intpd <- which( pds == '2016_01')
df2$order_pd_num_2[df2$order_pd_num > intpd] <- 1
df2$gname[df2$gname > 0] <- 1
attgt.2 <- att_gt(yname = .yname, ## "order_sum_percap_mean", ## "Y",  
                  tname = "order_pd_num_2",
                  idname = "seller_mem_no",
                  gname = "gname",
                  # xformla = ~age_mean + married_y_prop + sex_f_prop,
                  data = df2 , 
                  panel = F)  ## NOT BALANCED



#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
l.x <- l.m
##
pds <- sort(unique(l.x$matchsample$order_pd))
npds <- length(pds)
cohort.pds <- sort(unique(l.x$matchsample$cohort_pd))
intpd <- which( pds == cohort.pds[1] )
# pds <- sort(unique(l.y$matchsample$order_pd_num))
# cohort.pds <- sort(unique(l.y$matchsample$cohort_pd_num))
ms <- l.x$matchsample

ms <- ms[! is.na(ms$match_pd), ]


##====================
## BSTS Timseries Setup
##--------------------
## Aggregate into timeseries dataframe
tsdf <- ms %>% ## simdf %>%
  dplyr::filter( ! is.na(match_id)) %>%
  group_by(order_pd_num, gname, group) %>%   ##match_pd
  dplyr::summarize(
    n_in_pd = n(),
    actors = paste(unique(seller_mem_no), collapse = '|'),
    ##
    y_cnt_percap_mean = sum(order_cnt,na.rm=T)/n(),
    y_sum_percap_mean = sum(order_sum,na.rm=T)/n(),
    ##
    y_sum_mean = mean(order_sum, na.rm=T),
    y_sum_sum = sum(order_sum, na.rm=T),
    y_sum_sd = sd(order_sum, na.rm=T),
    y_sum_min = min(order_sum, na.rm=T),
    y_sum_max = max(order_sum, na.rm=T),
    ##
    y_cnt_mean = mean(order_cnt, na.rm=T),
    y_cnt_sum = sum(order_cnt, na.rm=T),
    y_cnt_sd = sd(order_cnt, na.rm=T),
    y_cnt_min = min(order_cnt, na.rm=T),
    y_cnt_max = max(order_cnt, na.rm=T),
    ##
    age_mean = mean(age, na.rm=T),
    age_cat_mode = mode(age_cat),
    married_y_prop = sum(married=='Y',na.rm=T)/n(),
    sex_f_prop = sum(sex=='F',na.rm=T)/n()
    # ##
    # x1_mean = mean(x1, na.rm=T),
    # x2_mean = mean(x2, na.rm=T),
    # x3_mean = mean(x3, na.rm=T),
    # ##
    # c1_mean = mean(c1, na.rm=T),
    # c2_mean = mean(c2, na.rm=T),
    # c3_mean = mean(c3, na.rm=T),
    # #
    # b1_mean = mean(b1, na.rm=T),
    # b2_mean = mean(b2, na.rm=T),
    # b3_mean = mean(b3, na.rm=T),
    # #
    # u_mean = mean(u, na.rm=T),
    # v_mean = mean(v, na.rm=T)
  )
tsdf$.id <- 1:nrow(tsdf)

## MAKE WIDE TIMESERIES FOR treatment,control groups in n periods
val.cols <- c('y_cnt_percap_mean','y_sum_percap_mean',
              'y_sum_mean','y_sum_sum','y_sum_sd','y_sum_min','y_sum_max',
              'y_cnt_mean','y_cnt_sum','y_cnt_sd','y_cnt_min','y_cnt_max',
              'age_mean','age_cat_mode','married_y_prop','sex_f_prop')
ts <- unique(tsdf$order_pd_num)
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
dat <- tsdfw[,c('treatment_y_cnt_percap_mean','control_y_cnt_percap_mean','control_y_sum_sd',#'control_y_sum', 'control_y_min','control_y_max',
                'control_age_mean','control_married_y_prop','control_sex_f_prop'#,
                # 'treatment_c1_mean','treatment_c2_mean','treatment_c3_mean',
                # 'control_u_mean','control_v_mean'
)]
## Train on y pre-treatment but NA's post-treatment
y.pre.treat.NAs.post.treat <- c(dat$treatment_y_cnt_percap_mean[1:(intpd-1)], rep(NA,npds-intpd+1))
## Then use the post-treatment response for causal impact estimation
post.period.response <- dat$treatment_y_cnt_percap_mean[intpd:npds]
## Covariates (predictors) - Matrix for "formula = y ~ predictors" argument
predictors <- dat[, ! names(dat) %in% 'treatment_y_cnt_percap_mean'] ## remove response; convert to matrix
# ## Covariates (predictors) - Dataframe for "data" argument
# predictors <- as.matrix(predictors) 

## ADD temporal trend to covariates
predictors$covar_temp_trend <- (1:npds) + rnorm(npds, 0, 1.2)   ## * simlist$`0.rand.base`$b5

st.sp <- list()
st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post.treat)
st.sp <- AddTrig(st.sp, y.pre.treat.NAs.post.treat, period = 12, frequencies = c(1))

bsts.niter <- 1000
bsts1 <- bsts(y.pre.treat.NAs.post.treat ~ . ,
     state.specification = st.sp,
     data = predictors,
     niter = bsts.niter)

summary(bsts1)
plot(bsts1)





##
plot(bsts1)
PlotBstsComponents(bsts1) ## , main=sprintf('BSTS Components: %s: %s',effect.type,paste(state.comps,collapse = ' + '))
PlotBstsState(bsts1)
PlotBstsResiduals(bsts1)
PlotBstsPredictionErrors(bsts1)
# PlotBstsForecastDistribution(bsts.model, main=sprintf('BSTS Forecast Dist: %s: %s',effect.type,paste(state.comps,collapse = ' + ')))
PlotBstsSize(bsts1)
##

# ## BSTS model for Dynamic Regression
# bsts.model <- bsts(y.pre.treat.NAs.post.treat,
#                    state.specification = st.sp,
#                    niter = 5000)

## Use BSTS prediction of counterfactual to estimate CausalImpact
impact_amount <- CausalImpact(bsts.model=bsts1,
                              post.period.response = post.period.response,
                              alpha=0.05, model.args = list(niter = bsts.niter))
# ##
# summary(impact_amount)
# summary(impact_amount$model$bsts.model)
# plot(impact_amount)

# summary(impact_amount)
# png(filename=sprintf('single_intervention_BSTS_CausalImpact_plot_%s_%s_%s.png',
#                         key,effect.type,sim.id))
(p.bsts.impact.all <- plot(impact_amount, c('original','pointwise','cumulative')) ) # pointwise','cumulative

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################





























ccattgt <- att_gt(yname = "Y", ## "Y",
                  tname = "period",
                  idname = "id",
                  gname = "G",
                  xformla = ~X,
                  data = matchdata
)
















##
##
sp <- reset.sim()

set.seed(1814)

time.periods <- 4

## cannot reproduce results at 
## https://cran.r-project.org/web/packages/did/vignettes/did-basics.html
## sp$te.e <- 1:time.periods
sp$te.e <- (1:time.periods) - 1
 
dta <- build_sim_dataset(sp)

nrow(dta)

head(dta)

example_attgt <- att_gt(yname = "Y",
                        tname = "period",
                        idname = "id",
                        gname = "G",
                        xformla = ~X,
                        data = dta
)

summary(example_attgt)

ggdid(example_attgt)


##------------------------
## AGGREGATION
##------------------------


agg.simple <- aggte(example_attgt, type='simple')
summary(agg.simple)

## DYNAMIC EFFECTS AND EVENT STUDIES
agg.es <- aggte(example_attgt, type = "dynamic")
summary(agg.es)
ggdid(agg.es)

## Group-Specific Effects
agg.gs <- aggte(example_attgt, type = "group")
summary(agg.gs)
ggdid(agg.gs)


## Calendar Time Effects
agg.ct <- aggte(example_attgt, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)



##-------------------------
## Alternative Control Groups
##-------------------------
## Not yet treated
example_attgt_altcontrol <- att_gt(yname = "Y",
                                   tname = "period",
                                   idname = "id",
                                   gname = "G",
                                   xformla = ~X,
                                   data = dta,
                                   control_group = "notyettreated"          
)
summary(example_attgt_altcontrol)




##------------------------------
## EXAMPLE DATA 
##  - Minimum Wage effect on county-level teen employment rates
##------------------------------
data(mpdta)
head(mpdta)

# estimate group-time average treatment effects without covariates
mw.attgt <- att_gt(yname = "lemp",
                   gname = "first.treat",
                   idname = "countyreal",
                   tname = "year",
                   xformla = ~1,
                   data = mpdta,
)

# summarize the results
summary(mw.attgt)

ggdid(mw.attgt, ylim = c(-.3,.3))

# aggregate the group-time average treatment effects
mw.dyn <- aggte(mw.attgt, type = "dynamic")
summary(mw.dyn)
ggdid(mw.dyn, ylim = c(-.3,.3))


## Balance across group exposure levels
## (avoid confounding dynamics and selective treatment timing among different groups)
mw.dyn.balance <- aggte(mw.attgt, type = "dynamic", balance_e=1)
summary(mw.dyn.balance)














