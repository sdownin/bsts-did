library(DBI)
library(RMySQL)
library(dplyr)
library(zoo)
library(xts)
library(forecast)
library(Boom)
library(BoomSpikeSlab)
library(bsts)
library(CausalImpact)
library(lubridate)

############   Seller analysis  #################
### p1: Sellers who adopted the social network
### p2: Sellers who did not adopt the social network
### p1: treatment, p0: control, apply CausalImpact library

## Load data
data<- read.csv("E:/My Drive/Research/Collaboration/With_SDowning/BSTS/R_scripts/purchase_seller.csv")

data_in <- data[data$s_in_nw=='1',]  ## purchase data which sellers joined the network
data_out <- data[data$s_in_nw=='0',]  ## purchase data which sellers did not join the network                 

count(distinct(data_in, seller_mem_no))  ## 55
count(distinct(data_out, seller_mem_no))  ## 292

## Aggregate purchase data to daily, US$ 1000
p1 <- data_in %>% group_by(order_date)%>%summarise(sales_m1 = sum(rev_krw_sum)/1000000)
p1$order_date <- as.Date(p1$order_date)

p0 <- data_out %>% group_by(order_date)%>%summarise(sales_m0 = sum(rev_krw_sum)/1000000)
p0$order_date <- as.Date(p0$order_date)

## Convert daily data to weekly
p1$week <- floor_date(p1$order_date, "week")
p1_w <- p1 %>% group_by(week)%>%summarise(sales_w=sum(sales_m1))

p0$week <- floor_date(p0$order_date, "week")
p0_w <- p0 %>% group_by(week)%>%summarise(sales_w=sum(sales_m0))

## Limit the time period from 2012-12-30 to 2018-12-31
p1_w <- p1_w[53:366,]
p0_w <- p0_w[53:366,]

## Convert the dataframe to a time series
p1_w_ts <- xts(p1_w$sales_w, p1_w$week)
p0_w_ts <- xts(p0_w$sales_w, p0_w$week)

plot.xts(p1_w_ts)
plot.xts(p0_w_ts)

## Define pre- and post-intervention periods
pre.period <- as.Date(c("2012-12-30", "2016-01-17"))
post.period <- as.Date(c("2016-01-24", "2018-12-30"))

data <- cbind(p1_w_ts, p0_w_ts)

## CausalImpact analysis using the default model in the package
impact_s <- CausalImpact(data, pre.period, post.period, alpha=0.05, 
                         model.args=list(nseasons = 52,season.duration = 1,niter = 10000))
plot(impact_s)
impact_s$model$bsts.model$state.specification
impact_s$model$bsts.model$prior


####################################################################
########### CausalImpact analysis using the custom bsts model

y.post <- p1_w_ts
x1 <- p0_w_ts

post.period.response <- p1_w$sales_w[161:314]
pre.period.response <- p1_w$sales_w[1:160]
y0 = (p1_w$sales_w[1]-mean(pre.period.response, na.rm=T))/sd(pre.period.response, na.rm = T)
      
y.post[161:314,] <- NA

#Replicate baseline with manual spec
sigma_epsilon = c(1,0.01) #Sets average variance and dof hyperparameter for observed
sigma_delta   = c(0.01,32 ) #Sets average variance and dof hyperparameter for trend (state)
sigma_gamma   = c(0.01,0.01 ) #Sets average variance and dof hyperparameter for seasonality (state)

levelsdprior      =  SdPrior(sigma_delta[1], sample.size = sigma_delta[2] , fixed = FALSE, upper.limit =1)
seasonsdprior     =  SdPrior(sigma_gamma[1], sample.size = sigma_gamma[2] , fixed = FALSE, upper.limit =1)
priorobs          =  SdPrior(sigma_epsilon[1], sample.size =sigma_epsilon[2] ,fixed = FALSE, upper.limit =1.2)

priorvar              = 1
seasonvar             = 1
initialvalueprior     = NormalPrior(y0, priorvar, fixed = FALSE)
initialseasonprior    = NormalPrior(0, priorvar, fixed = FALSE)

st.sp <- list()

st.sp <- AddLocalLevel(st.sp, y.post, sigma.prior =levelsdprior, initial.state.prior =initialvalueprior)
st.sp <- AddSeasonal(st.sp, y.post, nseasons =52, season.duration = 1, sigma.prior=seasonsdprior,initial.state.prior=initialseasonprior)

AddHolidayEffect <- function(y, dates, effect) {
  ## Adds a holiday effect to simulated data.
  ## Args:
  ##   y: A zoo time series, with Dates for indices.
  ##   dates: The dates of the holidays.
  ##   effect: A vector of holiday effects of odd length.  The central effect is
  ##     the main holiday, with a symmetric influence window on either side.
  ## Returns:
  ##   y, with the holiday effects added.
  time <- dates - (length(effect) - 1) / 2
  for (i in 1:length(effect)) {
    y[time] <- y[time] + effect[i]
    time <- time + 1
  }
  return(y)
}

cny.day <- NamedHoliday("NewYearsDay") ##" " should be chosen from a fixed list
cny.day.effect <- c(.3, 2, .2)
cny.day.dates <- as.Date(c("2013-01-27", "2014-01-19", "2015-02-08", "2016-01-24", "2017-01-15", "2018-02-04"))
y.post <- AddHolidayEffect(y.post, cny.day.dates, cny.day.effect)

ctg.day <- NamedHoliday("Thanksgiving")
ctg.day.effect <- c(.3, 2, .2)
ctg.day.dates <- as.Date(c("2013-09-08", "2014-08-31", "2015-09-13", "2016-09-04", "2017-09-17", "2018-09-09"))
y.post <- AddHolidayEffect(y.post, ctg.day.dates, ctg.day.effect)

holiday.list <- list(cny.day, ctg.day)
number.of.holidays <- length(holiday.list)

st.sp <- AddRegressionHoliday(st.sp, y.post, holiday.list = holiday.list)


bsts_s1 <- bsts(y.post ~ x1, st.sp, niter = 5000, expected.model.size=1, ping=0)
plot.bsts(bsts_s1)
impact_s1 <- CausalImpact(bsts.model = bsts_s1,
                       post.period.response = post.period.response,
                       alpha=0.05, model.args=list(niter = 10000))
plot(impact_s1)

plot(impact_s$model$bsts.model, "coefficients")
plot(impact_s1$model$bsts.model, "coefficients")

impact_s$model$bsts.model$state.specification
impact_s1$model$bsts.model$state.specification
impact_s1$model$bsts.model$prior

summary(impact_s)
summary(impact_s1)
######## Default model options
# impact_s$model$bsts.model$state.specification
# [[1]]
# $name
# [1] "trend"
# 
# $sigma.prior
# $prior.guess
# [1] 0.01
# 
# $prior.df
# [1] 32
# 
# $initial.value
# [1] 0.01
# 
# $fixed
# [1] FALSE
# 
# $upper.limit
# [1] 1
# 
# attr(,"class")
# [1] "SdPrior"         "DiffDoubleModel" "DoubleModel"     "Prior"          
# 
# $initial.state.prior
# $mu
# 1 
# -0.004314896 
# 
# $sigma
# [1] 1
# 
# $initial.value
# 1 
# -0.004314896 
# 
# $fixed
# [1] FALSE
# 
# attr(,"class")
# [1] "NormalPrior"     "DiffDoubleModel" "DoubleModel"     "Prior"          
# 
# $size
# [1] 1
# 
# attr(,"class")
# [1] "LocalLevel" "StateModel"
# 
# [[2]]
# $name
# [1] "seasonal.52.1"
# 
# $nseasons
# [1] 52
# 
# $season.duration
# [1] 1
# 
# $sigma.prior
# $prior.guess
# [1] 0.01
# 
# $prior.df
# [1] 0.01
# 
# $initial.value
# [1] 0.01
# 
# $fixed
# [1] FALSE
# 
# $upper.limit
# [1] 1
# 
# attr(,"class")
# [1] "SdPrior"         "DiffDoubleModel" "DoubleModel"     "Prior"          
# 
# $initial.state.prior
# $mu
# [1] 0
# 
# $sigma
# [1] 1
# 
# $initial.value
# [1] 0
# 
# $fixed
# [1] FALSE
# 
# attr(,"class")
# [1] "NormalPrior"     "DiffDoubleModel" "DoubleModel"     "Prior"          
# 
# $size
# [1] 51
# 
# attr(,"class")
# [1] "Seasonal"   "StateModel"
# 
# > impact_s$model$bsts.model$prior
# $prior.inclusion.probabilities
# [1] 1 1
# 
# $mu
# [1] 0 0
# 
# $sigma.guess
# [1] 0.4472136
# 
# $prior.df
# [1] 50
# 
# $sigma.upper.limit
# [1] 1.2
# 
# $max.flips
# [1] -1
# 
# $siginv
# (Intercept)       p0_w_ts
# (Intercept)  0.0100000000 -0.0007509018
# p0_w_ts     -0.0007509018  0.0098560645
# 
# attr(,"class")
# [1] "SpikeSlabPrior"     "SpikeSlabPriorBase" "Prior"             

