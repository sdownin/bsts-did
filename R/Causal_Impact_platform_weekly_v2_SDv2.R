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

############   Overall platform causal effect analysis  ################# 
## Load data
# data_0<- read.csv("E:/My Drive/Research/Collaboration/With_SDowning/BSTS/R_scripts/purchase_hw_clean_v2.csv")
data_0<- read.csv("C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\BSTS\\purchase_hw_clean_v2.csv")

## Aggregate purchase data to daily
p <- data_0 %>% group_by(order_date)%>%summarise(sales_m = sum(rev_krw_sum)/1000000)
p$order_date <- as.Date(p$order_date)

## Convert daily data to weekly
p$week <- floor_date(p$order_date, "week")
p_w <- p %>% group_by(week)%>%summarise(sales_w=sum(sales_m))

## Convert the dataframe to time series and cut to 2018
p_w_ts <- xts(p_w$sales_w, p_w$week)
p_w_ts <- window(p_w_ts, end="2018-12-30")

plot.xts(p_w_ts, major.ticks="weeks", main="Cconma weekly sales", type="b")

### Scenario 1: no counterfactual, build a BSTS model for the counterfactual
## Preintervention: 1-160, post-intervention 161-314
## causal impact analysis
# pre.period <- as.Date(c("2012-12-30", "2016-01-17"))
# post.period <- as.Date(c("2016-01-24", "2018-12-30"))

## Step (1): 
y.post <- p_w_ts
y.post[161:314] <- NA
post.period.response <- p_w$sales_w[161:314]
pre.period.response <- p_w$sales_w[1:160]
y0 =(p_w$sales_w[1]-mean(pre.period.response, na.rm=T))/sd(pre.period.response, na.rm = T)

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


#build a bsts model
#Initiate empty state space configuration list
st.sp <- list()

## Step (2): Add a state-space: local trend
st.sp <- AddLocalLevel(st.sp, y.post, sigma.prior =levelsdprior, initial.state.prior =initialvalueprior)

## run bsts model based on state.specification found using pre-intervention data
## bsts1: local trend only

bsts1 <- bsts(y.post, state.specification = st.sp,
              niter = 5000)
plot(bsts1)
plot(bsts1, 'components')

## Step (3): Add a state-space: seasonality (weekly)
# Saturday & Sunday are the lowest
st.sp <- AddSeasonal(st.sp, y.post, nseasons =52, season.duration = 1, 
                     sigma.prior=seasonsdprior,initial.state.prior=initialseasonprior)
                    

## bsts2: local trend + seasonality
bsts2 <- bsts(y.post, state.specification = st.sp,
              niter = 5000)
plot(bsts2)
plot(bsts2, 'components') 

## Step (3): Add state-space: named holidays
######### Add specific holidays
# Shopping increases 10-14 days before a holiday
# Chinese new year and Chinese thanksgiving holidays according to lunar calendar 
# are the main holidays which influence the platform
# Monday-Friday
# Chinese new year: 2013-02-10, 2014-01-31, 2015-02-19, 2016-02-08, 2017-01-28, 2018-02-16
# Mid-Autumn holiday: 2013-09-19, 2014-09-08, 2015-09-27, 2016-09-15, 2017-10-04, 2018-09-24 

## Add Chinese new years and Mid-autumn holiday
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

#### end of adding Chinese new year and thanksgiving holidays

## bsts3: local trend + seasonality + named holidays
bsts3 <- bsts(y.post, state.specification = st.sp, prior=priorobs, ping=0,
              niter = 5000)
plot(bsts3)

## causalimpact analysis
impact1 <-CausalImpact(bsts.model=bsts3,
                      post.period.response=post.period.response,
                      alpha=0.05, model.args=list(niter = 10000))

plot(impact1)

plot(bsts3)
plot.bsts(bsts3, "component")
summary(impact1)
summary(impact1,"report")

##############
### Below is the different way to add specific holidays. But, in this case, I couldn't add to causalimpact function.


# # ## MCMC convergence diagnostics: This does not work in this format
# r <- residuals(bsts3)
# par(mfrow = c(1,2))  ## show the plot split into two
# qqdist(r)   ## A bit of departure in the upper tail
# AcfDist(r)   ## posterior distribution of the autocorrelation function using a set of side-by-side boxplots.
# summary(bsts3, burn=1000)
# plot.bsts(bsts3, "component")


### Another way to add holidays
## Add ranges of holidays, Chinese 
chinese.new.year <- DateRangeHoliday("ChineseNewYearandfullmoon", 
                     start = as.Date(c("2013-01-27",
                                       "2014-01-19",
                                       "2015-02-01",
                                       "2016-01-24",
                                       "2017-01-15",
                                       "2018-02-04")
                                       ),
                                       
                     end=as.Date(c("2013-01-27",
                                   "2014-01-19",
                                   "2015-02-08",
                                    "2016-01-31",
                                    "2017-01-15",
                                    "2018-02-04")
                                  ))

st.sp <- AddRandomWalkHoliday(st.sp,
                              y.post,
                              chinese.new.year,
                              time0 = NULL, 
                              sigma.prior = NULL,
                              initial.state.prior = NULL,
                              sdy = sd(as.numeric(y.post), na.rm = TRUE))

## bsts4: local trend + seasonal + Chinese new year (CNY) (random walk)
bsts4 <- bsts(y.post,  ## This specifies including regression (all columns in pre)http://127.0.0.1:16997/graphics/plot_zoom_png?width=1664&height=948
                       state.specification = st.sp,
                       niter = 5000)
plot(bsts4)

### use range holidays for Chinese Thanksgiving (CTG) holidays
chinese.fullmoon <- DateRangeHoliday("ChineseFullmoon", 
                                     start = as.Date(c("2013-09-01",
                                                       "2014-08-24",
                                                       "2015-09-06",
                                                       "2016-08-28",
                                                       "2017-09-17",
                                                       "2018-09-09")),
                                     
                                     end=as.Date(c("2013-09-08",
                                                   "2014-08-31",
                                                   "2015-09-13",
                                                   "2016-09-04",
                                                   "2017-09-24",
                                                   "2018-09-16")))

st.sp <- AddRandomWalkHoliday(st.sp,
                              y.post,
                              chinese.fullmoon,
                              time0 = NULL, 
                              sigma.prior = NULL,
                              initial.state.prior = NULL,
                              sdy = sd(as.numeric(y.post), na.rm = TRUE))

## bsts5: local trend + seasonal + Chinese new year (random walk) + Chinese thanksgiving (randomwalk)
bsts5 <- bsts(y.post,  
              state.specification = st.sp,
              niter = 5000)
plot(bsts5)

## causal impact analysis
impact2 <-CausalImpact(bsts.model=bsts5,
                      post.period.response=post.period.response,
                      alpha=0.05, model.args=list(niter = 5000))

plot(impact2)

#### other things to check
## compare multiple BSTS models
CompareBstsModels(list(trend=bsts1,
                       "trend and seasonal" = bsts2,
                       "trend, seasonal, and named holidays"=bsts3))
                      # "trend, seasonal, and CNY (randomwalk)"=bsts4,
                       #"trend, seasonal, CNY, and CTG (randomwalk)"= bsts5))





