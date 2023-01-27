############################ This is the same as the results for the file, Causal_Impact_sellers_2013-2017.R"

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

## Load data
data_purchase<- read.csv("E:/My Drive/Research/Collaboration/With_SDowning/BSTS/R_scripts/data_purchase.csv")   ## purchase data
data_mem <- read.csv("E:/My Drive/Research/Collaboration/With_SDowning/BSTS/R_scripts/mem_clean.csv")           ## member info

## Limit the data_purchase (i.e. purchase data) period from 2013 to 2017
data_purchase <-data_purchase %>% filter(order_date > "2012-12-31")
data_purchase <-data_purchase %>% filter(order_date < "2017-12-31")

############# Basic information of the data
## count the number of sellers in the purchase data, 485 
count(distinct(data_purchase, seller_mem_no))

## count the number of buyers in the purchase data,118,504 
count(distinct(data_purchase, mem_no))

################################################################################
############# Data preparation for sellers analysis#########################

list_seller <- unique(data_purchase$seller_mem_no) #485

## create the data of sellers whose member information is available, 334 
data_mem_seller <- data_mem[data_mem$mem_no %in% list_seller, ]

## keep the sellers who registered before 2016-01-25, 315
data_mem_seller <- data_mem_seller %>% filter(reg_date < "2016-01-25")
data_mem_seller_in <- data_mem_seller[data_mem_seller$in_nw=='1',]
data_mem_seller_out <- data_mem_seller[data_mem_seller$in_nw=='0',]

list_in <- unique(data_mem_seller_in$mem_no)  ##52
list_out <- unique(data_mem_seller_out$mem_no)  ## 263

## purchase data limited to sellers whose member information is avaiable.
data_in <- data_purchase[data_purchase$seller_mem_no %in% list_in,]
data_out <- data_purchase[data_purchase$seller_mem_no %in% list_out,]

## Aggregate purchase data to daily, US$1000
p1 <- data_in %>% group_by(order_date)%>%summarise(sales_d = sum(rev_krw_sum)/1000000)
p1$order_date <- as.Date(p1$order_date)

p0 <- data_out %>% group_by(order_date)%>%summarise(sales_d = sum(rev_krw_sum)/1000000)
p0$order_date <- as.Date(p0$order_date)

## Convert daily data to weekly
p1$week <- floor_date(p1$order_date, "week")
p1_w <- p1 %>% group_by(week)%>%summarise(sales_w=sum(sales_d))

p0$week <- floor_date(p0$order_date, "week")
p0_w <- p0 %>% group_by(week)%>%summarise(sales_w=sum(sales_d))

## Convert the dataframe to xts time series
p1_w_ts <- xts(p1_w$sales_w, p1_w$week)
p0_w_ts <- xts(p0_w$sales_w, p0_w$week)

events <- xts(letters[1],as.Date(c("2016-01-24")))
dev.off()

plot.xts(p1_w_ts, major.ticks="quarters", ylim=c(0,600), ylab="US$1000 (weekly)")
lines(p0_w_ts, main="Sellers who did not join the network", col="blue")
addEventLines(events, main="Intervention", col="red")
addLegend('topleft', 
       legend.names = c("Join", "Not join"), 
       col = c("black", "blue"),
       lty = c(1, 1))


################################################################################
############ Overall platform analysis: Custom BSTS Model for Causal Inference #

### Scenario 1: NO counterfactual, build a BSTS model for the counterfactual
## Pre-intervention: 1-160, post-intervention 161-261
## causal impact analysis
# pre.period <- as.Date(c("2012-12-30", "2016-01-17"))
# post.period <- as.Date(c("2016-01-24", "2017-12-30"))

################################################################################
################### With default bsts model using a covariate

## Define pre- and post-intervention periods
pre.period <- as.Date(c("2012-12-30", "2016-01-17"))
post.period <- as.Date(c("2016-01-24", "2017-12-24"))

data <- cbind(p1_w_ts, p0_w_ts)

## CausalImpact analysis using the default model in the package
impact0 <- CausalImpact(data, pre.period, post.period, alpha=0.05, 
                         model.args=list(nseasons = 52,
                                         season.duration = 1,niter = 5000,
                                         standardize.data=FALSE))

plot(impact0)
impact0$model$bsts.model$prior
impact0$model$bsts.model$state.specification

################# With custom bsts model 
y.post <- p1_w_ts
y.post[161:261] <- NA
x1 <- p0_w_ts


post.period.response <- p1_w$sales_w[161:261]
pre.period.response <- p1_w$sales_w[1:160]
sdy=sd(pre.period.response)   #54.156
mu=mean(pre.period.response)

sigma.guess=sdy*0.001

## bsts option: prior -- this is only used when there is no regressor
# sigma_epsilon = c(1,0.01) #Sets average variance and dof hyperparameter for observed
# priorobs          =  SdPrior(sigma_epsilon[1]*sdy, sample.size =sigma_epsilon[2] ,fixed = FALSE, upper.limit =sdy)

#build a bsts model
#Initiate empty state space configuration list
st.sp <- list()

## Step (2): Add a state-space: local trend
st.sp <- AddLocalLevel(st.sp, y.post, 
                       sigma.prior =SdPrior(sigma.guess=sigma.guess, 
                                            sample.size = 32, 
                                            initial.value = sigma.guess,
                                            fixed = FALSE, upper.limit = sdy),
                       
                       initial.state.prior =NormalPrior(mu=p1_w$sales_w[1],
                                                        sigma=sdy,
                                                        initial.value=p1_w$sales_w[1],
                                                        fixed=FALSE))

## Step (3): Add a state-space: seasonality (weekly)
# Saturday & Sunday are the lowest
st.sp <- AddSeasonal(st.sp, y.post, nseasons =52, season.duration = 1, 
                     sigma.prior =SdPrior(sigma.guess=sigma.guess, 
                                          sample.size = 0.01, 
                                          initial.value = sigma.guess,
                                          fixed = FALSE, upper.limit = sdy),
                     
                     initial.state.prior =NormalPrior(mu=0,
                                                      sigma=sdy,
                                                      initial.value=0,
                                                      fixed=FALSE))

## Step (3): Add state-space: named holidays: This causes more problems....
## state-space: AddRegressionHoliday, AddHierarchicalRegressionHoliday,  AddRandomwalkHoliday
## The first two need at least three types of holidays
## state models: Namedholiday, DateRangeHoliday, 
## It is better to have a covariate which captures the holiday effects

######### Add specific holidays
# Shopping increases 10-14 days before a holiday
# Chinese new year and Chinese thanksgiving holidays according to lunar calendar 
# are the main holidays which influence the platform
# Monday-Friday
# Lunar New Year (lny): 2013-02-10, 2014-01-31, 2015-02-19, 2016-02-08, 2017-01-28, 2018-02-16
# Lunar Mid-Autumn (lma): 2013-09-19, 2014-09-08, 2015-09-27, 2016-09-15, 2017-10-04, 2018-09-24 

### Or another way to add specific holidays
holiday.name = "main_holiday_lny_lma"
start.date = as.Date (c("2013-01-20", "2014-01-12", "2015-02-01", "2016-01-17", "2017-01-08",
                        "2013-09-01", "2014-08-24", "2015-09-06", "2016-08-28", "2017-09-17"))

end.date = as.Date (c("2013-01-27", "2014-01-19", "2015-02-08", "2016-01-24", "2017-01-15",
                      "2013-09-08", "2014-08-31", "2015-09-13", "2016-09-04", "2017-09-24"))

holidays=DateRangeHoliday(holiday.name, start.date, end.date)
st.sp <- AddRegressionHoliday(st.sp, y.post, holiday.list=list(holidays),
                              prior=NormalPrior(mu=0,sigma=sdy,
                                                initial.value=0,
                                                fixed=FALSE))


################################################################################
## Custom bsts
## bsts2: local trend + seasonality 
bsts2 <- bsts(y.post ~ x1, 
              state.specification = st.sp, 
             # prior=priorobs,
             # expected.model.size=1, 
              niter = 10000)
             # data=p_c_all)
plot(bsts2)
plot.bsts(bsts2, "component")

## causalimpact analysis
impact1 <-CausalImpact(bsts.model=bsts2,
                      post.period.response=post.period.response,
                      alpha=0.05, model.args=list(niter = 10000))

plot(impact1)


############################ This is the same as the results for the file, Causal_Impact_sellers_2013-2017.R"

############### Individual seller

s13474 <- data_in[data_in$seller_mem_no=='13474',]

p1 <- s13474 %>% group_by(order_date)%>%summarise(sales_d = sum(rev_krw_sum)/1000000)
p1$order_date <- as.Date(p1$order_date)

## Convert daily data to weekly
p1$week <- floor_date(p1$order_date, "week")
p1_w <- p1 %>% group_by(week)%>%summarise(sales_w=sum(sales_d))


## Convert the dataframe to xts time series
p1_w_ts <- xts(p1_w$sales_w, p1_w$week)
plot.xts(p1_w_ts)
