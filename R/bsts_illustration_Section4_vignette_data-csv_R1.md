---
title: "A Bayesian Counterfactual Approach to Dynamic Causal Inference: R Code and Tutorial"
author: 
  - "Author 1^[Information redacted for review]"
  - "Author 2^[Information redacted for review]"
email: "<email redacted for review>"
date: "January 28, 2023"
output:
  pdf_document:
    keep_md: true
    toc: true
---



\newpage

This vignette presents `R` code and a tutorial example of  Bayesian Counterfactual 
Analysis (BCA). This approach includes two main steps: (1) building a Bayesian 
structural time series (BSTS) model in [Section 1](#Section1), and (2) using the BSTS model as a
counterfactual for causal inference in [Section 2](#Section2).

* [**Section 1**](#Section1) covers the complexity of the Bayesian approach. The basic structure of a BSTS model is local trend (state space) + seasonality (state space) + regression. Each component has several ways to be defined, and prior settings are slightly different. Thus, we aim to show a step-by-step process of adding different components in [Section 1](#Section1) accompanied by checking resulting BSTS model diagnostics, setting up Bayesian priors, and adding regression component (Spike and slab priors included). This part requires the `bsts` package in R: [https://cran.r-project.org/web/packages/bsts/bsts.pdf](https://cran.r-project.org/web/packages/bsts/bsts.pdf). 
* [**Section 2**](#Section2) shows how to assess causal effects using the BSTS model from [Section 1](#Section1) as a counterfactual, and is a straightforward process using the `CausalImpact` R package: [https://google.github.io/CausalImpact/CausalImpact.html](https://google.github.io/CausalImpact/CausalImpact.html)  


# 0. Simulated Data

First, load the data simulated from a data generating process (DGP) with the 
following features: 

- local trend involving a random walk of small step sizes (i.e., Gaussian noise standard deviation `sd=0.01` in the local level stochastic process), 
- weekly seasonality with a frequency of one per cycle (i.e., 52 weeks = 1 year), and 
- positive linear trend (i.e., upward linear trajectory over time, on average)

This data set includes the outcome series `y_observed`, and nine 
covariate series, `cov1-cov9`, from which the BSTS models that have regression
components will select their predictors. These covariates were created to have
some type of relevant structure (e.g., noisy linear trend or moderate correlation
with lagged mean of the outcome series) to add to the BSTS model via regression.

Note that the width limit prevents all columns from displaying, but the data 
`tibble` object summary displays the total number of rows and truncated
series with their corresponding data classes.

```r
## load BSTS library
library(bsts)
library(CausalImpact)
library(tibble) ## load a data table class with convenience functions
## LOAD FROM FILE
filename <-'data_bsts_vignette_S4__sim-sd0.01_quadratic_agg_bsts_df1-NOctrl.csv'
df1 <- as_tibble(read.csv(filename)) ##convert table to print labels
print(df1)
```

```
## # A tibble: 520 x 10
##    y_observed    cov1      cov2    cov3   cov4   cov5    cov6     cov7     cov8
##         <dbl>   <dbl>     <dbl>   <dbl>  <dbl>  <dbl>   <dbl>    <dbl>    <dbl>
##  1     0.0859 0.00199 -0.00108   0.0198 0.0490 0.0527 0.00989 -0.207    0.229  
##  2     0.276  0.0156  -0.000129  0.0343 0.0497 0.0484 0.0107   0.126    0.0874 
##  3    -0.127  0.0331   0.0118    0.0336 0.0468 0.0496 0.00943 -0.342    0.0130 
##  4     0.757  0.0343   0.00424  -0.0229 0.0500 0.0446 0.00992 -0.286   -0.00126
##  5     0.695  0.0494   0.00786   0.145  0.0537 0.0506 0.00913  0.153   -0.359  
##  6     1.45   0.0719   0.00344   0.149  0.0434 0.0492 0.0107   0.0614  -0.285  
##  7     1.54   0.0770   0.00629   0.282  0.0520 0.0518 0.00991 -0.00741 -0.268  
##  8     1.75   0.0823   0.00924   0.329  0.0488 0.0565 0.00996 -0.120    0.415  
##  9     2.43   0.0824  -0.00421   0.343  0.0496 0.0501 0.00954 -0.216   -0.321  
## 10     1.92   0.0962  -0.00722   0.488  0.0473 0.0488 0.00981  0.291   -0.203  
## # ... with 510 more rows, and 1 more variable: cov9 <dbl>
```
We can plot this data to visualize the post-intervention change in the 
outcome series `y_observed`.

```r
npds <- nrow(df1)
intpd <- round(npds * 0.6)  # time period to split series (60% pre-, 40% post-)
## PLOT
plot(df1$y_observed, main='Simulated Outcome Time Series (Y)',ylab='Y',xlab='Time')
abline(v=intpd, lty=2)
legend('topleft',legend=c('observed','intervention'), pch=c(1,NA),lty=c(NA,2))
```


\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/outcome_plot-1} 
\newline
The question of interest for research and practice is whether there
is a significant effect caused by the intervention at time point
`t=312`. We will apply BSTS for counterfactual causal inference to 
answer this question in a manner robust to the shape of the onset and 
decay structure of the dynamic treatment effect in the post-intervention window. 

For this purpose, we first build a BSTS model based on observed data before the intervention, 
and the prediction of the BSTS model after the intervention will serve as a counterfactual to assess the causal effects of the intervention. 

Thus, we create the pre-intervention outcome series `y.pre.treat.NAs.post`
to pass into the `bsts()` function. 
Here we replace the post-intervention periods with `NA`'s that 
represent missing values within the `R` language. This is because the 
BSTS model is trained on the pre-intervention outcome data. 
We also define covariates as predictors to add a regression component
to the BSTS model ([Section 1.4](#Section1-4)). 

```r
## Indices of pre-intervention and post-intervention windows
pre.idx <- 1:(intpd-1)
post.idx <- intpd:npds
## Outcome (response) series
y <- df1$y_observed
## Train on y pre-treatment but NA's post-treatment
y.pre.treat.NAs.post <- c( y[pre.idx], rep(NA,length(post.idx)) )
## Then use the post-treatment response for causal impact estimation
post.period.response <- y[post.idx]
## Covariates (predictors)  data
predictors <- df1[ , ! names(df1) %in% 'y_observed']
```


# 1. Bayesian Structural Time Series (BSTS) Modeling {#Section1}

## 1.1. State Space Specification: Local Trend Only

The simplest BSTS model has only the local trend in the state space. 
There are several state spaces in the `bsts` package for defining local trend including local level, 
local linear trend, student local linear trend, and semilocal 
(or generalized local) linear trend.

In this section, we start with the simplest local trend, local level 
(see [Section 1.6](#Section1-6) for other types of local trend).

First, configure the model's state space by adding components to the 
state space list `st.sp`. Here `y.pre.treat.NAs.post` is passed into the state
space functions (after the state space list object `st.sp`)
as the outcome series being modeled. 

```r
## Initiate empty state space configuration list
st.sp <- list()
## Add local level to trend component of state space
st.sp <- AddLocalLevel(st.sp, y.pre.treat.NAs.post)
## Set BSTS MCMC iterations
bsts.niter <- 200 ## suggest on the order of 1k-10k at least
```
Then, fit the BSTS model by calling the `bsts()` function to run MCMC estimation 
for `bsts.niter` draws from the stationary distribution of the Markov chain, which  
(assuming convergence was reached) specifies the posterior predictive distribution.
This will be addressed below in [Section 1.6](#Section1-6)).

```r
## Fit BSTS model
bsts.fit <- bsts(y.pre.treat.NAs.post,
                 state.specification = st.sp,
                 niter = bsts.niter,
                 seed = rand.seed)
```

```
## =-=-=-=-= Iteration 0 Sat Jan 28 14:30:24 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 20 Sat Jan 28 14:30:24 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 40 Sat Jan 28 14:30:24 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 60 Sat Jan 28 14:30:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 80 Sat Jan 28 14:30:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 100 Sat Jan 28 14:30:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 120 Sat Jan 28 14:30:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 140 Sat Jan 28 14:30:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 160 Sat Jan 28 14:30:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 180 Sat Jan 28 14:30:25 2023
##  =-=-=-=-=
```
Visualizing the model's behavior can help build intuition about 
what the model is doing and, ideally, about how to improve its fit to the 
characteristics of the data set. To inspect the contributions of the model 
components, we write a function that will report the model's summary, 
plot the state space components, and if there is a regression in the model, 
it will plot the coefficients, predictors, and distribution of model sizes.

```r
## Define a function to summarize and visualize the BSTS model components 
getBstsSummaryPlots <- function(bsts.fit) {
  par(mfrow=c(1,1)) ## set plotting params to one image in plot window
  print(summary(bsts.fit))
  plot(bsts.fit, 'components')
  if (bsts.fit$has.regression) {
    plot(bsts.fit, 'coefficients', main='Coefficients')
    plot(bsts.fit, 'size', main='Model Size Distribution')
    plot(bsts.fit, 'predictors', main='Predictors')
  }
  plot(bsts.fit, 'state', 
       main='Observations on Fitted State Space with Forecast')
}
## Call our newly defined function
getBstsSummaryPlots(bsts.fit)
```

```
## $residual.sd
## [1] 0.3768937
## 
## $prediction.sd
## [1] 0.5398463
## 
## $rsquare
## [1] 0.9411521
## 
## $relative.gof
## [1] NA
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-5-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-5-2} \end{center}
There appears to be a seasonal component to the time series, which is currently
being captured by the trend. However, this would be a poor model to use as a 
counterfactual because it does not capture how seasonality continues
after the training window. This is visible in the lack of seasonal pattern and
quickly widening forecast intervals in the plot of the trend component 
of the state space after the intervention period, t=312. 

Next, we plot the contributions to the model's state space against
the observed outcome (black dots in figure below) to help with visual 
understanding of this model's prediction, or fitted values (blue lines), 
and later for comparison with other models introduced below. This uses 
a custom function `plotBstsStateComps()` we created for this purpose, which
can be found in the `R/bsts_helper_functions.R` script in the
online repository. 

```r
check <- plotBstsStateComps(bsts.fit, intpd=intpd, return.val = T) 
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-6-1} \end{center}
This model appears to match the data with the trend only. 
However, this is merely accurate for short (e.g., one-step ahead) predictions 
within the training window. In terms of capturing the structure of the time 
series, this model is underspecified--while at the same time, 
it is essentially overfitted to the training data.  

Although the predicted values (blue lines) seem accurate as one-step ahead, 
in-sample predictions, the model would lose the seasonal information 
outside the training window (post intervention), leaving the counterfactual series
devoid of seasonality when used for causal inference computations
([Section 2](#Section2)).

Therefore, it is crucial to capture seasonality in the outcome time series 
and incorporate it into the BSTS model that will be used as the 
counterfactual for causal inference with the `CasualImpact()` function in
[Section 2](#Section2). 

## 1.2. Seasonality:  State Space = Local Trend + Seasonality

To account for seasonality in the time series data generating process, 
we add a seasonal component to the state space defined by the number of seasons and the duration of each season. For instance, we use 52 periods with duration =1 to produce the state contribution of seasonal cycles (e.g., weekly periods each year). As a comparison, for daily data, 
it is recommended to use 7 for a 
day-of-week component with the duration of each season of 1.

```r
## Initiate empty state space configuration list
st.sp2 <- list()
## Add local level to trend component of state space
st.sp2 <- AddLocalLevel(st.sp2, y.pre.treat.NAs.post)
## Add seasonality to state space (52 weekly periods = 1 yearly cycle) 
st.sp2 <- AddSeasonal(st.sp2, y.pre.treat.NAs.post, nseasons = 52)
## Fit BSTS model and get summary plots
bsts.fit2 <- bsts(y.pre.treat.NAs.post ,
                 state.specification = st.sp2, 
                 niter = bsts.niter,
                 seed = rand.seed, ping=0)
getBstsSummaryPlots(bsts.fit2)
```

```
## $residual.sd
## [1] 0.4995245
## 
## $prediction.sd
## [1] 0.581718
## 
## $rsquare
## [1] 0.8966271
## 
## $relative.gof
## [1] NA
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-7-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-7-2} \end{center}

```r
check2 <- plotBstsStateComps(bsts.fit2, intpd=intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-7-3} \end{center}
Following BSTS (and other forecasting) conventions, 
it is common to use a model performance metric, 
such as the one-step ahead prediction (forecast) error within sample
during the pre-intervention (or training) window. We might expect that the augmented
state space should improve the model performance (e.g., in terms of 
decreased MAE). However, the MAE of this model is now 0.569.

Here we have a seasonal component separate from the trend (which includes 
only the local level at this point). Now the state space of the fitted model 
exhibits the seasonal structure of the time series into the post-intervention 
forecast window. 

However, besides seasonality, the observed data (Y) also appears to have
a positive linear trend, which is not currently captured by the 
local level (default without slope or drift) of the state space trend. 
This is apparent in trend component plot where the upward linear 
trajectory before the intervention turns into a constant expectation after the 
intervention -- as a simple random walk where the expectation at `t` is 
the trend at `t-1` (also just the local level without slope or drift) 
plus Gaussian noise. Therefore, this seasonal model still does not capture the linear trend in the 
DGP and the fit could also be improved with more information about the outcome from 
covariate series with regression. But before moving on to addressing 
regression in BSTS models, we first detail the process of parameter tuning. 

A crucial step in any machine learning (ML) analysis is hyperparameter tuning,
which applies to our investigation of dynamic causal inference using a Bayesian
counterfactual approach drawing upon the ML tradition. Therefore,
adjusting the values of prior distributions in BCA plays an important role in 
achieving adequate model fit for complex time series structures. 
The BSTS model should, of course, be able to mimic the observed outcome series 
pre-intervention and project that state space forward into the post-intervention
window while also characterizing the uncertainty via credible intervals. 

Additionally, in terms of research design and process, 
prior parameter tuning represents an important step for researchers 
to incorporate information from theory and context into their empirical
analyses--coherently within a Bayesian framework. 

## 1.3. Bayesian Priors:  Parameter Tuning

Setting up Bayesian priors is notorious for its complication. 
The conventional wisdom is to exploit information of the observed outcome. 
There are two main types of priors to set up for local level and seasonality,  priors for the variance of the state innovation errors and the initial value of the state at time 0, provided by Boom package.  The former is defined by the option of `sigma.prior` created by `SdPrior` which describes the prior distribution for the standard deviation of the random walk increment (e.g., inverse Gamma priors). The latter one is `initial.state.prior` defined by `NormalPrior` to describe the prior distribution of the initial state vector (e.g., Gaussian distribution). If not defined, the default setting will be applied.
Examples include:
```
sigma.prior = SdPrior(sigma.guess= , 
                      sample.size = , 
                      initial.value = ,
                      fixed = , 
                      upper.limit = ),
                       
initial.state.prior = NormalPrior(mu= ,
                                  sigma= ,
                                  initial.value= ,
                                  fixed= ,
                                  upper.limit= )
```
Most of the particulars of prior elicitation we necessarily bound outside the 
scope of this tutorial vignette. Beginning at the assumption that the 
researcher has theoretical and/or contextual knowledge that can inform  
prior distribution attributes, we focus next on illustrating how to adjust 
the prior distribution settings for the `bsts()` function. 

In order to simplify the logic for specific values, we set parameter settings
for prior variances in terms of multiples of the empirical 
standard deviation of the observed outcome in the pre-intervention
window, `y[pre.idx]` (which is the same as `y.pre.treat.NAs.post`
within the `pre.idx` interval). 

```r
## Get outcome (Y) standard deviation to use for adjusting priors below
y.sd <- sd(y[pre.idx], na.rm = T)
## Initiate empty state space configuration list
st.sp3 <- list()
## Add local level to trend component of state space
st.sp3 <- AddLocalLevel(st.sp3, y.pre.treat.NAs.post, 
  initial.state.prior = NormalPrior(mu=1, sigma = 0.01 * y.sd, fixed = FALSE),
  initial.y = y.pre.treat.NAs.post[1]
)
## Add seasonality to state space (52 weekly periods = 1 yearly cycle) 
st.sp3 <- AddSeasonal(st.sp3, y.pre.treat.NAs.post, 
  nseasons = 52, 
  season.duration = 1, 
  sigma.prior = SdPrior(
    sigma.guess = 0.01 * y.sd,      ##try a low weight of observed Y SD
    sample.size = round(0.1 * (intpd-1)), ## portion of sample size (pre-int.)
    upper.limit= 1.5 * y.sd ## use SD of observed Y to set upper limit
  ),
  initial.state.prior = NormalPrior(
    mu = 0,   
    sigma = 1.0 * y.sd  ## increase initial state prior sigma
  )
)
## Fit BSTS model:
bsts.fit3 <- bsts(y.pre.treat.NAs.post,
                 state.specification = st.sp3,
                 niter = bsts.niter,
                 seed = rand.seed, ping=0)
getBstsSummaryPlots(bsts.fit3)
```

```
## $residual.sd
## [1] 0.4976422
## 
## $prediction.sd
## [1] 0.5903305
## 
## $rsquare
## [1] 0.8974047
## 
## $relative.gof
## [1] NA
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-8-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-8-2} \end{center}

```r
check3 <- plotBstsStateComps(bsts.fit3, intpd=intpd, return.val = T) 
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-8-3} \end{center}
The MAE of this model is now 0.567.

This model does not appear to improve (in terms of MAE),
but the model now has seasonal structure, in addition to the local level
captured in the trend component, which collectively would inform a better
counterfactual post-intervention. So, essentially, the addition of seasonal
structure reduces a portion of the overfitting against current data
in order to incorporate a temporal structure (seasonality) 
that can continue into the post-intervention forecast.
This, in turn, will be used as the counterfactual for causal inference. 

Beyond seasonality, another way to incorporate structure from 
the training window into the post-intervention window
is including a regression component in BSTS model. Assuming the 
covariates included in the regression carry at least some information
about the outcome series before the intervention (and that relationship
does not change for the covaritate after the intervention), then including
regression in the BSTS model should generally improve the model's performance 
as a counterfactual (i.e., as a representation of the pre-intervention outcome 
time series structure projected forward into the post-intervention window).
In using BCA for causal inference, the researcher can include 
covariates with information about (i.e., correlation with) the outcome 
of interest and therefore improve the representatives of the 
counterfactual series.

## 1.4. Regression:  State Space = Local Trend + Seasonality + Regression {#Section1-4}

### 1.4.1. Regression `formula` in BSTS function

This adds a static regression component (i.e., time-invariant betas) with
covariate series in the `predictors` dataframe, which is used as
the `data` input for the `bsts()` function. When we do not specify this, 
the default setting is 1, which entails the model choosing the one
covariate that best explains the outcome, Y.  

In this model we set `expected.model.size=3` to allow the model to 
incorporate information from, on average, three covariate series in the
regression component. The default setting for this parameter is 
`expected.model.size=1`.


```r
## Initiate empty state space configuration list
st.sp4 <- list()
## Add local level to trend component of state space
st.sp4 <- AddLocalLevel(st.sp4, y.pre.treat.NAs.post,
  initial.state.prior = NormalPrior(mu=1, sigma = .01 * y.sd, fixed = FALSE),
  initial.y = y.pre.treat.NAs.post[1]
)
## Add seasonality to state space (52 weekly periods = 1 yearly cycle) 
st.sp4 <- AddSeasonal(st.sp4, y.pre.treat.NAs.post, 
  nseasons = 52, 
  season.duration = 1, 
  sigma.prior = SdPrior(
    sigma.guess = 0.01 * y.sd,      ## use SD of observed y
    sample.size = round(0.1 * (intpd-1)),  ## use portion of pre-int window
    upper.limit= 1.5 * y.sd ## use SD of observed y to set upper limit
  ),
  initial.state.prior = NormalPrior(
    mu= 0,  
    sigma= 1.0 * y.sd  
  )
)
## Fit BSTS model
bsts.fit4 <- bsts(y.pre.treat.NAs.post ~ . ,
                 state.specification = st.sp4,
                 data = predictors,  ## covariates
                 niter = bsts.niter, 
                 expected.model.size=3, ##expected number of covariates in model
                 seed = rand.seed, ping=0)
getBstsSummaryPlots(bsts.fit4)
```

```
## $residual.sd
## [1] 0.4991445
## 
## $prediction.sd
## [1] 0.5838873
## 
## $rsquare
## [1] 0.8967843
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  1.0000  1.0000  0.8982  1.0000  2.0000 
## 
## $coefficients
##                      mean           sd    mean.inc      sd.inc    inc.prob
## cov5         2.2931603036  4.142207451  8.14803767  3.62970151 0.281437126
## cov4         1.6665100232  3.416565314  7.32387300  3.13316944 0.227544910
## cov6         7.4448125579 19.815602382 42.87185163 27.47426336 0.173652695
## cov3         0.0348657237  0.115328734  0.30645136  0.18657603 0.113772455
## cov2        -0.0207215441  0.358426788 -0.38449976  1.58379003 0.053892216
## cov1         0.0008115065  0.014431605  0.02258693  0.07952607 0.035928144
## cov9        -0.0003017008  0.003898833 -0.05038403  0.00000000 0.005988024
## cov7         0.0005707255  0.007375399  0.09531116  0.00000000 0.005988024
## cov8         0.0000000000  0.000000000  0.00000000  0.00000000 0.000000000
## (Intercept)  0.0000000000  0.000000000  0.00000000  0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-9-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-9-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-9-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-9-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-9-5} \end{center}

```r
check4 <- plotBstsStateComps(bsts.fit4, intpd=intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-9-6} \end{center}
The MAE of this model is now 0.567.

### 1.4.2.  Variable Selection: Spike-and-Slab Priors

Having proper variables in regression is a crucial step in building a counterfactual time series
that sufficiently approximates (i.e., correlates to a sufficient degree with) 
the outcomes series in the pre-intervention window. 

One of the advantages of BSTS models is that variable selection among covariates 
is incorporated into the model within a Bayesian framework via 
spike-and-slab priors. 

There are two main ways to modify the the number of covariate series included
in the BSTS model, which is also referred to as the "model size," or 
size of the regression component of the model:

1. Pass `expected.model.size` argument into `bsts()` directly
2. Create a `SpikeSlabPrior` object to pass into `bsts()` by:
    1. Using `expected.model.size` argument of `SpikeSlabPrior()`
    2. Specifying inclusion probability for each covariate in `SpikeSlabPrior()`

First, the simplest method is to set the `expected.model.size` attribute 
(of the `SpikeSlabPrior`) by passing the `expected.model.size` argument 
into the `bsts()` function directly.

```r
st.sp5 <- st.sp4  ## copy previous state space config to new model
#=================
# Option 1. Specify expected model size (internally creates SpikeSlabPrior)
#----------------
#   This uniformly sets prior (spike) probabilities = (1 / expected.model.size)
bsts.fit5 <- bsts(y.pre.treat.NAs.post ~ .,
                 state.specification = st.sp5,
                 data = predictors,
                 expected.model.size = 6,  ## argument passed to SpikeSlabPrior
                 niter = bsts.niter,
                 seed = rand.seed, ping=0)
summary(bsts.fit5)
```

```
## $residual.sd
## [1] 0.5058727
## 
## $prediction.sd
## [1] 0.5647206
## 
## $rsquare
## [1] 0.893983
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.000   1.000   1.000   1.503   2.000   4.000 
## 
## $coefficients
##                    mean          sd    mean.inc     sd.inc   inc.prob
## cov3         0.52459194  0.42410042  0.71808897  0.3265812 0.73053892
## cov2         0.63995379  1.92175105  2.37493961  3.1175355 0.26946108
## cov1        -0.01104471  0.07750426 -0.06360227  0.1793049 0.17365269
## cov6         4.25873450 15.21192722 37.43203484 28.6935093 0.11377246
## cov4         0.71120179  3.37298566  6.25108940  8.2684509 0.11377246
## cov5         0.43971761  2.65780633  4.89552271  7.7708737 0.08982036
## cov8         0.00305984  0.02941997  0.25549664  0.1211827 0.01197605
## cov9         0.00000000  0.00000000  0.00000000  0.0000000 0.00000000
## cov7         0.00000000  0.00000000  0.00000000  0.0000000 0.00000000
## (Intercept)  0.00000000  0.00000000  0.00000000  0.0000000 0.00000000
```

```r
par(mfrow=c(2,1))
plot(bsts.fit5, 'coefficients', main='Expected Size = 6')
plot(bsts.fit4, 'coefficients', main='Expected Size = 3')
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-10-1} \end{center}
Second, the user can create a `SpikeSlabPrior` object to pass into the `bsts()`
function. One way to do this is by using the `expected.model.size` argument 
of `SpikeSlabPrior()`. This allows you to set the expected model size in a 
simple manner when all covariates are assumed to have equal, naive inclusion
probabilities. This is just like the first method above but the intermediate
step of creating a `SpikeSlabPrior` object enables the researcher to adjust 
the expected model size and all other regression prior parameter settings 
in the same step. 

```r
st.sp6 <- st.sp5  ## copy previous state space config to new model
#=================
# Option 2. Create SpikeSlabPrior object and pass it into BSTS function
#-----------------
###---------------
#  2A. Set expected.model.size in SpikeSlabPrior (with other options)
priorObjA <- SpikeSlabPrior(
  x = model.matrix(y.pre.treat.NAs.post ~ ., data=predictors),
  y = y.pre.treat.NAs.post,
  expected.model.size = 3, ## size in SpikeSlabPrior
  expected.r2 = .9,
  prior.df = .1,## weight in terms of observation count
  prior.information.weight = .1
) 

# Run the bsts model with the SpikeSlabPrior object ‘priorObjA’
bsts.fit6 <- bsts(y.pre.treat.NAs.post ~ ., 
                  state.specification = st.sp6,
                  data = predictors,
                  niter = bsts.niter,
                  prior= priorObjA, ## spike-and-slab priors defined above
                  seed = rand.seed, ping = 0)
getBstsSummaryPlots(bsts.fit6)
```

```
## $residual.sd
## [1] 0.4966362
## 
## $prediction.sd
## [1] 0.5649397
## 
## $rsquare
## [1] 0.8978191
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.519   2.000   5.000 
## 
## $coefficients
##                      mean          sd   mean.inc      sd.inc    inc.prob
## cov3         0.6077633084 0.331656729  0.6748783  0.27684645 0.900552486
## cov1        -0.0380157296 0.141094742 -0.1495836  0.25014258 0.254143646
## cov2         0.4454722951 1.555109728  3.5056733  2.92761489 0.127071823
## cov4         0.3216657789 1.402144377  4.1586790  3.16816989 0.077348066
## cov5         0.1830146101 0.857191510  2.7604704  2.06441777 0.066298343
## cov6         0.3776655725 3.605753802  8.5446836 15.95523295 0.044198895
## cov8         0.0050506581 0.034144330  0.1828338  0.10877323 0.027624309
## cov7        -0.0024266850 0.020703780 -0.1464100  0.08341039 0.016574586
## cov9        -0.0007169068 0.009644995 -0.1297601  0.00000000 0.005524862
## (Intercept)  0.0000000000 0.000000000  0.0000000  0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-11-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-11-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-11-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-11-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-11-5} \end{center}

```r
check6 <- plotBstsStateComps(bsts.fit6, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-11-6} \end{center}
The MAE of this model is now 0.549.

Finally, the other way to create a `SpikeSlabPrior` object to pass
into the `bsts()`function is by specifying prior inclusion probabilities
(and estimates) for each covariates series in `predictors`. This is useful
when the researcher has prior knowledge that certain covariates are 
more likely to be important contributors to the counterfactual series
and should therefore be included in the model in more iterations 
(i.e., more MCMC draws from the posterior predictive distribution). 

```r
st.sp7 <- st.sp6  ## copy previous state space config to new model
###---------------
# 2B. Specify prior.inclusion.probabilities for each covariate where 
# length of prior spikes and means vectors should equal the number of covariates
prior.spikes <- c(.7,.8,.9,.1,.1,.1,.1,.1,.1,.1) ## +1 column for intercept
prior.means <- c(.2,.2,.3,0,0,0,0,0,0,0) ## +1 column for intercept 
# Directly set prior spike probabilities for each covariate in the prior object
priorObjB <- SpikeSlabPrior(x = model.matrix(y.pre.treat.NAs.post ~ ., 
                                             data=predictors),
                            y = y.pre.treat.NAs.post,
                            expected.r2 = .9,
                            prior.df = .1,
                            prior.information.weight = .1, ##  ~ obs. count
                            prior.inclusion.probabilities = prior.spikes,
                            optional.coefficient.estimate = prior.means)
# Run the bsts model with the SpikeSlabPrior object ‘priorObjB’
bsts.fit7 <- bsts(y.pre.treat.NAs.post ~ ., 
                  state.specification = st.sp7,
                  data = predictors,
                  niter = bsts.niter,
                  prior= priorObjB, ## spike-and-slab prior defined above
                  seed = rand.seed, ping = 0)
getBstsSummaryPlots(bsts.fit7)
```

```
## $residual.sd
## [1] 0.5024276
## 
## $prediction.sd
## [1] 0.5703497
## 
## $rsquare
## [1] 0.8954221
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   2.000   2.005   3.000   4.000 
## 
## $coefficients
##                      mean          sd    mean.inc     sd.inc    inc.prob
## cov1         0.3912775424 0.179951921  0.40461655  0.1675259 0.967032967
## cov2         0.9596318979 2.904413440  1.58775459  3.6057335 0.604395604
## cov3         0.3256057355 0.488124900  0.92594131  0.3461832 0.351648352
## cov6         0.5926489797 6.153024610 13.48276429 27.9342410 0.043956044
## cov4         0.3127586955 1.902085012  9.48701376  5.1502845 0.032967033
## cov9        -0.0003870325 0.005221354 -0.07043992  0.0000000 0.005494505
## cov8         0.0000000000 0.000000000  0.00000000  0.0000000 0.000000000
## cov7         0.0000000000 0.000000000  0.00000000  0.0000000 0.000000000
## cov5         0.0000000000 0.000000000  0.00000000  0.0000000 0.000000000
## (Intercept)  0.0000000000 0.000000000  0.00000000  0.0000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-12-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-12-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-12-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-12-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-12-5} \end{center}

```r
check7 <- plotBstsStateComps(bsts.fit7, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-12-6} \end{center}
The MAE of this model is now 0.558.

## 1.5. Local Trends Comparison: Level vs. Linear, Short(er)- vs. Long(er)-Term

Besides the default state space containing only a local level component, 
we can evaluate different state space configurations, such as those
including a local linear trend and semi-local linear trend.

### 1.5.0. Local Level (Trend without Slope)

This was the simple model introduced above. As a state space without
a linear trend, this model will serve as a baseline for comparison 
with other state space configurations that contain a type of linear trend. 
Given that we know our simulation DGP includes a small positive linear trend,
this should accentuate performance differences between different types of
local trends (specifically, those with vs. without a slope) to make evident the 
differences between the types of linear trend components that can be 
added to the BSTS model's state space. 

### 1.5.1. Local Linear Trend (Trend with Level and Slope)

The local linear trend is useful for short term predictions of relatively 
well-behaved time series. 

```r
## Initiate empty state space configuration list
st.sp8 <- list()
## Add local linear trend component of state space
st.sp8 <- AddLocalLinearTrend(st.sp8, y.pre.treat.NAs.post, 
  level.sigma.prior=SdPrior(sigma.guess=.01 * y.sd, 
                            sample.size=round( 0.1 * (intpd-1)), 
                            upper.limit = 1.5 * y.sd,
                            fixed = F), 
  slope.sigma.prior=SdPrior(sigma.guess=.0001 * y.sd, 
                            sample.size=round( 0.1 * (intpd-1)), 
                            upper.limit = 1.5 * y.sd,
                            fixed = F), 
  initial.level.prior=NormalPrior(mu=1, sigma=.001 * y.sd, fixed = F), 
  initial.slope.prior=NormalPrior(mu=0, sigma=.001 * y.sd, fixed = F)
)
## Add seasonality to state space
st.sp8 <- AddSeasonal(st.sp8, y.pre.treat.NAs.post, 
  nseasons = 52, 
  season.duration = 1, 
  sigma.prior = SdPrior(
    sigma.guess = 0.01 * y.sd,  
    sample.size = round( 0.1 * (intpd-1)),    
    upper.limit = 1.5 * y.sd
  ), 
  initial.state.prior = NormalPrior(
    mu = 0,   
    sigma = 1 * y.sd  
  )
)
## Fit BSTS model
bsts.fit8 <- bsts(y.pre.treat.NAs.post ~ . ,
                 state.specification = st.sp8,
                 data = predictors, 
                 niter = bsts.niter, 
                 expected.model.size=3, ## instead of:  prior = priorObjB
                 seed=rand.seed, ping = 0)
getBstsSummaryPlots(bsts.fit8)
```

```
## $residual.sd
## [1] 0.5180204
## 
## $prediction.sd
## [1] 0.5520123
## 
## $rsquare
## [1] 0.8888302
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   2.000   1.591   2.000   3.000 
## 
## $coefficients
##                     mean          sd   mean.inc      sd.inc    inc.prob
## cov3         1.299304388  0.44157176  1.2993044  0.44157176 1.000000000
## cov2         0.857570517  2.05563453  5.1111203  1.84711061 0.167785235
## cov6         5.427879387 14.33667496 35.1632186 17.01788019 0.154362416
## cov5         0.991884190  2.59849354  7.0376545  2.30013119 0.140939597
## cov4         0.624424940  2.26523180  8.4581196  1.74891741 0.073825503
## cov1         0.010175736  0.05156756  0.2526974  0.07191086 0.040268456
## cov9        -0.001203903  0.01469551 -0.1793815  0.00000000 0.006711409
## cov8         0.002250653  0.02747272  0.3353472  0.00000000 0.006711409
## cov7         0.000000000  0.00000000  0.0000000  0.00000000 0.000000000
## (Intercept)  0.000000000  0.00000000  0.0000000  0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-13-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-13-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-13-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-13-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-13-5} \end{center}

```r
check8 <- plotBstsStateComps(bsts.fit8, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-13-6} \end{center}
The MAE of this model is now 0.545.

### 1.5.2. Semilocal Linear Trend (Trend with Level and Slope with AR(1) Drift)

The semilocal linear trend component in the BSTS model state space is 
useful for longer-term forecasts wherein the linear trend persists longer, 
so tighter prediction bounds are reasonable over a longer post-intervention 
window. 

```r
## Initiate empty state space configuration list
st.sp9 <- list()
## Add semilocal linear trend component of state space
st.sp9 <- AddSemilocalLinearTrend(st.sp9, y.pre.treat.NAs.post, 
  level.sigma.prior=SdPrior(
    sigma.guess=.001 * y.sd, 
    sample.size=round( 0.1 * (intpd-1)), 
    upper.limit = 1.5 * y.sd, 
    fixed = F),
  slope.mean.prior=NormalPrior(
    mu = 0, 
    sigma = .01 * y.sd), 
  slope.ar1.prior=Ar1CoefficientPrior(
    mu = 0, 
    sigma = .001 * y.sd, 
    force.stationary = FALSE, 
    force.positive = FALSE),
  slope.sigma.prior=SdPrior(
    sigma.guess=.0001 * y.sd, 
    sample.size=round( 0.1 * (intpd-1))),
  initial.level.prior=NormalPrior(
    mu = 1, ##adjusted
    sigma = .01 * y.sd, ##adjusted
    fixed = F),
  initial.slope.prior=NormalPrior(
    mu = 0, 
    sigma = .001 * y.sd, 
    fixed = F)
)
## Add seasonality to state space
st.sp9 <- AddSeasonal(st.sp9, y.pre.treat.NAs.post, 
  nseasons = 52, 
  season.duration = 1, 
  sigma.prior = SdPrior(
    sigma.guess = .01 * y.sd,   
    sample.size = round( 0.1 * (intpd-1)),    
    upper.limit= 1.5 * y.sd
  ), 
  initial.state.prior = NormalPrior(
    mu= 0,   
    sigma = 1 * y.sd,  
  )
)
## Fit BSTS model
bsts.fit9 <- bsts(y.pre.treat.NAs.post ~ . ,
                 state.specification = st.sp9,
                 data = predictors,  
                 niter = bsts.niter, 
                 expected.model.size=3,
                 seed=rand.seed, ping=0)
getBstsSummaryPlots(bsts.fit9)
```

```
## $residual.sd
## [1] 0.5192231
## 
## $prediction.sd
## [1] 0.5232126
## 
## $rsquare
## [1] 0.8883134
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   2.000   2.000   2.000   2.455   3.000   4.000 
## 
## $coefficients
##                     mean           sd     mean.inc     sd.inc    inc.prob
## cov3        1.269078e+00 0.3591098920  1.269077904 0.35910989 1.000000000
## cov1        9.507625e-01 0.0942565084  0.950762503 0.09425651 1.000000000
## cov4        1.018083e+00 2.1985885807  5.057575679 1.88183390 0.201298701
## cov5        7.691033e-01 2.1548895746  5.640090770 2.58393079 0.136363636
## cov2        2.011587e-01 0.8344036541  3.097843910 1.36578511 0.064935065
## cov6        8.595982e-01 4.2054778398 18.911161112 7.26099871 0.045454545
## cov9        4.434925e-05 0.0005503598  0.006829785 0.00000000 0.006493506
## cov8        0.000000e+00 0.0000000000  0.000000000 0.00000000 0.000000000
## cov7        0.000000e+00 0.0000000000  0.000000000 0.00000000 0.000000000
## (Intercept) 0.000000e+00 0.0000000000  0.000000000 0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-14-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-14-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-14-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-14-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-14-5} \end{center}

```r
check9 <- plotBstsStateComps(bsts.fit9, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-14-6} \end{center}
The MAE of this model is now 0.516.

Zooming into the first 104 periods shows how the BSTS model learns the structure
of the time series through the first couple of cycles that train the 
the seasonality and regression structures into the fitted model (via
recursive least squares in the Kalman filter). 


```r
par(mfrow=c(1,2), mar=c(4,4,.5,1))
plotBstsStateComps(bsts.fit9, intpd=intpd, pd.ids=1:104, title.show=F)##2 cycles
plotBstsStateComps(bsts.fit9, intpd=intpd, pd.ids=1:520, title.show=F)##10 cycles
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-15-1} \end{center}

## 1.6. MCMC Diagnostics: Checking Convergence and Fit {#Section1-6}

Draws from the posterior distribution of the Markov chain(s)--once
converged at the stationary distribution--are used to estimate the 
model prediction at each time period through 
Markov chain Monte Carlo (MCMC) estimation. However, when the Markov chains
have not yet converged (usually due to insufficient MCMC burn-in period
or sampled iterations), then the iterations cannot be treated as draws from
the stationary distribution to sample the posterior predictive distribution. 
In such case, the estimates may be unstable and/or biased, so users must assess 
convergence of the MCMC chains. 

There are several ways to check the convergence 
of MCMC including the trace plots, Heidelberger-Welch tests, and Geweke 
diagnostic tests. Additionally, the researcher may evalute residual
autocorrelation using the autocorrelation function (ACF) plots, among other tests.
Similarly, there is need to assess how the estimated model fits the observed data, 
in terms of how well new samples simulated from the fitted model create 
simulated distribution(s) that resemble the observed value(s). For this purpose, 
we examine MCMC convergence trace figures and plot model fit 
checks for the posterior predictive distribution of the outcome
(Y, top row panels) and the standardized residuals (bottom row panels). Here the
residuals are the 1-step prediction errors during the pre-intervention (training)
window. 

To check the convergence of the model estimation procedure and inspect 
the model's state space component contributions, we wrap a custom 
function created for this analysis, `bstsPostPredChecks()`, 
inside a formatting function, `getBstsDiagnosticPlots()`. Calling this on
several candidate BSTS models `bsts.fit`, `bsts.fit2`,
`bsts.fit8`, `bsts.fit9`, allows for comparison. 

```r
## Define a function to call all of our diagnostics and plotting
getBstsDiagnosticPlots <- function(bsts.fit, return.val=F) {
  par(mfrow=c(1,1))
  checks <- bstsPostPredChecks(bsts.fit, save.plot = F, return.val = T)
  cat('MCMC Convergence Diagnostics:\n')
  cat(checks$summary)
  if(return.val)
    return(checks)
}
```
Following recent empirical progress in applying Bayesian structural time series
for dynamic causal inference, we adapt the diagnostics panel below from 
Menchetti & Bojinov (2022) to cover the following scope of diagnostic checks: 

- Column 1 (left) addresses MCMC convergence along two dimensions 
(posterior predictive distribution chain, and standardized residual chain) 
- Columns 2-3 (middle, right) address model fit or problems, 
conditional on convergence established in the first column. 

Regarding MCMC convergence, the first column shows the trace plots
of the Markov chain for the posterior predictive distribution of the outcome (Y). 
For model fit, the observed and predicted values of Y are compared (top-center), 
where tighter alignment indicates better fit. 

Three MCMC convergence diagnostics are used to assess convergence of the Markov
chains for each of the two traces (posterior predictive distribution and
strandardized residuals). 

1. *Geweke's Convergence Test* (see `coda::geweke.diag()`) assess convergence of the MCMC chain with the Geweke Z-score based on equality of means of the first and last portions of the chain (defaults: 0.1 and 0.5).  If the chain has reached the stationary distribution, the two measures are equal and the low Z-score would not reject the null hypothesis of convergence.
2. *Heidelberg & Welch's Stationary Test* (see `coda::heidel.diag()`) evaluates if the chain samples come from the same stationary distribution according to the Cramer-von Mises (CMV) statistic, for judging the goodness of fit of a cumulative distribution function, under the null hypothesis that the sampled values come from a stationary distribution.
3. *Heidelberg & Welch's Halfwidth Ratio Test* (see `coda::heidel.diag()`) calcuates a 95% interval for the mean of the Markov chain, using the portion of the chain that passed the H&W stationary test (in the CMV criterion above). If the ratio of interval halfwidth-to-mean is not less than some tolerance (default `tol=0.1`), then the sample chain length is not long enough to estimate the mean with sufficient accuracy, suggesting the need for more MCMC iterations. 

Additionally, it is necessary to assess model fit and residual assumptions.
We use the maximum value of Y as an auxiliary statistic for computing 
a Bayesian p-val (Menchetti & Bojinov, 2022), which represents the proportion 
of the maximum Y values (from MCMC iterations) that are larger than the 
observed largest value of Y. Smaller values of Bayesian p-val indicate worse 
fit, which offers an analogous interpretation to the frequentist null 
hypothesis statistical test (NHST) under a null assumption of a well fitting model. 
Values of Bayesian p-val >= alpha (conventionally 0.05) would therefore not reject the 
null of good fit, whereas Bayesian p-val < alpha heuristically suggests
problem(s) of fit to be addressed (e.g., by increasing MCMC iterations, 
or respecifying the state space).
\newpage

```r
## Baseline (Trend = local level)
getBstsDiagnosticPlots(bsts.fit)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-17-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=0.20, p=0.84)
## H&W Stationarity: PASS (CMV=0.05, p=0.86)
## H&W Halfwidth: FAIL (hw/mean=0.17 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=7.39, p=0.00)
## H&W Stationarity: PASS (CMV=11605.93, p=0.71)
## H&W Halfwidth: FAIL (hw/mean=1.30 < eps=0.10)
```
\newpage

```r
## Seasonality Baseline
getBstsDiagnosticPlots(bsts.fit2)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-18-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=1.37, p=0.17)
## H&W Stationarity: PASS (CMV=0.07, p=0.76)
## H&W Halfwidth: PASS (hw/mean=0.03 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=3.15, p=0.00)
## H&W Stationarity: PASS (CMV=0.08, p=0.71)
## H&W Halfwidth: PASS (hw/mean=0.04 < eps=0.10)
```
\newpage

```r
## Regression, Seasonality, Trend = Local linear
getBstsDiagnosticPlots(bsts.fit8)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-19-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=1.72, p=0.09)
## H&W Stationarity: PASS (CMV=0.17, p=0.35)
## H&W Halfwidth: PASS (hw/mean=0.02 < eps=0.10)
## Std. Residuals:-----
## Geweke: PASS (z=0.25, p=0.80)
## H&W Stationarity: PASS (CMV=0.11, p=0.56)
## H&W Halfwidth: FAIL (hw/mean=0.10 < eps=0.10)
```
\newpage

```r
## Regression, Seasonality, Trend = Semilocal linear
getBstsDiagnosticPlots(bsts.fit9)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-20-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: FAIL (z=3.34, p=0.00)
## H&W Stationarity: PASS (CMV=0.36, p=0.09)
## H&W Halfwidth: PASS (hw/mean=0.00 < eps=0.10)
## Std. Residuals:-----
## Geweke: PASS (z=1.50, p=0.13)
## H&W Stationarity: PASS (CMV=0.13, p=0.46)
## H&W Halfwidth: FAIL (hw/mean=0.42 < eps=0.10)
```

Each transition between diagnostic panels contrasted here
(from `bsts.fit` to `bsts.fit2`, from `bsts.fit2` to `bsts.fit8`, 
and `bsts.fit8` to `bsts.fit9` ) should exhibit (assuming adequate priors
and enough MCMC iterations) improvement in the following visually 
discernible ways, which also relate to quantifiable fit improvement 
(i.e., decreased MAE). 

1. Adding seasonality in `bsts.fit2` begins to address the baseline model's
overly narrow prediction distribution 
2. Changing trend to local linear trend, and including regression, in `bsts.fit8` scales the predicted SD to the observed SD; however, the predictions are shifted (overestimated) and the standardized one-step
ahead prediction errors exhibit non-normality and autocorrelation. 
3. Changing trend to semilocal linear trend (and still including seasonality + 
regression) in `bsts.fit9` generally shows multifaceted improvements: (i) fixes the bias in the predicted distribution (realigned with observed distribution), and addresses the (ii) non-normality and (iii) autocorrelation problems in the residual chain. 

If the trace plots and convergence tests show that convergence is not reached 
for a given number of MCMC iterations, it is advised to 
repeat the analysis with a larger value of `bsts.niter`. If convergence is 
reached but issues of model fit or residual assumption violations persist, 
then respecifying the state space for the BSTS model may be necessary 
(hint: revisit parameter tuning, check alternative trend components, etc.).

With all diagnostic checks passed by the fitted BSTS model--that is, 
equipped with a quantitative and qualitative diagnosis of convergence and 
fit -- now the BSTS model is ready to serve as the counterfactual
for dynamic causal inference.


# 2. Dynamic Causal Inference: BSTS as Counterfactual {#Section2}

The preferred fitted bsts model `bsts.fit9` (i.e., the model with the most 
suitable combination of performance [low MAE] and diagnostics 
[convergence and fit]) can now be passed into the `CausalImpact()`function. 
This will use `bsts.fit9`, which was fitted on the pre-intervention 
data, in order to forecast (i.e., predict outside of the pre-intervention
sample time window) the counterfactual series (i.e., the hypothetical 
outcome series if there had been no exposure to treatment).  The pointwise 
causal impact of the intervention is then computed as the 
difference between the predicted (counterfactual untreated outcome) and 
observed (actual treated outcome). The cumulative impact estimate is the 
cumulative sum of the pointwise impact estimates. 

In this section, we demonstrate using a custom BSTS model in the 
`CausalImpact()` function. This allows for more modeling flexibility
than the default BSTS model created by `CausalImpact()` if 
no BSTS model is set for the `bsts.model` argument. Despite the 
simplicity or convenience of this default, we advise researchers 
to build custom BSTS models not only to take advantage of the flexibility of 
configuring the state space components added to the model but also to take
ownership of (and provide better explanations of) the modeling choices that are 
central to research designs involving BCA.

Here we compare the causal inference results for two different 
counterfactuals (i.e., two different BSTS models) that differ in terms of 
their trend component: 

* `bsts.fit8` has a local linear trend (suited for shorter forecasts)
* `bsts.fit9` has a semilocal linear trend (suited for longer forecasts)

Below we echo the summary of causal inference results for comparing these
candidate models.

```r
## Causal impact estimation: fitted BSTS model forecasts the counterfactual
impact8 <- CausalImpact(bsts.model = bsts.fit8,
                        post.period.response = post.period.response,
                        alpha=0.05, model.args=list(niter = bsts.niter))
impact9 <- CausalImpact(bsts.model = bsts.fit9,
                        post.period.response = post.period.response,
                        alpha=0.05, model.args=list(niter = bsts.niter))
summary(impact8)
```

```
## Posterior inference {CausalImpact}
## 
##                          Average       Cumulative    
## Actual                   7.6           1584.1        
## Prediction (s.d.)        6.6 (0.37)    1389.8 (77.31)
## 95% CI                   [6, 7.4]      [1256, 1543.6]
##                                                      
## Absolute effect (s.d.)   0.93 (0.37)   194.29 (77.31)
## 95% CI                   [0.19, 1.6]   [40.58, 327.8]
##                                                      
## Relative effect (s.d.)   14% (5.6%)    14% (5.6%)    
## 95% CI                   [2.9%, 24%]   [2.9%, 24%]   
## 
## Posterior tail-area probability p:   0.01333
## Posterior prob. of a causal effect:  98.667%
## 
## For more details, type: summary(impact, "report")
```

```r
summary(impact9)
```

```
## Posterior inference {CausalImpact}
## 
##                          Average       Cumulative      
## Actual                   7.6           1584.1          
## Prediction (s.d.)        6.7 (0.13)    1391.7 (26.15)  
## 95% CI                   [6.4, 6.9]    [1328.1, 1440.9]
##                                                        
## Absolute effect (s.d.)   0.92 (0.13)   192.47 (26.15)  
## 95% CI                   [0.69, 1.2]   [143.24, 256.1] 
##                                                        
## Relative effect (s.d.)   14% (1.9%)    14% (1.9%)      
## 95% CI                   [10%, 18%]    [10%, 18%]      
## 
## Posterior tail-area probability p:   0.00645
## Posterior prob. of a causal effect:  99.35484%
## 
## For more details, type: summary(impact, "report")
```
The `CausalImpact` package includes a convenience function to print a 
helpful (i.e., plain English) explanation of these Bayesian counterfactual 
analysis results from the BSTS model. We use this by calling the report
summary for a given model:
```
summary(impact8, 'report'); summary(impact9, 'report')
```
In paragraph formatting, the reports would look as follows:

**Model 8: Local Linear Trend**

> <br /><br />  <br /><br />  During the post-intervention period, the response variable had an average value of approx. 7.58. By contrast, in the absence of an intervention, we would have expected an average response of 6.65. The 95% interval of this counterfactual prediction is [6.01, 7.39]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 0.93 with a 95% interval of [0.19, 1.57]. For a discussion of the significance of this effect, see below.<br /><br />  <br /><br />  Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.58K. By contrast, had the intervention not taken place, we would have expected a sum of 1.39K. The 95% interval of this prediction is [1.26K, 1.54K].<br /><br />  <br /><br />  The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +14%. The 95% interval of this percentage is [+3%, +24%].<br /><br />  <br /><br />  This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (0.93) to the original goal of the underlying intervention.<br /><br />  <br /><br />  The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.013). This means the causal effect can be considered statistically significant.

**Model 9: Semilocal Linear Trend [Preferred]**

> <br /><br />  <br /><br />  During the post-intervention period, the response variable had an average value of approx. 7.58. By contrast, in the absence of an intervention, we would have expected an average response of 6.66. The 95% interval of this counterfactual prediction is [6.35, 6.89]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 0.92 with a 95% interval of [0.69, 1.23]. For a discussion of the significance of this effect, see below.<br /><br />  <br /><br />  Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.58K. By contrast, had the intervention not taken place, we would have expected a sum of 1.39K. The 95% interval of this prediction is [1.33K, 1.44K].<br /><br />  <br /><br />  The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +14%. The 95% interval of this percentage is [+10%, +18%].<br /><br />  <br /><br />  This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (0.92) to the original goal of the underlying intervention.<br /><br />  <br /><br />  The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.006). This means the causal effect can be considered statistically significant.

Then plotting the causal impact figures (original, pointwise, cumulative) is
useful for visually assessing the onset and decay structure of the 
dynamic treatment effect. In the simulation scenario, the DGP treatment
effect has an inverse U-shaped trajectory post intervention.
Therefore the onset and decay structure of the dynamic causal effect 
creates an s-shaped cumulative effect curve, which becomes cumulatively 
significant after a short span and remains so through the post-intervention 
window. 

The 95% Bayesian credible intervals of the model 8's 
(trend = local linear) posterior predictive distribution widen more quickly 
over time and display substantially more uncertainty in
the post-intervention window than the intervals of model 9
(trend = semilocal linear).
\newpage

```r
plot(impact8)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-22-1} \end{center}
\newpage

```r
plot(impact9)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_data-csv_R1_files/figure-latex/unnamed-chunk-23-1} \end{center}
Finally, we can print the results table of pointwise estimates of the 
average treatment effect on the treated (ATT), along 
with 95% Bayesian credible intervals (i.e., 95% posterior distribution interval).
In this example, we take the `$series` attribute of the `impact9` object
and add columns for the event time and significance codes ('*') where
the pointwise interval does not contain zero. 

Due to the time series length, we truncate the output within the interval from 
1 period before the intervention to 20 periods after the intervention. 

```r
## Show Pointwise ATT Estimates from BCA
restbl <- as_tibble(impact9$series)
restbl$event.time <- (1:npds) - intpd
restbl$`(sig.)` <- apply(restbl[,c('point.effect.lower','point.effect.upper')],
                       1, function(x){ ifelse(prod(x)>0,'*',"") })
rescols <- c('event.time','point.effect','(sig.)',
             'point.effect.lower', 'point.effect.upper')
knitr::kable(restbl[(intpd-1):(intpd+20), rescols], n = 22, digits = 4, 
             caption = 'Pointwise ATT Estimates from BCA')
```



Table: Pointwise ATT Estimates from BCA

| event.time| point.effect|(sig.) | point.effect.lower| point.effect.upper|
|----------:|------------:|:------|------------------:|------------------:|
|         -1|       0.3058|       |            -0.6211|             1.4109|
|          0|       0.0770|       |            -1.1952|             1.1243|
|          1|       0.2629|       |            -0.9424|             1.4688|
|          2|       0.1418|       |            -1.0907|             1.0672|
|          3|       0.0436|       |            -0.9956|             0.9392|
|          4|      -0.1031|       |            -1.0524|             0.8915|
|          5|       0.6123|       |            -0.5151|             1.6384|
|          6|       0.8806|       |            -0.0790|             1.9538|
|          7|      -0.3677|       |            -1.4743|             0.9219|
|          8|       1.0350|       |            -0.1934|             2.0893|
|          9|      -0.9399|       |            -2.0360|             0.0503|
|         10|       1.1860|       |            -0.1230|             2.3082|
|         11|      -0.5096|       |            -1.7559|             0.5981|
|         12|       0.0475|       |            -1.0300|             1.2055|
|         13|       1.3457|*      |             0.2401|             2.6005|
|         14|      -0.6073|       |            -1.5695|             0.5451|
|         15|       0.9648|       |            -0.1134|             1.8282|
|         16|       1.7068|*      |             0.6699|             2.8705|
|         17|       0.9508|       |            -0.0895|             1.9549|
|         18|       1.1626|*      |             0.0138|             2.3163|
|         19|      -0.2871|       |            -1.4504|             0.7035|
|         20|      -0.0223|       |            -1.0101|             1.0367|

The same process can then be followed to return the cumulative effect estimates 
for each post-intervention period along with their credible intervals. 

```r
## Show Cumulative ATT Estimates from BCA
cumutbl <- as_tibble(impact9$series)
cumutbl$event.time <- (1:npds) - intpd
cumutbl$`(sig.)` <- apply(cumutbl[,c('cum.effect.lower','cum.effect.upper')],
                       1, function(x){ ifelse(prod(x)>0,'*','') })
cumucols <- c('event.time','cum.effect','(sig.)',
             'cum.effect.lower', 'cum.effect.upper')
knitr::kable(cumutbl[(intpd-1):(intpd+20), cumucols], n = 22, digits = 4,
             caption = 'Cumulative ATT Estimates from BCA')
```



Table: Cumulative ATT Estimates from BCA

| event.time| cum.effect|(sig.) | cum.effect.lower| cum.effect.upper|
|----------:|----------:|:------|----------------:|----------------:|
|         -1|     0.0000|       |           0.0000|           0.0000|
|          0|     0.0770|       |          -1.1952|           1.1243|
|          1|     0.3399|       |          -1.2413|           1.8424|
|          2|     0.4817|       |          -1.4603|           2.4550|
|          3|     0.5253|       |          -1.8019|           2.6432|
|          4|     0.4223|       |          -2.1715|           2.7244|
|          5|     1.0346|       |          -1.3487|           3.7741|
|          6|     1.9152|       |          -0.9311|           4.9450|
|          7|     1.5474|       |          -1.2846|           4.9799|
|          8|     2.5825|       |          -0.7785|           5.9150|
|          9|     1.6426|       |          -2.2611|           5.0320|
|         10|     2.8286|       |          -1.6283|           6.4596|
|         11|     2.3190|       |          -2.2491|           6.1403|
|         12|     2.3665|       |          -2.0723|           6.1448|
|         13|     3.7122|       |          -1.1916|           7.6024|
|         14|     3.1049|       |          -1.6536|           7.4760|
|         15|     4.0698|       |          -0.8468|           8.7444|
|         16|     5.7766|*      |           0.3393|          10.8347|
|         17|     6.7273|*      |           0.5512|          11.6827|
|         18|     7.8899|*      |           1.7157|          13.1500|
|         19|     7.6028|*      |           1.4678|          13.1419|
|         20|     7.5805|*      |           1.1615|          13.6034|
These results of BCA for dynamic causal inference are valuable because they
enable researchers and practitioners to move beyond just asking whether 
or not a causal effect existed. Instead, this approach begins to enable 
addressing temporal questions such as, when does a causal effect become 
significant? Or when does it switch to become no longer significant? Also, for 
how long does the causal effect remain significant? And what is the 
cumulative effect within a given post-intervention time interval?

This concludes the vignette illustrating how to use BSTS as the predictive model
in a Bayesian counterfactual approach to dynamic causal inference. 

Please note that references not yet included will be added to a future version of 
this document. 
