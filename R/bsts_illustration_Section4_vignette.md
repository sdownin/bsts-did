---
title: "A Bayesian Counterfactual Approach to Dynamic Causal Inference: R Code and Tutorial"
author: 
  - "Author 1^[Information redacted for review]"
  - "Author 2^[Information redacted for review]"
email: "<email redacted for review>"
date: "January 27, 2023"
output:
  pdf_document:
    keep_md: true
    toc: true
---




\newpage

This vignette presents `R` code and a tutorial example of  Bayesian Counterfactual 
Analysis (BCA). This approach includes two main steps: (1) building a Bayesian 
structural time series (BSTS) model 
[Section 1](#Section1), and (2) using the BSTS model as a
counterfactual for causal inference [Section 2](#Section2).

* [**Section 1**](#Section1) covers the complexity of the Bayesian approach. The basic structure of a BSTS model is local trend (state space) + seasonality (state space) + regression. Each component has several ways to be defined, and prior settings are slightly different. Thus, we aim to show a step-by-step process of adding different components in [Section 1](#Section1) accompanied by checking resulting BSTS model diagnostics, setting up Bayesian priors,and adding regression component (Spike and slab priors included). This part requires the `bsts` package in R: [https://cran.r-project.org/web/packages/bsts/bsts.pdf](https://cran.r-project.org/web/packages/bsts/bsts.pdf). 
* [**Section 2**](#Section2) shows how to assess causal effects using the BSTS model from [Section 1](#Section1) as a counterfactual, and is a straightforward process using the `CausalImpact` R package: [https://google.github.io/CausalImpact/CausalImpact.html](https://google.github.io/CausalImpact/CausalImpact.html)  


# 0. Simulated Data

First, load the data simulated from a data generating process (DGP) with the 
following features: 

- local trend involving a random walk of small step sizes (i.e., Gaussian noise standard deviation `sd=0.01` in the local level stochastic process), 
- weekly seasonality with a frequency of one per cycle (i.e., 52 weeks = 1 year), and 
- positive linear trend (i.e., upward linear trajectory over time, on average)

This data set includes the outcome series `y_observed`, and nine 
covariate series, `c1-c9`, from which the BSTS models that have regression
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
df1 <- as_tibble(read.csv(dat1$bsts.df1.filepath)) ##convert table, print labels
print(df1)
```

```
## # A tibble: 520 x 10
##    y_observed c1_mean   c2_mean  c3_mean  c1_sd  c2_sd   c3_sd  c1_skew  c2_skew
##         <dbl>   <dbl>     <dbl>    <dbl>  <dbl>  <dbl>   <dbl>    <dbl>    <dbl>
##  1     0.0234  0.0111 -0.00930   0.0187  0.0542 0.0490 0.0102   0.00423  2.64e-1
##  2     0.134   0.0272 -0.00735   0.00747 0.0453 0.0466 0.00905 -0.248    2.81e-1
##  3     0.815   0.0288  0.00206   0.0683  0.0503 0.0479 0.00817 -0.234   -3.39e-1
##  4     0.578   0.0519  0.00210   0.129   0.0411 0.0443 0.0113   0.317   -1.70e-1
##  5     0.278   0.0565 -0.00104   0.123   0.0441 0.0438 0.00974 -0.0773   2.48e-1
##  6    -0.367   0.0491 -0.000560  0.0670  0.0424 0.0539 0.0106  -0.353    9.81e-2
##  7     1.45    0.0728 -0.00903  -0.0777  0.0530 0.0526 0.00988  0.0278  -2.66e-1
##  8     1.85    0.0861 -0.00824   0.286   0.0557 0.0449 0.00913  0.194   -4.21e-3
##  9     1.90    0.0834  0.000128  0.396   0.0543 0.0433 0.0104  -0.331   -4.16e-2
## 10     3.13    0.106   0.00466   0.402   0.0441 0.0459 0.00964 -0.135    1.86e-4
## # ... with 510 more rows, and 1 more variable: c3_skew <dbl>
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


\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/outcome_plot-1} 
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
## INDICES OF PRE-INTERVENTION and POST-INTERVENTION WINDOWS
pre.idx <- 1:(intpd-1)
post.idx <- intpd:npds
## Outcome (response) series
y <- df1$y_observed
## Train on y pre-treatment but NA's post-treatment
y.pre.treat.NAs.post <- c( y[pre.idx], rep(NA,length(post.idx)) )
## Then use the post-treatment response for causal impact estimation
post.period.response <- y[post.idx]
## Covariates (predictors) - data for the "formula = y ~ predictors" argument
predictors <- df1[ , ! names(df1) %in% 'y_observed']
```


# 1. Bayesian Structural Time Series (BSTS) Modeling {#Section1}

## 1.1. State Space Specification: Local Trend Only

The simplest BSTS model is with only the local trend in the state space. 
There are several state spaces in the `bsts` package for defining local trend including local level, 
local linear trend, student local linear trend, and generalized local linear trend.

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
bsts.niter <- 5000 ## suggest on the order of 1k-10k at least
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
## =-=-=-=-= Iteration 0 Fri Jan 27 22:35:14 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 500 Fri Jan 27 22:35:16 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 1000 Fri Jan 27 22:35:18 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 1500 Fri Jan 27 22:35:19 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 2000 Fri Jan 27 22:35:21 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 2500 Fri Jan 27 22:35:23 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 3000 Fri Jan 27 22:35:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 3500 Fri Jan 27 22:35:26 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 4000 Fri Jan 27 22:35:28 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 4500 Fri Jan 27 22:35:30 2023
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
## [1] 0.4046933
## 
## $prediction.sd
## [1] 0.5837429
## 
## $rsquare
## [1] 0.9499213
## 
## $relative.gof
## [1] NA
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-5-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-5-2} \end{center}
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



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-6-1} \end{center}
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
st.sp2 <- AddSeasonal(st.sp2, y.pre.treat.NAs.post, 
                      nseasons = 52, season.duration = 1)
## Fit BSTS model and get summary plots
bsts.fit2 <- bsts(y.pre.treat.NAs.post ,
                 state.specification = st.sp2, 
                 niter = bsts.niter,
                 seed = rand.seed, ping=0)
getBstsSummaryPlots(bsts.fit2)
```

```
## $residual.sd
## [1] 0.5438631
## 
## $prediction.sd
## [1] 0.6354091
## 
## $rsquare
## [1] 0.9095559
## 
## $relative.gof
## [1] NA
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-7-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-7-2} \end{center}

```r
check2 <- plotBstsStateComps(bsts.fit2, intpd=intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-7-3} \end{center}
Following BSTS (and other forecasting) conventions, 
it is common to use a model performance metric, 
such as the one-step ahead prediction (forecast) error within sample
during the pre-intervention (or training) window. We might expect that the augmented
state space should improve the model performance (e.g., in terms of 
decreased MAE). However, the MAE of this model is now 0.599.

Here we have a seasonal component separate from the trend (which includes 
one the local level at this point). Now the state space of the fitted model 
exhibits the seasonal structure of the time series into the post-intervention 
forecast window. 

However, besides seasonality, the observed data (Y) also appears to have
a positive linear trend, which is not currently captured by the 
local level (default without slope or drift) of the state space trend. 
This is apparent in trend component plot where the upward linear 
trajectory before the intervention turns into a constant expectation after the 
intervention -- as a simple random walk where the expetation at `t` is 
the trend at `t-1` (also just the local level without slope or drift) 
plus Gaussian noise. Therefore, this seasonal model still does not capture the linear trend in the 
DGP and the fit could also be improved with more information about the outcome from 
covariate series with regression. But before moving on to addressing 
regression in BSTS models, we first detail the process of parameter tuning. 

A crucial step in any machine learning (ML) analysis is hyperparameter tuning,
which applies to our investigation of dynamic causal inference using a Bayesian
counterfactual approach drawing upon the ML tradition. Therefore,
adjusting the values of prior distributions in the Bayesian counterfactual
approach plays an important role in achieving adequate model fit for 
complex time series structures. The BSTS model should, of course, be able to 
mimic the observed outcome series pre-intervention and project that 
state space forward into the post-intervention window while also 
characterizing the uncertainty via credible intervals. 

Additionally, in terms of research design and process, 
prior parameter tuning represents an important step for researchers 
to incorporate information from theory and context into their empirical
analyses--coherently within a Bayesian framework. 

## 1.3. Bayesian Priors:  Parameter Tuning

Setting up Bayesian priors is notorious for its complication. 
The conventional wisdom is to exploit information of the observed outcome. 
There are two main types of priors to set up for local level and seasonality,  priors for the variance of the state innovation errors and the initial value of the state at time 0, provided by Boom package.  
The former is defined by the option of `sigma.prior` created by `SdPrior` which describes the prior distribution for the standard deviation of the random walk increment (e.g., inverse Gamma priors). 
The latter one is `initial.state.prior` defined by `NormalPrior` to describe the prior distribution of the initial state vector (e.g., Gaussian distribution). If not defined, the default setting will be applied.
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
window, `y[pre.idx]` (which is the same as `y.pre.treat.NAs.post` vector
within the `pre.idx` interval). 

```r
## Get outcome (Y) standard deviation to use for adjusting priors below
y.sd <- sd(y[pre.idx], na.rm = T)
## Initiate empty state space configuration list
st.sp3 <- list()
## Add local level to trend component of state space
st.sp3 <- AddLocalLevel(st.sp3, y.pre.treat.NAs.post, 
  initial.state.prior = NormalPrior(mu=1, sigma = 0.25 * y.sd, fixed = FALSE),
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
## [1] 0.5490828
## 
## $prediction.sd
## [1] 0.6385154
## 
## $rsquare
## [1] 0.9078115
## 
## $relative.gof
## [1] NA
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-8-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-8-2} \end{center}

```r
check3 <- plotBstsStateComps(bsts.fit3, intpd=intpd, return.val = T) 
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-8-3} \end{center}
The MAE of this model is now 0.602.

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
about the outcome series pre-intervention (and that relationship
does not change for the covaritate post-intervention), then including 
regression in the BSTS model should generally improve the model's performance 
as a counterfactual (i.e., as a representation of the pre-intervention trend 
of the outcome series projected forward into the post-intervention window).
In using BSTS for causal inference, the researcher can include 
covariates with information about (i.e., correlation with) the outcome 
of interest and therefore improve the representativeness of the 
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
  initial.state.prior = NormalPrior(mu=1, sigma = .25 * y.sd, fixed = FALSE),
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
## [1] 0.5550138
## 
## $prediction.sd
## [1] 0.6227451
## 
## $rsquare
## [1] 0.9058091
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.0000  1.0000  0.6444  1.0000  4.0000 
## 
## $coefficients
##                      mean           sd     mean.inc      sd.inc    inc.prob
## c3_mean      3.959457e-01  0.496903347   0.89961323  0.32837778 0.440128799
## c3_sd       -2.888089e+00 13.958671675 -46.44309312 33.36794863 0.062185550
## c1_mean      1.624523e-02  0.086719868   0.29353657  0.23379885 0.055343127
## c2_mean     -6.864345e-02  0.664991915  -1.67200646  2.85100482 0.041054538
## c2_sd        1.052319e-01  1.137725551   5.80997063  6.22398670 0.018112296
## c1_sd        4.332151e-02  0.760389836   2.79564370  5.47720625 0.015496076
## c3_skew     -6.621589e-04  0.011578760  -0.11750956  0.10212357 0.005634937
## c2_skew     -6.286286e-05  0.005542671  -0.01644029  0.09052538 0.003823707
## c1_skew      4.148803e-05  0.004101215   0.01585800  0.08180309 0.002616221
## (Intercept)  0.000000e+00  0.000000000   0.00000000  0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-9-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-9-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-9-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-9-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-9-5} \end{center}

```r
check4 <- plotBstsStateComps(bsts.fit4, intpd=intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-9-6} \end{center}
The MAE of this model is now 0.594.

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
## [1] 0.5719478
## 
## $prediction.sd
## [1] 0.5907259
## 
## $rsquare
## [1] 0.8999738
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   2.000   3.000   2.683   3.000   5.000 
## 
## $coefficients
##                      mean           sd      mean.inc      sd.inc    inc.prob
## c1_mean      1.073194e+00  0.171866160   1.073194212  0.17186616 1.000000000
## c3_mean      1.111433e+00  0.364257725   1.131857984  0.33465551 0.981954225
## c3_sd       -2.866113e+01 37.838101591 -64.923325576 29.80501549 0.441461268
## c2_mean     -1.122123e-01  0.854776353  -1.053496896  2.42678038 0.106514085
## c2_sd        8.944287e-02  1.542325912   1.373068888  5.91397162 0.065140845
## c1_sd       -3.217632e-02  1.343161233  -0.580195273  5.69702220 0.055457746
## c3_skew     -2.205028e-03  0.022955217  -0.139161756  0.12077799 0.015845070
## c2_skew      4.335044e-04  0.011152362   0.049246097  0.11106646 0.008802817
## c1_skew      6.430438e-05  0.008182274   0.008594092  0.09708127 0.007482394
## (Intercept)  0.000000e+00  0.000000000   0.000000000  0.00000000 0.000000000
```

```r
par(mfrow=c(2,1))
plot(bsts.fit5, 'coefficients', main='Expected Size = 6')
plot(bsts.fit4, 'coefficients', main='Expected Size = 3')
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-10-1} \end{center}
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
## [1] 0.5706167
## 
## $prediction.sd
## [1] 0.5942611
## 
## $rsquare
## [1] 0.9004388
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   2.000   2.000   2.437   3.000   5.000 
## 
## $coefficients
##                      mean           sd     mean.inc     sd.inc    inc.prob
## c1_mean      1.092807e+00  0.194486998   1.09280704  0.1944870 1.000000000
## c3_mean      9.744147e-01  0.457457470   1.08244526  0.3398228 0.900197628
## c3_sd       -1.526289e+01 27.625428099 -51.40115629 26.6925291 0.296936759
## c1_sd        1.436006e-01  1.708107007   1.93765107  6.0094601 0.074110672
## c2_sd        2.759625e-01  1.803810308   4.04744947  5.7157373 0.068181818
## c2_mean     -2.205036e-02  0.592119459  -0.39848157  2.4977012 0.055335968
## c3_skew     -2.183906e-03  0.022161250  -0.09822725  0.1137186 0.022233202
## c1_skew     -3.153322e-04  0.013244727  -0.02363823  0.1143535 0.013339921
## c2_skew      1.954095e-04  0.009540232   0.02825064  0.1153687 0.006916996
## (Intercept)  0.000000e+00  0.000000000   0.00000000  0.0000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-11-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-11-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-11-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-11-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-11-5} \end{center}

```r
check6 <- plotBstsStateComps(bsts.fit6, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-11-6} \end{center}
The MAE of this model is now 0.572.

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
# 2B. Specify prior.inclusion.probabilities for each covariate
#   where length of prior spikes and means vectors should equal the number of
#    covariates
#   (i.e., this example with vector lengths = 11 implies
#      ‘data’ has 11 columns) <-- this includes intercept in example below
prior.spikes <- c(.7,.8,.9,.1,.1,.1,.1,.1,.1,.1) ## +1 column for intercept (?)
prior.means <- c(.2,.2,.3,0,0,0,0,0,0,0) ## +1 column for intercept (?)
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
## [1] 0.5677897
## 
## $prediction.sd
## [1] 0.6103475
## 
## $rsquare
## [1] 0.9014229
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   2.000   2.000   2.209   3.000   5.000 
## 
## $coefficients
##                      mean           sd     mean.inc      sd.inc    inc.prob
## c1_mean      9.902431e-01  0.322728499   0.99854701  0.31101662 0.991683992
## c2_mean     -1.301793e+00  2.444385586  -2.15732042  2.83864470 0.603430353
## c3_mean      5.160681e-01  0.562835518   1.01421356  0.34225162 0.508835759
## c3_sd       -3.936292e+00 16.839217896 -56.94304506 33.01076318 0.069126819
## c2_sd        1.014585e-02  0.648327050   0.62969734  5.15168790 0.016112266
## c1_sd        1.865224e-02  0.638489250   1.56030036  5.75492510 0.011954262
## c1_skew     -7.981994e-05  0.005744201  -0.01706373  0.08720861 0.004677755
## c2_skew     -6.787409e-05  0.003373802  -0.03264744  0.07666825 0.002079002
## c3_skew     -5.725580e-05  0.004220951  -0.03672005  0.12292725 0.001559252
## (Intercept)  0.000000e+00  0.000000000   0.00000000  0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-12-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-12-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-12-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-12-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-12-5} \end{center}

```r
check7 <- plotBstsStateComps(bsts.fit7, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-12-6} \end{center}
The MAE of this model is now 0.587.

## 1.5. Local Trends Comparison: Level vs. Linear, Short(er)-Term vs. Long(er)-Term

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

## CALL OUR CUSTOM BSTS SUMMARY
getBstsSummaryPlots(bsts.fit8)
```

```
## $residual.sd
## [1] 0.575704
## 
## $prediction.sd
## [1] 0.5901336
## 
## $rsquare
## [1] 0.8986556
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   2.000   2.000   2.016   2.000   5.000 
## 
## $coefficients
##                      mean          sd      mean.inc      sd.inc    inc.prob
## c3_mean      1.252420e+00 0.264843483   1.252420057  0.26484348 1.000000000
## c1_mean      9.131656e-01 0.332629549   0.957548547  0.27112884 0.953649441
## c2_mean      4.942052e-03 0.353408561   0.223523656  2.38052440 0.022109750
## c2_sd        1.510250e-02 0.418068018   1.453712611  3.88787215 0.010388918
## c1_sd        8.963142e-03 0.313677494   0.862759875  2.99360207 0.010388918
## c3_sd       -1.037051e-01 2.091544865 -11.797243124 19.25594350 0.008790623
## c3_skew     -3.328163e-04 0.007933837  -0.089242303  0.09811793 0.003729355
## c1_skew      1.247191e-05 0.005736277   0.003344255  0.09740321 0.003729355
## c2_skew      3.821980e-05 0.008137571   0.011956426  0.14979194 0.003196590
## (Intercept)  0.000000e+00 0.000000000   0.000000000  0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-13-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-13-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-13-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-13-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-13-5} \end{center}

```r
check8 <- plotBstsStateComps(bsts.fit8, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-13-6} \end{center}
The MAE of this model is now 0.549.

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
## [1] 0.5837771
## 
## $prediction.sd
## [1] 0.5970417
## 
## $rsquare
## [1] 0.8957934
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   2.000   2.000   2.114   2.000   4.000 
## 
## $coefficients
##                      mean          sd      mean.inc      sd.inc    inc.prob
## c1_mean      1.184172e+00 0.104132663   1.184171947  0.10413266 1.000000000
## c3_mean      1.085464e+00 0.234456155   1.086337364  0.23251819 0.999195818
## c2_mean     -1.237295e-01 0.623613164  -2.270960912  1.50615436 0.054483313
## c3_sd       -5.013669e-01 4.190292303 -19.036631712 17.77973061 0.026336952
## c1_sd       -2.348517e-02 0.363988602  -1.769927612  2.64531598 0.013268999
## c2_sd       -1.257029e-02 0.333527907  -0.992454407  2.81680994 0.012665862
## c3_skew     -6.232189e-04 0.011367683  -0.123995626  0.10412154 0.005026136
## c1_skew      9.546639e-06 0.002394752   0.005276109  0.05944458 0.001809409
## c2_skew      3.980692e-05 0.002655263   0.028285662  0.07008219 0.001407318
## (Intercept)  0.000000e+00 0.000000000   0.000000000  0.00000000 0.000000000
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-14-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-14-2} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-14-3} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-14-4} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-14-5} \end{center}

```r
check9 <- plotBstsStateComps(bsts.fit9, intpd = intpd, return.val = T)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-14-6} \end{center}
The MAE of this model is now 0.553.

Zooming into the first 100 periods shows how the BSTS model learns the structure
of the time series through the first couple of cycles that train the 
the seasonality and regression structures into the fitted model (via
recursive least squares in the Kalman filter). 


```r
par(mfrow=c(1,2), mar=c(4,4,.5,1))
plotBstsStateComps(bsts.fit9, intpd=intpd, pd.ids=1:104, title.show=F)##2 cycles
plotBstsStateComps(bsts.fit9, intpd=intpd, pd.ids=1:520, title.show=F)##10 cycles
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-15-1} \end{center}

## 1.6. MCMC Diagnostics: Checking Convergence and Fit {#Section1-6}

Draws from the posterior distribution of the Markov chain(s)--once
converged at the stationary distribution--are used to estimate the 
model prediction at each time period though 
Markov chain Monte Carlo (MCMC) estimation. However, when the Markov chains
have not yet converged (usually due to insufficient MCMC burn-in period
or sampled iterations), then the iterations cannot be treated as draws from
the stationary distribution to sample the posterior predictive distribution. 
Thus, the estimates may be unstable and/or biased, so users must assess 
convergence of the MCMC chains. 

There are several ways to check the convergence 
of MCMC including the trace plots, Heidelberger-Welch tests, and Geweke 
diagnostic tests. Additionally, the researcher may evalue autocorrelation 
using the autocorrelation function (ACF) plots, among other tests.

Similarly, there is need to assess how the estimated model fits the observed data, 
in terms of how well new samples simulated from the fitted model create 
simulated distribution(s) that resemble the observed value(s). 

For this purpose, we examine MCMC convergence trace figures and plot model fit 
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
###
## Define a funciton to call all of our plotting, 
##    diagnostics and converngence checks functions
###
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
for dyamic causal inference, we adapt the diagnostics panel below from 
Menchetti & Bojinov (2022) to cover the following scope of diagnostic checks: 

- Column 1 (left) addresses MCMC convergence along two dimensions 
(posterior predictive distribution chain, and standardized residual chain) 
- Columns 2-3 (middle, right) address model fit or problems, 
conditional on convergence established in the first column. 

Regarding MCMC convergence, the first column shows the trace plots
of the Markov chain (i.e., mean of ATT for each draw) for the 
posterior predictive distribution of the outcome (Y). 
For model fit, the observed and predicted values of Y are compared (top-center), 
where tighter alignment indicates better fit. 

Three MCMC convergence diagnostics are used to assess convergence of the Markov
chains for each of the two traces (posterior predictive distribution and
strandardized residuals). 

1. *Geweke's Convergence Test* (see `coda::geweke.diag()`) assess convergence of the MCMC chain with the Geweke Z-score based on equality of means of the fist and last portions of the chain (defaults: 0.1 and 0.5).  If the chain has reached the stationary distribution, the two measures are equal and the low Z-score would not reject the null hypothesis of convergence.
2. *Heidelberg & Welch's Stationary Test* (see `coda::heidel.diag()`) evaluates the null hypothesis that the samples come from the same stationary distribution according to the Cramer-von-Mises (CMV) statistic.
3. *Heidelberg & Welch's Halfwidth Ratio Test* (see `coda::heidel.diag()`) calcuates a 95% interval for the mean of the Markov chain, using the portion of the chain that possed the H&W stationary test above. If the ratio of interval halfwidth-to-mean is not less than some tolerance (default `tol=0.1`), then the sample chain length is not long enough to estimate the mean with sufficient accuracy, suggesting the need for more MCMC iterations. 

Additionally, it is necessary to assess model fit and residual assumptions.
We use the maximum value of Y as an auxiliary statistic for computing 
a Bayesian p-val (Menchetti & Bojinov, 2022), which represents the proportion 
of the maximum Y values (from MCMC iterations) that are larger than the 
observed largest value of Y. Smaller values of Bayesian p-val indicate worse 
fit, which offers a similar interpretation to the frequentist null 
hypothesis statistical test (NHST) under a null assumption of a well fitting model. 
Values of Bayesian p-val >= alpha (conventionally 0.05) would therefore not reject the 
null of good fit, whereas Bayesian p-val < alpha indicates problem(s) of fit to 
be addressed (e.g., by increasing MCMC iterations, or respecifying the state space).
\newpage

```r
## Baseline (Trend = local level)
getBstsDiagnosticPlots(bsts.fit)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-17-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=0.29, p=0.77)
## H&W Stationarity: PASS (CMV=0.05, p=0.86)
## H&W Halfwidth: PASS (hw/mean=0.03 < eps=0.10)
## Std. Residuals:-----
## Geweke: PASS (z=0.44, p=0.66)
## H&W Stationarity: PASS (CMV=0.20, p=0.26)
## H&W Halfwidth: PASS (hw/mean=0.01 < eps=0.10)
```
\newpage

```r
## Seasonality Baseline
getBstsDiagnosticPlots(bsts.fit2)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-18-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=0.09, p=0.93)
## H&W Stationarity: PASS (CMV=0.06, p=0.81)
## H&W Halfwidth: PASS (hw/mean=0.01 < eps=0.10)
## Std. Residuals:-----
## Geweke: PASS (z=1.44, p=0.15)
## H&W Stationarity: PASS (CMV=0.09, p=0.63)
## H&W Halfwidth: PASS (hw/mean=0.05 < eps=0.10)
```
\newpage

```r
## Regression, Seasonality, Trend = Local linear
getBstsDiagnosticPlots(bsts.fit8)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-19-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=1.69, p=0.09)
## H&W Stationarity: PASS (CMV=0.10, p=0.59)
## H&W Halfwidth: FAIL (hw/mean=0.20 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=2.07, p=0.04)
## H&W Stationarity: PASS (CMV=0.17, p=0.34)
## H&W Halfwidth: FAIL (hw/mean=0.50 < eps=0.10)
```
\newpage

```r
## Regression, Seasonality, Trend = Semilocal linear
getBstsDiagnosticPlots(bsts.fit9)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-20-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: FAIL (z=7.22, p=0.00)
## H&W Stationarity: PASS (CMV=0.11, p=0.52)
## H&W Halfwidth: PASS (hw/mean=0.02 < eps=0.10)
## Std. Residuals:-----
## Geweke: PASS (z=0.90, p=0.37)
## H&W Stationarity: PASS (CMV=0.17, p=0.33)
## H&W Halfwidth: FAIL (hw/mean=0.23 < eps=0.10)
```

Each transition between diagnostic panels contrasted here
(from `bsts.fit` to `bsts.fit2`, from `bsts.fit2` to `bsts.fit8`, 
and `bsts.fit8` to `bsts.fit9` ) should exhibit (assuming adequate priors
and enough MCMC iterations) incremental improvement in the following visually 
discernible ways, which also relate to quantifiable fit improvement 
(i.e., decreased MAE). 

1. Adding seasonality in `bsts.fit2` begins to address the baseline model's
overly narrow prediction distribution 
2. Changing trend to local linear trend, and including regression, in `bsts.fit8` scales the predicted SD to the observed SD; however, the predictions are shifted (overestimated) and the standardized one-step
ahead prediction errors exhibit non-normality and autocorrelation. 
3. Changing trend to semilocal linear trend (and still including regression) in 
`bsts.fit9` should (a) fix the bias in the predicted distribution
(realigned with observed distribution), and address the (b) non-normality
and (c) autocorrelation problems in the residual chain.

If the trace plots and convergence tests show that convergence is not reached 
for a given number of MCMC iterations, it is advised to 
repeat the analysis with a larger value of `bsts.niter`. If convergence is 
reached but issues of model fit or residual assumption violations persist, 
then respecifying the state space for the BSTS model may be necessary 
(hint: revisit parameter tuning, check alternative trend components, etc.).

With all diagnostic checks passed by the fitted BSTS model--that is, 
equipped with a quantitative and qualitative diagnosis of convergence and good
fit -- now the BSTS model is ready to serve as the counterfactual
for dynamic causal inference.


# 2. Dynamic Causal Inference: BSTS as Counterfactual {#Section2}

The preferred fitted bsts model `bsts.fit9` (i.e., the model with the least 
pre-intervention, one-step-ahead cumulative absolute prediction error) 
can now be passed into the `CausalImpact()`function. 
This will use `bsts.fit9`, which was fitted on the pre-intervention 
data, in order to forecast (i.e., predict outside of the pre-intervention
sample time window) the counterfactual series (i.e., the hypothetical 
outcome series if there had been no exposure to treatment).  The pointwise 
causal impact of the intervention is then computed as the 
difference between the predicted (counterfactual untreated outcome) and 
observed (actual treated outcome). The cumulative impact estimate is the 
cumulative sum of the pointwise impact estimates. 

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
##                          Average        Cumulative      
## Actual                   8.1            1685.5          
## Prediction (s.d.)        7.5 (0.36)     1566.8 (75.63)  
## 95% CI                   [6.8, 8.2]     [1418.9, 1715.5]
##                                                         
## Absolute effect (s.d.)   0.57 (0.36)    118.71 (75.63)  
## 95% CI                   [-0.14, 1.3]   [-29.92, 266.7] 
##                                                         
## Relative effect (s.d.)   7.6% (4.8%)    7.6% (4.8%)     
## 95% CI                   [-1.9%, 17%]   [-1.9%, 17%]    
## 
## Posterior tail-area probability p:   0.05806
## Posterior prob. of a causal effect:  94%
## 
## For more details, type: summary(impact, "report")
```

```r
summary(impact9)
```

```
## Posterior inference {CausalImpact}
## 
##                          Average        Cumulative      
## Actual                   8.1            1685.5          
## Prediction (s.d.)        7.4 (0.13)     1554.6 (27.69)  
## 95% CI                   [7.2, 7.7]     [1503.9, 1614.0]
##                                                         
## Absolute effect (s.d.)   0.63 (0.13)    130.97 (27.69)  
## 95% CI                   [0.34, 0.87]   [71.53, 181.59] 
##                                                         
## Relative effect (s.d.)   8.4% (1.8%)    8.4% (1.8%)     
## 95% CI                   [4.6%, 12%]    [4.6%, 12%]     
## 
## Posterior tail-area probability p:   8e-04
## Posterior prob. of a causal effect:  99.9196%
## 
## For more details, type: summary(impact, "report")
```
The `CausalImpact` package includes a convenience function to print a 
helpful (i.e., plain English) explanation of these Bayesian counterfactual 
analysis results from the BSTS model.

```r
summary(impact9, 'report')
```

```
## Analysis report {CausalImpact}
## 
## 
## During the post-intervention period, the response variable had an average value of approx. 8.06. By contrast, in the absence of an intervention, we would have expected an average response of 7.44. The 95% interval of this counterfactual prediction is [7.20, 7.72]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 0.63 with a 95% interval of [0.34, 0.87]. For a discussion of the significance of this effect, see below.
## 
## Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.69K. By contrast, had the intervention not taken place, we would have expected a sum of 1.55K. The 95% interval of this prediction is [1.50K, 1.61K].
## 
## The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +8%. The 95% interval of this percentage is [+5%, +12%].
## 
## This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (0.63) to the original goal of the underlying intervention.
## 
## The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.001). This means the causal effect can be considered statistically significant.
```

In paragraph formatting, this report would look as follows:

> <br /><br />  <br /><br />  During the post-intervention period, the response variable had an average value of approx. 8.06. By contrast, in the absence of an intervention, we would have expected an average response of 7.44. The 95% interval of this counterfactual prediction is [7.20, 7.72]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 0.63 with a 95% interval of [0.34, 0.87]. For a discussion of the significance of this effect, see below.<br /><br />  <br /><br />  Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.69K. By contrast, had the intervention not taken place, we would have expected a sum of 1.55K. The 95% interval of this prediction is [1.50K, 1.61K].<br /><br />  <br /><br />  The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +8%. The 95% interval of this percentage is [+5%, +12%].<br /><br />  <br /><br />  This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (0.63) to the original goal of the underlying intervention.<br /><br />  <br /><br />  The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.001). This means the causal effect can be considered statistically significant.

Then plotting the causal impact figures (original, pointwise, cumulative) is
useful for visually assessing the onset and decay structure of the 
dynamic treatment effect. In simulation scenario, the DGP has an inverse 
U-shape trajectory post intervention. Therefore the onset and
decay structure of the dynamic causal effect creates an s-shaped cumulative 
effect curve, which becomes cumulatively significant early and remains so 
through the post-intervention window. 

The 95% Bayesian credible intervals of the first model's 
(trend = local linear) posterior predictive distribution widen more quickly 
over time and display substantially more uncertainty in
the post-intervention window than the intervals of second model
(trend = semilocal linear). 

```r
plot(impact8); plot(impact9)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-23-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-23-2} \end{center}
Finally, we can print the results table of pointwise estimates of the ATT, along 
with 95% Bayesian credible intervals from the posterior predictive distribution.
In this example, we take the `$series` attribute of the `impact9` object
and add columns for the event time and significance codes ('*') where
the pointwise interval does not contain zero. 

Due to the time series length, we truncate the output within the interval from 
1 period before the intervention to 20 periods after the intervention. 

```r
## Show Pointwise ATT Estimates from BSTS CausalImpact
restbl <- as_tibble(impact9$series)
restbl$event.time <- (1:npds) - intpd
restbl$`(sig.)` <- apply(restbl[,c('point.effect.lower','point.effect.upper')],
                       1, function(x){ ifelse(prod(x)>0,'*',"") })
rescols <- c('event.time','point.effect','(sig.)',
             'point.effect.lower', 'point.effect.upper')
restbl <- restbl[,rescols]
knitr::kable(restbl[(intpd-1):(intpd+20), ], n = 22, digits = 4, 
             caption = 'Pointwise ATT Estimates from BCA')
```



Table: Pointwise ATT Estimates from BCA

| event.time| point.effect|(sig.) | point.effect.lower| point.effect.upper|
|----------:|------------:|:------|------------------:|------------------:|
|         -1|      -0.1441|       |            -1.3979|             1.0720|
|          0|       0.0557|       |            -1.2462|             1.3260|
|          1|       0.1830|       |            -1.0773|             1.3887|
|          2|      -0.6502|       |            -1.8756|             0.5536|
|          3|       0.4652|       |            -0.7692|             1.7104|
|          4|       0.4094|       |            -0.8611|             1.6591|
|          5|       0.8922|       |            -0.3939|             2.1296|
|          6|       0.6106|       |            -0.6108|             1.8742|
|          7|      -0.0847|       |            -1.3427|             1.1791|
|          8|       0.4701|       |            -0.7707|             1.7085|
|          9|      -1.5660|*      |            -2.8069|            -0.3191|
|         10|       0.7427|       |            -0.4937|             1.9800|
|         11|       0.0550|       |            -1.1644|             1.2922|
|         12|       1.5387|*      |             0.3053|             2.7847|
|         13|       0.5285|       |            -0.7050|             1.7992|
|         14|       1.6737|*      |             0.4627|             2.9168|
|         15|      -1.2059|       |            -2.4699|             0.0420|
|         16|       1.3630|*      |             0.0749|             2.6498|
|         17|      -0.3194|       |            -1.5840|             0.9252|
|         18|       1.4112|*      |             0.1530|             2.6694|
|         19|       1.2892|*      |             0.0516|             2.5791|
|         20|       0.4588|       |            -0.7708|             1.6635|

The same process can then be followed to return the cumulative effect estimates 
for each post-intervention period along with their credible intervals. 

```r
## Show Cumulative ATT Estimates from BSTS CausalImpact
cumutbl <- as_tibble(impact9$series)
cumutbl$event.time <- (1:npds) - intpd
cumutbl$`(sig.)` <- apply(cumutbl[,c('cum.effect.lower','cum.effect.upper')],
                       1, function(x){ ifelse(prod(x)>0,'*','') })
cumucols <- c('event.time','cum.effect','(sig.)',
             'cum.effect.lower', 'cum.effect.upper')
cumutbl <- cumutbl[,cumucols]
knitr::kable(cumutbl[(intpd-1):(intpd+20), ], n = 22, digits = 4,
             caption = 'Cumulative ATT Estimates from BCA')
```



Table: Cumulative ATT Estimates from BCA

| event.time| cum.effect|(sig.) | cum.effect.lower| cum.effect.upper|
|----------:|----------:|:------|----------------:|----------------:|
|         -1|     0.0000|       |           0.0000|           0.0000|
|          0|     0.0557|       |          -1.2462|           1.3260|
|          1|     0.2387|       |          -1.5820|           1.9912|
|          2|    -0.4115|       |          -2.6000|           1.7769|
|          3|     0.0538|       |          -2.4711|           2.5264|
|          4|     0.4631|       |          -2.4530|           3.2917|
|          5|     1.3554|       |          -1.8802|           4.5109|
|          6|     1.9660|       |          -1.4659|           5.3508|
|          7|     1.8812|       |          -1.8344|           5.5689|
|          8|     2.3513|       |          -1.5060|           6.2508|
|          9|     0.7854|       |          -3.4430|           4.9100|
|         10|     1.5281|       |          -2.9237|           5.8303|
|         11|     1.5831|       |          -3.0610|           6.1036|
|         12|     3.1217|       |          -1.8053|           7.9545|
|         13|     3.6502|       |          -1.4631|           8.6746|
|         14|     5.3239|       |          -0.0205|          10.4937|
|         15|     4.1180|       |          -1.3628|           9.5098|
|         16|     5.4810|       |          -0.2301|          11.0989|
|         17|     5.1616|       |          -0.6995|          10.9791|
|         18|     6.5728|*      |           0.3602|          12.6706|
|         19|     7.8620|*      |           1.4332|          14.1705|
|         20|     8.3208|*      |           1.7723|          14.6759|
This concludes the vignette illustrating how to use BSTS as the predictive model
in a Bayesian counterfactual approach to dynamic causal inference. 

Please note that references not yet included will be added to a future version of 
this document. 
