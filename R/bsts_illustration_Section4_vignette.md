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

* [**Section 1**](#Section1) includes complexity of the Bayesian approach. The basic structure of a BSTS model is local trend (state space) + seasonality (state space) + regression. Each component has several ways to define, and prior settings are slightly different. Thus, we aim to show step-by-step process of adding different components in Section 1 accompanied by checking resulting BSTS models, setting up Bayesian priors,and adding regression component (Spike and slab priors included). This part requires the `bsts` package in R: [https://cran.r-project.org/web/packages/bsts/bsts.pdf](https://cran.r-project.org/web/packages/bsts/bsts.pdf). 
* [**Section 2**](#Section2) shows how to assess causal effects using the BSTS model from [Section 1] as a counterfactual, and is a straightforward process using the `CausalImpact` R package: [https://google.github.io/CausalImpact/CausalImpact.html](https://google.github.io/CausalImpact/CausalImpact.html)  


# 0. Simulated Data

First load the data simulated from a data generating process (DGP) with 

- a random walk with small step size in the trend (i.e., local level standard deviation [SD] = 0.01), 
- weekly seasonality with frequency of 1 per cycle (i.e., 52 weeks = 1 year), and 
- a positive linear trend over time.  

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
## =-=-=-=-= Iteration 0 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 20 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 40 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 60 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 80 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 100 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 120 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 140 Fri Jan 27 13:13:11 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 160 Fri Jan 27 13:13:12 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 180 Fri Jan 27 13:13:12 2023
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
## [1] 0.4124693
## 
## $prediction.sd
## [1] 0.5836572
## 
## $rsquare
## [1] 0.9479783
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

## 1.2. State-Space Specification: Local Trend + Seasonality

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
## [1] 0.543113
## 
## $prediction.sd
## [1] 0.6358756
## 
## $rsquare
## [1] 0.9098052
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
decreased MAE). However, the MAE of this model is now 0.598.

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

## 1.3. Bayesian Priors: Parameter Tuning

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
## [1] 0.5454359
## 
## $prediction.sd
## [1] 0.6385696
## 
## $rsquare
## [1] 0.909032
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

## 1.4. Covariate Information: Local Trend + Seasonality + Regression {#Section1-4}

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
## [1] 0.5548387
## 
## $prediction.sd
## [1] 0.6220116
## 
## $rsquare
## [1] 0.9058686
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.453   2.000   4.000 
## 
## $coefficients
##                      mean           sd     mean.inc     sd.inc    inc.prob
## c3_mean      4.013109e-01  0.376634370   0.65118366  0.2583713 0.616279070
## c1_mean      1.961679e-01  0.244927500   0.43257540  0.1721372 0.453488372
## c3_sd       -1.431977e+01 27.961814708 -58.64286655 24.4604120 0.244186047
## c2_mean     -1.471514e-01  0.763529036  -1.94692676  2.1240017 0.075581395
## c2_skew     -9.612169e-04  0.011391969  -0.05510977  0.0813619 0.017441860
## c1_sd       -2.283948e-02  0.651971382  -1.30946373  5.8151543 0.017441860
## c1_skew     -9.094739e-04  0.011905296  -0.07821475  0.1102000 0.011627907
## c2_sd        4.803638e-02  0.478399111   4.13112877  2.3238923 0.011627907
## c3_skew     -6.608201e-04  0.008666574  -0.11366106  0.0000000 0.005813953
## (Intercept)  0.000000e+00  0.000000000   0.00000000  0.0000000 0.000000000
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
The MAE of this model is now 0.59.

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
                 expected.model.size = 5,  ## argument passed to SpikeSlabPrior
                 niter = bsts.niter,
                 seed = rand.seed, ping=0)
summary(bsts.fit5)
```

```
## $residual.sd
## [1] 0.5573911
## 
## $prediction.sd
## [1] 0.6160306
## 
## $rsquare
## [1] 0.9050005
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.339   2.000   4.000 
## 
## $coefficients
##                      mean           sd     mean.inc      sd.inc    inc.prob
## c3_mean      0.6033311103  0.364283082   0.71886260  0.27322195 0.839285714
## c2_mean     -0.5248657079  1.767352286  -3.39143996  3.27940860 0.154761905
## c1_mean      0.0170138166  0.075944256   0.11909672  0.17085359 0.142857143
## c3_sd       -4.7466645798 20.576144027 -46.90821467 48.14222723 0.101190476
## c1_sd        0.1337004404  1.413421691   3.20881057  6.64020513 0.041666667
## c2_sd       -0.0218834521  1.170746667  -0.61273666  6.73388374 0.035714286
## c2_skew     -0.0001437257  0.002266009  -0.01207296  0.02386355 0.011904762
## c3_skew     -0.0012508593  0.016212990  -0.21014436  0.00000000 0.005952381
## c1_skew     -0.0005309230  0.006881548  -0.08919506  0.00000000 0.005952381
## (Intercept)  0.0000000000  0.000000000   0.00000000  0.00000000 0.000000000
```

```r
par(mfrow=c(2,1))
plot(bsts.fit5, 'coefficients', main='Expected Size = 5')
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
## [1] 0.5567967
## 
## $prediction.sd
## [1] 0.6162523
## 
## $rsquare
## [1] 0.905203
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.616   2.000   4.000 
## 
## $coefficients
##                      mean          sd     mean.inc      sd.inc   inc.prob
## c3_mean       0.546877379  0.37699480   0.74653103  0.21020723 0.73255814
## c1_mean       0.157981705  0.22348512   0.39380947  0.17709187 0.40116279
## c3_sd       -13.716912169 25.13947982 -50.19806155 21.86946402 0.27325581
## c2_sd         0.008857213  0.94842750   0.15234406  4.13116104 0.05813953
## c1_sd         0.073637540  1.14342382   1.58320712  5.40434737 0.04651163
## c2_mean      -0.040145345  0.71386408  -0.86312491  3.41131007 0.04651163
## c3_skew      -0.001124105  0.02374404  -0.04833652  0.17056570 0.02325581
## c1_skew      -0.002191650  0.01983895  -0.09424094  0.10424911 0.02325581
## c2_skew       0.001588694  0.01485798   0.13662768  0.02915062 0.01162791
## (Intercept)   0.000000000  0.00000000   0.00000000  0.00000000 0.00000000
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
The MAE of this model is now 0.586.

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
## [1] 0.5576675
## 
## $prediction.sd
## [1] 0.6286836
## 
## $rsquare
## [1] 0.9049063
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   2.000   1.788   2.000   4.000 
## 
## $coefficients
##                      mean           sd     mean.inc     sd.inc    inc.prob
## c2_mean     -2.2701261483  2.429915377  -2.96862650  2.3760633 0.764705882
## c1_mean      0.1424793147  0.188994330   0.26912759  0.1826422 0.529411765
## c3_mean      0.3089015471  0.458994248   0.84698811  0.3470606 0.364705882
## c3_sd       -6.5231955141 20.842354105 -65.23195514 22.7891799 0.100000000
## c1_sd        0.0232107708  0.615842454   1.97291552  7.5101305 0.011764706
## c2_skew      0.0004595991  0.005992439   0.07813185  0.0000000 0.005882353
## c1_skew     -0.0003123743  0.004072863  -0.05310363  0.0000000 0.005882353
## c2_sd        0.0186878593  0.243659874   3.17693608  0.0000000 0.005882353
## c3_skew      0.0000000000  0.000000000   0.00000000  0.0000000 0.000000000
## (Intercept)  0.0000000000  0.000000000   0.00000000  0.0000000 0.000000000
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
The MAE of this model is now 0.594.

## 1.5. Local Trends Comparison: Level vs. Linear Short-Term and Long(er)-Term

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
## [1] 0.5781393
## 
## $prediction.sd
## [1] 0.5976851
## 
## $rsquare
## [1] 0.8977964
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.215   1.000   2.000 
## 
## $coefficients
##                    mean         sd  mean.inc     sd.inc   inc.prob
## c3_mean     1.632767907 0.21602370 1.6327679 0.21602370 1.00000000
## c1_mean     0.034348221 0.08974432 0.2297037 0.09513022 0.14953271
## c2_mean     0.005417359 0.20378840 0.1932191 1.46514315 0.02803738
## c2_sd       0.043841506 0.31922227 2.3455206 0.06657041 0.01869159
## c1_sd       0.045826181 0.41939099 2.4517007 2.61672408 0.01869159
## c3_skew     0.000000000 0.00000000 0.0000000 0.00000000 0.00000000
## c2_skew     0.000000000 0.00000000 0.0000000 0.00000000 0.00000000
## c1_skew     0.000000000 0.00000000 0.0000000 0.00000000 0.00000000
## c3_sd       0.000000000 0.00000000 0.0000000 0.00000000 0.00000000
## (Intercept) 0.000000000 0.00000000 0.0000000 0.00000000 0.00000000
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
The MAE of this model is now 0.571.

### 1.5.2. Semilocal Linear Trend (Trend with Level and Slope with AR(1) Drift)

The semilocal linear trend component in the BSTS model state space is 
useful for longer-term forecasts wherein the linear trend persists longer, 
so tighter prediction bounds are reasonable over longer post-intervention 
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
    sigma = .001 * y.sd), 
  slope.ar1.prior=Ar1CoefficientPrior(
    mu = 0, 
    sigma = .001 * y.sd, 
    force.stationary = FALSE, 
    force.positive = FALSE),
  slope.sigma.prior=SdPrior(
    sigma.guess=.0001 * y.sd, 
    sample.size=round( 0.1 * (intpd-1))),
  initial.level.prior=NormalPrior(
    mu = 0,
    sigma = .001 * y.sd, 
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
## [1] 0.5835679
## 
## $prediction.sd
## [1] 0.5958242
## 
## $rsquare
## [1] 0.8958681
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   2.000   2.000   2.000   2.444   3.000   5.000 
## 
## $coefficients
##                      mean           sd     mean.inc      sd.inc   inc.prob
## c3_mean      1.3930235012  0.341652999   1.39302350  0.34165300 1.00000000
## c1_mean      1.1124745582  0.115471456   1.11247456  0.11547146 1.00000000
## c2_sd        2.4588586164  5.966762581  14.52263370  5.92124359 0.16931217
## c1_sd        1.7246746925  5.175021270  12.53705834  7.77266052 0.13756614
## c3_sd       -4.3073955093 17.720799082 -50.88110945 37.52462977 0.08465608
## c2_mean     -0.0481851007  0.468002096  -1.51783067  2.35766164 0.03174603
## c3_skew     -0.0009113768  0.009605044  -0.08612511  0.05163534 0.01058201
## c2_skew     -0.0006658868  0.006515292  -0.06292630  0.01203254 0.01058201
## c1_skew      0.0000000000  0.000000000   0.00000000  0.00000000 0.00000000
## (Intercept)  0.0000000000  0.000000000   0.00000000  0.00000000 0.00000000
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
The MAE of this model is now 0.556.

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

To check the convergence of the model estimation procedure and inspect 
the model's state space component contributions, we wrap a custom 
function created for this analysis, `bstsPostPredChecks()`, 
inside a formatting function, `getBstsDiagnosticPlots()`. Calling this on
several candidate BSTS models `bsts.fit`, `bsts.fit2`,
`bsts.fit8`, `bsts.fit9`, allows for comparison. 

Draws from the posterior distribution of the Markov chain(s)--once
converged at the stationary distribution--are used to estimate the 
model prediction at each time period though 
Markov chain Monte Carlo (MCMC) estimation. However, when the Markov chains
have not yet converged (usually due to insufficient MCMC burn-in period
or sampled iterations), then the iterations cannot be treated as draws from
the stationary distribution to sample the posterior predictive distribution. 
Thus, the estimates may be unstable and/or biased, so users must assess 
convergence of the MCMC chains. 

```
TODO: There are several ways to check the convergence of MCMC including the trace plots, Heidelberger-Welch tests, and Geweke diagnostic tests. We may evalue autocorrelation using ACF/PACF autocorrelation plots and Durbin-Watson tests. Precision was estimated from mean absolute 1-step prediction errors and range. 
```

Similarly, there is need to assess how the fitted model fits the observed data, 
in terms of how well new samples simulated from the fitted model create 
simulated distribution(s) that resemble the observed value(s). For this purpose, 
we examine MCMC convergence trace figures and plot model fit 
checks for the posterior predictive distribution of the outcome
(Y, top row panels) and the standardized residuals (bottom row panels). 

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
In the first column, the trace of the Markov chain (mean of draws per period) 
is shown for the posterior predictive distribution of the outcome (Y). 
For model fit, the observed and predicted values of Y are compared (top-center), 
where tighter alignment indicates better fit. 

Additionally, using the maximum value of Y as an auxiliary statistic, 
the Bayesian p-val presents the proportion 
of maximum Y values (from MCMC) that are larger than the observed largest 
value of Y. Smaller values of Bayesian p-val indicate worse fit, which offers a 
similar interpretation to the frequentist null hypothesis statistical test 
(NHST) under a null assumption of a well fitting model. 
Values of Bayesian p-val >= alpha (conventionally 0.05) would therefore not reject the 
null of good fit, whereas Bayesian p-val < alpha indicates problem(s) of fit to 
be addressed (e.g., by increasing MCMC iterations, or respecifying the state space).

```
##TODO: ADD explantion about THE MCMC CONVERGENCE CHECKS: Geweke, H&W stationarity CMV and halfwidth-mean ratio check
```
\newpage

```r
## Call the new diagnostics function on the candidate models:
## Baseline (Trend = local level)
getBstsDiagnosticPlots(bsts.fit)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-17-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=0.72, p=0.47)
## H&W Stationarity: PASS (CMV=0.05, p=0.87)
## H&W Halfwidth: FAIL (hw/mean=0.19 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=4.92, p=0.00)
## H&W Stationarity: PASS (CMV=7304.28, p=0.68)
## H&W Halfwidth: FAIL (hw/mean=1.09 < eps=0.10)
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
## Geweke: PASS (z=1.15, p=0.25)
## H&W Stationarity: PASS (CMV=0.06, p=0.81)
## H&W Halfwidth: PASS (hw/mean=0.03 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=3.32, p=0.00)
## H&W Stationarity: PASS (CMV=0.40, p=0.07)
## H&W Halfwidth: PASS (hw/mean=0.07 < eps=0.10)
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
## Geweke: FAIL (z=5.62, p=0.00)
## H&W Stationarity: PASS (CMV=0.33, p=0.11)
## H&W Halfwidth: PASS (hw/mean=0.02 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=7.40, p=0.00)
## H&W Stationarity: PASS (CMV=0.24, p=0.21)
## H&W Halfwidth: FAIL (hw/mean=0.13 < eps=0.10)
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
## Geweke: FAIL (z=2.35, p=0.02)
## H&W Stationarity: FAIL (CMV=0.77, p=0.01)
## H&W Halfwidth: NA (hw/mean=NA < eps=0.10)
## Std. Residuals:-----
## Geweke: PASS (z=0.54, p=0.59)
## H&W Stationarity: PASS (CMV=0.16, p=0.35)
## H&W Halfwidth: FAIL (hw/mean=0.22 < eps=0.10)
```


Note how each transition between diagnostic panels contrasted here
(from `bsts.fit` to `bsts.fit2`, from `bsts.fit2` to `bsts.fit8`, 
and `bsts.fit8` to `bsts.fit9` ) exhibits 
incremental improvement in the following visually discernible ways that also 
relate to quantifiable fit improvement (i.e., decreased MAE)

1. Adding seasonality in `bsts.fit2` begins to address the baseline model's
overly narrow prediction distribution 
2. Changing trend to local linear trend, and including regression, in `bsts.fit8` scales the predicted SD to the observed SD; however, the predictions are shifted (overestimated) and the standardized one-step
ahead prediction errors now exhibit non-normality and worsened autocorrelation. 
3. Changing trend to semilocal linear trend (and still including regression) in 
`bsts.fit9` simultaneously (a) fixes the bias in the predicted distribution
(realigned with observed distribution) and fixes the (b) non-normality
and (c) autocorrelation problems in the standardized one-step ahead prediction 
(i.e., forecast) error. 

With all convergence checks passed, and model diagnostics that can be 
quantitatively and qualitatively characterized as well fitting, 
we now have a BSTS model ready to serve as the counterfactual
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
## Prediction (s.d.)        7.4 (0.42)     1545.4 (87.14)  
## 95% CI                   [6.7, 8.3]     [1395.0, 1734.2]
##                                                         
## Absolute effect (s.d.)   0.67 (0.42)    140.13 (87.14)  
## 95% CI                   [-0.23, 1.4]   [-48.63, 290.5] 
##                                                         
## Relative effect (s.d.)   9.1% (5.6%)    9.1% (5.6%)     
## 95% CI                   [-3.1%, 19%]   [-3.1%, 19%]    
## 
## Posterior tail-area probability p:   0.05556
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
## Prediction (s.d.)        7.4 (0.13)     1555.7 (26.39)  
## 95% CI                   [7.2, 7.7]     [1506.3, 1607.4]
##                                                         
## Absolute effect (s.d.)   0.62 (0.13)    129.84 (26.39)  
## 95% CI                   [0.37, 0.86]   [78.14, 179.22] 
##                                                         
## Relative effect (s.d.)   8.3% (1.7%)    8.3% (1.7%)     
## 95% CI                   [5%, 12%]      [5%, 12%]       
## 
## Posterior tail-area probability p:   0.00526
## Posterior prob. of a causal effect:  99.47368%
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
## During the post-intervention period, the response variable had an average value of approx. 8.06. By contrast, in the absence of an intervention, we would have expected an average response of 7.44. The 95% interval of this counterfactual prediction is [7.21, 7.69]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 0.62 with a 95% interval of [0.37, 0.86]. For a discussion of the significance of this effect, see below.
## 
## Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.69K. By contrast, had the intervention not taken place, we would have expected a sum of 1.56K. The 95% interval of this prediction is [1.51K, 1.61K].
## 
## The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +8%. The 95% interval of this percentage is [+5%, +12%].
## 
## This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (0.62) to the original goal of the underlying intervention.
## 
## The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.005). This means the causal effect can be considered statistically significant.
```

In paragraph formatting, this report would look as follows:

> <br /><br />  <br /><br />  During the post-intervention period, the response variable had an average value of approx. 8.06. By contrast, in the absence of an intervention, we would have expected an average response of 7.44. The 95% interval of this counterfactual prediction is [7.21, 7.69]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 0.62 with a 95% interval of [0.37, 0.86]. For a discussion of the significance of this effect, see below.<br /><br />  <br /><br />  Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.69K. By contrast, had the intervention not taken place, we would have expected a sum of 1.56K. The 95% interval of this prediction is [1.51K, 1.61K].<br /><br />  <br /><br />  The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +8%. The 95% interval of this percentage is [+5%, +12%].<br /><br />  <br /><br />  This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (0.62) to the original goal of the underlying intervention.<br /><br />  <br /><br />  The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.005). This means the causal effect can be considered statistically significant.

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
             caption = 'Pointwise ATT Estimates of BCA')
```



Table: Pointwise ATT Estimates of BCA

| event.time| point.effect|(sig.) | point.effect.lower| point.effect.upper|
|----------:|------------:|:------|------------------:|------------------:|
|         -1|      -0.1914|       |            -1.5509|             1.2417|
|          0|       0.0898|       |            -1.0910|             1.4696|
|          1|       0.1666|       |            -0.9813|             1.5302|
|          2|      -0.6711|       |            -2.0022|             0.4851|
|          3|       0.4764|       |            -0.7808|             1.6049|
|          4|       0.4124|       |            -0.6920|             1.4517|
|          5|       0.8555|       |            -0.2296|             1.9719|
|          6|       0.5269|       |            -0.7821|             1.6037|
|          7|      -0.1411|       |            -1.2078|             1.1117|
|          8|       0.4512|       |            -0.7097|             1.4765|
|          9|      -1.6034|*      |            -2.9728|            -0.5743|
|         10|       0.8738|       |            -0.3568|             2.2070|
|         11|      -0.0412|       |            -1.2502|             1.0413|
|         12|       1.5437|*      |             0.3139|             2.9735|
|         13|       0.4274|       |            -0.6358|             1.4826|
|         14|       1.6559|*      |             0.4273|             2.8519|
|         15|      -1.3028|       |            -2.5969|             0.0492|
|         16|       1.4405|*      |             0.2926|             2.6916|
|         17|      -0.3101|       |            -1.4632|             0.9144|
|         18|       1.3983|*      |             0.2112|             2.6643|
|         19|       1.1887|       |            -0.0541|             2.2337|
|         20|       0.4133|       |            -0.8465|             1.4949|

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
             caption = 'Cumulative ATT Estimates of BCA')
```



Table: Cumulative ATT Estimates of BCA

| event.time| cum.effect|(sig.) | cum.effect.lower| cum.effect.upper|
|----------:|----------:|:------|----------------:|----------------:|
|         -1|     0.0000|       |           0.0000|           0.0000|
|          0|     0.0898|       |          -1.0910|           1.4696|
|          1|     0.2564|       |          -1.4664|           2.0464|
|          2|    -0.4147|       |          -2.4908|           1.4399|
|          3|     0.0617|       |          -2.6170|           2.2375|
|          4|     0.4741|       |          -2.2714|           3.1735|
|          5|     1.3296|       |          -2.0633|           4.3641|
|          6|     1.8565|       |          -1.6728|           5.3980|
|          7|     1.7154|       |          -2.2455|           5.6566|
|          8|     2.1666|       |          -2.0284|           5.9700|
|          9|     0.5632|       |          -4.4479|           5.0215|
|         10|     1.4370|       |          -3.1854|           6.2051|
|         11|     1.3958|       |          -3.4520|           5.9871|
|         12|     2.9396|       |          -1.6974|           7.9136|
|         13|     3.3670|       |          -1.4587|           7.7134|
|         14|     5.0229|       |          -0.0394|           9.5602|
|         15|     3.7201|       |          -1.3513|           8.5238|
|         16|     5.1607|       |          -0.4821|          11.1823|
|         17|     4.8505|       |          -1.0223|          11.4223|
|         18|     6.2489|       |          -0.6371|          13.3267|
|         19|     7.4376|*      |           0.6457|          14.6214|
|         20|     7.8509|*      |           0.9611|          15.5195|
This concludes the vignette illustrating how to use BSTS as the predictive model
in a Bayesian counterfactual approach to dynamic causal inference. 




