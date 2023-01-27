---
title: "A Bayesian Counterfactual Approach to Dynamic Causal Inference: R Code and Tutorial"
author: 
  - "Anon.1^[redacted for review]"
  - "Anon.2^[redacted for review]"
email: ""
date: "January 27, 2023"
output:
  pdf_document:
    keep_md: true
    toc: true
---




\newpage

This is a tutorial and illustration for (1) building a Bayesian structural 
time series (BSTS) model [Section 1](#Section1) and (2) using the BSTS model as a
counterfactual for causal inference [Section 2](#Section2).

*  [**Section 1**](#Section1) includes complexity of the Bayesian approach. The basic structure of a BSTS model is local trend (state space) + seasonality (state space) + regression. Each component has several ways to define, and prior settings are slightly different. Thus, we aim to show step-by-step process of adding different components in Section 1 accompanied by checking resulting BSTS models, setting up Bayesian priors,and adding regression component (Spike and Slab priors included). This part requires bsts package in R. 
* [**Section 2**](#Section2) shows how to assess causal effects using the BSTS model from [Section 1] as a counterfactual, and is a straightforward process using `CausalImpact` packgage in R. 


# 0. Simulated Data

First load the data simulated from a data generating process with 
a low prior standard deviation weight (`dgp.prior.sd.weight.low=0.01`).
This data set includes the outcome series `y_observed`, and nine 
covariate series, `c1-c9`, from which the BSTS models that have regression
components will select their predictors

Note that the width limit prevents all columns from displaying, but the data 
`tibble` object summary displays the total number of rows and truncated
series with their corresponding data classes.


```r
## load BSTS library
library(bsts)
library(CausalImpact)
library(tibble) ## load a data table class with extra capabilities
## LOAD FROM FILE
df1 <- as_tibble(read.csv(dat1$bsts.df1.filepath)) ## convert dataframe to tibble
print(df1)
```

```
## # A tibble: 520 x 10
##    y_treatment  c1_mean   c2_mean  c3_mean  c1_sd  c2_sd   c3_sd c1_skew c2_skew
##          <dbl>    <dbl>     <dbl>    <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>
##  1      0.0879 -0.00141 -0.00333   0.0211  0.0430 0.0476 0.0110  -0.455   0.390 
##  2     -0.158   0.00148 -0.000158  0.0192  0.0397 0.0449 0.0115   0.255   0.220 
##  3      0.139   0.0382   0.00620  -0.0293  0.0419 0.0550 0.0111   0.0851  0.582 
##  4      0.725   0.0337   0.00613   0.00798 0.0464 0.0447 0.00874 -0.233  -0.0608
##  5      0.889   0.0697   0.00332   0.147   0.0572 0.0465 0.0114   0.588  -0.203 
##  6      1.59    0.0565   0.00410   0.185   0.0570 0.0479 0.00962  0.175  -0.0117
##  7      1.35    0.0915  -0.00824   0.312   0.0387 0.0530 0.0101   0.250   0.916 
##  8      2.19    0.0746  -0.00404   0.280   0.0487 0.0549 0.00971 -0.143  -0.780 
##  9      3.23    0.105    0.00477   0.401   0.0546 0.0392 0.00910  0.525  -1.01  
## 10      3.91    0.101    0.00310   0.662   0.0438 0.0597 0.00986 -0.110  -0.0513
## # ... with 510 more rows, and 1 more variable: c3_skew <dbl>
```
We can plot this data to visualize the post-intervention change in the 
outcome series `y_treatment`.

```r
npds <- nrow(df1)
intpd <- round(npds * 0.6)  # time period to split series (60% pre-, 40% post-)
## PLOT
plot(df1$y_treatment, main='Simulated Outcome Time Series (Y)',ylab='Y',xlab='Time')
abline(v=intpd, lty=2)
legend('topleft',legend=c('observed','intervention'), pch=c(1,NA),lty=c(NA,2))
```


\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/outcome_plot-1} 
\newline

The question of interest for research and practice is whether there
is a significant effect caused by the intervention at time point
`t=312`. We will apply BSTS for counterfactual causal inference to 
answer this question in a manner robust to the shape of the causal cycle
curve post-intervention. 

For this purpose, we first build a BSTS model based on observed data before the intervention, 
and the prediction of the BSTS model after the intervention will serve as a counterfactual to assess the causal effects of the intervention. 

Thus, we create the pre-intervention outcome series `y.pre.treat.NAs.post`
to pass into the `bsts()` function. 
Here we replace the post-intervention periods with `NA`'s that 
represent missing values within the `R` language. This is because the 
BSTS model is trained on the pre-intervention outcome (y) data. 
We also define covariates as predictors to add regression later (section 1.4)

```r
## INDICES OF PRE-INTERVENTION and POST-INTERVENTION WINDOWS
pre.idx <- 1:(intpd-1)
post.idx <- intpd:npds
## Outcome (response) series
y <- df1$y_treatment
## Train on y pre-treatment but NA's post-treatment
y.pre.treat.NAs.post <- c( y[pre.idx], rep(NA,length(post.idx)) )
## Then use the post-treatment response for causal impact estimation
post.period.response <- y[post.idx]
## Covariates (predictors) - data for the "formula = y ~ predictors" argument
predictors <- df1[ , ! names(df1) %in% 'y_treatment']
```


# 1. Bayesian Structural Time Series (BSTS) Modeling {#Section1}

## 1.1. State Space Specification: Local Trend Only

The simplest BSTS model is with only the local level in the state space. 
There are several state spaces in the bsts package for defining local trend including local level, 
local linear trend, student loca linear trend, and generalized local linear trend.

In this section, we start with the simplest local trend, local level 
(see section 1.6. for other types of local trend).

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
bsts.niter <- 150 ## suggest on the order of 1k-10k at least
```
Then, fit the BSTS model by calling the `bsts()` function to run MCMC estimation 
for `bsts.niter` draws from the stationary distribution of the Markov chain that 
specifies the posterior distribution (assuming convergence was reached),
which will be addressed below [Section 1.6](#Section1-6))

```r
## Fit BSTS model
bsts.fit <- bsts(y.pre.treat.NAs.post,
                 state.specification = st.sp,
                 niter = bsts.niter,
                 seed = rand.seed)
```

```
## =-=-=-=-= Iteration 0 Fri Jan 27 08:56:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 15 Fri Jan 27 08:56:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 30 Fri Jan 27 08:56:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 45 Fri Jan 27 08:56:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 60 Fri Jan 27 08:56:25 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 75 Fri Jan 27 08:56:26 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 90 Fri Jan 27 08:56:26 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 105 Fri Jan 27 08:56:26 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 120 Fri Jan 27 08:56:26 2023
##  =-=-=-=-=
## =-=-=-=-= Iteration 135 Fri Jan 27 08:56:26 2023
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
## [1] 0.3113475
## 
## $prediction.sd
## [1] 0.5404869
## 
## $rsquare
## [1] 0.9614202
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
project online repository. 

```r
check <- plotBstsStateComps(bsts.fit, intpd=intpd, return.val = T) 
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-6-1} \end{center}
This model appears to match the data with the trend only. 
However, this is only accurate for short (e.g., one-step ahead) predictions 
within the training window. In terms of capturing the structure of the time 
series, this model is underspecified and, essentially, overfitted to the 
training data. While the predicted values (blue lines) seem accurate
as one-step ahead, in-sample predictions, the model would lose the seasonal 
information outside the training window (post intervention), 
leaving the counterfactual series
devoid of seasonality when used for causal inference computations (see below).

Therefore, it is crucial to capture seasonality in the outcome time series 
and incorporate it into the BSTS model that will be used as the 
counterfactual for causal inference with the `CasualImpact()` function in
[Section 2](#2.). 

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
## [1] 0.5153044
## 
## $prediction.sd
## [1] 0.5753828
## 
## $rsquare
## [1] 0.8943189
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
decreased MAE). However, the MAE of this model is now 0.567.

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
covariate series. But before moving on to addressing regression in BSTS
models, we first detail the process of parameter tuning. 

A crucial step in machine learning analysis is 
hyperparameter tuning, and in our ivestigation of dynamic causal inference,
adjusting the values of Bayesian priors plays an important part in achieving 
adequate model fit for different time series structures that should 
be able to mimic the obsered outcome series pre-intervention.  
Therefore, this is an important step for researchers 
to incorporate information from theory and context into their empirical
analyses--within a Bayesian framework. 

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
                                  upper.limit= ))
```
Most of the particulars of prior elicitation we necessarily bound outside the 
scope of this tutorial vignette. Beginning at the assumption that the 
researcher has theoretical and/or contextual knowledge that can inform  
prior distribution attributes, we focus next on illustrating how to adjust 
the prior distribution settings for the `bsts()` function. 

In order to simplify the logic for specific values, we set parameter settings
for prior variances in terms of multiples of the empirical 
standard deviation of the observed outcome (y) in the pre-intervention
window. 

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
    sigma.guess = 0.01 * y.sd,      ##try a low weight of empirical y SD
    sample.size = round(0.1 * (intpd-1)), ##prop. empirical sample size (pre-int.)
    upper.limit= 1.5 * y.sd ##try a multiple of empirical y SD to set upper limit
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
## [1] 0.5157383
## 
## $prediction.sd
## [1] 0.5793489
## 
## $rsquare
## [1] 0.8941409
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
The MAE of this model is now 0.569.

This model does not appear to improve (in terms of MAE),
but the model now has seasonal structure, in addition to the local level
captured in the trend component, which collectively would inform a better
counterfactual post-intervention. So, essentially, the addition of seasonal
structure reduces a portion of the overfitting against current data
in order to incorporate a temporal structure (seasonality) 
that can continue into the post-intervention forecast.
THis, in turn, will be used as the counterfactual for causal inference. 

Beyond seasonality, another way to incorporate structure from 
the training window into the post-intervention window
is including a regression component in the state space. Assuming the 
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

## 1.4. Regression: Adding Covariates to State Space with Trend + Seasonality

### 1.4.1. Regression `formula` in BSTS function

This adds a static regression component (i.e., time-invariant betas) with
covariate series in the `predictors` dataframe, which is used as
the `data` input for the `bsts()` function. 

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
## [1] 0.5198743
## 
## $prediction.sd
## [1] 0.5419406
## 
## $rsquare
## [1] 0.8924362
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.228   1.000   3.000 
## 
## $coefficients
##                      mean          sd     mean.inc      sd.inc   inc.prob
## c3_mean      0.9702949528 0.244941401   0.97029495  0.24494140 1.00000000
## c1_sd       -0.1899909468 0.860120223  -3.09413828  1.85729796 0.06140351
## c1_mean     -0.0055067622 0.029164132  -0.10462848  0.08232379 0.05263158
## c3_sd       -0.4564851245 2.940762141 -10.40786084 10.71646624 0.04385965
## c2_mean      0.0040755896 0.245088506   0.09292344  1.29869610 0.04385965
## c2_sd        0.0433002187 0.797490355   2.46811246  7.73935070 0.01754386
## c2_skew      0.0007500031 0.008007842   0.08550035  0.00000000 0.00877193
## c3_skew      0.0000000000 0.000000000   0.00000000  0.00000000 0.00000000
## c1_skew      0.0000000000 0.000000000   0.00000000  0.00000000 0.00000000
## (Intercept)  0.0000000000 0.000000000   0.00000000  0.00000000 0.00000000
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
The MAE of this model is now 0.536.

### 1.4.2.  Variable Selection: Spike-and-Slab Priors

Having proper variables in regression is a crucial step in building a counterfactual time series
that sufficiently approximates (i.e., correlates to a sufficient degree with) 
the outcomes series in the pre-intervention window. 

One of the advantages in BSTS model is that variable selection among covariates is incorporated into the model within a Bayesian framework via spike-and-slab priors. 

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
##TODO
st.sp5 <- st.sp4
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
## [1] 0.5174937
## 
## $prediction.sd
## [1] 0.5611402
## 
## $rsquare
## [1] 0.893419
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    0.00    1.00    1.00    1.21    2.00    4.00 
## 
## $coefficients
##                      mean          sd    mean.inc      sd.inc   inc.prob
## c3_mean      0.5265388462 0.557109935  0.96397112  0.37945581 0.54621849
## c2_mean     -0.8658933995 1.865998649 -3.22004108  2.32999890 0.26890756
## c1_sd       -1.1093586220 3.745155224 -8.80091173  6.75964888 0.12605042
## c1_mean      0.0099162798 0.049549677  0.07866915  0.12205946 0.12605042
## c3_sd        1.1745230546 5.335747539 19.96689193 11.05220509 0.05882353
## c2_sd        0.2150401946 2.000977680  4.26496386  8.58865425 0.05042017
## c3_skew      0.0002055999 0.008852245  0.01223319  0.09461751 0.01680672
## c1_skew      0.0015905270 0.012286315  0.09463636  0.01419442 0.01680672
## c2_skew      0.0000000000 0.000000000  0.00000000  0.00000000 0.00000000
## (Intercept)  0.0000000000 0.000000000  0.00000000  0.00000000 0.00000000
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
st.sp6 <- st.sp5
#=================
# Option 2. Create SpikeSlabPrior object and pass it into BSTS function
#-----------------
##----------------
#  2A. Set expected.model.size in SpikeSlabPrior (with other options)
priorObjA <- SpikeSlabPrior(
  x = model.matrix(y.pre.treat.NAs.post ~ ., data=predictors),
  y = y.pre.treat.NAs.post,
  expected.model.size = 3, ## size in SpikeSlabPrior
  expected.r2 = .9,
  prior.df = .01,## obs. count
  prior.information.weight = .01
) 

# Run the bsts model with the SpikeSlabPrior object ‘priorObjA’
bsts.fit6 <- bsts(y.pre.treat.NAs.post ~ ., 
                  state.specification = st.sp6,
                  data = predictors,
                  niter = bsts.niter,
                  prior= priorObjA,
                  seed = rand.seed, ping = 0)
getBstsSummaryPlots(bsts.fit6)
```

```
## $residual.sd
## [1] 0.5184913
## 
## $prediction.sd
## [1] 0.537981
## 
## $rsquare
## [1] 0.8930077
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.077   1.000   2.000 
## 
## $coefficients
##                      mean          sd    mean.inc    sd.inc    inc.prob
## c3_mean      1.1006345237 0.239127636  1.10063452 0.2391276 1.000000000
## c1_sd       -0.0552438459 0.703233643 -2.15450999 4.6796200 0.025641026
## c2_sd        0.0441181395 0.420302210  2.58091116 2.7198131 0.017094017
## c2_skew     -0.0002672301 0.002890535 -0.03126592 0.0000000 0.008547009
## c1_skew     -0.0003821741 0.004133845 -0.04471437 0.0000000 0.008547009
## c3_sd        0.0480584575 0.519831699  5.62283953 0.0000000 0.008547009
## c2_mean      0.0029317074 0.031711264  0.34300977 0.0000000 0.008547009
## c3_skew      0.0000000000 0.000000000  0.00000000 0.0000000 0.000000000
## c1_mean      0.0000000000 0.000000000  0.00000000 0.0000000 0.000000000
## (Intercept)  0.0000000000 0.000000000  0.00000000 0.0000000 0.000000000
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
The MAE of this model is now 0.532.

Finally, the other way to create a `SpikeSlabPrior` object to pass
into the `bsts()`function is by specifying prior inclusion probabilities
(and estimates) for each covariates series in `predictors`. This is useful
when the researcher has prior knowledge that certain covariates are 
more likely to be important contributors to the counterfactual series
and should therefore be included in the model in more iterations 
(i.e., more MCMC draws from the posterior predictive distribution). 

```r
st.sp7 <- st.sp6
##----------------
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
                  prior= priorObjB,
                  seed = rand.seed, ping = 0)
getBstsSummaryPlots(bsts.fit7)
```

```
## $residual.sd
## [1] 0.519857
## 
## $prediction.sd
## [1] 0.5415484
## 
## $rsquare
## [1] 0.8924434
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   2.000   2.000   2.068   2.750   4.000 
## 
## $coefficients
##                      mean          sd    mean.inc     sd.inc    inc.prob
## c3_mean      0.9284980629 0.571035717  1.19089969  0.3214343 0.779661017
## c2_mean     -2.2882099600 3.229544122 -3.37510969  3.4252087 0.677966102
## c1_mean      0.0663305610 0.119678852  0.12041548  0.1398923 0.550847458
## c1_sd       -0.0649599689 0.694239820 -2.55509211  4.3187275 0.025423729
## c3_sd        0.1580007554 1.730980080  9.32204457 13.4056301 0.016949153
## c2_skew      0.0002896412 0.003146309  0.03417766  0.0000000 0.008474576
## c2_sd        0.0148515222 0.161328826  1.75247962  0.0000000 0.008474576
## c3_skew      0.0000000000 0.000000000  0.00000000  0.0000000 0.000000000
## c1_skew      0.0000000000 0.000000000  0.00000000  0.0000000 0.000000000
## (Intercept)  0.0000000000 0.000000000  0.00000000  0.0000000 0.000000000
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
The MAE of this model is now 0.539.

## 1.5. Compare Local Trend Models 

Besides the default state space containing only a local level component, 
we can evaluate different state space configurations, such as those
including a local linear trend and semi-local linear trend.

### 1.5.0. Local Level (Trend without Slope)

This was the simple model already introduced above. As a state space without
a linear trend, this model will serve as a baseline for comparison 
with other state space configurations that contain a type of linear trend. 

### 1.5.1. Local linear trend (Trend with Level and Slope)

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
    sigma.guess = 0.01 * y.sd,    ## guess
    sample.size = round( 0.1 * (intpd-1)),    ## default
    upper.limit = 1.5 * y.sd
  ), 
  initial.state.prior = NormalPrior(
    mu = 0,   ##guess
    sigma = 1 * y.sd  ## guess
  )
)
## Fit BSTS model
bsts.fit8 <- bsts(y.pre.treat.NAs.post ~ . ,
                 state.specification = st.sp8,
                 data = predictors, 
                 niter = bsts.niter, 
                 expected.model.size=3, ## instead of prior = priorObjB
                 seed=rand.seed, ping = 0)

## CALL OUR CUSTOM BSTS SUMMARY
getBstsSummaryPlots(bsts.fit8)
```

```
## $residual.sd
## [1] 0.5283619
## 
## $prediction.sd
## [1] 0.5227257
## 
## $rsquare
## [1] 0.8888953
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   1.000   1.143   1.000   2.000 
## 
## $coefficients
##                      mean          sd    mean.inc     sd.inc   inc.prob
## c3_mean      1.7009710911 0.212855283  1.70097109 0.21285528 1.00000000
## c1_mean      0.0063478220 0.028245948  0.09627530 0.06302103 0.06593407
## c2_mean     -0.0755205215 0.603263145 -2.29078915 2.96050992 0.03296703
## c3_skew     -0.0015133685 0.014436615 -0.13771653 0.00000000 0.01098901
## c2_skew     -0.0001251624 0.001193973 -0.01138978 0.00000000 0.01098901
## c2_sd        0.0236920613 0.226007861  2.15597758 0.00000000 0.01098901
## c1_sd       -0.0071127173 0.067850998 -0.64725727 0.00000000 0.01098901
## c1_skew      0.0000000000 0.000000000  0.00000000 0.00000000 0.00000000
## c3_sd        0.0000000000 0.000000000  0.00000000 0.00000000 0.00000000
## (Intercept)  0.0000000000 0.000000000  0.00000000 0.00000000 0.00000000
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
The MAE of this model is now 0.516.

### 1.5.2. Semilocal linear trend (Trend with Level and slope with AR(1) drift)

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
    sigma.guess = .01 * y.sd,    ## guess
    sample.size = round( 0.1 * (intpd-1)),    ## default
    upper.limit= 1.5 * y.sd
  ), 
  initial.state.prior = NormalPrior(
    mu= 0,   
    sigma = 1 * y.sd,  
  )
)
## Fit BSTS model: formula “y.pre ~ .” indicates inclusion of regression  
## component (all covariates in df1[,-1] dataframe)
bsts.fit9 <- bsts(y.pre.treat.NAs.post ~ . ,
                 state.specification = st.sp9,
                 data = predictors,  ## -1 excludes the outcome in 1st column 
                 niter = bsts.niter, 
                 # prior = priorObjB,
                 expected.model.size=3,
                 seed=rand.seed, ping=0)
## Call summary plots and visualize state components together
getBstsSummaryPlots(bsts.fit9)
```

```
## $residual.sd
## [1] 0.5381207
## 
## $prediction.sd
## [1] 0.5131102
## 
## $rsquare
## [1] 0.8847532
## 
## $relative.gof
## [1] NA
## 
## $size
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   2.000   2.000   2.000   2.352   3.000   4.000 
## 
## $coefficients
##                      mean           sd    mean.inc      sd.inc    inc.prob
## c3_mean      1.7460597765  0.307125958  1.74605978  0.30712596 1.000000000
## c1_mean      0.8003134707  0.088194062  0.80031347  0.08819406 1.000000000
## c2_sd        3.0333586134  6.497127689 16.29025922  3.08687572 0.186206897
## c3_sd        5.4182602120 20.696614684 78.56477307 21.64097941 0.068965517
## c2_mean     -0.0502105148  0.491324114 -1.04007495  2.14300331 0.048275862
## c1_sd       -0.0767289382  1.078772423 -2.78142401  6.76974584 0.027586207
## c1_skew     -0.0005952569  0.007125428 -0.04315613  0.06031255 0.013793103
## c3_skew      0.0003626064  0.004366360  0.05257793  0.00000000 0.006896552
## c2_skew      0.0000000000  0.000000000  0.00000000  0.00000000 0.000000000
## (Intercept)  0.0000000000  0.000000000  0.00000000  0.00000000 0.000000000
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
The MAE of this model is now 0.508.

Zooming into the first 100 periods shows how the BSTS model learns the structure
of the time series through the first couple of cycles that train the 
the seasonality and regression structures into the fitted model (via
recursive least squares in the Kalman filter). 


```r
par(mfrow=c(1,2), mar=c(4,4,.5,1))
plotBstsStateComps(bsts.fit9, pd.ids=1:104, title.show=F) ## 2 cycles (52-pds)
plotBstsStateComps(bsts.fit9, pd.ids=1:520, title.show=F) ## 10 cycles (52-pds)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-15-1} \end{center}

## 1.6. Check MCMC Diagnostics: Convergence and Fit {#Section1-6}

To check the convergence of the model estimation procedure and inspect 
the model's state space component contributions, we wrap a custom 
function created for this analysis, `bstsPostPredChecks()`, 
inside a formatting function, `getBstsDiagnosticPlots()`. Calling this on
several candidate BSTS models `bsts.fit`,`bsts.fit2`,`bsts.fit8`,`bsts.fit9`,
allows for comparison. 

Draws from the posterior distribution of the Markov chain(s)--once
converged at the stationary distribution--are used to estimate the 
model prediction at each time period though 
Markov chain Monte Carlo (MCMC) estimation. However, when the Markov chains
have not yet converged (usually due to insufficient MCMC burn-in period
or sampled iterations), then the iterations cannot be treated as draws from
the stationary distribution to sample the posterior predicted distribution. 
Thus, the estimates may be unstable and/or biased, so users must assess 
convergence of the MCMC chains. 

> There are several ways to check the convergence of the MCMC including the trace plots, Heidelberger-Welch tests, and Geweke diagnostic tests. We may evalue autocorrelation using ACF/PACF autocorrelation plots and Durbin-Watson tests. Precision was estimated from mean absolute 1-step prediction errors and range. 

Similarly, there is need to assess how the fitted model fits the observed data, 
in terms of how well new samples simulated from the fitted model create 
simulated distribution(s) that resemble the observed value(s). For this purpose, 
we examine MCMC convergence trace figures and plot model fit 
checks for the posterior predictive distribution of the outcome
(Y, top row panels) and the standardized residuals (bottom row panels). 

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

>##TODO: ADD paragraph about THE MCMC CONVERGENCE CHECKS: Geweke, H&W stationarity CMV and halfwidth-mean ratio check


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
## Call the new diagnostics function on the candidate models:
## Baseline (Trend = local level)
getBstsDiagnosticPlots(bsts.fit)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-16-1} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=1.05, p=0.30)
## H&W Stationarity: PASS (CMV=0.16, p=0.36)
## H&W Halfwidth: FAIL (hw/mean=0.19 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=27.93, p=0.00)
## H&W Stationarity: PASS (CMV=4190.77, p=0.63)
## H&W Halfwidth: FAIL (hw/mean=1.97 < eps=0.10)
```

```r
## Seasonality Baseline
getBstsDiagnosticPlots(bsts.fit2)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-16-2} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=1.03, p=0.31)
## H&W Stationarity: PASS (CMV=0.06, p=0.84)
## H&W Halfwidth: PASS (hw/mean=0.04 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=14.68, p=0.00)
## H&W Stationarity: PASS (CMV=35.50, p=0.06)
## H&W Halfwidth: FAIL (hw/mean=0.16 < eps=0.10)
```

```r
## Regression, Seasonality, Trend = Local linear
getBstsDiagnosticPlots(bsts.fit8)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-16-3} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=0.70, p=0.48)
## H&W Stationarity: PASS (CMV=0.17, p=0.35)
## H&W Halfwidth: PASS (hw/mean=0.03 < eps=0.10)
## Std. Residuals:-----
## Geweke: FAIL (z=6.92, p=0.00)
## H&W Stationarity: PASS (CMV=0.41, p=0.07)
## H&W Halfwidth: FAIL (hw/mean=0.14 < eps=0.10)
```

```r
## Regression, Seasonality, Trend = Semilocal linear
getBstsDiagnosticPlots(bsts.fit9)
```



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-16-4} \end{center}

```
## MCMC Convergence Diagnostics:
## 
## Posterior Predictive:-----
## Geweke: PASS (z=0.31, p=0.76)
## H&W Stationarity: PASS (CMV=0.10, p=0.59)
## H&W Halfwidth: PASS (hw/mean=0.01 < eps=0.10)
## Std. Residuals:-----
## Geweke: PASS (z=1.18, p=0.24)
## H&W Stationarity: PASS (CMV=0.34, p=0.11)
## H&W Halfwidth: FAIL (hw/mean=14.82 < eps=0.10)
```
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
##                          Average       Cumulative      
## Actual                   7.7           1617.5          
## Prediction (s.d.)        6.7 (0.36)    1403.6 (74.78)  
## 95% CI                   [6.1, 7.5]    [1285.0, 1571.6]
##                                                        
## Absolute effect (s.d.)   1 (0.36)      214 (74.78)     
## 95% CI                   [0.22, 1.6]   [45.97, 332.6]  
##                                                        
## Relative effect (s.d.)   15% (5.3%)    15% (5.3%)      
## 95% CI                   [3.3%, 24%]   [3.3%, 24%]     
## 
## Posterior tail-area probability p:   0.03261
## Posterior prob. of a causal effect:  96.739%
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
## Actual                   7.7           1617.5          
## Prediction (s.d.)        6.6 (0.11)    1370.3 (22.42)  
## 95% CI                   [6.3, 6.8]    [1326.9, 1411.9]
##                                                        
## Absolute effect (s.d.)   1.2 (0.11)    247.3 (22.42)   
## 95% CI                   [0.98, 1.4]   [205.60, 290.6] 
##                                                        
## Relative effect (s.d.)   18% (1.6%)    18% (1.6%)      
## 95% CI                   [15%, 21%]    [15%, 21%]      
## 
## Posterior tail-area probability p:   0.00685
## Posterior prob. of a causal effect:  99.31507%
## 
## For more details, type: summary(impact, "report")
```
The `CausalImpact` package includes a convenience function to print a 
helpful (i.e., plain English) explanation of these Bayesian counterfactual 
analysis results from the BSTS model.

```r
cat(summary(impact8, 'report'))
```

```
## Analysis report {CausalImpact}
## 
## 
## During the post-intervention period, the response variable had an average value of approx. 7.74. By contrast, in the absence of an intervention, we would have expected an average response of 6.72. The 95% interval of this counterfactual prediction is [6.15, 7.52]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 1.02 with a 95% interval of [0.22, 1.59]. For a discussion of the significance of this effect, see below.
## 
## Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.62K. By contrast, had the intervention not taken place, we would have expected a sum of 1.40K. The 95% interval of this prediction is [1.28K, 1.57K].
## 
## The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +15%. The 95% interval of this percentage is [+3%, +24%].
## 
## This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (1.02) to the original goal of the underlying intervention.
## 
## The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.033). This means the causal effect can be considered statistically significant.
```

```r
cat(summary(impact9, 'report'))
```

```
## Analysis report {CausalImpact}
## 
## 
## During the post-intervention period, the response variable had an average value of approx. 7.74. By contrast, in the absence of an intervention, we would have expected an average response of 6.56. The 95% interval of this counterfactual prediction is [6.35, 6.76]. Subtracting this prediction from the observed response yields an estimate of the causal effect the intervention had on the response variable. This effect is 1.18 with a 95% interval of [0.98, 1.39]. For a discussion of the significance of this effect, see below.
## 
## Summing up the individual data points during the post-intervention period (which can only sometimes be meaningfully interpreted), the response variable had an overall value of 1.62K. By contrast, had the intervention not taken place, we would have expected a sum of 1.37K. The 95% interval of this prediction is [1.33K, 1.41K].
## 
## The above results are given in terms of absolute numbers. In relative terms, the response variable showed an increase of +18%. The 95% interval of this percentage is [+15%, +21%].
## 
## This means that the positive effect observed during the intervention period is statistically significant and unlikely to be due to random fluctuations. It should be noted, however, that the question of whether this increase also bears substantive significance can only be answered by comparing the absolute effect (1.18) to the original goal of the underlying intervention.
## 
## The probability of obtaining this effect by chance is very small (Bayesian one-sided tail-area probability p = 0.007). This means the causal effect can be considered statistically significant.
```
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



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-19-1} \end{center}



\begin{center}\includegraphics[width=1\linewidth]{bsts_illustration_Section4_vignette_files/figure-latex/unnamed-chunk-19-2} \end{center}
Finally, we can print the results table of pointwise estimates of the ATT, along 
with 95% Bayesian credible intervals from the posterior predictive distribution.
In this example, we take the `$series` attribute of the `impact9` object
and add columns for the event time and significance codes ('*') where
the pointwise interval does not contain zero. Due to the time series length,
we truncate the output within the interval from 1 period before the intervention 
to 15 periods after the intervention. 

```r
## Show Pointwise ATT Estimates from BSTS CausalImpact
restbl <- as_tibble(impact9$series)
restbl$event.time <- (1:npds) - intpd
restbl$signif <- apply(restbl[,c('point.effect.lower','point.effect.upper')],
                       1, function(x){ ifelse(prod(x)>0,'*',NA) })
rescols <- c('event.time','point.effect',
             'point.effect.lower', 'point.effect.upper','signif')
restbl <- restbl[,rescols]
# print(restbl[(intpd-1):(intpd+15),], n = 17)
knitr::kable(print(restbl[(intpd-1):(intpd+15),], n = 17))
```

```
## # A tibble: 17 x 5
##    event.time point.effect point.effect.lower point.effect.upper signif
##         <dbl>        <dbl>              <dbl>              <dbl> <chr> 
##  1         -1      -0.0695             -1.23               1.01  <NA>  
##  2          0      -0.0575             -1.10               0.949 <NA>  
##  3          1      -0.0381             -0.997              1.17  <NA>  
##  4          2      -0.0169             -1.02               1.01  <NA>  
##  5          3       0.223              -0.761              1.68  <NA>  
##  6          4      -0.0570             -1.24               1.14  <NA>  
##  7          5       0.253              -0.748              1.39  <NA>  
##  8          6      -0.604              -1.87               0.394 <NA>  
##  9          7      -0.159              -1.21               0.912 <NA>  
## 10          8      -0.230              -1.29               0.838 <NA>  
## 11          9       0.752              -0.465              1.81  <NA>  
## 12         10       0.615              -0.600              1.51  <NA>  
## 13         11       0.129              -0.987              1.44  <NA>  
## 14         12      -0.509              -1.44               0.540 <NA>  
## 15         13       1.36                0.372              2.42  *     
## 16         14       1.36                0.329              2.45  *     
## 17         15       0.296              -1.04               1.42  <NA>
```



| event.time| point.effect| point.effect.lower| point.effect.upper|signif |
|----------:|------------:|------------------:|------------------:|:------|
|         -1|   -0.0694561|         -1.2283218|          1.0142982|NA     |
|          0|   -0.0575052|         -1.0985913|          0.9485404|NA     |
|          1|   -0.0381200|         -0.9969796|          1.1650916|NA     |
|          2|   -0.0169254|         -1.0190324|          1.0145411|NA     |
|          3|    0.2228731|         -0.7605816|          1.6805695|NA     |
|          4|   -0.0570220|         -1.2402293|          1.1360634|NA     |
|          5|    0.2529257|         -0.7476692|          1.3893319|NA     |
|          6|   -0.6042383|         -1.8657878|          0.3940220|NA     |
|          7|   -0.1585199|         -1.2060858|          0.9119843|NA     |
|          8|   -0.2300821|         -1.2857383|          0.8381722|NA     |
|          9|    0.7518433|         -0.4650095|          1.8146801|NA     |
|         10|    0.6151029|         -0.6001188|          1.5083301|NA     |
|         11|    0.1291817|         -0.9871554|          1.4415724|NA     |
|         12|   -0.5086227|         -1.4416702|          0.5396444|NA     |
|         13|    1.3565324|          0.3718393|          2.4201996|*      |
|         14|    1.3573883|          0.3290755|          2.4501451|*      |
|         15|    0.2961276|         -1.0404902|          1.4204023|NA     |


```r
knitr::kable(print(restbl[(intpd-1):(intpd+15),], n = 17))
```

```
## # A tibble: 17 x 5
##    event.time point.effect point.effect.lower point.effect.upper signif
##         <dbl>        <dbl>              <dbl>              <dbl> <chr> 
##  1         -1      -0.0695             -1.23               1.01  <NA>  
##  2          0      -0.0575             -1.10               0.949 <NA>  
##  3          1      -0.0381             -0.997              1.17  <NA>  
##  4          2      -0.0169             -1.02               1.01  <NA>  
##  5          3       0.223              -0.761              1.68  <NA>  
##  6          4      -0.0570             -1.24               1.14  <NA>  
##  7          5       0.253              -0.748              1.39  <NA>  
##  8          6      -0.604              -1.87               0.394 <NA>  
##  9          7      -0.159              -1.21               0.912 <NA>  
## 10          8      -0.230              -1.29               0.838 <NA>  
## 11          9       0.752              -0.465              1.81  <NA>  
## 12         10       0.615              -0.600              1.51  <NA>  
## 13         11       0.129              -0.987              1.44  <NA>  
## 14         12      -0.509              -1.44               0.540 <NA>  
## 15         13       1.36                0.372              2.42  *     
## 16         14       1.36                0.329              2.45  *     
## 17         15       0.296              -1.04               1.42  <NA>
```



| event.time| point.effect| point.effect.lower| point.effect.upper|signif |
|----------:|------------:|------------------:|------------------:|:------|
|         -1|   -0.0694561|         -1.2283218|          1.0142982|NA     |
|          0|   -0.0575052|         -1.0985913|          0.9485404|NA     |
|          1|   -0.0381200|         -0.9969796|          1.1650916|NA     |
|          2|   -0.0169254|         -1.0190324|          1.0145411|NA     |
|          3|    0.2228731|         -0.7605816|          1.6805695|NA     |
|          4|   -0.0570220|         -1.2402293|          1.1360634|NA     |
|          5|    0.2529257|         -0.7476692|          1.3893319|NA     |
|          6|   -0.6042383|         -1.8657878|          0.3940220|NA     |
|          7|   -0.1585199|         -1.2060858|          0.9119843|NA     |
|          8|   -0.2300821|         -1.2857383|          0.8381722|NA     |
|          9|    0.7518433|         -0.4650095|          1.8146801|NA     |
|         10|    0.6151029|         -0.6001188|          1.5083301|NA     |
|         11|    0.1291817|         -0.9871554|          1.4415724|NA     |
|         12|   -0.5086227|         -1.4416702|          0.5396444|NA     |
|         13|    1.3565324|          0.3718393|          2.4201996|*      |
|         14|    1.3573883|          0.3290755|          2.4501451|*      |
|         15|    0.2961276|         -1.0404902|          1.4204023|NA     |


The same process can be followed to return the cumulative estimates.

```r
## Show Cumulative ATT Estimates from BSTS CausalImpact
cumutbl <- as_tibble(impact9$series)
cumutbl$event.time <- (1:npds) - intpd
cumutbl$signif <- apply(cumutbl[,c('cum.effect.lower','cum.effect.upper')],
                       1, function(x){ ifelse(prod(x)>0,'*','') })
cumucols <- c('event.time','cum.effect',
             'cum.effect.lower', 'cum.effect.upper','signif')
cumutbl <- cumutbl[,cumucols]
print(cumutbl[(intpd-1):(intpd+15),], n = 17)
```

```
## # A tibble: 17 x 5
##    event.time cum.effect cum.effect.lower cum.effect.upper signif
##         <dbl>      <dbl>            <dbl>            <dbl> <chr> 
##  1         -1     0                  0               0     ""    
##  2          0    -0.0575            -1.10            0.949 ""    
##  3          1    -0.0956            -1.43            1.70  ""    
##  4          2    -0.113             -1.61            2.07  ""    
##  5          3     0.110             -1.75            2.60  ""    
##  6          4     0.0533            -2.42            3.06  ""    
##  7          5     0.306             -2.60            3.25  ""    
##  8          6    -0.298             -3.10            2.59  ""    
##  9          7    -0.457             -3.66            2.88  ""    
## 10          8    -0.687             -4.08            2.79  ""    
## 11          9     0.0652            -3.64            3.67  ""    
## 12         10     0.680             -3.57            4.02  ""    
## 13         11     0.810             -3.29            4.30  ""    
## 14         12     0.301             -3.99            3.80  ""    
## 15         13     1.66              -2.98            5.34  ""    
## 16         14     3.01              -1.96            7.03  ""    
## 17         15     3.31              -1.89            7.81  ""
```
This concludes the vignette illustrating how to use BSTS as the predictive model
in a Bayesian counterfactual approach to dynamic causal inference. 


