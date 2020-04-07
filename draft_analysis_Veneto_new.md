COVID19 - Forecast and predictions using a time dependent SEIR model for
the Veneto Region
================
Paolo Girardi
07 Aprile, 2020

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This
work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative
Commons Attribution-NonCommercial 4.0 International License</a>.

# Disclaimer

  - We want to investigate the evolution of the coronavirus pandemic in
    the Veneto Region from a statistical perspective using aggregated
    data.

  - Our point of view is that of surveillance with the goal of detecting
    important changes in the underlying (random) process as soon as
    possible after it has occured.

  - We use data provided by Italian Civil Protection Department and the
    analysis was restricted to the Veneto Region

  - This document is in a draft mode, and it is continuously updated.

  - The layout of the draft must definitely be improved.

## The COVID dataset

The present analysis started from the dataset on COVID19 updated in
<https://github.com/pcm-dpc/COVID-19>, database provided by the Italian
Civil Protection.

# Software

Install packages `dygraphs`, `xts` and `EpiDynamics` if not available

``` r
checkpackage <- function(package) {
  if (!package %in% installed.packages()) install.packages(package)
}
checkpackage("dygraphs")
checkpackage("xts")
checkpackage("EpiDynamics")
checkpackage("webshot")
checkpackage("bsts")
checkpackage("ggplot2")
checkpackage("knitr")
checkpackage("splines")
```

and load them.

``` r
library(dygraphs)
library(xts)
```

    ## Loading required package: zoo

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

``` r
library(EpiDynamics)
library(webshot)
library(bsts)
```

    ## Loading required package: BoomSpikeSlab

    ## Loading required package: Boom

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'Boom'

    ## The following object is masked from 'package:stats':
    ## 
    ##     rWishart

    ## 
    ## Attaching package: 'BoomSpikeSlab'

    ## The following object is masked from 'package:stats':
    ## 
    ##     knots

    ## 
    ## Attaching package: 'bsts'

    ## The following object is masked from 'package:BoomSpikeSlab':
    ## 
    ##     SuggestBurn

``` r
library(EpiDynamics)
library(ggplot2)
library(knitr)
library(splines)
```

# Source of the data

Download the data from

<https://github.com/pcm-dpc/COVID-19/>

# Results

## Load dataset

``` r
rm(list=ls())
###import italian updated dataset 
dat_csv<-read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv",header=T)
#we restricted the analysis to the Veneto Region
dat_csv<-dat_csv[dat_csv$codice_regione==5,]
days<-dim(dat_csv)[1]
dat_csv$t<-1:days
# The total number of epidemic day is
days
```

    ## [1] 44

Several outcomes can be potentially monitored, that is

``` r
names(dat_csv[,-c(1:6,18:20)])
```

    ##  [1] "ricoverati_con_sintomi"     "terapia_intensiva"         
    ##  [3] "totale_ospedalizzati"       "isolamento_domiciliare"    
    ##  [5] "totale_positivi"            "variazione_totale_positivi"
    ##  [7] "nuovi_positivi"             "dimessi_guariti"           
    ##  [9] "deceduti"                   "totale_casi"               
    ## [11] "tamponi"

It is worth noting that some outcomes present negative counts in some
regions. It looks like some of these negative counts are redesignations.
Outcomes presenting negative values cannot be analyzed using the
proposed model.

Then we extract the
timeseries.

``` r
myDateTimeStr1 <- paste(substr(dat_csv$data,1,10),substr(dat_csv$data,12,19))
myPOSIXct1 <- as.POSIXct(myDateTimeStr1, format="%Y-%m-%d %H:%M:%S")
days_dy<-as.Date(myPOSIXct1)
dat_csv_dy<-xts(dat_csv[,-c(1:6,18:20)], order.by = days_dy, frequency = 7)
```

``` r
p <- dygraph(dat_csv_dy,main=paste("Veneto Region",sep =""),xlab="Day",height=400,width=800) 
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
The time series of COVID-19 in the Veneto Region was influenced by the
presence of a hotspot in the small village of Vo
Euganeo.

``` r
plot(days_dy[-1],diff(dat_csv$totale_casi), ylab="New cases",xlab="Date",type="l")
```

![](draft_analysis_Veneto_new_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
\#\#\# The S(E)IR model (to be revised)

With the aim of predicting the future number of COVID19 cases on the
basis of the actual data, we used a SEIR model applied to the COVID19
epidemic to the Veneto Region

We will consider the classical [SIR
model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology)
\[@Kermack1927\].

The model divides a population of hosts into three classes: susceptible,
infected, recovered. The model describes how the portion of the
population in each of these classes changes with time. Births are
modeled as flows from “elsewhere” into the susceptible class; deaths are
modeled as flows from the \(S\), \(I\), or \(R\) compartment into
“elsewhere”. If \(S\), \(I\), and \(R\) refer to the numbers of
individuals in each compartment, then these **state variables** change
according to the following system of differential equations:
\[\begin{aligned}
\frac{d}{dt}S(t) &= B(t)-\lambda\,S(t)-\mu\,S(t)\\
\frac{d}{dt}I(t) &= \lambda\,S(t)-\gamma\,I(t)-\mu\,I(t)\\
\frac{d}{dt}R(t) &= \gamma\,I(t)-\mu\,R(t).\\
\end{aligned}\] Here, \(B\) is the crude birth rate (births per unit
time), \(\mu\) is the death rate and \(\gamma\) is the recovery rate.
We’ll assume that the force of infection, \(\beta\), for a constant
population \(N\) \[\lambda = \beta\,\frac{I}{N},\] so that the risk of
infection a susceptible faces is proportional to the *prevalence* (the
fraction of the population that is infected). This is known as the
assumption of frequency-dependent transmission.

# The reproduction number of COVID19.

The number of infected individuals \(I\) at time \(t\) is approximately
\[I(t)\;\approx\;I_0\,e^{(R_0-1)\,(\gamma+\mu)\,t}\] where \(I_0\) is
the (small) number of infectives at time \(0\), \(\frac{1}{\gamma}\) is
the infectious period, and \(\frac{1}{\mu}\) is the host lifespan.

\(R_0\) is the reproduction number
(<https://en.wikipedia.org/wiki/Basic_reproduction_number>) and
indicates how contagious an infectious disease is.

Taking logs of both sides, we get

\[\log{I}(t)\;\approx\;\log{I_0}+(R_0-1)\,(\gamma+\mu)\,t,\] which
implies that a semi-log plot of \(I\) vs \(t\) should be approximately
linear with a slope proportional to \(R_0\) and the recovery rate.

``` r
dat_csv_dy$log_totale_positivi<-log(dat_csv_dy$totale_positivi)
p <- dygraph(dat_csv_dy$log_totale_positivi,main=paste("Veneto Region"),ylab="Log Infected case",xlab="Day",height=400,width=800) 
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

We estimate the \(R_0\) parameter in the linear model.

\[
\log(I(t))= \beta_0 + \beta_1  t +e_t.
\]

The estimated slope coefficient \(\hat\beta_1\) is used to estimate
\(R_0\) as in the following formula:

\[
\widehat\beta_1=(\widehat{R_0}-1)\,(\gamma+\mu)
\] The parameter \(\mu\)\<\<\(\gamma\) and it can not be considered. As
consequence, R0 can be estimated as follows \[
\hat{R_0}=1+\frac{\hat{\beta_1}}{\gamma}
\] Respect to the SIR model, \(R_0\) can be estimated as follows: \[
\hat{R_0}=\frac{\beta}{\gamma}
\] And this was we can retrive the value of \(\beta\) in the SEIR model
by means of \[
\hat{R_0}=\frac{\hat{\beta}}{\gamma}=1+\frac{\hat{\beta_1}}{\gamma}\\
\hat{\beta}=(1+\frac{\hat{\beta_1}}{\gamma})*{\gamma}={\gamma}+\hat{\beta}_1
\] where \(\beta_1\) is the slope coefficient. \\

The incubation period for the coronavirus is in mean 5.1 days with a
range from 2-14 days. Please see
<https://www.worldometers.info/coronavirus/coronavirus-incubation-period/>.
However, the incubation period is used for epidemic diseases that causes
the immediate home isolation of infected subjects. The duration of the
diseases is about 2 weeks.

However, in the calculation of R0 we considered an infectious period of
18 days on the basis of recent publications
(<https://www.nature.com/articles/s41421-020-0148-0>).

We calculate several R0 values, each one based on a mobile window of 5
days, that can be sufficient to estimate a local trend, in order to
assess if the R0 trend is decreasing (how is expected to be). In this
way, the R0 for the first and the last two days of observation is
impossibile to estimate.

``` r
#calculate r0 based with a mobile window of 5 days
#vector for beta and standard deviation
duration<-18
beta1<-NULL
sd_beta1<-NULL
#for cycle for R0 estimates from days-2 to days+2
for (i in 3:(days-2)){
fit <- lm(log(totale_positivi)~t,data=dat_csv[(i-2):(i+2),])
beta1<-c(beta1,coef(fit)[2])
sd_beta1<-c(sd_beta1,coef(summary(fit))[2,2])
}

label<-as.Date(substr(dat_csv$data,1,10))[3:(days-2)]


mean  <- 1+(beta1*duration)
lower <- 1+((beta1-1.96*sd_beta1)*duration)
upper <- 1+((beta1+1.96*sd_beta1)*duration)

df <- data.frame(label, mean, lower, upper)


fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Date") + ylab("R0 Mean (95% CI)") +
  theme_bw() 
print(fp)
```

![](draft_analysis_Veneto_new_files/figure-gfm/R0%20trend-1.png)<!-- -->

The R0 shows a decreasing trend in the last period. We use the estimated
trend between R0 and time to calculate the future R0 value for the next
14 days. We predict beta (and R0) for the next 14 days by means of a
linear regressione model, assuming a Log-Normal distribution for the the
\(\beta1\) and as consequence for the \(\beta\) for the regression,
rescaled of \(1/\gamma\). We implement function to deal with log-normal
distribution, in order to recover the parameter \(\mu\) and \(\sigma\)
belonging normal distribution starting from this one:
\[\beta_{LogN}(t) \sim Log N (\beta_1(t)+\frac{1}{\gamma}, \sigma_{\beta_1(t)}^2)\]
In order to deal with estimations of \(\hat\beta_1(t)\) that we have
estimated before, we adopted an inverse transformation as reported here
<https://ocw.mit.edu/courses/civil-and-environmental-engineering/1-151-probability-and-statistics-in-engineering-spring-2005/lecture-notes/briefnts8_relnl.pdf>.
The following formulas provide us the mean and variance of the normal
distribution

\[\beta_{N}(t) \sim N (\mu_{\beta_N},\sigma_{\beta_N}^2)\] derived from
the log-normal distribution \(\beta_{LogN}(t)\).
\[\mu_{\beta_N} =e^{\mu_{\beta_{LogN}}+0.5*\sigma_{\beta_{LogN}}^2}\]

\[\sigma_{\beta_N}^2 =-2\log(\mu_{\beta_{LogN}})+\log(\mu_{\beta_{LogN}}^2+\sigma_{\beta_{LogN}}^2)\]
This transformation permits us to perform an inference based on the
classical regression model for the predict forward of the future values
of \(\beta_{LogN}(t)\) with a direct transformation.

``` r
# define direct transf function
M_from_Norm<-function(mu,sigma) exp(mu+0.5*sigma^2)
SD_from_Norm<-function(mu,sigma)  sqrt(exp(2*mu+sigma^2)*(exp(sigma^2)-1))
# define inverse transf function
M_from_LNorm<-function(mu,sigma) 2*log(mu)-0.5*log(mu^2+sigma^2)
SD_from_LNorm<-function(mu,sigma) sqrt( -2*log(mu)+log(mu^2+sigma^2))
```

``` r
#start from the day 7, since vo hotspot can influence the trend
time<-3:(days-2)
weekend<-rep(c(1,2,3,4,5,6,7),ceiling(days/7))[3:(days-2)]
data<-data.frame(time)
beta_norm<-M_from_LNorm(beta1+1/duration,sd_beta1)
sd_beta_norm<-SD_from_LNorm(beta1+1/duration,sd_beta1)
sd_norm_mean<-sqrt(mean(SD_from_LNorm(beta1,sd_beta1)^2))

beta.model<-glm(beta_norm~time,weights = 1/(sd_beta_norm)^2,family=gaussian,data)
summary(beta.model)
```

    ## 
    ## Call:
    ## glm(formula = beta_norm ~ time, family = gaussian, data = data, 
    ##     weights = 1/(sd_beta_norm)^2)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -7.7317  -1.2874   0.9356   2.9659   7.3324  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -1.022301   0.051449  -19.87   <2e-16 ***
    ## time        -0.035022   0.001541  -22.73   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 14.0055)
    ## 
    ##     Null deviance: 7766.68  on 39  degrees of freedom
    ## Residual deviance:  532.21  on 38  degrees of freedom
    ## AIC: -43.745
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
anova(beta.model,test="Chisq")
```

    ## Analysis of Deviance Table
    ## 
    ## Model: gaussian, link: identity
    ## 
    ## Response: beta_norm
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##      Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
    ## NULL                    39     7766.7              
    ## time  1   7234.5        38      532.2 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#there is an effect of the week, as supposed 
#we can use
#summary(gamlss(formula = beta_norm ~ time + bs(weekend, 3),Lp=1) )
forecast=14
# add 'fit', 'lwr', and 'upr' columns to dataframe (generated by predict)
pre<-predict(beta.model,type='response',newdata=data.frame(time=1:(days+forecast)),se.fit=TRUE)
date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+forecast-1+dim(dat_csv)[1],1)
predict <- data.frame(beta=c(rep(NA,2),beta_norm,rep(NA,forecast+2)),time=date,fit=pre$fit,lwr=pre$fit-1*1.96*pre$se.fit,upr=pre$fit+1*1.96*pre$se.fit)
beta.predict<-predict 
r0.predict<-beta.predict
r0.predict[,c(1,3:5)]<-M_from_Norm(r0.predict[,c(1,3:5)],sd_norm_mean)*duration
# plot the points (actual observations), regression line, and confidence interval
p <- ggplot(r0.predict, aes(date,beta))
p <- p + geom_point() +labs(x="Date",y="R0 value") 
p <- p + geom_line(aes(date,fit))
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3)
p
```

    ## Warning: Removed 18 rows containing missing values (geom_point).

![](draft_analysis_Veneto_new_files/figure-gfm/R0%20forecast-1.png)<!-- -->

R0 passes from a value of NA in the initial phase to an estimated value
of 0.85 at the ending of the 14-days forecast.  
We use the library(EpiDynamics) and the function SEIR() to implement a
SEIR model:

\[\begin{aligned}
\frac{d}{dt}S(t) &= \mu (N-S)-\beta \frac{SI}{N}-\nu S\\  
\frac{d}{dt}E(t) &= \beta \frac{SI}{N}-(\mu+\sigma) E    \\
\frac{d}{dt}I(t) &=\sigma E- (\mu+\gamma)I \\
\frac{d}{dt}R(t) &= \gamma\,I(t)-\mu\,R(t)+\nu S.\\
\end{aligned}\]

<img src="https://upload.wikimedia.org/wikipedia/commons/3/3d/SEIR.PNG"/>

Respect to the previous formulation, the compartment of exposed people
(E) was inserted between the susceptible and infected compartments.  
(<https://en.wikipedia.org/wiki/Bayesian_structural_time_series>). The
parameter \(\nu\) is the rate of vaccination.

We want to make a short term forecast (14 days).

We made a forecast by means of a SEIR model fixing a series of initial
status:  
*S:N, the size of the Veneto Region population  
*E: The number of exposed people, but supposed to be at least 3 times
the infected (3 x I\_start)  
*I\_start: initial number of COVID-19 cases  
*R\_start: initial number of recovered

and parameters:  
\-beta: gamma + slope  
\-gamma= 1/duration of diseases (duration=21 days)  
\-sigma: the coronavirus transmission rate  
\-mu0: the overall mortality rate

``` r
# initial number of infectus
I_start<-dat_csv$totale_positivi[dim(dat_csv)[1]]; I_start
```

    ## [1] 9965

``` r
# initial number of recovered, based on the proportion of discharged from the health services
prop<-dat_csv$dimessi_guariti[dim(dat_csv)[1]]/dat_csv$totale_ospedalizzati[dim(dat_csv)[1]]
R_start<-prop*dat_csv$totale_casi[dim(dat_csv)[1]]; R_start
```

    ## [1] 8041.111

``` r
# Veneto population
N=4905864
#  duration of COVID19 diseases 
duration<-18
#mortality rate 
mu0<-1/(82*365.25) # 1/lifespan
```

We try to estimates sigma (from exposed to infected) and the number of
exposed population on the basis of the last 5 days of observed number of
total infected people.

We suppose that the number of total infected \(I(t)\) can be the
realization of a Poisson random variable \[
I(t) \sim Poisson (\theta_t)
\] where \(\theta_t\) is the unknown mean of the random process at the
time \(t\).

The Poisson distribution has a fixed variance that is equal to the mean;
in the next part we try to estimates the parameters considering a
Negative Binomial distribution which allows for overdispersion.

However, the parameter \(\theta_t\) can be supposed to be a realization
of a SEIR model with initial parameters
\(\theta=(\mu,\beta,\sigma, \gamma)\) and initial status
\(S=(S_{start},E_{start},I_{start},R_{start})\).

The parameter \(\sigma\) is our parameter of interest; we estimate
\(\sigma\) fixing \(E_{start}= 3*I_{start}\) minimizing the log
likelihood of the described Poisson distribution (or Negative Binomial
one).

``` r
# number of the last considered days for calibration
last<-10
#I_start
I_1<-dat_csv$totale_positivi[(days-last)]
#R_start
R_1<-prop*dat_csv$totale_casi[(days-last)]
#loglikelihood function
LogLikfun <- function (initials,parameters,obs) {
n<-length(obs)
N=4905864
seir_fit <- SEIR(pars = parameters, init = initials, time = 0:last)
#Poisson
sum(dpois(x=obs[-1],lambda = seir_fit$results$I[-1]*N,log=TRUE))
#Neg Binomial
#sum(dnbinom(x=obs[-1],mu=seir_fit$results$I[-1]*N,size=N,log=TRUE))
#SSE
#-sum((obs[-1]-seir_fit$results$I[-1]*N)^2)
}
### logit and its inverse for sigma 
logit <- function (p) log(p/(1-p))    # the logit transform
expit <- function (x) 1/(1+exp(-x))   # inverse logit

f1<-function(par){
#max 0.20 according to https://www.nejm.org/doi/10.1056/NEJMoa2001316 is 1/5.2
sigma<-expit(par)
f_exp<-3*I_1/N
parameters <- c(mu = mu0, beta =(as.numeric(beta1[(days-last-2)])+1/duration), sigma = sigma, gamma =1/duration)
initials <- c(S = (1-(f_exp+I_1/N+R_start/N)), E = f_exp, I = I_1/N, R = R_1/N)
-LogLikfun(initials,parameters,dat_csv$totale_positivi[c((days-last):days)])
}
est<-optim(fn=f1,par=logit(0.1))
```

    ## Warning in optim(fn = f1, par = logit(0.1)): one-dimensional optimization by Nelder-Mead is unreliable:
    ## use "Brent" or optimize() directly

``` r
est
```

    ## $par
    ## [1] -3.269444
    ## 
    ## $value
    ## [1] 55.93137
    ## 
    ## $counts
    ## function gradient 
    ##       30       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

``` r
expit(est$par)
```

    ## [1] 0.03663443

``` r
f_exp<-3*I_1/N
sigma<-expit(est$par)

parameters <- c(mu = mu0, beta =(as.numeric(beta1[(days-last-2)])+1/duration), sigma = sigma, gamma = 1/duration)
initials <- c(S = (1-(f_exp+I_1/N+R_1/N)), E = f_exp, I = I_1/N, R = R_1/N)
pro<-SEIR(pars = parameters, init = initials, time = 0:last)
#predicted vs observed
cbind(pro$results$I*N,dat_csv$totale_positivi[c((days-last):days)])
```

    ##           [,1] [,2]
    ##  [1,] 6913.000 6913
    ##  [2,] 7278.078 7251
    ##  [3,] 7623.586 7564
    ##  [4,] 7951.848 7850
    ##  [5,] 8264.710 8224
    ##  [6,] 8564.425 8578
    ##  [7,] 8852.569 8861
    ##  [8,] 9130.654 9093
    ##  [9,] 9400.042 9409
    ## [10,] 9661.960 9722
    ## [11,] 9917.512 9965

The values of initial parameters are

``` r
#mu0
mu0;
```

    ## [1] 3.338842e-05

``` r
#sigma
sigma
```

    ## [1] 0.03663443

``` r
#gamma
1/duration
```

    ## [1] 0.05555556

``` r
#mthe unknowrn fraction of exposed people is
pro$results$E[last+1]
```

    ## [1] 0.004473631

``` r
#that is 
pro$results$E[last+1]/pro$results$I[last+1]
```

    ## [1] 2.212957

``` r
#times the infected
```

# Forecast 1 - Results based on beta coefficient trend

For the beta parameter, we perfomed a simulation on its trend by means
of a Bayesian Structural Time Series using the library bsts of R.  
We estimated a BSTS model specifing a local linear trend and a
seasonality component. The model predicted 1000 trajectories forward of
the beta coefficent.

We made a forecast forward of 14 days dropping the first 100 simulation
as burn-in.

``` r
# Bayesian Structural Time Series
#we assume beta ~ LN (mu_b,sigma_b)

#we drop the first 6 beta_norm
ss <- AddStudentLocalLinearTrend(list(), beta_norm[6:(days-4)])
ss <- AddSeasonal(ss, beta_norm[6:(days-4)], nseasons = 7)
model1 <- bsts(beta_norm[6:(days-4)],
               state.specification = ss,
               niter = 1000,seed=123)
```

    ## =-=-=-=-= Iteration 0 Tue Apr  7 22:41:39 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 100 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 200 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 300 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 400 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 500 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 600 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 700 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 800 Tue Apr  7 22:41:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 900 Tue Apr  7 22:41:40 2020 =-=-=-=-=

``` r
par(mfrow = c(1,1))
plot(model1, "components")
```

![](draft_analysis_Veneto_new_files/figure-gfm/bsts%20model%20fit%20and%20prediction-1.png)<!-- -->

``` r
#previsioni
pred1 <- predict(model1, horizon = 16, burn = 100)
plot(pred1 , ylab="Log Beta1 coefficient",main="Data and predictions",ylim=c(-3,3)*IQR(pred1$distribution)+median(pred1$distribution))
```

![](draft_analysis_Veneto_new_files/figure-gfm/bsts%20model%20fit%20and%20prediction-2.png)<!-- -->

``` r
# matrix of beta coefficients
coef<-M_from_Norm(pred1$distribution[,3:16],sd_norm_mean)
par(mfrow = c(1,1))
```

For each vector of simulated beta coefficients we perform a SEIR model.
We save the results and plot a credible interval at 90% (5-95%
percentile).

``` r
seir1_sim<-NULL
for(s in 1:dim(coef)[1]){
  # average number of single connections of an infected person
  # less contacts, less probability of new infections
  # we keep constant the other parameters
  forecast<-14
  seir1<-NULL
  for(i in 1:forecast){
    parameters <- c(mu = mu0, beta = matrix(coef[s,i]), sigma = sigma, gamma = 1/duration)
    if( i==1) initials <- c(S = 1-(pro$results$E[last+1]+I_start/N+R_start/N), E = pro$results$E[last+1], I = I_start/N, R = R_start/N)
    if( i>1) initials <- c(S = seir1_temp$results$S[2], E = seir1_temp$results$E[2], I =seir1_temp$results$I[2], R = seir1_temp$results$R[2])
    seir1_temp <- SEIR(pars = parameters, init = initials, time = 0:1)
    seir1 <- rbind(seir1,SEIR(pars = parameters, init = initials, time = 0:1)$results[2,])
  }
  seir1_sim<-rbind(seir1_sim,cbind(rep(s,forecast),seir1))

}
seir1_sim[,2]<-rep(1:forecast,dim(coef)[1])
colnames(seir1_sim)[1]<-"sim"

### confidence limits
I_seir_med<-tapply(seir1_sim$I,seir1_sim$time,median)
I_seir_lwr<-tapply(seir1_sim$I,seir1_sim$time,quantile,p=0.05)
I_seir_upr<-tapply(seir1_sim$I,seir1_sim$time,quantile,p=0.95)

days.before<-date[1:days]
days.ahead<-date[(days+1):(days+forecast)]
step.ahead<-forecast+1
mu.lower<-c(dat_csv$totale_positivi,I_seir_lwr*N)
mu.upper<-c(dat_csv$totale_positivi,I_seir_upr*N)
mu.med<-xts(c(dat_csv$totale_positivi,I_seir_med*N),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 1 (Credible Interval ",100*0.90,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario%20plot%20-1.png)<!-- -->

At the end of the 2 weeks (2020-04-21): \*the number of infected is
(1.245618110^{4}).

\*the total number of COVID19 cases is expected to be (2.933575810^{4}).

The estimated (median) numbers of the current scenario by date are:

``` r
S_seir_med<-tapply(seir1_sim$S,seir1_sim$time,median)
E_seir_med<-tapply(seir1_sim$E,seir1_sim$time,median)
R_seir_med<-tapply(seir1_sim$R,seir1_sim$time,median)
forecast<-data.frame(S_seir_med,E_seir_med,I_seir_med,R_seir_med)*N
colnames(forecast)<-c("Susceptible","Exposed","Infected","Removed")
rownames(forecast)<-days.ahead
kable(forecast)
```

|            | Susceptible |  Exposed | Infected |   Removed |
| ---------- | ----------: | -------: | -------: | --------: |
| 2020-04-08 |     4865054 | 21999.69 | 10209.25 |  8601.230 |
| 2020-04-09 |     4864173 | 22073.55 | 10442.53 |  9174.590 |
| 2020-04-10 |     4863308 | 22128.51 | 10665.59 |  9760.604 |
| 2020-04-11 |     4862450 | 22172.22 | 10878.23 | 10358.705 |
| 2020-04-12 |     4861655 | 22155.60 | 11080.10 | 10968.280 |
| 2020-04-13 |     4860874 | 22120.01 | 11271.20 | 11588.757 |
| 2020-04-14 |     4860085 | 22125.19 | 11450.08 | 12219.548 |
| 2020-04-15 |     4859225 | 22153.30 | 11618.75 | 12860.288 |
| 2020-04-16 |     4858394 | 22179.72 | 11778.86 | 13509.741 |
| 2020-04-17 |     4857573 | 22183.26 | 11935.74 | 14167.547 |
| 2020-04-18 |     4856823 | 22136.65 | 12080.58 | 14834.266 |
| 2020-04-19 |     4856076 | 22060.60 | 12213.26 | 15508.094 |
| 2020-04-20 |     4855392 | 21983.07 | 12339.49 | 16190.022 |
| 2020-04-21 |     4854658 | 21913.04 | 12456.18 | 16879.576 |

# Scenario 2 - Forecast based on the last beta coefficient

We estimate a second scenario using a Random Walk on the basis of the
last value of the beta.

``` r
ss<-NULL
set.seed(123)
#we assum log_beta ~ N()
ss <-   AddLocalLevel(, beta_norm[(days-9):(days-4)],sdy=sd_norm_mean)
model1 <- bsts(beta_norm[(days-9):(days-4)],
               state.specification = ss,
               niter = 1000,seed=123)
```

    ## =-=-=-=-= Iteration 0 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 100 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 200 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 300 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 400 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 500 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 600 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 700 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 800 Tue Apr  7 22:42:06 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 900 Tue Apr  7 22:42:06 2020 =-=-=-=-=

``` r
par(mfrow = c(1,1))
plot(model1, "components")
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario2%20sim%20and%20plot%20-1.png)<!-- -->

``` r
#previsioni
pred1 <- predict(model1, horizon = 16, burn = 100)
par(mfrow = c(1,1))
plot(pred1 , ylab="Beta1 coefficient",main="Data and predictions",ylim=c(-1,1)*IQR(pred1$distribution)+median(pred1$distribution))
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario2%20sim%20and%20plot%20-2.png)<!-- -->

``` r
# matrix of beta coefficients
coef<-M_from_Norm(pred1$distribution[,3:16],sd_norm_mean)
par(mfrow = c(1,1))

seir2_sim<-NULL
for(s in 1:dim(coef)[1]){
  # average number of single connections of an infected person
  # less contacts, less probability of new infections
  # we keep constant the other parameters
  forecast<-14
  seir2<-NULL
  for(i in 1:forecast){
    parameters <- c(mu = mu0, beta = (matrix(coef[s,i])+1/duration), sigma = sigma, gamma = 1/duration)
    if( i==1) initials <- c(S = 1-(pro$results$E[last+1]+I_start/N+R_start/N), E = pro$results$E[last+1], I = I_start/N, R = R_start/N)
    if( i>1) initials <- c(S = seir2_temp$results$S[2], E = seir2_temp$results$E[2], I =seir2_temp$results$I[2], R = seir2_temp$results$R[2])
    seir2_temp <- SEIR(pars = parameters, init = initials, time = 0:1)
    seir2 <- rbind(seir2,SEIR(pars = parameters, init = initials, time = 0:1)$results[2,])
  }
  seir2_sim<-rbind(seir2_sim,cbind(rep(s,forecast),seir2))

}
seir2_sim[,2]<-rep(1:forecast,dim(coef)[1])
colnames(seir2_sim)[1]<-"sim"

### confidence limits
I_seir_med<-tapply(seir2_sim$I,seir2_sim$time,median)
I_seir_lwr<-tapply(seir2_sim$I,seir2_sim$time,quantile,p=0.05)
I_seir_upr<-tapply(seir2_sim$I,seir2_sim$time,quantile,p=0.95)

days.before<-date[1:days]
days.ahead<-date[(days+1):(days+forecast)]
step.ahead<-forecast+1
mu.lower<-c(dat_csv$totale_positivi,I_seir_lwr*N)
mu.upper<-c(dat_csv$totale_positivi,I_seir_upr*N)
mu.med<-xts(c(dat_csv$totale_positivi,I_seir_med*N),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 2 (Credible Interval ",100*0.90,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario2%20sim%20and%20plot%20-3.png)<!-- -->

At the end of the 2 weeks (2020-04-21): \*the number of infected is
(1.429743510^{4}).

\*the total number of COVID19 cases is expected to be (3.165758310^{4}).

The estimated (median) numbers of the second scenario by date are:

``` r
S_seir_med<-tapply(seir2_sim$S,seir2_sim$time,median)
E_seir_med<-tapply(seir2_sim$E,seir2_sim$time,median)
R_seir_med<-tapply(seir2_sim$R,seir2_sim$time,median)
forecast<-data.frame(S_seir_med,E_seir_med,I_seir_med,R_seir_med)*N
colnames(forecast)<-c("Susceptible","Exposed","Infected","Removed")
rownames(forecast)<-days.ahead
kable(forecast)
```

|            | Susceptible |  Exposed | Infected |   Removed |
| ---------- | ----------: | -------: | -------: | --------: |
| 2020-04-08 |     4864491 | 22552.82 | 10219.11 |  8601.503 |
| 2020-04-09 |     4863036 | 23170.57 | 10481.25 |  9176.215 |
| 2020-04-10 |     4861541 | 23805.09 | 10751.63 |  9765.699 |
| 2020-04-11 |     4860011 | 24454.56 | 11030.34 | 10370.432 |
| 2020-04-12 |     4858439 | 25115.54 | 11316.72 | 10990.854 |
| 2020-04-13 |     4856819 | 25803.33 | 11612.79 | 11627.331 |
| 2020-04-14 |     4855171 | 26497.09 | 11916.55 | 12280.484 |
| 2020-04-15 |     4853470 | 27209.86 | 12229.34 | 12950.867 |
| 2020-04-16 |     4851724 | 27951.18 | 12549.71 | 13638.795 |
| 2020-04-17 |     4849945 | 28691.67 | 12881.69 | 14344.923 |
| 2020-04-18 |     4848098 | 29483.74 | 13222.97 | 15069.349 |
| 2020-04-19 |     4846211 | 30258.16 | 13574.34 | 15813.020 |
| 2020-04-20 |     4844284 | 31068.22 | 13930.41 | 16575.953 |
| 2020-04-21 |     4842295 | 31906.17 | 14297.43 | 17360.148 |

# Number of admitted in intensive care

The number of person admitted in intensive care repart is a fraction of
the COVID 19 case. The trend can be estimated by means a binomial
regression model.  
We report the number of intensive care patients on the basis of the
estimated COVID19 cases of the previous scenario 1 and
2.

``` r
fit_int_care<-glm(cbind(terapia_intensiva,totale_positivi-terapia_intensiva)~t+I(t^2),family=binomial,data=dat_csv[])




forecast=14
# add 'fit', 'lwr', and 'upr' columns to dataframe (generated by predict)
pre<-predict(fit_int_care,type='response',newdata=data.frame(t=1:(days+forecast)),se.fit=TRUE)
date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+forecast-1+dim(dat_csv)[1],1)
plot(date,pre$fit,type="l",main="Proportion of admitted COVID 19 in intensive care",ylim=range(c(dat_csv$terapia_intensiva/dat_csv$totale_positivi,pre$fit)))
points(date[1:days],dat_csv$terapia_intensiva/dat_csv$totale_positivi)
lines(date,pre$fit+1.96*pre$se.fit,col=3,lty=2)
lines(date,pre$fit-1.96*pre$se.fit,col=3,lty=2)
legend("topright",c("fitted","95% Confidence Interval"),lty=1:2,col=c(1,3))
```

![](draft_analysis_Veneto_new_files/figure-gfm/intensive%20care%20numbers-1.png)<!-- -->

``` r
days.before<-date[1:days]
days.ahead<-date[(days+1):(days+forecast)]
step.ahead<-forecast+1


#Scenario 1

I_seir_med<-tapply(seir1_sim$I,seir1_sim$time,median)
I_seir_lwr<-tapply(seir1_sim$I,seir1_sim$time,quantile,p=0.05)
I_seir_upr<-tapply(seir1_sim$I,seir1_sim$time,quantile,p=0.75)

mu.lower<-c(dat_csv$terapia_intensiva,I_seir_med*N*((pre$fit-1.645*pre$se.fit)[(days+1):(days+forecast)]))
mu.upper<-c(dat_csv$terapia_intensiva,I_seir_med*N*((pre$fit+1.645*pre$se.fit)[(days+1):(days+forecast)]))
mu.med<-xts(c(dat_csv$terapia_intensiva,I_seir_med*N*((pre$fit)[(days+1):(days+forecast)])),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 1 (Confidence Interval ",100*0.90,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/intensive%20care%20numbers-2.png)<!-- -->

``` r
#Numbers
intensive_care<-data.frame(mu.med,mu.lower,mu.upper)[(days+1):(days+forecast),]
kable(intensive_care)
```

|            |   mu.med | mu.lower | mu.upper |
| ---------- | -------: | -------: | -------: |
| 2020-04-08 | 299.3744 | 284.2009 | 314.5478 |
| 2020-04-09 | 288.6311 | 272.0970 | 305.1652 |
| 2020-04-10 | 277.3239 | 259.4568 | 295.1910 |
| 2020-04-11 | 265.5659 | 246.4211 | 284.7108 |
| 2020-04-12 | 253.4633 | 233.1200 | 273.8067 |
| 2020-04-13 | 241.1258 | 219.6829 | 262.5687 |
| 2020-04-14 | 228.6275 | 206.2032 | 251.0518 |
| 2020-04-15 | 216.1069 | 192.8287 | 239.3851 |
| 2020-04-16 | 203.6784 | 179.6800 | 227.6767 |
| 2020-04-17 | 191.4987 | 166.9097 | 216.0878 |
| 2020-04-18 | 179.4824 | 154.4538 | 204.5110 |
| 2020-04-19 | 167.6972 | 142.3822 | 193.0123 |
| 2020-04-20 | 156.2768 | 130.8152 | 181.7384 |
| 2020-04-21 | 145.2207 | 119.7558 | 170.6857 |

``` r
#Scenario 2

I_seir_med<-tapply(seir2_sim$I,seir2_sim$time,median)
I_seir_lwr<-tapply(seir2_sim$I,seir2_sim$time,quantile,p=0.25)
I_seir_upr<-tapply(seir2_sim$I,seir2_sim$time,quantile,p=0.75)

mu.lower<-c(dat_csv$terapia_intensiva,I_seir_med*N*((pre$fit-1.645*pre$se.fit)[(days+1):(days+forecast)]))
mu.upper<-c(dat_csv$terapia_intensiva,I_seir_med*N*((pre$fit+1.645*pre$se.fit)[(days+1):(days+forecast)]))
mu.med<-xts(c(dat_csv$terapia_intensiva,I_seir_med*N*((pre$fit)[(days+1):(days+forecast)])),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 2 (Confidence Interval ",100*0.9,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/intensive%20care%20numbers-3.png)<!-- -->

``` r
#Numbers
intensive_care<-data.frame(mu.med,mu.lower,mu.upper)[(days+1):(days+forecast),]
kable(intensive_care)
```

|            |   mu.med | mu.lower | mu.upper |
| ---------- | -------: | -------: | -------: |
| 2020-04-08 | 299.6635 | 284.4754 | 314.8516 |
| 2020-04-09 | 289.7012 | 273.1058 | 306.2967 |
| 2020-04-10 | 279.5610 | 261.5498 | 297.5722 |
| 2020-04-11 | 269.2793 | 249.8668 | 288.6919 |
| 2020-04-12 | 258.8762 | 238.0984 | 279.6540 |
| 2020-04-13 | 248.4336 | 226.3408 | 270.5264 |
| 2020-04-14 | 237.9416 | 214.6038 | 261.2794 |
| 2020-04-15 | 227.4638 | 202.9623 | 251.9653 |
| 2020-04-16 | 217.0078 | 191.4389 | 242.5766 |
| 2020-04-17 | 206.6757 | 180.1379 | 233.2136 |
| 2020-04-18 | 196.4550 | 169.0597 | 223.8504 |
| 2020-04-19 | 186.3859 | 158.2497 | 214.5222 |
| 2020-04-20 | 176.4254 | 147.6811 | 205.1698 |
| 2020-04-21 | 166.6871 | 137.4579 | 195.9162 |
