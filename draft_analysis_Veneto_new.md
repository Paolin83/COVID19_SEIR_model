COVID19 - Forecast and predictions using a time dependent SEIR model for
the Veneto Region
================
Paolo Girardi
28 Marzo, 2020

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

    ## [1] 34

Several outcomes can be potentially monitored, that is

``` r
names(dat_csv[,-c(1:2,17:19)])
```

    ##  [1] "codice_regione"              "denominazione_regione"      
    ##  [3] "lat"                         "long"                       
    ##  [5] "ricoverati_con_sintomi"      "terapia_intensiva"          
    ##  [7] "totale_ospedalizzati"        "isolamento_domiciliare"     
    ##  [9] "totale_attualmente_positivi" "nuovi_attualmente_positivi" 
    ## [11] "dimessi_guariti"             "deceduti"                   
    ## [13] "totale_casi"                 "tamponi"

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
dat_csv_dy<-xts(dat_csv[,-c(1:7,17:19)], order.by = days_dy, frequency = 7)
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
linear with a slope proportional to \(R_0\) and the recovery
rate.

``` r
dat_csv_dy$log_totale_attualmente_positivi<-log(dat_csv_dy$totale_attualmente_positivi)
p <- dygraph(dat_csv_dy$log_totale_attualmente_positivi,main=paste("Veneto Region"),ylab="Log Infected case",xlab="Day",height=400,width=800) 
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
beta_vec<-NULL
sd_vec<-NULL
#for cycle for R0 estimates from days-2 to days+2
for (i in 3:(days-2)){
fit <- lm(log(totale_attualmente_positivi)~t,data=dat_csv[(i-2):(i+2),])
beta_vec<-c(beta_vec,coef(fit)[2])
sd_vec<-c(sd_vec,coef(summary(fit))[2,2])
}

label<-as.Date(substr(dat_csv$data,1,10))[3:(days-2)]


mean  <- 1+(beta_vec*duration)
lower <- 1+((beta_vec-1.96*sd_vec)*duration)
upper <- 1+((beta_vec+1.96*sd_vec)*duration)

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
linear regressione model, assuming a Normal distribution for the beta
(the slope).

``` r
#start from the day 7, since vo hotspot can influence the trend
time<-3:(days-2)
weekend<-rep(c(1,2,3,4,5,6,7),ceiling(days/7))[3:(days-2)]
data<-data.frame(time)
beta.model<-glm(beta_vec~time,weights = 1/sd_vec,family=gaussian,data)
summary(beta.model)
```

    ## 
    ## Call:
    ## glm(formula = beta_vec ~ time, family = gaussian, data = data, 
    ##     weights = 1/sd_vec)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.89244  -0.18529   0.06661   0.25146   1.20164  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.2428613  0.0198899  12.210 9.88e-13 ***
    ## time        -0.0056055  0.0007546  -7.428 4.33e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.1981073)
    ## 
    ##     Null deviance: 16.478  on 29  degrees of freedom
    ## Residual deviance:  5.547  on 28  degrees of freedom
    ## AIC: -99.168
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
anova(beta.model,test="Chisq")
```

    ## Analysis of Deviance Table
    ## 
    ## Model: gaussian, link: identity
    ## 
    ## Response: beta_vec
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##      Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
    ## NULL                    29     16.478             
    ## time  1   10.931        28      5.547  1.1e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#there is an effect of the week, as supposed 
#we can use
#summary(gamlss(formula = beta_vec ~ time + bs(weekend, 3),      weights = 1/sd_vec,Lp=1) )
forecast=14
# add 'fit', 'lwr', and 'upr' columns to dataframe (generated by predict)
pre<-predict(beta.model,type='response',newdata=data.frame(time=1:(days+forecast)),se.fit=TRUE)
date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+forecast-1+dim(dat_csv)[1],1)
predict <- data.frame(beta_vec=c(rep(NA,2),beta_vec,rep(NA,forecast+2)),time=date,fit=pre$fit,lwr=pre$fit-1*1.96*pre$se.fit,upr=pre$fit+1*1.96*pre$se.fit)
beta.predict<-predict 
r0.predict<-beta.predict
r0.predict[,c(1,3:5)]<-r0.predict[,c(1,3:5)]*duration+1
# plot the points (actual observations), regression line, and confidence interval
p <- ggplot(r0.predict, aes(date,beta_vec))
p <- p + geom_point() +labs(x="Date",y="R0 value") 
p <- p + geom_line(aes(date,fit))
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3)
p
```

    ## Warning: Removed 18 rows containing missing values (geom_point).

![](draft_analysis_Veneto_new_files/figure-gfm/R0%20forecast-1.png)<!-- -->

R0 passes from a value of NA in the initial phase to an estimated value
of 0.53 at the ending of the 14-days forecast.  
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
I_start<-dat_csv$totale_attualmente_positivi[dim(dat_csv)[1]]; I_start
```

    ## [1] 6913

``` r
# initial number of recovered, based on the proportion of discharged from the health services
prop<-dat_csv$dimessi_guariti[dim(dat_csv)[1]]/dat_csv$totale_ospedalizzati[dim(dat_csv)[1]]
R_start<-prop*dat_csv$totale_casi[dim(dat_csv)[1]]; R_start
```

    ## [1] 2729.453

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
I_1<-dat_csv$totale_attualmente_positivi[(days-last)]
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
parameters <- c(mu = mu0, beta =(as.numeric(beta_vec[(days-last-2)])+1/duration), sigma = sigma, gamma =1/duration)
initials <- c(S = (1-(f_exp+I_1/N+R_start/N)), E = f_exp, I = I_1/N, R = R_1/N)
-LogLikfun(initials,parameters,dat_csv$totale_attualmente_positivi[c((days-last):days)])
}
est<-optim(fn=f1,par=logit(0.1))
```

    ## Warning in optim(fn = f1, par = logit(0.1)): one-dimensional optimization by Nelder-Mead is unreliable:
    ## use "Brent" or optimize() directly

``` r
est
```

    ## $par
    ## [1] -2.569079
    ## 
    ## $value
    ## [1] 64.9143
    ## 
    ## $counts
    ## function gradient 
    ##       26       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

``` r
expit(est$par)
```

    ## [1] 0.07115515

``` r
f_exp<-3*I_1/N
sigma<-expit(est$par)

parameters <- c(mu = mu0, beta =(as.numeric(beta_vec[(days-last-2)])+1/duration), sigma = sigma, gamma = 1/duration)
initials <- c(S = (1-(f_exp+I_1/N+R_1/N)), E = f_exp, I = I_1/N, R = R_1/N)
pro<-SEIR(pars = parameters, init = initials, time = 0:last)
#predicted vs observed
cbind(pro$results$I*N,dat_csv$totale_attualmente_positivi[c((days-last):days)])
```

    ##           [,1] [,2]
    ##  [1,] 2953.000 2953
    ##  [2,] 3404.142 3169
    ##  [3,] 3828.700 3677
    ##  [4,] 4233.272 4214
    ##  [5,] 4622.866 4644
    ##  [6,] 5003.372 4986
    ##  [7,] 5378.620 5351
    ##  [8,] 5752.207 5745
    ##  [9,] 6127.319 6140
    ## [10,] 6506.759 6648
    ## [11,] 6893.026 6913

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

    ## [1] 0.07115515

``` r
#gamma
1/duration
```

    ## [1] 0.05555556

``` r
#mthe unknowrn fraction of exposed people is
pro$results$E[last+1]
```

    ## [1] 0.002216225

``` r
#that is 
pro$results$E[last+1]/pro$results$I[last+1]
```

    ## [1] 1.577318

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
#we drop the first 6 beta_vec
ss <- AddStudentLocalLinearTrend(list(), beta_vec[6:(days-4)])
ss <- AddSeasonal(ss, beta_vec[6:(days-4)], nseasons = 7)
model1 <- bsts(beta_vec[6:(days-4)],
               state.specification = ss,
               niter = 1000,seed=123)
```

    ## =-=-=-=-= Iteration 0 Sat Mar 28 19:59:01 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 100 Sat Mar 28 19:59:01 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 200 Sat Mar 28 19:59:01 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 300 Sat Mar 28 19:59:01 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 400 Sat Mar 28 19:59:01 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 500 Sat Mar 28 19:59:02 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 600 Sat Mar 28 19:59:02 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 700 Sat Mar 28 19:59:02 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 800 Sat Mar 28 19:59:02 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 900 Sat Mar 28 19:59:02 2020 =-=-=-=-=

``` r
par(mfrow = c(1,1))
plot(model1, "components")
```

![](draft_analysis_Veneto_new_files/figure-gfm/bsts%20model%20fit%20and%20prediction-1.png)<!-- -->

``` r
#previsioni
pred1 <- predict(model1, horizon = 16, burn = 100)
par(mfrow = c(1,1))
plot(pred1 , ylab="Beta1 coefficient",main="Data and predictions",ylim=c(-1,1)*IQR(pred1$distribution)+median(pred1$distribution))
```

![](draft_analysis_Veneto_new_files/figure-gfm/bsts%20model%20fit%20and%20prediction-2.png)<!-- -->

``` r
# matrix of beta coefficients
coef<-(pred1$distribution[,3:16])
par(mfrow = c(1,1))
```

For each vector of simulated beta coefficients we perform a SEIR model.
We save the results and plot a credible interval at 50% (25-75%
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
    parameters <- c(mu = mu0, beta = (matrix(coef[s,i])+1/duration), sigma = sigma, gamma = 1/duration)
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
I_seir_lwr<-tapply(seir1_sim$I,seir1_sim$time,quantile,p=0.25)
I_seir_upr<-tapply(seir1_sim$I,seir1_sim$time,quantile,p=0.75)

days.before<-date[1:days]
days.ahead<-date[(days+1):(days+forecast)]
step.ahead<-forecast+1
mu.lower<-c(dat_csv$totale_attualmente_positivi,I_seir_lwr*N)
mu.upper<-c(dat_csv$totale_attualmente_positivi,I_seir_upr*N)
mu.med<-xts(c(dat_csv$totale_attualmente_positivi,I_seir_med*N),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 1 (Credible Interval ",100*0.50,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario%20plot%20-1.png)<!-- -->

At the end of the 2 weeks (2020-04-11): \*the number of infected is
(7716.9771735).

\*the total number of COVID19 cases is expected to be (1.676656210^{4}).

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

|            | Susceptible |   Exposed | Infected |  Removed |
| ---------- | ----------: | --------: | -------: | -------: |
| 2020-03-29 |     4884938 | 10522.841 | 7279.709 | 3123.600 |
| 2020-03-30 |     4884605 | 10119.805 | 7600.780 | 3536.820 |
| 2020-03-31 |     4884378 |  9636.958 | 7873.291 | 3966.524 |
| 2020-04-01 |     4884136 |  9226.677 | 8101.142 | 4410.050 |
| 2020-04-02 |     4883891 |  8842.366 | 8288.881 | 4865.159 |
| 2020-04-03 |     4883757 |  8368.189 | 8440.188 | 5329.339 |
| 2020-04-04 |     4883755 |  7789.699 | 8527.879 | 5799.920 |
| 2020-04-05 |     4884023 |  7052.824 | 8573.461 | 6274.868 |
| 2020-04-06 |     4884224 |  6214.268 | 8581.084 | 6752.419 |
| 2020-04-07 |     4884801 |  5288.862 | 8503.189 | 7223.839 |
| 2020-04-08 |     4885336 |  4481.485 | 8403.109 | 7697.492 |
| 2020-04-09 |     4885684 |  3813.974 | 8244.713 | 8163.199 |
| 2020-04-10 |     4886376 |  2849.263 | 8009.862 | 8614.065 |
| 2020-04-11 |     4887183 |  2013.219 | 7716.977 | 9049.585 |

# Scenario 2 - Forecast based on the last beta coefficient

We estimate a second scenario using a Random Walk on the basis of the
last value of the beta.

``` r
ss<-NULL
set.seed(123)
ss <-   AddLocalLevel(list(), beta_vec[(days-9):(days-4)],sigma=SdPrior(sd_vec[(days-4)]))
model1 <- bsts(beta_vec[(days-9):(days-4)],
               state.specification = ss,
               niter = 1000,seed=123)
```

    ## =-=-=-=-= Iteration 0 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 100 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 200 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 300 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 400 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 500 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 600 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 700 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 800 Sat Mar 28 19:59:30 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 900 Sat Mar 28 19:59:30 2020 =-=-=-=-=

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
coef<-(pred1$distribution[,3:16])
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
mu.lower<-c(dat_csv$totale_attualmente_positivi,I_seir_lwr*N)
mu.upper<-c(dat_csv$totale_attualmente_positivi,I_seir_upr*N)
mu.med<-xts(c(dat_csv$totale_attualmente_positivi,I_seir_med*N),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 2 (Credible Interval ",100*0.90,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario2%20sim%20and%20plot%20-3.png)<!-- -->

At the end of the 2 weeks (2020-04-11): \*the number of infected is
(1.187671610^{4}).

\*the total number of COVID19 cases is expected to be (2.192486610^{4}).

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

|            | Susceptible |  Exposed |  Infected |   Removed |
| ---------- | ----------: | -------: | --------: | --------: |
| 2020-03-29 |     4884489 | 10955.98 |  7294.703 |  3124.016 |
| 2020-03-30 |     4883565 | 11094.88 |  7663.527 |  3539.414 |
| 2020-03-31 |     4882630 | 11239.46 |  8021.658 |  3974.969 |
| 2020-04-01 |     4881644 | 11415.17 |  8372.768 |  4430.226 |
| 2020-04-02 |     4880611 | 11636.92 |  8716.818 |  4904.953 |
| 2020-04-03 |     4879476 | 11918.01 |  9060.186 |  5398.463 |
| 2020-04-04 |     4878328 | 12191.38 |  9407.108 |  5911.317 |
| 2020-04-05 |     4877214 | 12463.44 |  9758.104 |  6442.782 |
| 2020-04-06 |     4875994 | 12761.10 | 10106.900 |  6995.636 |
| 2020-04-07 |     4874776 | 13081.17 | 10454.135 |  7566.648 |
| 2020-04-08 |     4873535 | 13445.11 | 10802.896 |  8158.068 |
| 2020-04-09 |     4872165 | 13766.50 | 11157.608 |  8767.290 |
| 2020-04-10 |     4870844 | 14080.50 | 11504.195 |  9398.654 |
| 2020-04-11 |     4869491 | 14487.79 | 11876.716 | 10048.150 |

# Number of admitted in intensive care

The number of person admitted in intensive care repart is a fraction of
the COVID 19 case. The trend can be estimated by means a binomial
regression model.  
We report the number of intensive care patients on the basis of the
estimated COVID19 cases of the previous scenario 1 and
2.

``` r
fit_int_care<-glm(cbind(terapia_intensiva,totale_attualmente_positivi-terapia_intensiva)~t,family=binomial,data=dat_csv)




forecast=14
# add 'fit', 'lwr', and 'upr' columns to dataframe (generated by predict)
pre<-predict(fit_int_care,type='response',newdata=data.frame(t=1:(days+forecast)),se.fit=TRUE)
date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+forecast-1+dim(dat_csv)[1],1)
plot(date,pre$fit,type="l",main="Proportion of admitted COVID 19 in intensive care",ylim=range(dat_csv$terapia_intensiva/dat_csv$totale_attualmente_positivi))
points(date[1:days],dat_csv$terapia_intensiva/dat_csv$totale_attualmente_positivi)
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
| 2020-03-29 | 374.5412 | 359.2434 | 389.8391 |
| 2020-03-30 | 384.4619 | 367.5746 | 401.3493 |
| 2020-03-31 | 391.5204 | 373.0753 | 409.9654 |
| 2020-04-01 | 396.0414 | 376.0868 | 415.9961 |
| 2020-04-02 | 398.3640 | 376.9592 | 419.7688 |
| 2020-04-03 | 398.7674 | 375.9823 | 421.5526 |
| 2020-04-04 | 396.0826 | 372.0809 | 420.0843 |
| 2020-04-05 | 391.4460 | 366.3541 | 416.5378 |
| 2020-04-06 | 385.1435 | 359.0919 | 411.1952 |
| 2020-04-07 | 375.1640 | 348.4466 | 401.8814 |
| 2020-04-08 | 364.4453 | 337.1782 | 391.7123 |
| 2020-04-09 | 351.4917 | 323.9184 | 379.0650 |
| 2020-04-10 | 335.6649 | 308.1075 | 363.2223 |
| 2020-04-11 | 317.8806 | 290.6161 | 345.1452 |

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
| 2020-03-29 | 375.3126 | 359.9833 | 390.6420 |
| 2020-03-30 | 387.6358 | 370.6090 | 404.6626 |
| 2020-03-31 | 398.8983 | 380.1057 | 417.6910 |
| 2020-04-01 | 409.3204 | 388.6967 | 429.9442 |
| 2020-04-02 | 418.9307 | 396.4208 | 441.4405 |
| 2020-04-03 | 428.0600 | 403.6011 | 452.5189 |
| 2020-04-04 | 436.9189 | 410.4426 | 463.3951 |
| 2020-04-05 | 445.5342 | 416.9753 | 474.0932 |
| 2020-04-06 | 453.6265 | 422.9426 | 484.3104 |
| 2020-04-07 | 461.2405 | 428.3932 | 494.0878 |
| 2020-04-08 | 468.5247 | 433.4706 | 503.5788 |
| 2020-04-09 | 475.6753 | 438.3602 | 512.9904 |
| 2020-04-10 | 482.1001 | 442.5206 | 521.6795 |
| 2020-04-11 | 489.2301 | 447.2690 | 531.1913 |
