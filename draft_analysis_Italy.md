COVID19 - Forecast and predictions using a time dependent SEIR model
================
Paolo Girardi
15 Marzo, 2020

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This
work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative
Commons Attribution-NonCommercial 4.0 International License</a>.

# Disclaimer

  - We want to investigate the evolution of the coronavirus pandemic in
    Italy from a statistical perspective using aggregated data.

  - Our point of view is that of surveillance with the goal of detecting
    important changes in the underlying (random) process as soon as
    possible after it has occured.

  - We use data provided by Italian Civil Protection Department

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
```

and load them.

    ## Loading required package: zoo

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

    ## phantomjs has been installed to /Users/Paolo/Library/Application Support/PhantomJS

# Source of the data

Download the data from

<https://github.com/pcm-dpc/COVID-19/>

# Results

## Load dataset

    ## [1] 21

Several outcomes can be potentially monitored, that is

    ##  [1] "ricoverati_con_sintomi"      "terapia_intensiva"          
    ##  [3] "totale_ospedalizzati"        "isolamento_domiciliare"     
    ##  [5] "totale_attualmente_positivi" "nuovi_attualmente_positivi" 
    ##  [7] "dimessi_guariti"             "deceduti"                   
    ##  [9] "totale_casi"                 "tamponi"

It is worth noting that some outcomes present negative counts in some
regions. It looks like some of these negative counts are redesignations.
Outcomes presenting negative values cannot be analyzed using the
proposed model.

Then we extract the timeseries.

![](draft_analysis_Italy_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### The S(E)IR model (to be revised)

With the aim of predicting the future number of COVID19 cases on the
basis of the actual data, we used a SEIR model applied to the COVID19
epidemic in Italy.

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
We’ll assume that the force of infection, \(\lambda\), for a constant
population \(N\) \[\lambda = \beta\,\frac{I}{N},\] so that the risk of
infection a susceptible faces is proportional to the *prevalence* (the
fraction of the population that is infected). This is known as the
assumption of frequency-dependent transmission.

# The reproduction number of COVID19.

The number of infected individuals \(I\) at time \(t\) is approximately
\[I(t)\;\approx\;I_0\,e^{R_0\,(\gamma+\mu)\,t}\] where \(I_0\) is the
(small) number of infectives at time \(0\), \(\frac{1}{\gamma}\) is the
infectious period, and \(\frac{1}{\mu}\) is the host lifespan.

\(R_0\) is the reproduction number
(<https://en.wikipedia.org/wiki/Basic_reproduction_number>) and
indicates how contagious an infectious disease is.

Taking logs of both sides, we get

\[\log{I}(t)\;\approx\;\log{I_0}+(R_0)\,(\gamma+\mu)\,t,\] which implies
that a semi-log plot of \(I\) vs \(t\) should be approximately linear
with a slope proportional to \(R_0\) and the recovery rate.

![](draft_analysis_Italy_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

We estimate the \(R_0\) parameter in the linear model.

\[
\log(I(t))= \alpha + \beta  t +e_t
\]

The estimated slope coefficient \(\hat\beta\) is used to estimate
\(R_0\) as in the following formula:

\[\widehat\beta=(\widehat{R_0})\,(\gamma+\mu)\] The parameter
\(\mu\)\<\<\(\gamma\) and it can not be considered. As consequence, R0
can be estimated as follows \[\hat{R_0}=\frac{\hat{\beta}}{\gamma}
\]

The incubation period \(1/ \gamma\) for the coronavirus is in mean 5.1
days with a range from 2-14 days. Please see
<https://www.worldometers.info/coronavirus/coronavirus-incubation-period/>.
However, the incubation period is used for epidemic diseases that causes
the immediate home isolation of infected subjects.

In the calculation we considered an “incubation period” of 14 days for
two reasons:  
1\) the majority of cases is asymptomatic, contagiousness is greater
than 5, maybe 14. A minority (who made the swab) will have a duration of
about 5 days between the start of contagiousness and swab; 2) 14 days is
the worst scenario because in this period low impact COVID19 symptoms
can be confused with the concomitant FLU epidemic
(<https://www.webmd.com/lung/news/20200310/know-the-symptoms-of-covid19>).

We calculate several R0 values, each one based on a mobile window of 5
days, that can be sufficient to estimate a local trend, in order to
assess if the R0 trend is decreasing (how is expected to be).

![](draft_analysis_Italy_files/figure-gfm/R0%20trend-1.png)<!-- -->

The R0 shows a decreasing trend in the last period. We use the estimated
trend between R0 and time to calculate the future R0 value for the next
14 days. We predict beta (and R0) for the next 14 days by means of a
linear regressione model, assuming a Log-normal distribution for the
beta (the slope) forcing its value to be greater than 0.

    ## Warning: Removed 18 rows containing missing values (geom_point).

![](draft_analysis_Italy_files/figure-gfm/R0%20forecast-1.png)<!-- -->

R0 passes from a value of 4.57 in the initial phase to an estimated
value of 1.29 at the ending of the 14-days forecast.

We want to make a short term forecast (14 days) with 3 scenario, based
on the number of exposed people:

\-Scenario 1: 10 exposed people for each COVID-19 case (no home
restrictions made or even no effects)

\-Scenario 2: 5 exposed people for each COVID-19 case (-50% exposed
people)

\-Scenario 3: 3 exposed people for each COVID-19 case (-70% exposed
people)

We made a forecast by means of a SEIR model fixing a series of initial
status:  
\- S0=N, the size of Italian population  
\- E= f \* I0 (with f a fixed factor of the previous scenario)  
\- I0: initial number of COVID-19 cases  
\- R0: initial number of recovered

and parameters:  
\- beta: the quantity connected to R0 is considered to vary according
the previous estimation  
\- gamma= 1/duration (rate of infection duration of COVID-19, 14 days)  
\- sigma0: the coronavirus transmission rate (half of flu epidemic)  
\- mu0: the overall mortality rate

``` r
# initial number of infectus
I0<-dat_csv$totale_attualmente_positivi[dim(dat_csv)[1]]; I0
```

    ## [1] 20603

``` r
# initial number of recovered
R0<-dat_csv$dimessi_guariti[dim(dat_csv)[1]]; R0
```

    ## [1] 2335

``` r
# italian poulation
N=60480000
# duration of COVID19 
duration<-14
#sigma0 is the coronavirus transmission rate fixed to 5%  (half of flu epidemic)
sigma0<-0.05
#mortality rate 
mu0<-1/(82*365.25) # 1/lifespan
```

We use the library(EpiDynamics) and the function SEIR() to implement a
SEIR
model:

<img src="http://www.public.asu.edu/~hnesse/classes/seireqn.png"/>  
<img src="https://upload.wikimedia.org/wikipedia/commons/3/3d/SEIR.PNG"/>

where the parameter beta here is time dependent, as estimated before by
the gamma regression model.

``` r
library(EpiDynamics)

# average number of single connections of an infected person
# less contacts, less probability of new infections
# we keep constant the other parameters
forecast<-14
seir1<-seir2<-seir3<-NULL
for(i in 1:forecast){
parameters <- c(mu = mu0, beta = beta.predict$fit[days+i], sigma = sigma0, gamma = 1/duration)
f1<-10
if( i==1) initials <- c(S = 0.95, E = (f1*I0/N), I = I0/N, R = R0/N)
if( i>1) initials <- c(S = seir1_temp$results$S[2], E = seir1_temp$results$E[2], I =seir1_temp$results$I[2], R = seir1_temp$results$R[2])
seir1_temp <- SEIR(pars = parameters, init = initials, time = 0:1)
seir1 <- rbind(seir1,SEIR(pars = parameters, init = initials, time = 0:1)$results[2,])
f2<-5
if( i==1) initials <- c(S = 0.95, E = (f2*I0/N), I = I0/N, R = R0/N)
if( i>1) initials <- c(S = seir2_temp$results$S[2], E = seir2_temp$results$E[2], I =seir2_temp$results$I[2], R = seir2_temp$results$R[2])
seir2_temp <- SEIR(pars = parameters, init = initials, time = 0:1)
seir2 <- rbind(seir2,SEIR(pars = parameters, init = initials, time = 0:1)$results[2,])
f3<-3
if( i==1) initials <- c(S = 0.95, E = (f3*I0/N), I = I0/N, R = R0/N)
if( i>1) initials <- c(S = seir3_temp$results$S[2], E = seir3_temp$results$E[2], I =seir3_temp$results$I[2], R = seir3_temp$results$R[2])
seir3_temp <- SEIR(pars = parameters, init = initials, time = 0:1)
seir3 <- rbind(seir3,SEIR(pars = parameters, init = initials, time = 0:1)$results[2,])
}

date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+forecast-1+dim(dat_csv)[1],1)
plot(date,c(dat_csv$totale_attualmente_positivi,seir1$I*N),type="l",ylab="Cases",xlab="time",main="Infected")
lines(date,c(dat_csv$totale_attualmente_positivi,seir2$I*N),col=2)
lines(date,c(dat_csv$totale_attualmente_positivi,seir3$I*N),col=3)
lines(date[1:dim(dat_csv)[1]],dat_csv$totale_attualmente_positivi,lwd=2)
legend("topleft",c("first scenario - Exp=10*I","second scenario Exp=5*I","third scenario Exp=3*I"),lty=1,col=1:3)
```

![](draft_analysis_Italy_files/figure-gfm/scenario%20plot-1.png)<!-- -->

The 3 scenarios show different numbers. If we consider the second
scenario, at the end of the 2 weeks (2020-03-29) the number of infected
is (4.978022710^{4}).

In the next plot the cumulative number of infected.  
At the end of the 2 weeks (2020-03-29) the total number of COVID19 cases
is expected to be
(8.106112910^{4}).

``` r
plot(date,c(dat_csv$totale_casi,(seir1$I+seir1$R)*N),type="l",ylab="Cases",xlab="time",main="Cumulative Infected")
lines(date,c(dat_csv$totale_casi,(seir2$I+seir2$R)*N),col=2)
lines(date,c(dat_csv$totale_casi,(seir3$I+seir3$R)*N),col=3)
lines(date[1:dim(dat_csv)[1]],(dat_csv$totale_casi),lwd=2)
legend("topleft",c("Scenario 1 -Exp=10*I","Scenario 2 - Exp=5*I","Scenario 3- Exp=3*I"),lty=1,col=1:3)
```

![](draft_analysis_Italy_files/figure-gfm/cumulative%20plot-1.png)<!-- -->

``` r
plot(date[(days+1):(days+forecast)],seir2$E*N,type="l",ylab="Cases",xlab="Date",main="14 days - Forecast Scenario 2",ylim=c(0,max(seir2$E*N,seir2$I*N,seir2$R*N )))
lines(date[(days+1):(days+forecast)],seir2$I*N,col=2)
lines(date[(days+1):(days+forecast)],seir2$R*N,col=3)
legend("topleft",c("Exposed","Infected","Recovered"),lty=1,col=1:4)
```

![](draft_analysis_Italy_files/figure-gfm/forecast%20plot-1.png)<!-- -->
