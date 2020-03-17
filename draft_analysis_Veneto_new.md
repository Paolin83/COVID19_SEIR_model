COVID19 - Forecast and predictions using a time dependent SEIR model -
Veneto Region
================
Paolo Girardi
17 Marzo, 2020

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
    analysis was restricted to the Veneto Rego

  - This document is in a draft mode, and it is continuously updated.

  - The layout of the draft must definitely be improved.

\*NB: set the file output format to

\#output:html\_document:  
df\_print: paged  
pdf\_document:  
toc: yes

which performs the same analysis enabling Javascript Pictures.

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
```

and load them.

    ## Loading required package: zoo

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

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

# Source of the data

Download the data from

<https://github.com/pcm-dpc/COVID-19/>

# Results

## Load dataset

    ## [1] 23

Several outcomes can be potentially monitored, that is

    ##  [1] "codice_regione"              "denominazione_regione"      
    ##  [3] "lat"                         "long"                       
    ##  [5] "ricoverati_con_sintomi"      "terapia_intensiva"          
    ##  [7] "totale_ospedalizzati"        "isolamento_domiciliare"     
    ##  [9] "totale_attualmente_positivi" "nuovi_attualmente_positivi" 
    ## [11] "deceduti"                    "totale_casi"                
    ## [13] "tamponi"                     "t"

It is worth noting that some outcomes present negative counts in some
regions. It looks like some of these negative counts are redesignations.
Outcomes presenting negative values cannot be analyzed using the
proposed model.

Then we extract the
timeseries.

![](draft_analysis_Veneto_new_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### The S(E)IR model (to be revised)

With the aim of predicting the future number of COVID19 cases on the
basis of the actual data, we used a SEIR model applied to the COVID19
epidemic to the Veneto Region in Italy

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
with a slope proportional to \(R_0\) and the recovery
rate.

![](draft_analysis_Veneto_new_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

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
assess if the R0 trend is decreasing (how is expected to be). In this
way, the R0 for the first and the last two days of observation is
impossibile to estimate.

![](draft_analysis_Veneto_new_files/figure-gfm/R0%20trend-1.png)<!-- -->

The R0 shows a decreasing trend in the last period. We use the estimated
trend between R0 and time to calculate the future R0 value for the next
14 days. We predict beta (and R0) for the next 14 days by means of a
linear regressione model, assuming a Log-normal distribution for the
beta (the slope) and forcing its value to be greater than
0.

    ## Warning: Removed 18 rows containing missing values (geom_point).

![](draft_analysis_Veneto_new_files/figure-gfm/R0%20forecast-1.png)<!-- -->

R0 passes from a value of 5.64 in the initial phase to an estimated
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
\- S0=N, the size of Veneto Region population  
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

    ## [1] 2488

``` r
# initial number of recovered
R0<-dat_csv$dimessi_guariti[dim(dat_csv)[1]]; R0
```

    ## [1] 136

``` r
# Veneto Region population
N=4905854
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
Respect to the previous formulation, the compartment of exposed people
(E) was inserted between the susceptible and infected compartments.  
(<https://en.wikipedia.org/wiki/Bayesian_structural_time_series>)

We perfom simulation on the trend of beta coefficient by means of a
Bayesian Structural Time Series using the library bsts of R.  
We estimate a BSTS model specifing a local linear trend ans 1000
simulation.

We made a forecast forward of 14 days dropping the first 100 simulation
as burn-in.

``` r
# Bayesian Structural Time Series
ss <- AddLocalLinearTrend(list(), log(beta_vec))
model1 <- bsts(log(beta_vec),
               state.specification = ss,
               niter = 1000)
```

    ## =-=-=-=-= Iteration 0 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 100 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 200 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 300 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 400 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 500 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 600 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 700 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 800 Tue Mar 17 18:27:40 2020 =-=-=-=-=
    ## =-=-=-=-= Iteration 900 Tue Mar 17 18:27:41 2020 =-=-=-=-=

``` r
par(mfrow = c(1,1))
plot(model1, "components", ylab="Log beta coefficient", xlab="days")
```

![](draft_analysis_Veneto_new_files/figure-gfm/bsts%20model%20fit%20and%20prediction-1.png)<!-- -->

``` r
#previsioni
pred1 <- predict(model1, horizon = 16, burn = 100)
par(mfrow = c(1,1))
plot(pred1 , ylab="Log beta coefficient",main="Data and predictions")
```

![](draft_analysis_Veneto_new_files/figure-gfm/bsts%20model%20fit%20and%20prediction-2.png)<!-- -->

``` r
# matrix of beta coefficients
coef<-(exp(pred1$distribution)[,3:16])
par(mfrow = c(1,1))
```

For each vectors of simulated beta coefficients we perform a SEIR model
(one for each scenario).  
We save the results and plot a credible interval at 50% (25-75%
percentile).

``` r
seir1_sim<-seir2_sim<-seir3_sim<-NULL
for(s in 1:dim(coef)[1]){
  # average number of single connections of an infected person
  # less contacts, less probability of new infections
  # we keep constant the other parameters
  forecast<-14
  seir1<-seir2<-seir3<-NULL
  for(i in 1:forecast){
    parameters <- c(mu = mu0, beta = matrix(coef[s,i]), sigma = sigma0, gamma = 1/duration)
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
  seir1_sim<-rbind(seir1_sim,cbind(rep(s,forecast),seir1))
  seir2_sim<-rbind(seir2_sim,cbind(rep(s,forecast),seir2))
  seir3_sim<-rbind(seir3_sim,cbind(rep(s,forecast),seir3))
}
seir1_sim[,2]<-seir2_sim[,2]<-seir3_sim[,2]<-rep(1:forecast,dim(coef)[1])
colnames(seir1_sim)[1]<-colnames(seir2_sim)[1]<-colnames(seir3_sim)[1]<-"sim"

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

``` r
### confidence limits - scenario 2
I_seir_med<-tapply(seir2_sim$I,seir2_sim$time,median)
I_seir_lwr<-tapply(seir2_sim$I,seir2_sim$time,quantile,p=0.25)
I_seir_upr<-tapply(seir2_sim$I,seir2_sim$time,quantile,p=0.75)

mu.lower<-c(dat_csv$totale_attualmente_positivi,I_seir_lwr*N)
mu.upper<-c(dat_csv$totale_attualmente_positivi,I_seir_upr*N)
mu.med<-xts(c(dat_csv$totale_attualmente_positivi,I_seir_med*N),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 2  (Credible Interval ",100*0.50,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario%20plot%20-2.png)<!-- -->

``` r
### confidence limits - scenario 3
I_seir_med<-tapply(seir3_sim$I,seir3_sim$time,median)
I_seir_lwr<-tapply(seir3_sim$I,seir3_sim$time,quantile,p=0.25)
I_seir_upr<-tapply(seir3_sim$I,seir3_sim$time,quantile,p=0.75)


mu.lower<-c(dat_csv$totale_attualmente_positivi,I_seir_lwr*N)
mu.upper<-c(dat_csv$totale_attualmente_positivi,I_seir_upr*N)
mu.med<-xts(c(dat_csv$totale_attualmente_positivi,I_seir_med*N),order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Veneto Region: Scenario 3  (Credible Interval ",100*0.50,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p
```

![](draft_analysis_Veneto_new_files/figure-gfm/scenario%20plot%20-3.png)<!-- -->

The 3 scenarios show different numbers. If we consider the second
scenario, at the end of the 2 weeks (2020-03-31) the number of infected
is (5608.1314227).

In the next plot the cumulative number of infected.  
At the end of the 2 weeks (2020-03-31) the total number of COVID19 cases
is expected to be (1.019809810^{4}).
