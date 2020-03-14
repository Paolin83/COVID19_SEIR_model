COVID19 - Forecast and predictions using a time dependent SEIR model
================
PG
3/13/2020

## The COVID dataset

The present analysis started from the dataset on COVID19 updated in
<https://github.com/pcm-dpc/COVID-19>, database provided by the Italian
Civil Protection.

``` r
rm(list=ls())
###import italian dataset updated 13 March 2020
dat_csv<-read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv",header=T)
days<-dim(dat_csv)[1]
dat_csv$t<-1:days
#database contains  the following variables
names(dat_csv)
```

    ##  [1] "data"                        "stato"                      
    ##  [3] "ricoverati_con_sintomi"      "terapia_intensiva"          
    ##  [5] "totale_ospedalizzati"        "isolamento_domiciliare"     
    ##  [7] "totale_attualmente_positivi" "nuovi_attualmente_positivi" 
    ##  [9] "dimessi_guariti"             "deceduti"                   
    ## [11] "totale_casi"                 "tamponi"                    
    ## [13] "t"

``` r
# The total number of epidemic day is
days
```

    ## [1] 20

``` r
plot(as.Date(substr(dat_csv$data,1,10)),dat_csv$totale_attualmente_positivi,ylab=" Infected",xlab="Date",type="l")
points(as.Date(substr(dat_csv$data,1,10)),dat_csv$totale_attualmente_positivi)
```

![](draft_analysis_Italy_files/figure-gfm/setup-1.png)<!-- -->

With the aim of predicting the future number of COVID19 cases on the
basis of the actual data, we used a SEIR model applied to the COVID19
epidemic in Italy. A first parameter to estimate is R0, the reproduction
number (<https://en.wikipedia.org/wiki/Basic_reproduction_number>). R0
indicates how contagious an infectious disease is. It is also referred
to as “the reproduction number” of COVID19.

The R0 value has been estimated by means of a linear regression model as
reported in <https://kingaa.github.io/clim-dis/parest/parest.html>.

The following plot presented the actual status.

``` r
library(tidyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

``` r
df <- dat_csv %>%
  select(data, ricoverati_con_sintomi, terapia_intensiva, 
         totale_ospedalizzati, isolamento_domiciliare, 
         totale_attualmente_positivi, nuovi_attualmente_positivi, 
         dimessi_guariti, deceduti, totale_casi) %>%
  gather(key = "variable", value = "value", -data)
head(df, 3)
```

    ##                  data               variable value
    ## 1 2020-02-24 18:00:00 ricoverati_con_sintomi   101
    ## 2 2020-02-25 18:00:00 ricoverati_con_sintomi   114
    ## 3 2020-02-26 18:00:00 ricoverati_con_sintomi   128

``` r
ggplot(df, aes(x = as.Date(data), y = value)) + 
  geom_line(aes(color = variable), size = 1) 
```

![](draft_analysis_Italy_files/figure-gfm/plots-1.png)<!-- -->

The plot shows an exponential grow of cases.\\

We estimate the R0 parameter by means of a linear model.

Y\_t= a + beta \* t +e\_t

where \(`Y_t`\) is the number of infected at the time t, while b is
beta, the slope of the regression line.

The slope coefficient is used to estimate R0 as in the following
formula:

R0=beta\*incubation period.

The incubation period for the coronavirus is in mean 5.1 days with a
range from 2-14 days. Please see
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

We calculate several R0 values, each one based on a different number of
days before the last day, in order to assess if the R0 trend is
decreasing (how is expected to be).

``` r
#calculate r0 based in the last XX observation 
#vector for beta and standard deviation
beta_vec<-NULL
sd_vec<-NULL
#for cycle for all days to only the last 5.
for (i in 1:(days-4)){
fit <- lm(log(totale_attualmente_positivi)~t,data=dat_csv[i:days,])
beta_vec<-c(beta_vec,coef(fit)[2])
sd_vec<-c(sd_vec,coef(summary(fit))[2,2])
}


label <-days:5

mean  <- (beta_vec*14)
lower <- ((beta_vec-1.96*sd_vec)*14)
upper <- ((beta_vec+1.96*sd_vec)*14)

df <- data.frame(label, mean, lower, upper)[order(label),]
df$label <- factor(df$label, levels=rev(df$label))

library(ggplot2)
fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Last days used for the calculation") + ylab("R0 Mean (95% CI)") +
  theme_bw() 
print(fp)
```

![](draft_analysis_Italy_files/figure-gfm/R0%20trend-1.png)<!-- -->

The R0 shows a decreasing trend in the last period. We use the estimated
trend between R0 and time to calculate the future R0 value for the next
14 days. We predict beta (and R0) for the next 14 days assuming a Gamma
distribution for the beta (the slope) forcing its value to be greater
than 0. The trend was not monotonic, we use a simple splines to increase
the fitting of the model to the data.

``` r
library(splines)
time<-1:(days-4)
# splines with only 3 basis, weighted for standard deviation of estimates.
beta.model<-glm(beta_vec~bs(time,3),weights = 1/sd_vec,family=Gamma)
forecast=14
# add 'fit', 'lwr', and 'upr' columns to dataframe (generated by predict)
pre<-predict(beta.model,type='response',newdata=data.frame(time=1:(days+forecast)),se.fit=TRUE)
```

    ## Warning in bs(time, degree = 3L, knots = numeric(0), Boundary.knots =
    ## c(1L, : some 'x' values beyond boundary knots may cause ill-conditioned
    ## bases

``` r
date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+forecast-1+dim(dat_csv)[1],1)
beta.predict <- data.frame(beta_vec=c(beta_vec,rep(NA,forecast+4)),time=date,fit=pre$fit,lwr=pre$fit-1*1.96*pre$se.fit,upr=pre$fit+1*1.96*pre$se.fit)

r0.predict<-beta.predict
r0.predict[,c(1,3:5)]<-r0.predict[,c(1,3:5)]*14
# plot the points (actual observations), regression line, and confidence interval
p <- ggplot(r0.predict, aes(date,beta_vec))
p <- p + geom_point() +labs(x="Date",y="R0 value") 
p <- p + geom_line(aes(date,fit))
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3)
p
```

    ## Warning: Removed 18 rows containing missing values (geom_point).

![](draft_analysis_Italy_files/figure-gfm/R0%20forecast-1.png)<!-- -->

R0 passes from a value of 3.16 in the initial phase to an estimated
value of 1.87 at the ending of the 14-days forecast.

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

    ## [1] 17750

``` r
# initial number of recovered
R0<-dat_csv$dimessi_guariti[dim(dat_csv)[1]]; R0
```

    ## [1] 1966

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
scenario, at the end of the 2 weeks (2020-03-28) the number of infected
is (4.581521310^{4}).

In the next plot the cumulative number of infected.  
At the end of the 2 weeks (2020-03-28) the total number of COVID19 cases
is expected to be
(7.344144410^{4}).

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
