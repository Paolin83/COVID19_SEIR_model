COVID19 - Forecast analysis
================
PG
3/11/2020

## The COVID dataset

The present analysis used the dataset on COVID19 updated in
<https://github.com/pcm-dpc/COVID-19>. We used a SEIR model to estimate
the prediction of this epidemic in Italy. We estimated the R0 parameter
by means of a simple linear regression ad in
<https://kingaa.github.io/clim-dis/parest/parest.html>

``` r
plot(dat_csv$data,dat_csv$totale_attualmente_positivi,ylab="Total Covid cases",xlab="Date")
```

![](draft_analysis_Italy_files/figure-gfm/plot%20data-1.png)<!-- -->

``` r
mean(diff(log(dat_csv$totale_attualmente_positivi)))
```

    ## [1] 0.2389459

The grow is exponential (rate of increase of about 25%). We estimate the
R0 parameter by means of a linear model. R0 indicates how contagious an
infectious disease is. It’s also referred to as the reproduction number.

\#estimate r0 \#see
<https://kingaa.github.io/clim-dis/parest/parest.html> \#calculate r0
based in the last XX observation

``` r
beta_vec<-NULL
sd_vec<-NULL
for (i in 1:(days-6)){
fit <- lm(log(totale_attualmente_positivi)~t,data=dat_csv[i:days,])
beta_vec<-c(beta_vec,coef(fit)[2])
sd_vec<-c(sd_vec,coef(summary(fit))[2,2])
}

label <-days:7
mean  <- (beta_vec*14+1)
lower <- ((beta_vec-1.96*sd_vec)*14+1)
upper <- ((beta_vec+1.96*sd_vec)*14+1)

df <- data.frame(label, mean, lower, upper)
df$label <- factor(df$label, levels=rev(df$label))

library(ggplot2)
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

``` r
fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Last days used for the calculation") + ylab("R0 Mean (95% CI)") +
  theme_bw()  # use a white background
print(fp)
```

![](draft_analysis_Italy_files/figure-gfm/R0%20trend-1.png)<!-- -->

  - <img src="https://latex.codecogs.com/gif.latex?log(Y_t) = a + b*t + e_t" />

where \(`Y_t`\) is the cumulative number of infected at the time t,
while b is beta, the slope of the regression line.  
We calculate the R0 values based on a different number of days before
the last. The R0 shows a decreasing trend, but it is stable in the last
12 days.

The slope b indicate the rate of exponetial increase, used for the R0
calculation.  
R code:  
fit1 \<-
    lm(log(totale\_attualmente\_positivi)~t,data=dat\_csv)

``` r
head(dat_csv)
```

    ##                  data stato ricoverati_con_sintomi terapia_intensiva
    ## 1 2020-02-24 18:00:00   ITA                    101                26
    ## 2 2020-02-25 18:00:00   ITA                    114                35
    ## 3 2020-02-26 18:00:00   ITA                    128                36
    ## 4 2020-02-27 18:00:00   ITA                    248                56
    ## 5 2020-02-28 18:00:00   ITA                    345                64
    ## 6 2020-02-29 18:00:00   ITA                    401               105
    ##   totale_ospedalizzati isolamento_domiciliare totale_attualmente_positivi
    ## 1                  127                     94                         221
    ## 2                  150                    162                         311
    ## 3                  164                    221                         385
    ## 4                  304                    284                         588
    ## 5                  409                    412                         821
    ## 6                  506                    543                        1049
    ##   nuovi_attualmente_positivi dimessi_guariti deceduti totale_casi tamponi
    ## 1                        221               1        7         229    4324
    ## 2                         90               1       10         322    8623
    ## 3                         74               3       12         400    9587
    ## 4                        203              45       17         650   12014
    ## 5                        233              46       21         888   15695
    ## 6                        228              50       29        1128   18661
    ##   t
    ## 1 1
    ## 2 2
    ## 3 3
    ## 4 4
    ## 5 5
    ## 6 6

``` r
fit1 <- lm(log(totale_attualmente_positivi)~t,data=dat_csv)
summary(fit1)
```

    ## 
    ## Call:
    ## lm(formula = log(totale_attualmente_positivi) ~ t, data = dat_csv)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.28204 -0.15480  0.03257  0.10196  0.27824 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 5.446063   0.079891   68.17  < 2e-16 ***
    ## t           0.234140   0.007381   31.72 7.12e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1625 on 16 degrees of freedom
    ## Multiple R-squared:  0.9844, Adjusted R-squared:  0.9834 
    ## F-statistic:  1006 on 1 and 16 DF,  p-value: 7.115e-16

The prediction is good (\(R^2\) near 1). However the tendency is for a
reduction of the
slope.

``` r
plot(dat_csv$t,log(dat_csv$totale_attualmente_positivi),ylab="log cases",xlab="time")
abline(coef(summary(fit1))[,1])
```

![](draft_analysis_Italy_files/figure-gfm/model%20plot-1.png)<!-- -->
The slope coefficient estimated in the linear regression model can be
used to estimate R0.

R0=1+b\*incubation period.

The incubation period for the coronavirus is in mean 5.1 days with a
range from 2-14 days. Please see
<https://www.worldometers.info/coronavirus/coronavirus-incubation-period/>.
In the calculation we considered an incubation period of 14 days for two
reasons: 1) the majority of cases is asymptomatic, contagiousness is
greater than 5, maybe 14. A minority (who made the swab) will have a
duration of about 5 days between the start of contagiousness and swab;
2) 14 days is the worst scenario to consider.

``` r
slope <-coef(summary(fit1))[2,1]; slope
```

    ## [1] 0.2341403

``` r
slope.se <- coef(summary(fit1))[2,2]; slope.se
```

    ## [1] 0.007380637

``` r
### R0 estimates and 95% IC 
R_0=slope*14+1;R_0
```

    ## [1] 4.277964

``` r
(slope+c(-1,1)*1.96*slope.se)*14+1
```

    ## [1] 4.075439 4.480489

We want to make a short term forecast (14 days) with 3 scenario:  
\-Scenario 1: 10 exposed people for each COVID-19 case and beta the same
(no restrictions made or even no effects) -Scenario 2: 10 exposed people
for each COVID-19 case and beta reduced of 50% (-50% exposed people)
-Scenario 3: 5 exposed people for each COVID-19 case and beta reduced of
50% (-50% both exposed people and COVID19 contagious power)

We fix a series of initial parameters: -I0: initial number of COVID-19
cases  
\- R0: initial number of recovered  
\- beta: the quantity connected to R0  
\- N: Italian population  
\- duration: infection duration of COVID-19  
\- sigma0: the coronavirus transmission rate (half of flu epidemic)  
\- mu0: the overall mortality rate

``` r
# initial number of infectus
I0<-dat_csv$totale_attualmente_positivi[dim(dat_csv)[1]]; I0
```

    ## [1] 12839

``` r
# initial number of recovered
R0<-dat_csv$dimessi_guariti[dim(dat_csv)[1]]; R0
```

    ## [1] 1258

``` r
#beta 
beta0<-R_0/(14)
# italian poulation
N=60480000
# duration of COVID19 
duration<-14
#sigma0 is the coronavirus transmission rate fixed to 5%  (half of flu epidemic)
sigma0<-0.05
#mortality rate 
mu0<-1/(82*365.25) # 1/lifespan
```

We use the library(EpiDynamics)

``` r
library(EpiDynamics)

# average number of single connections of an infected person
# less contacts, less probability of new infections
# we keep constant the other parameters
forecast<-14
parameters <- c(mu = mu0, beta = beta0, sigma = sigma0, gamma = 1/duration)
f1<-10
initials <- c(S = 0.95, E = (f1*I0/N), I = I0/N, R = R0/N)
seir1 <- SEIR(pars = parameters, init = initials, time = 0:forecast)
parameters <- c(mu = mu0, beta = beta0, sigma = sigma0, gamma = 1/duration)
f2<-5
initials <- c(S = 0.95, E = (f2*I0/N), I = I0/N, R = R0/N)
seir2 <- SEIR(pars = parameters, init = initials, time = 0:forecast)
parameters <- c(mu = mu0, beta = beta0*1/2, sigma = sigma0, gamma = 1/duration)
f3<-5
initials <- c(S = 0.95, E = (f3*I0/N), I = I0/N, R = R0/N)
seir3 <- SEIR(pars = parameters, init = initials, time = 0:forecast)


date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+forecast-1+dim(dat_csv)[1],1)
plot(date,c(dat_csv$totale_attualmente_positivi,seir1$results$I[-1]*N),type="l",ylab="Cases",xlab="time",main="Infected")
lines(date,c(dat_csv$totale_attualmente_positivi,seir2$results$I[-1]*N),col=2)
lines(date,c(dat_csv$totale_attualmente_positivi,seir3$results$I[-1]*N),col=3)
lines(date[1:dim(dat_csv)[1]],dat_csv$totale_attualmente_positivi,lwd=2)
legend("topleft",c("first scenario","second scenario","third scenario"),lty=1,col=1:3)
```

![](draft_analysis_Italy_files/figure-gfm/scenario%20plot-1.png)<!-- -->
The 3 scenarios show how measures of restriction can help to reduce the
number of infected.  
At the end of the 2 weeks (2020-03-26) the number of infected is
(4.391582110^{4}). In the next plot the cumulative number of
infected.

``` r
plot(date,c(dat_csv$totale_casi,(seir1$results$I[-1]+seir1$results$R[-1])*N),type="l",ylab="Cases",xlab="time",main="Cumulative Infected")
lines(date,c(dat_csv$totale_casi,(seir2$results$I[-1]+seir2$results$R[-1])*N),col=2)
lines(date,c(dat_csv$totale_casi,(seir3$results$I[-1]+seir3$results$R[-1])*N),col=3)
lines(date[1:dim(dat_csv)[1]],(dat_csv$totale_casi),lwd=2)
legend("topleft",c("first scenario","second scenario","third scenario"),lty=1,col=1:3)
```

![](draft_analysis_Italy_files/figure-gfm/cumulative%20plot-1.png)<!-- -->

``` r
PlotMods(seir2)
```

![](draft_analysis_Italy_files/figure-gfm/cumulative%20plot-2.png)<!-- -->
