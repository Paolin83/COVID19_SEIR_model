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

fit1 \<-
    lm(log(totale\_casi)~t,data=dat\_csv)

``` r
show(dat_csv)
```

    ##                   data stato ricoverati_con_sintomi terapia_intensiva
    ## 1  2020-02-24 18:00:00   ITA                    101                26
    ## 2  2020-02-25 18:00:00   ITA                    114                35
    ## 3  2020-02-26 18:00:00   ITA                    128                36
    ## 4  2020-02-27 18:00:00   ITA                    248                56
    ## 5  2020-02-28 18:00:00   ITA                    345                64
    ## 6  2020-02-29 18:00:00   ITA                    401               105
    ## 7  2020-03-01 18:00:00   ITA                    639               140
    ## 8  2020-03-02 18:00:00   ITA                    742               166
    ## 9  2020-03-03 18:00:00   ITA                   1034               229
    ## 10 2020-03-04 18:00:00   ITA                   1346               295
    ## 11 2020-03-05 18:00:00   ITA                   1790               351
    ## 12 2020-03-06 18:00:00   ITA                   2394               462
    ## 13 2020-03-07 18:00:00   ITA                   2651               567
    ## 14 2020-03-08 18:00:00   ITA                   3557               650
    ## 15 2020-03-09 18:00:00   ITA                   4316               733
    ## 16 2020-03-10 18:00:00   ITA                   5038               877
    ##    totale_ospedalizzati isolamento_domiciliare totale_attualmente_positivi
    ## 1                   127                     94                         221
    ## 2                   150                    162                         311
    ## 3                   164                    221                         385
    ## 4                   304                    284                         588
    ## 5                   409                    412                         821
    ## 6                   506                    543                        1049
    ## 7                   779                    798                        1577
    ## 8                   908                    927                        1835
    ## 9                  1263                   1000                        2263
    ## 10                 1641                   1065                        2706
    ## 11                 2141                   1155                        3296
    ## 12                 2856                   1060                        3916
    ## 13                 3218                   1843                        5061
    ## 14                 4207                   2180                        6387
    ## 15                 5049                   2936                        7985
    ## 16                 5915                   2599                        8514
    ##    nuovi_attualmente_positivi dimessi_guariti deceduti totale_casi tamponi
    ## 1                         221               1        7         229    4324
    ## 2                          90               1       10         322    8623
    ## 3                          74               3       12         400    9587
    ## 4                         203              45       17         650   12014
    ## 5                         233              46       21         888   15695
    ## 6                         228              50       29        1128   18661
    ## 7                         528              83       34        1694   21127
    ## 8                         258             149       52        2036   23345
    ## 9                         428             160       79        2502   25856
    ## 10                        443             276      107        3089   29837
    ## 11                        590             414      148        3858   32362
    ## 12                        620             523      197        4636   36359
    ## 13                       1145             589      233        5883   42062
    ## 14                       1326             622      366        7375   49937
    ## 15                       1598             724      463        9172   53826
    ## 16                        529            1004      631       10149   60761
    ##     t
    ## 1   1
    ## 2   2
    ## 3   3
    ## 4   4
    ## 5   5
    ## 6   6
    ## 7   7
    ## 8   8
    ## 9   9
    ## 10 10
    ## 11 11
    ## 12 12
    ## 13 13
    ## 14 14
    ## 15 15
    ## 16 16

``` r
fit1 <- lm(log(totale_casi)~t,data=dat_csv)
summary(fit1)
```

    ## 
    ## Call:
    ## lm(formula = log(totale_casi) ~ t, data = dat_csv)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.23731 -0.10073  0.02064  0.09894  0.24917 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 5.414859   0.079782   67.87  < 2e-16 ***
    ## t           0.252974   0.008251   30.66 3.09e-14 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1521 on 14 degrees of freedom
    ## Multiple R-squared:  0.9853, Adjusted R-squared:  0.9843 
    ## F-statistic:   940 on 1 and 14 DF,  p-value: 3.089e-14

Estimates

You can also embed plots, for example:

``` r
plot(dat_csv$t,log(dat_csv$totale_casi),ylab="log cases",xlab="time")
abline(coef(summary(fit1))[,1])
```

![](draft_analysis_files/figure-gfm/model%20plot-1.png)<!-- --> The
slope coefficient estimated in the linear regression model can be used
to estimate R0, the number of

``` r
slope <-coef(summary(fit1))[2,1]; slope
```

    ## [1] 0.2529736

``` r
slope.se <- coef(summary(fit1))[2,2]; slope.se
```

    ## [1] 0.008250869

``` r
### R0 estimates and 95%IC 
R_0=slope*14+1;R_0
```

    ## [1] 4.54163

``` r
(slope+c(-1,1)*1.96*slope.se)*14+1
```

    ## [1] 4.315226 4.768034

We want to make a short term forecast (14 days) with 3 scenario. We fix
a series of initial parameters.

``` r
# initial number of infectus
I0<-max(dat_csv$totale_casi); I0
```

    ## [1] 10149

``` r
# initial number of recovered
R0<-max(dat_csv$dimessi_guariti); R0
```

    ## [1] 1004

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
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

``` r
parameters <- c(mu = mu0, beta = beta0, sigma = sigma0, gamma = 1/duration)
# average number of single connections of an infected person
# less contacts, less probability of new infections
# we keep constant the other parameters
f1<-10
initials <- c(S = 0.95, E = (f1*I0/N), I = I0/N, R = R0/N)
seir1 <- SEIR(pars = parameters, init = initials, time = 0:14)

f2<-5
initials <- c(S = 0.95, E = (f2*I0/N), I = I0/N, R = R0/N)
seir2 <- SEIR(pars = parameters, init = initials, time = 0:14)

f3<-3
initials <- c(S = 0.95, E = (f3*I0/N), I = I0/N, R = R0/N)
seir3 <- SEIR(pars = parameters, init = initials, time = 0:14)



plot(c(dat_csv$totale_casi,seir1$results$I[-1]*N),type="l",ylab="Number of infectus",xlab="time")
lines(c(dat_csv$totale_casi,seir2$results$I[-1]*N),col=2)
lines(c(dat_csv$totale_casi,seir3$results$I[-1]*N),col=3)
legend("topleft",c("first scenario","second scenario","third scenario"),lty=1,col=1:3)
```

![](draft_analysis_files/figure-gfm/first%20scenario%20plot-1.png)<!-- -->
