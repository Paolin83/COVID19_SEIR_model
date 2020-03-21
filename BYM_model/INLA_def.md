COVID19 - Forecast and predictions using a BYM model in Italy
================
Paolo Girardi
21 Marzo, 2020

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This
work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative
Commons Attribution-NonCommercial 4.0 International License</a>.

# Disclaimer

  - We want to investigate the evolution of the coronavirus pandemic in
    Italy

  - Our point of view is that of surveillance with the goal of detecting
    important changes in the underlying (random) process as soon as
    possible after it has occured.

  - We use data provided by Italian Civil Protection Department and the
    analysis was restricted to the Lombardy Region

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
checkpackage("ggplot2")
checkpackage("spdep")
checkpackage("maptools")
checkpackage("INLA")
checkpackage("rgdal")
```

and load them.

    ## Loading required package: zoo

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

    ## Loading required package: sp

    ## Loading required package: spData

    ## To access larger datasets in this package, install the spDataLarge
    ## package with: `install.packages('spDataLarge',
    ## repos='https://nowosad.github.io/drat/', type='source')`

    ## Loading required package: sf

    ## Linking to GEOS 3.6.1, GDAL 2.1.3, PROJ 4.9.3

    ## Checking rgeos availability: TRUE

    ## Loading required package: Matrix

    ## This is INLA_18.07.12 built 2018-07-12 11:07:12 UTC.
    ## See www.r-inla.org/contact-us for how to get help.
    ## To enable PARDISO sparse library; see inla.pardiso()

    ## rgdal: version: 1.4-4, (SVN revision 833)
    ##  Geospatial Data Abstraction Library extensions to R successfully loaded
    ##  Loaded GDAL runtime: GDAL 2.1.3, released 2017/20/01
    ##  Path to GDAL shared files: /Library/Frameworks/R.framework/Versions/3.6/Resources/library/rgdal/gdal
    ##  GDAL binary built with GEOS: FALSE 
    ##  Loaded PROJ.4 runtime: Rel. 4.9.3, 15 August 2016, [PJ_VERSION: 493]
    ##  Path to PROJ.4 shared files: /Library/Frameworks/R.framework/Versions/3.6/Resources/library/rgdal/proj
    ##  Linking to sp version: 1.3-1

# Datasets loading

    ## Warning: readShapePoly is deprecated; use rgdal::readOGR or sf::st_read

# Spatial analysis with R-INLA and BYM

We modelled COVID 19 cases by means of a BYM (Besag, York and Molli'e)
model using an Integrated Nested Laplace Approximation (INLA).  

For the \(i\)-th Nuts-3 Region (Italian province), the observed number
of COVID19, , was modelled as follows:  
\[
y_{i} \sim Poisson (\lambda_{i})\\
\] with i=1,…, 12. We modelled the incidence of COVID-19 by means of a
BYM model including the population size as offset variable as:  
\[
    log(\frac{\lambda_{i}}{N_{i}})=\alpha+\mu_i+\nu_i,\\
\] where \(\alpha\) is the intercept, \(\mu_i\) and \(\nu_i\) are two
area specific effects with normal distribution modelled using an
intrinsic conditional autoregressive structure (iCAR) and \(N_i\) is the
population size of each NUTS-3 Region.

The parameter were estimated by INLA (Integrated Nested Laplace
Approximation) and R software.

![](INLA_def_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> The number
of observed cases \(y_i\) is very different among NUTS-3 regions.

![](INLA_def_files/figure-gfm/unnamed-chunk-6-1.png)<!-- --> A correct
comparison is the made by incidence cases \(y_i/N_i\), here riported for
x1000 inhabitants. \#Spatial model  
We estimate a BYM model fixing prior distributions for iid and besag
components as follows:

We create new variables that are required by INLA procedure.

    ## 
    ## Call:
    ## c("inla(formula = formula.bym, family = \"poisson\", data = dat_csv, ",  "    E = pop, control.compute = list(dic = T))")
    ## 
    ## Time used:
    ##  Pre-processing    Running inla Post-processing           Total 
    ##          2.3680          1.1186          0.1320          3.6187 
    ## 
    ## Fixed effects:
    ##               mean     sd 0.025quant 0.5quant 0.975quant    mode    kld
    ## (Intercept) -9.544 0.0688    -9.6881   -9.544    -9.4001 -9.5439 0.0023
    ## 
    ## Random effects:
    ## Name   Model
    ##  ID   BYM model 
    ## 
    ## Model hyperparameters:
    ##                                         mean      sd 0.025quant 0.5quant
    ## Precision for ID (iid component)     96.9836 98.0821     6.6029  67.8352
    ## Precision for ID (spatial component)  0.5221  0.2109     0.2103   0.4902
    ##                                      0.975quant    mode
    ## Precision for ID (iid component)        356.824 18.0201
    ## Precision for ID (spatial component)      1.026  0.4255
    ## 
    ## Expected number of effective parameters(std dev): 12.10(0.013)
    ## Number of equivalent replicates : 25.79 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 23766.41
    ## Deviance Information Criterion (DIC, saturated) ....: 22397.55
    ## Effective number of parameters .....................: 12.11
    ## 
    ## Marginal log-Likelihood:  -11929.22 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ##  [1] "Bergamo"               "Brescia"              
    ##  [3] "Como"                  "Cremona"              
    ##  [5] "Lecco"                 "Lodi"                 
    ##  [7] "Mantova"               "Milano"               
    ##  [9] "Monza e della Brianza" "Pavia"                
    ## [11] "Sondrio"               "Varese"

The image reported the IRR=exp(\(\mu_i\)), the increase of the Incidence
of Covid-19 in each NUTS-3 Region respect to the overall mean.

\#Temporal model  
We consider the temporal aspect estimating a Bayesian RW2 model to model
the temporal trend of COVID-19 in the reported temporal window. The
observed number of COVID-19 a the time \(t\), , is modelled as follows  
\[
y_{t} \sim Poisson (\lambda_{t})
\] with t=1,…,26. \\end{center} where the quantity
\(\frac{\lambda_{t}}\), is modelled by a Random Walk of order 2 variable
\[
log({\lambda_{it}})=\alpha+\gamma_t+\phi_t, 
\] where \(\alpha\) is the intercept and \(\gamma_t \sim RW2\) are the
coefficients related to the random walk process while
\(\phi_t \sim N(0, \tau_t^{-1} )\) are the temporal specific random
errors. ![](INLA_def_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

    ## 
    ## Call:
    ## c("inla(formula = formula_t, family = \"poisson\", data = dat_csv, ",  "    control.compute = list(dic = T))")
    ## 
    ## Time used:
    ##  Pre-processing    Running inla Post-processing           Total 
    ##          1.9732          8.0011          0.3042         10.2785 
    ## 
    ## Fixed effects:
    ##               mean     sd 0.025quant 0.5quant 0.975quant   mode kld
    ## (Intercept) 1.7204 0.3223     1.0866   1.7208     2.3519 1.7215   0
    ## 
    ## Random effects:
    ## Name   Model
    ##  t   RW2 model 
    ## t2   IID model 
    ## 
    ## Model hyperparameters:
    ##                      mean        sd 0.025quant 0.5quant 0.975quant    mode
    ## Precision for t  5844.389 6512.5356   173.8706 3690.171  23333.545 326.415
    ## Precision for t2    1.883    0.6434     0.9523    1.771      3.448   1.572
    ## 
    ## Expected number of effective parameters(std dev): 27.16(0.2187)
    ## Number of equivalent replicates : 11.49 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 24847.07
    ## Deviance Information Criterion (DIC, saturated) ....: 23478.21
    ## Effective number of parameters .....................: 25.80
    ## 
    ## Marginal log-Likelihood:  -12564.51 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->
We now consider a BYM model with spatial and temporal components and
their complete interactions (Type IV). The observed number of COVID-19,
, was modelled as: \[
y_{it} \sim Poisson (\lambda_{it})\\
\] with i=1,..,12, and t=1,…,26. We modelled the ratio
\(\frac{\lambda_{it}}{N_{i}}\), the Incidence of COVID-19, with a BYM
model including temporal covariates as follows  
\[
log(\frac{\lambda_{it}}{N_{i}})=\alpha+\mu_i+\nu_i+\gamma_t+\phi_t+\delta_{it},\\
\]

where \(\alpha\) is the intercept and:  
\*space: \(\mu_i\) and \(\nu_i\) are two area specific effects with
normal distribution modelled using an intrinsic conditional
autoregressive structure (iCAR);

\*time: \(\gamma_t \sim RW2\) are the coefficients related to the random
walk process, while \(\phi_t \sim N(0, \tau_t^{-1} )\) are the temporal
specific random errors;

\*space-time: \(\delta_{it}\) are the coefficients related space-time
interactions that in the its type 4 formulation (please see Spatial and
Spatio-temporal Bayesian Models with R - INLA, Blangiardo and Cameletti)
is made by the Kronecker product of time and space indexes.

    ## 
    ## Call:
    ## c("inla(formula = formula.intIV, family = \"poisson\", data = dat_csv, ",  "    E = pop, control.compute = list(dic = T))")
    ## 
    ## Time used:
    ##  Pre-processing    Running inla Post-processing           Total 
    ##          8.4928        132.1726          0.1547        140.8201 
    ## 
    ## Fixed effects:
    ##               mean     sd 0.025quant 0.5quant 0.975quant    mode   kld
    ## (Intercept) -4.316 0.4163    -5.1577  -4.3063    -3.5228 -4.2831 1e-04
    ## 
    ## Random effects:
    ## Name   Model
    ##  ID   BYM model 
    ## t   RW2 model 
    ## t2   IID model 
    ## ID2   Besags ICAR model 
    ## 
    ## Model hyperparameters:
    ##                                           mean        sd 0.025quant
    ## Precision for ID (iid component)     1849.0197 1852.3643    120.821
    ## Precision for ID (spatial component) 1865.1563 1837.7823    126.591
    ## Precision for t                      5477.0060 5471.8849    434.543
    ## Precision for t2                        0.2659    0.1702      0.084
    ## Precision for ID2                       0.0106    0.0015      0.008
    ##                                       0.5quant 0.975quant      mode
    ## Precision for ID (iid component)     1298.0208  6.747e+03  326.5889
    ## Precision for ID (spatial component) 1322.9358  6.725e+03  345.9891
    ## Precision for t                      3873.2893  2.010e+04 1218.5206
    ## Precision for t2                        0.2200  7.139e-01    0.1597
    ## Precision for ID2                       0.0105  1.370e-02    0.0103
    ## 
    ## Expected number of effective parameters(std dev): 237.38(3.30)
    ## Number of equivalent replicates : 1.314 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 2023.75
    ## Deviance Information Criterion (DIC, saturated) ....: 654.90
    ## Effective number of parameters .....................: 230.93
    ## 
    ## Marginal log-Likelihood:  -2470.14 
    ## Posterior marginals for linear predictor and fitted values computed

The time and space component mean estimates can be extracted here

When can represent the
IRR(=exp(\(\mu_i+\frac{1}{T} \sum_{i=t}^{T} \delta_{it}\))) for each
NUTS-3 regions.

![](INLA_def_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Trend for each
province.

![](INLA_def_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#predictions 3 days forward of new
CODID 19 cases

![](INLA_def_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

We modelled the number of cases.

In the rw2 I set constr = FALSE is set to FALSE and that, for this
reason, the intercept is not included in the linear predictor.

    ## 
    ## Call:
    ## c("inla(formula = formula.intIVn, family = \"poisson\", data = dat_csv_n, ",  "    control.compute = list(dic = T), control.predictor = list(link = 1))" )
    ## 
    ## Time used:
    ##  Pre-processing    Running inla Post-processing           Total 
    ##          3.3433        196.7985          0.2997        200.4416 
    ## 
    ## Fixed effects:
    ##                  mean     sd 0.025quant 0.5quant 0.975quant   mode kld
    ## (Intercept)    0.6981 0.3519     0.0070   0.6981     1.3882 0.6982   0
    ## log(hubei + 1) 0.0852 0.0311     0.0242   0.0852     0.1463 0.0851   0
    ## log(pop)       0.0554 0.0264     0.0035   0.0554     0.1072 0.0554   0
    ## 
    ## Random effects:
    ## Name   Model
    ##  t   RW2 model 
    ## t2   IID model 
    ## over   IID model 
    ## ID   BYM model 
    ## ID2   Besags ICAR model 
    ## 
    ## Model hyperparameters:
    ##                                          mean        sd 0.025quant
    ## Precision for t                      2667.057 8286.3146     61.891
    ## Precision for t2                        2.180    0.8377      0.957
    ## Precision for over                      4.274    0.5635      3.272
    ## Precision for ID (iid component)     1849.233 1830.8839    127.562
    ## Precision for ID (spatial component) 1846.158 1829.6823    126.816
    ## Precision for ID2                      27.832   26.2410      3.783
    ##                                      0.5quant 0.975quant    mode
    ## Precision for t                       855.705  16477.224 141.525
    ## Precision for t2                        2.045      4.197   1.792
    ## Precision for over                      4.239      5.486   4.169
    ## Precision for ID (iid component)     1308.799   6688.380 348.477
    ## Precision for ID (spatial component) 1305.704   6683.763 345.977
    ## Precision for ID2                      20.299     96.913  10.140
    ## 
    ## Expected number of effective parameters(std dev): 230.65(3.613)
    ## Number of equivalent replicates : 1.353 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 1948.46
    ## Deviance Information Criterion (DIC, saturated) ....: 579.61
    ## Effective number of parameters .....................: 223.65
    ## 
    ## Marginal log-Likelihood:  -2075.68 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->
