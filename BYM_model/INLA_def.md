COVID19 - Forecast and predictions using a BYM model in Italy
================
Paolo Girardi
30 Marzo, 2020

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

    ## Warning in dat_csv$totale_casi[dat_csv$data != "2020-02-24 18:00:00"] <-
    ## unlist(tapply(dat_csv$totale_casi, : il numero di elementi da sostituire
    ## non è un multiplo della lunghezza di sostituzione

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
    ##          3.0132          3.1908          0.3932          6.5972 
    ## 
    ## Fixed effects:
    ##                mean     sd 0.025quant 0.5quant 0.975quant    mode   kld
    ## (Intercept) -9.1279 0.0637    -9.2625  -9.1279    -8.9932 -9.1278 0.001
    ## 
    ## Random effects:
    ## Name   Model
    ##  ID   BYM model 
    ## 
    ## Model hyperparameters:
    ##                                         mean      sd 0.025quant 0.5quant
    ## Precision for ID (iid component)     94.8509 94.5801     6.5082  66.8721
    ## Precision for ID (spatial component)  0.9509  0.3915     0.3777   0.8897
    ##                                      0.975quant    mode
    ## Precision for ID (iid component)        346.818 17.8407
    ## Precision for ID (spatial component)      1.892  0.7669
    ## 
    ## Expected number of effective parameters(std dev): 12.22(0.0093)
    ## Number of equivalent replicates : 35.34 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 37088.56
    ## Deviance Information Criterion (DIC, saturated) ....: 34925.78
    ## Effective number of parameters .....................: 12.24
    ## 
    ## Marginal log-Likelihood:  -18598.49 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-8-1.png)<!-- --> The image
reported the IRR=exp(\(\mu_i\)), the increase of the Incidence of
Covid-19 in each NUTS-3 Region respect to the overall mean.

\#Temporal model  
We consider the temporal aspect estimating a Bayesian RW2 model to model
the temporal trend of COVID-19 in the reported temporal window. The
observed number of COVID-19 a the time \(t\), , is modelled as follows  
\[
y_{t} \sim Poisson (\lambda_{t})
\] with t=1,…,36. \\end{center} where the quantity
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
    ##          4.1249          8.1949          0.1426         12.4624 
    ## 
    ## Fixed effects:
    ##               mean     sd 0.025quant 0.5quant 0.975quant   mode kld
    ## (Intercept) 2.1492 0.2671     1.6248   2.1492     2.6732 2.1492   0
    ## 
    ## Random effects:
    ## Name   Model
    ##  t   RW2 model 
    ## t2   IID model 
    ## 
    ## Model hyperparameters:
    ##                    mean      sd 0.025quant 0.5quant 0.975quant   mode
    ## Precision for t  558.96 357.841     128.46   477.12    1475.74 323.84
    ## Precision for t2  28.19   7.976      15.43    27.25      46.49  25.47
    ## 
    ## Expected number of effective parameters(std dev): 36.24(0.2968)
    ## Number of equivalent replicates : 11.92 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 50638.53
    ## Deviance Information Criterion (DIC, saturated) ....: 48475.76
    ## Effective number of parameters .....................: 35.23
    ## 
    ## Marginal log-Likelihood:  -25513.77 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->
We now consider a BYM model with spatial and temporal components and
their complete interactions (Type IV). The observed number of COVID-19,
, was modelled as: \[
y_{it} \sim Poisson (\lambda_{it})\\
\] with i=1,..,12, and t=1,…,36. We modelled the ratio
\(\frac{\lambda_{it}}{N_{i}}\), the Incidence of COVID-19, with a BYM
model including temporal covariates as follows  
\[
log(\frac{\lambda_{it}}{N_{i}})=\alpha+\mu_i+\nu_i+\gamma_t+\phi_t+\delta_{it},\\
\]

where \(\alpha\) is the intercept and:  
\*space: \(\mu_i\) and \(\nu_i\) are two area specific effects with
normal distribution modelled using an intrinsic conditional
autoregressive structure (iCAR);

\*time: \(\gamma_t\) are the coefficients related to the time process,
while \(\phi_t \sim N(0, \tau_t^{-1} )\) are the temporal specific
random errors. Time process can be autoregressive or a random walk
process (order 1 or 2);

\*space-time: \(\delta_{it}\) are the coefficients related space-time
interactions.

There are many type of spatio-temporal interactions (please see Spatial
and Spatio-temporal Bayesian Models with R - INLA, Blangiardo and
Cameletti).

We explored four spatio-temporal interactions: -Type I: coefficients
\(\delta_{it}\) do not have a spatial or/and temporal structure
(\(\delta_{it}\sim N(0, 1/\tau_{\delta})\)). -TypeII: coefficients
\(\delta_{it}\) have a structured temporal main effect and the
unstructured spatial effect.  
\-TypeIII: coefficients \(\delta_{it}\) have a spatial structure (as
Besag), but are temporally unstructured. -TypeVI: coefficients
\(\delta_{it}\) have both a spatial and temporal structure. It is the
most complex type of space-time interaction.

    ## [1] 3073.61

    ## [1] 3087.757

    ## [1] 3072.782

    ## [1] 3072.805

The time and space component mean estimates can be extracted here

When can represent the
IRR(=exp(\(\mu_i+\frac{1}{T} \sum_{i=t}^{T} \delta_{it}\))) for each
NUTS-3 regions.

![](INLA_def_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Trend for each
province.

![](INLA_def_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

    ##           used  (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
    ## Ncells 2753205 147.1    4581944 244.8         NA  4581944 244.8
    ## Vcells 7162042  54.7   12255594  93.6      16384 10146329  77.5

![](INLA_def_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#predictions 7 days forward of new
CODID 19 cases

We modelled the number of cases.

In the rw2 I set constr = FALSE is set to FALSE and that, for this
reason, the intercept is not included in the linear predictor.

    ## 
    ## Call:
    ## c("inla(formula = formula.int_n, family = \"poisson\", data = dat_csv_n, ",  "    control.compute = list(dic = T), control.predictor = list(link = link))" )
    ## 
    ## Time used:
    ##  Pre-processing    Running inla Post-processing           Total 
    ##          3.1972         85.4394          0.7757         89.4124 
    ## 
    ## Fixed effects:
    ##               mean     sd 0.025quant 0.5quant 0.975quant   mode kld
    ## (Intercept) 0.7587 0.2926     0.1839   0.7587     1.3324 0.7589   0
    ## log(pop)    0.0635 0.0221     0.0201   0.0635     0.1067 0.0635   0
    ## hubei       0.0000 0.0000     0.0000   0.0000     0.0001 0.0000   0
    ## 
    ## Random effects:
    ## Name   Model
    ##  ID   BYM model 
    ## t   RW2 model 
    ## t2   IID model 
    ## over   IID model 
    ## 
    ## Model hyperparameters:
    ##                                           mean        sd 0.025quant
    ## Precision for ID (iid component)     1.147e+00 4.985e-01     0.4373
    ## Precision for ID (spatial component) 1.836e+03 1.825e+03   125.7771
    ## Precision for t                      5.572e+02 3.882e+02   109.3506
    ## Precision for t2                     1.836e+04 1.825e+04  1279.9079
    ## Precision for over                   5.613e-01 4.660e-02     0.4741
    ##                                       0.5quant 0.975quant      mode
    ## Precision for ID (iid component)     1.062e+00  2.362e+00    0.8980
    ## Precision for ID (spatial component) 1.296e+03  6.661e+03  343.8411
    ## Precision for t                      4.638e+02  1.561e+03  288.3685
    ## Precision for t2                     1.297e+04  6.655e+04 3524.8877
    ## Precision for over                   5.598e-01  6.575e-01    0.5574
    ## 
    ## Expected number of effective parameters(std dev): 412.69(1.509)
    ## Number of equivalent replicates : 1.047 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 3078.39
    ## Deviance Information Criterion (DIC, saturated) ....: 915.61
    ## Effective number of parameters .....................: 400.49
    ## 
    ## Marginal log-Likelihood:  -2353.42 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

\#Forecast for Lodi province

![](INLA_def_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->
