COVID19 - Forecast and predictions using a BYM model in Italy
================
Paolo Girardi
25 Marzo, 2020

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
    ##          2.4459          1.4682          0.2672          4.1814 
    ## 
    ## Fixed effects:
    ##                mean     sd 0.025quant 0.5quant 0.975quant    mode   kld
    ## (Intercept) -9.3276 0.0657    -9.4646  -9.3276    -9.1907 -9.3276 6e-04
    ## 
    ## Random effects:
    ## Name   Model
    ##  ID   BYM model 
    ## 
    ## Model hyperparameters:
    ##                                         mean      sd 0.025quant 0.5quant
    ## Precision for ID (iid component)     96.7557 96.7986     6.7748  68.1280
    ## Precision for ID (spatial component)  0.6252  0.2534     0.2512   0.5866
    ##                                      0.975quant    mode
    ## Precision for ID (iid component)        354.651 18.5563
    ## Precision for ID (spatial component)      1.231  0.5085
    ## 
    ## Expected number of effective parameters(std dev): 12.15(0.0103)
    ## Number of equivalent replicates : 29.62 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 30849.47
    ## Deviance Information Criterion (DIC, saturated) ....: 29172.81
    ## Effective number of parameters .....................: 12.16
    ## 
    ## Marginal log-Likelihood:  -15475.91 
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
\] with t=1,…,30. \\end{center} where the quantity
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
    ##          1.8884          6.7836          0.1213          8.7934 
    ## 
    ## Fixed effects:
    ##               mean     sd 0.025quant 0.5quant 0.975quant   mode kld
    ## (Intercept) 1.8219 0.2994     1.2335   1.8221     2.4086 1.8226   0
    ## 
    ## Random effects:
    ## Name   Model
    ##  t   RW2 model 
    ## t2   IID model 
    ## 
    ## Model hyperparameters:
    ##                      mean        sd 0.025quant 0.5quant 0.975quant    mode
    ## Precision for t  3752.477 7307.1653     203.54 1750.058  19800.319 507.374
    ## Precision for t2    2.262    0.7552       1.09    2.163      4.027   1.969
    ## 
    ## Expected number of effective parameters(std dev): 31.20(0.2224)
    ## Number of equivalent replicates : 11.54 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 32450.51
    ## Deviance Information Criterion (DIC, saturated) ....: 30773.85
    ## Effective number of parameters .....................: 29.81
    ## 
    ## Marginal log-Likelihood:  -16398.12 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->
We now consider a BYM model with spatial and temporal components and
their complete interactions (Type IV). The observed number of COVID-19,
, was modelled as: \[
y_{it} \sim Poisson (\lambda_{it})\\
\] with i=1,..,12, and t=1,…,30. We modelled the ratio
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

    ## [1] 2373.575

    ## [1] 2423.39

    ## [1] 2391.416

    ## [1] 2408.876

The time and space component mean estimates can be extracted here

When can represent the
IRR(=exp(\(\mu_i+\frac{1}{T} \sum_{i=t}^{T} \delta_{it}\))) for each
NUTS-3 regions.

![](INLA_def_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Trend for each
province.

![](INLA_def_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](INLA_def_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

    ##           used  (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
    ## Ncells 2743715 146.6    4582575 244.8         NA  4582575 244.8
    ## Vcells 6959683  53.1   12255594  93.6      16384 10146329  77.5

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
    ##          2.7898         67.6722          0.2082         70.6702 
    ## 
    ## Fixed effects:
    ##               mean     sd 0.025quant 0.5quant 0.975quant   mode kld
    ## (Intercept) 0.6673 0.3182     0.0424   0.6674     1.2913 0.6676   0
    ## log(pop)    0.0556 0.0239     0.0085   0.0556     0.1025 0.0556   0
    ## hubei       0.0001 0.0000     0.0000   0.0001     0.0001 0.0001   0
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
    ## Precision for ID (iid component)        0.7939    0.3218     0.3198
    ## Precision for ID (spatial component) 1890.1248 1849.3391   126.0569
    ## Precision for t                       401.7140  398.2981    50.1432
    ## Precision for t2                        3.5775    1.4267     1.5520
    ## Precision for over                      1.6947    0.1721     1.3782
    ##                                       0.5quant 0.975quant    mode
    ## Precision for ID (iid component)        0.7445      1.565   0.645
    ## Precision for ID (spatial component) 1345.0353   6783.182 343.784
    ## Precision for t                       285.7903   1455.365 134.653
    ## Precision for t2                        3.3258      7.069   2.873
    ## Precision for over                      1.6873      2.055   1.674
    ## 
    ## Expected number of effective parameters(std dev): 305.10(2.552)
    ## Number of equivalent replicates : 1.18 
    ## 
    ## Deviance Information Criterion (DIC) ...............: 2374.92
    ## Deviance Information Criterion (DIC, saturated) ....: 698.26
    ## Effective number of parameters .....................: 294.65
    ## 
    ## Marginal log-Likelihood:  -1716.57 
    ## Posterior marginals for linear predictor and fitted values computed

![](INLA_def_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

\#Forecast for Lodi province

![](INLA_def_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->
