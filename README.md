# misspi <a href='https://github.com/catstats/misspi'><img src='man/figures/logo_speed_final.png' align='right' height="150" width="130" /></a>
Missing Value Imputation in Parallel


  <!-- badges: start -->
  [![CRAN status](https://www.r-pkg.org/badges/version/misspi)](https://CRAN.R-project.org/package=misspi)
  [![](http://cranlogs.r-pkg.org/badges/grand-total/misspi?color=blue)](https://cran.r-project.org/package=misspi)
  [![](http://cranlogs.r-pkg.org/badges/last-month/misspi?color=red)](https://cran.r-project.org/package=misspi)
  [![CRAN downloads](https://cranlogs.r-pkg.org/badges/misspi)](https://CRAN.R-project.org/package=misspi)
  [![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
  [![DOI](https://zenodo.org/badge/706371189.svg)](https://zenodo.org/doi/10.5281/zenodo.12563060)
  <!-- badges: end --> 



## Install From R CRAN
install.packages("misspi")


## Tutorial
Please find a more detailed tutorial at catstats.github.io/misspi/


## Quick Start 
```r
data(toxicity, package = "misspi")
set.seed(0)
toxicity.miss <- missar(toxicity, 0.4, 0.2)
toxicity.impute <- misspi(toxicity.miss)
toxicity.impute
```


