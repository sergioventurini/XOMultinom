# XOMultinom

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/XOMultinom)](https://cran.r-project.org/package=XOMultinom)
[![](http://cranlogs.r-pkg.org/badges/grand-total/XOMultinom?color=blue)](https://cran.r-project.org/package=XOMultinom)
[![Travis build status](https://travis-ci.org/sergioventurini/XOMultinom.svg?branch=master)](https://travis-ci.org/sergioventurini/XOMultinom)

<!-- badges: end -->

## Overview

###### Current release: 0.6.4
###### R version required: at least 3.6.0
`R` package for computing the exact distributions of some functions of the
ordered multinomial counts.

## Installation

Since the package requires some code to be compiled, you need a working C++
compiler. To get it:

- On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
- On Mac, install Xcode from the app store.
- On Linux, `sudo apt-get install r-base-dev` or similar.

Then, the easiest way to get the package is to install it from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("sergioventurini/XOMultinom")
```

See the main help page of the package by executing `?XOMultinom` or run any of
the demos available in the package by executing the code
`demo(demo_name, package = "XOMultinom")`.

## Authors
Sergio Venturini, Department of Economic and Social Sciences, Università Cattolica del Sacro Cuore, Cremona, Italy

E-mail: sergio.venturini@unicatt.it

Marco Bonetti, Department of Social and Political Sciences, Bocconi University, Milan, Italy

E-mail: marco.bonetti@unibocconi.it
