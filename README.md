
# netReg <img src="https://rawgit.com/dirmeier/netReg/master/inst/sticker/sticker.png" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/dirmeier/netReg.svg?branch=master)](https://travis-ci.org/dirmeier/netReg)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dirmeier/netReg?branch=master&svg=true)](https://ci.appveyor.com/project/dirmeier/netReg)
[![codecov](https://codecov.io/gh/dirmeier/netReg/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/netReg)
[![bioc](https://bioconductor.org/shields/years-in-bioc/netReg.svg)](https://bioconductor.org/packages/release/bioc/html/netReg.html)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/netreg/README.html)

Network-penalized generalized linear models in `R` and `C++`.

## Introduction

`netReg` is an R/C++ implementation of a network-regularized linear regression model.
Graph prior knowledge, in the form of biological networks,
is being incorporated into the loss function of the linear model
which allows better estimation of regression coefficients. 

Depending on your installed libraries `netReg` uses `OpenBLAS` or `BLAS` and `Lapack` for
fast computation of matrix operations in an `Armadillo\RcppArmadillo` framework. 
We use `Dlib` in order to calculate the most optimal set of shrinkage parameters using k-fold cross-validation.

For instance, using R, you could fit a network-regularized model like that:

``` r
> X <- matrix(rnorm(10*5), 10)
> Y <- matrix(rnorm(10*5), 10)
> aff.mat <- abs(rWishart(1, 10, diag(5))[,,1])

> fit <- edgenet(X=X, Y=Y, G.X=aff.mat, 
                 lambda=1, psigx=1, family="gaussian")
> print(fit)

#>Call: edgenet.default(X = X, Y = Y, G.X = aff.mat, lambda = 1, psigx = 1, 
#>                      family = "gaussian")

#>Coefficients:
#>             [,1]       [,2]       [,3]      [,4]        [,5]
#>[1,]  0.015999757  0.0000000  0.1614923 0.1200417  0.11848110
#>[2,]  0.023446755 -0.2021715  0.0000000 0.4316659  0.00000000
#>[3,] -0.139861486  0.0000000  0.2917003 0.3362468  0.82524790
#>[4,]  0.000000000  0.4495721 -0.2507407 0.3454617 -0.01794482
#>[5,] -0.007317568  0.0000000  0.2590443 0.0000000 -0.10663651

#>Intercept:
            [,1]
#>[1,]  0.18825966
#>[2,] -0.36970048
#>[3,]  0.27555425
#>[4,]  0.04295871
#>[5,]  0.25898488

#>Parameters:
#>lambda psi_gx psi_gy 
#>     1      1      0 

#>Family:
#>[1] "gaussian"
```

## Installation

`netReg` comes as a stand alone `C++` command line tool shipped with [`bioconda`](https://anaconda.org/bioconda/netreg) as well as a [`Bioconductor`](https://bioconductor.org/packages/release/bioc/html/netReg.html) package.

### Installation of the R package

You can install and use `netReg` either as an `R` library from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/netReg.html),
or by downloading the [tarball](https://github.com/dirmeier/netReg/releases).

If you want to use the **recommended way** using Bioconductor just call:

```r
> source("https://bioconductor.org/biocLite.R")
> biocLite("netReg")
  
> library(netReg)
```
 
from the R-console. Installing the R package using the downloaded tarball works like this:

```bash
$ R CMD install <netReg-x.y.z.tar.gz>
```

I **do not** recommend using `devtools`, but preferring tarballed releases over the most recent commit.

### Installation of the command-line tool

You can install the `C++` command line tool using `conda`. For that you should download [Anaconda](https://www.continuum.io/downloads) and create a [virtual environment](https://conda.io/docs/using/envs.html).
Then install the tool using:

```sh
$ conda install -c bioconda netreg
  
$ netReg -h
```

You can find more instructions for (manual) installation [here](https://dirmeier.github.io/netReg/articles/netReg_commandline.html).

## Documentation

* For the R package: load the package using `library(netReg)`. We provide a vignette for the package that can be called using: `vignette("netReg")`. 
* For the command-line tool: You can also use the online [tutorial](https://dirmeier.github.io/netReg/articles/netReg_R.html).
  Have a look at the command line [tutorial](https://dirmeier.github.io/netReg/articles/netReg_commandline.html).

## Citation

It would be great if you cited `netReg` like this:

Simon Dirmeier, Christiane Fuchs, Nikola S Mueller, Fabian J Theis; 
netReg: network-regularized linear models for biological association studies, 
*Bioinformatics*, Volume 34, Issue 5, 1 March 2018, Pages 896–898, https://doi.org/10.1093/bioinformatics/btx677

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
