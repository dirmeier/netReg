
# netReg <img src="https://rawgit.com/dirmeier/netReg/master/inst/sticker/sticker.png" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/dirmeier/netReg.svg?branch=master)](https://travis-ci.org/dirmeier/netReg)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dirmeier/netReg?branch=master&svg=true)](https://ci.appveyor.com/project/dirmeier/netReg)
[![codecov](https://codecov.io/gh/dirmeier/netReg/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/netReg)
[![bioc](https://bioconductor.org/shields/years-in-bioc/netReg.svg)](https://bioconductor.org/packages/release/bioc/html/netReg.html)

Generalized linear models combined with network penalties. Now with `TensorFlow`.

## Introduction

`netReg` implements generalized linear models that use network penalties for regularization. Network regularization uses graphs or trees to incorporate information about interactions of covariables, or responses, into the loss function of a GLM. Ideally this allows better (i.e., lower variance) estimation of regression coefficients. 

From version `v1.9.0` on, we use `TensorFlow` instead of custom `C++` and `Dlib` which allowed deleting of major parts of the code base.
Some matrix operations are still done in `RcppArmadillo`. 

For instance, in R, you could fit a network-regularized model like that:

``` r
> X <- matrix(rnorm(10*5), 10)
> Y <- matrix(rnorm(10*5), 10)
> aff.mat <- abs(rWishart(1, 10, diag(5))[,,1])

> fit <- edgenet(X=X, Y=Y, G.X=aff.mat, 
                 lambda=1, psigx=1, family="gaussian")
> fit

#>Call: edgenet(X = X, Y = Y, G.X = aff.mat, lambda = 1, psigx = 1, \
#               family = gaussian)

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

You can install and use `netReg` either from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/netReg.html),
or by downloading the latest [tarball](https://github.com/dirmeier/netReg/releases).

If you want to use the **recommended way** using Bioconductor just call:

```r
> if (!requireNamespace("BiocManager", quietly=TRUE))
>   install.packages("BiocManager")
> BiocManager::install("netReg")
  
> library(netReg)
```
 
nstalling the R package using the downloaded tarball works like this:

```bash
$ R CMD install <netReg-x.y.z.tar.gz>
```

I **do not** recommend using `devtools`, but preferring tarballed releases over the most recent commit.

## Documentation

* For the R package: load the package using `library(netReg)`. We provide a vignette for the package that can be called using: `vignette("netReg")`. 
* For the command-line tool: You can also use the online [vignette](https://dirmeier.github.io/netReg/articles/netReg.html).

## Citation

It would be great if you cited `netReg` like this:

Simon Dirmeier, Christiane Fuchs, Nikola S Mueller, Fabian J Theis; 
netReg: network-regularized linear models for biological association studies, 
*Bioinformatics*, Volume 34, Issue 5, 1 March 2018, Pages 896–898, https://doi.org/10.1093/bioinformatics/btx677

## Author

Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
