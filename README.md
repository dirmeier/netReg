# netReg <img src="https://rawgit.com/dirmeier/netReg/master/inst/sticker/sticker.png" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Project Life](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://travis-ci.org/dirmeier/netReg.svg?branch=master)](https://travis-ci.org/dirmeier/netReg)
[![codecov](https://codecov.io/gh/dirmeier/netReg/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/netReg)
[![bioc](https://bioconductor.org/shields/years-in-bioc/netReg.svg)](https://bioconductor.org/packages/release/bioc/html/netReg.html)

Generalized linear regression models with network-regularization in `R`. Now with `TensorFlow`.

## About

Modelling dependencies using linear regression models is often complicated when the 
analysed data-sets are high-dimensional and less observations than variables 
are available (n << p). `netReg` implements generalized linear models 
that utilize network penalties for regularization. Network regularization uses graphs
or trees to incorporate information about interactions of covariables, 
or responses, into the loss function of a GLM. Ideally this allows better (i.e., lower variance)
estimation of regression coefficients. 

For instance, in `R`, you could fit a network-regularized model like that:

```r
> library(netReg)

> X <- matrix(rnorm(100 * 10), 100)
> Y <- matrix(rnorm(100 * 10), 100)

> G.X <- abs(rWishart(1, 10, diag(10))[,,1])
> G.Y <- abs(rWishart(1, 10, diag(10))[,,1])

> fit <- edgenet(X, Y, G.X, G.Y)

> summary(fit)

#>call:
#>edgenet(X = X, Y = Y, G.X = G.X, G.Y = G.Y)

#>parameters:
#>lambda  psigx  psigy 
#>     1      1      1 

#>family:  gaussian
#>link:  identity 

#>-> call coef(x) for coefficients
```

From version `v1.9.0` on, we use `TensorFlow`, instead of custom `C++` and `Dlib`, for
estimation of regression coefficients replacing a custom *cyclic coordinate descent*. This allowed deleting of major parts of the code base.
`netReg` still uses some `RcppArmadillo` for fast matrix algebra.

In order to estimate the optimal hyperparameters, i.e., the regularization parameters
of the network models, we use Powell's BOBYQA algorithm in a standard cross-validation framework.

For more details, please check out the respective vignettes of the single models.

## Installation

Before installing the package, make sure to have these Python dependencies installed:

* `tensorflow==1.14.0`,
* `tensorflow-probabiltiy==0.7.0`

The easiest way to install `tfprobability` and all dependencies from `R` is using:

```r
> install.packages(c("tensorflow", "tfprobability"))
> tfprobability::install_tfprobability(version = "0.7.0", tensorflow = "1.14.0",
                                       method="conda", envname="r-tensorflow")
```

That creates a `conda` environment (in case you use it) called `r-tensorflow` and installs the Python dependencies automatically.

You can install `netReg` either from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/netReg.html),
or using the latest GitHub [tag](https://github.com/dirmeier/netReg/tags).

If you want to use the recommended way using Bioconductor just call:

```r
> if (!requireNamespace("BiocManager", quietly=TRUE))
>   install.packages("BiocManager")
> BiocManager::install("netReg")
  
> library(netReg)
```
 
Installing the R package using the latest GitHub tag/release works like this:

```r
> devtools::install_github("dirmeier/netReg@v1.11.3")
```

The GitHub versions are usually ahead of Bioconductor.

## Documentation

* Load the package using `library(netReg)`. We provide vignettes for the package that can be called using: `vignette(package="netReg")`. 
* You can also use the online [vignette](https://dirmeier.github.io/netReg).

## Citation

If `netReg` was useful for you or your work, it would be great if you cited it like this:

```
@article{,
  title={netReg: network-regularized linear models for biological association studies},
  author={Dirmeier, Simon and Fuchs, Christiane and Mueller, Nikola S and Theis, Fabian J},
  journal={Bioinformatics},
  volume={34},
  number={5},
  pages={896--898},
  year={2017},
  publisher={Oxford University Press}
}

```

## Author

Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier @ gmx.de</a>
