# netReg <img src="https://rawgit.com/dirmeier/netReg/master/inst/sticker/sticker.png" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Project Life](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://travis-ci.org/dirmeier/netReg.svg?branch=master)](https://travis-ci.org/dirmeier/netReg)
[![codecov](https://codecov.io/gh/dirmeier/netReg/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/netReg)

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

* `tensorflow>=2.2.0`,
* `tensorflow-probabiltiy>=0.10.0`

The easiest way is probably to install `TensorFlow` from with `R`

```{r}
install.packages(c("tensorflow", "tfprobability"))
tensorflow::install_tensorflow(extra_packages = "tensorflow-probability")
```

That should do it.

If this does not work for you, try this approach on the command line:

```{bash}
conda create -n r-tensorflow python=3.7
source activate r-tensorflow
conda install tensorflow==2.2.0 tensorflow-probability==0.10.0
```

You can then install and use `netReg` by downloading the latest [release](https://github.com/dirmeier/netReg/releases)

```{r}
remotes::install_github("dirmeier/netReg@v1.12.0")
```

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

Simon Dirmeier <a href="mailto:simon.dirmeier @ web.de">simon.dirmeier @ web.de</a>
