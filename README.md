# netReg <img src="https://cdn.rawgit.com/dirmeier/netReg/_fig/netReg_hexsticker.svg" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/dirmeier/netReg.svg?branch=master)](https://travis-ci.org/dirmeier/netReg)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dirmeier/netReg?branch=master&svg=true)](https://ci.appveyor.com/project/dirmeier/netReg)
[![codecov](https://codecov.io/gh/dirmeier/netReg/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/netReg)

Network-penalized generalized linear models in R.

## Introduction

`netReg` is an R/C++ implementation of a network-regularized linear regression model.
It incorporates prior knowledge in the form of graphs into the model's likelihood and by that allows better estimation of regression coefficients.
The main routines for estimation of coefficients and shrinkage parameters are implemented in C++11. 
Depending on your installed libraries `netReg` uses `OpenBLAS` or `BLAS`, and `Lapack` for fast computation of matrix operations in an `RcppArmadillo` framework. We use `Dlib` in order to calculate the most optimal set of shrinkage parameters using k-fold cross-validation.

## Installation
 
Install `netReg` using:

```{R}
  source("https://bioconductor.org/biocLite.R")
  biocLite("netReg")
 ```
from the R-console.

Alternatively you can download the tarball from the latest [release](https://github.com/dirmeier/netReg/releases/tag/v1.0.0) 
and install it on the command line.

### Installation on Mac

In some cases it is required to install `gfortan` for Mac first (which is needed by `RcppArmadillo`). I that case run:

```sh
  curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
  sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Afterwards just install the package as described above.

## Usage

Load the package using `library(netReg)`. We provide a vignette for the package that can be called using: `vignette("netReg")`.
Basically that is all you have to know.

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
