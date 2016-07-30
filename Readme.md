<h1 align="center"> netReg </h1>

#TODO code cov

An R-library for network-regularized generalized linear regression. 

## Introduction

`netReg` is a generalization of the LASSO and as such suited to estimate high-dimensional linear dependencies of two sets of variables. We use graph prior knowledge, as for example given in co-expression or gene-regulatory networks in order to specifically shrink coefficient profiles as suggested by a graph-defined similarity measure. The main routines for estimation of coefficients and shrinkage parameters are implemented in C++11. 

Depending on your installed libraries `netReg` uses `Intel MKL`, `Armadillo`, `OpenBLAS`, `BLAS`, and `Lapack` for fast computation of matrix operations. 

We use `Dlib` in order to calculate the most optimal set of shrinkage parameters using k-fold cross-validation.

## Installation
 
Install `netReg` using:
```{r}
library(devtools)
install_github("rafstraumur/netreg") 
```
from the R-console.

## Usage

Load the package using `library(netReg)`. We provide a vignette for the package that can be called using: `vignette("netReg")`.
Basically that is all you have to know.

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
