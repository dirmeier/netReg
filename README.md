# netReg <img src="https://cdn.rawgit.com/dirmeier/netReg/7b8e31e0/_fig/sticker.svg" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/dirmeier/netReg.svg?branch=master)](https://travis-ci.org/dirmeier/netReg)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dirmeier/netReg?branch=master&svg=true)](https://ci.appveyor.com/project/dirmeier/netReg)
[![codecov](https://codecov.io/gh/dirmeier/netReg/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/netReg)
[![bioc](https://bioconductor.org/shields/years-in-bioc/netReg.svg)](https://bioconductor.org/packages/release/bioc/html/netReg.html)
[![conda](https://anaconda.org/bioconda/netreg/badges/installer/conda.svg)](https://anaconda.org/bioconda/netreg)

Network-penalized generalized linear models in R and C++.

## Introduction

`netReg` is an R/C++ implementation of a network-regularized linear regression model.
It incorporates prior knowledge in the form of graphs into the model's likelihood and by that allows better estimation of regression coefficients.
The main routines for estimation of coefficients and shrinkage parameters are implemented in `C++11`. 
Depending on your installed libraries `netReg` uses `OpenBLAS` or `BLAS`, and `Lapack` for fast computation of matrix operations in an `Armadillo\RcppArmadillo` framework. We use `Dlib` in order to calculate the most optimal set of shrinkage parameters using k-fold cross-validation.

`netReg` comes as a standalone `C++` tool shipped with [`bioconda`](https://anaconda.org/bioconda/netreg) as well as as [`Bioconductor`](https://bioconductor.org/packages/release/bioc/html/netReg.html) package. Check out the [landing page](https://dirmeier.github.io/netReg) for more information.

## Installation
 
1) If you want to use the `R` version of `netReg` call this:

```{R}
  source("https://bioconductor.org/biocLite.R")
  biocLite("netReg")
```
 
from the `R`-console. 

2) If you want to use the `C++` commandline tool you can do this using `conda`. For that you should download [Anaconda](https://www.continuum.io/downloads) and create a [virtual environment](https://conda.io/docs/using/envs.html).
Then install the tool using:

```{bash}
  conda install -c bioconda netreg
```

3) Alternatively you can download the `tarball` from the latest [release](https://github.com/dirmeier/netReg/releases/tag/v1.0.0) and install both.

**For the commandline tool you need *recent* `CMake`, `Armadillo`, `Boost`, `OpenBLAS\BLAS` and `Lapack` versions installed**. Optionally `OpenMP` is also supported.

```{bash}
  mkdir build && cd build
  cmake .. && make
```

If you want the tool to be installed at some place you would add:

```{bash}
  make install --prefix=/some/path
```

Installing the `R` using the downloaded `tarball` works like this:

```{bash}
  R CMD install <netReg-x.y.z.tar.gz>
```

### Installation on Mac

In some cases it is required to install `gfortan` for Mac first (which is needed by `Armadillo/RcppArmadillo`). I that case run:

```sh
  curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
  sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Afterwards just install the package as described above.

## Documentation

### R

Load the package using `library(netReg)`. We provide a vignette for the package that can be called using: `vignette("netReg")`. Basically that is all you have to know.

### C++

See the [tutorial](https://dirmeier.github.io/netReg/tutorial).


## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
