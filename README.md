
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
It incorporates prior knowledge in the form of graphs into the model's likelihood and by that allows better estimation of regression coefficients.
The main routines for estimation of coefficients and shrinkage parameters are implemented in `C++11`. 

Depending on your installed libraries `netReg` uses `OpenBLAS` or `BLAS`, and `Lapack` for fast computation of matrix operations in an `Armadillo\RcppArmadillo` framework. We use `Dlib` in order to calculate the most optimal set of shrinkage parameters using k-fold cross-validation.

`netReg` comes as a stand alone `C++` command line tool shipped with [`bioconda`](https://anaconda.org/bioconda/netreg) as well as a [`Bioconductor`](https://bioconductor.org/packages/release/bioc/html/netReg.html) package.

## Installation
 
You can install and use `netReg` either

* as an `R` library from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/netReg.html),
* or as a `C++` command line tool from [bioconda](https://anaconda.org/bioconda/netreg),
* or by downloading the [tarball](https://github.com/dirmeier/netReg/releases) and doing either of the previous options manually.

### Installation for R with Bioconductor

If you want to use the `R` version of `netReg` call:

```r
> source("https://bioconductor.org/biocLite.R")
> biocLite("netReg")
  
> library(netReg)
```
 
from the `R`-console. 

### Installation for command line with bioconda

If you want to use the `C++` command line tool you can do this using conda. 
For that you should download [Anaconda](https://www.continuum.io/downloads) and create a [virtual environment](https://conda.io/docs/using/envs.html).
Then install the tool using:

```sh
$ conda install -c bioconda netreg
  
$ netReg -h
```

### Manual installation

If you don't like package managers you can download the tarball of the latest [release](https://github.com/dirmeier/netReg/releases/tag/v1.0.0) and install both or either from the two.

#### Command line tool

The command line tool has the **following dependencies** which need to be installed:

* `CMake >= 3.6`,
* `Boost >= 1.6.x`,
* `Armadillo >= 7.800.3`,
* `OpenBLAS/BLAS` and `Lapack` (older versions should work),
* *optional*: `OpenMP` (older versions should work).

To install the command line tool manually:

```sh
$ mkdir build && cd build
$ cmake .. && make
  
$ ./netReg -h
```

If you want the tool to be installed in some specific folder you would also call:

```sh
$ make install --prefix=/some/path
```

#### R package

Installing the `R` using the downloaded tarball works like this:

```bash
$ R CMD install <netReg-x.y.z.tar.gz>
```

### Installation on Mac

In some cases it is required to install `gfortan` for Mac first (which is needed by `Armadillo/RcppArmadillo`). In that case run:

```sh
$ curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
$ sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Afterwards just install the package as described above.

## Documentation

There are two tutorials for either `R` or the `C++` command line tool available.
We are always glad to take questions, so feel free to write or open up an issue.

### R package

Load the package using `library(netReg)`. 
We provide a vignette for the package that can be called using: `vignette("netReg")`. You can also use the online [tutorial](https://dirmeier.github.io/netReg/articles/netReg_R.html).

### Command line tool

Have a look at the command line [tutorial](https://dirmeier.github.io/netReg/articles/netReg_commandline.html).

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
