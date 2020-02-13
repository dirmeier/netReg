# Installation

## Dependencies

Before installing the package, make sure to have these Python dependencies installed:

* `tensorflow==1.14.0`,
* `tensorflow-probabiltiy==0.5.0`

The easiest way is probably to install `TensorFlow` from with `R` and then call:

```{r}
tensorflow::install_tensorflow(extra_packages = "tensorflow-probability")
```

That creates a `conda` environment (in case you use it) called `r-tensorflow` and installs the Python dependencies automatically.

If this does not work for you, try this approach on the command line:

```{r}
conda create -n r-tensorflow python=3.6
source activate r-tensorflow
conda install tensorflow==1.14.0 tensorflow-probability==0.5.0
```

## netReg

You can install and use `netReg` either from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/netReg.html),
or by downloading the latest [tarball](https://github.com/dirmeier/netReg/releases).

If you want to use the **recommended way** using Bioconductor just call:

```{r}
> if (!requireNamespace("BiocManager", quietly=TRUE))
>   install.packages("BiocManager")
> BiocManager::install("netReg")

> library(netReg)
```

Installing the R package using the downloaded tarball works like this:

```{bash}
$ R CMD install <netReg-x.y.z.tar.gz>
```
