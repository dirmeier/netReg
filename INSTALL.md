# Installation

## Dependencies

Before installing the package, make sure to have these Python dependencies installed:

* `tensorflow>=2.2.0`,
* `tensorflow-probabiltiy>=0.10.0`

The easiest way is probably to install `TensorFlow` from with `R` and then call:

```{r}
tensorflow::install_tensorflow(extra_packages = "tensorflow-probability")
```

That creates a `conda` environment (in case you use it) called `r-tensorflow` and installs the Python dependencies automatically.

If this does not work for you, try this approach on the command line:

```{r}
conda create -n r-tensorflow python=3.7
source activate r-tensorflow
conda install tensorflow==2.2.0 tensorflow-probability==0.10.0
```

## netReg

You can install and use `netReg` by downloading the latest [tarball](https://github.com/dirmeier/netReg/releases) and then calling:

```{bash}
$ R CMD install <netReg-x.y.z.tar.gz>
```
