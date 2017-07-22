## Data

This file explains the available data sets.

* `X.tsv` is a `(n x p)` matrix of genetic markers where `n` is the number of samples (112) and `p` is the number of markers. The data has been taken from Brem <i>et al.</i>, Nature (2005).
* `Y.tsv` is a `(n x q)` matrix of expression values for `q` yeast genes. `n` is again the numer of samples (112).
* `GY.tsv` is a `(q x q)` adjacency matrix representing protein-protein interactions for the `q` response variables.
* `yeast.rda` is the three data sets compiled together for `R` usage.

## References

Data from `X.tsv` and `Y.tsv` has been compiled from:

* Rachel Brem *et al.*, Genetic interactions between polymorphisms that affect gene expression in yeast, Nature (2005).
* John D Storey *et al.*, Multiple locus linkage analysis of genomewide expression in yeast (2005), PLoS Biology.
* Wei Cheng *et al.*, Graph-regularized dual Lasso for robust eQTL mapping (2014), Bioinformatics.

`GY.tsv` is a *S. cerevisiae* protein-protein network taken and parsed from [BioGRID](https://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.150/BIOGRID-ORGANISM-3.4.150.tab.zip)
