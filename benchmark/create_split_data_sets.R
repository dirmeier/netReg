#!/usr/bin/env Rscript

ar <- commandArgs(TRUE)

fl.name <- ar[1]
X <- as.matrix(read.csv(fl.name, sep="\t", header=FALSE))

folds <- sample(rep(1:10, length.out=nrow(X)))

for (i in 1:10)
{
    test <- which(folds == i)
    train <- which(folds != i)
    X.train  <- X[train,]
    X.test  <- X[test,]
    out.train <- paste0(sub(".tsv", "", fl.name), "_train_", i, ".tsv")
    out.test <- paste0(sub(".tsv", "", fl.name), "_test_", i, ".tsv")

    write.table(X.train, out.train, quote=F, row.names=F, col.names=F)
    write.table(X.test, out.test, quote=F, row.names=F, col.names=F)
}
