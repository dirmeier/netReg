library(data.table)
library(igraph)
library(magrittr)
library(netReg)

load("data/yeast.rda")
X <- yeast$X[,sample(1:ncol(X), 500)]
Y <- yeast$Y
GY <- yeast$GY

response.genes <- colnames(Y)
graph.genes    <- rownames(GY)
genes.intersec <- intersect(response.genes, graph.genes)

Y  <-  Y[, response.genes %in% genes.intersec]
GY <- GY[graph.genes %in% genes.intersec,
         graph.genes %in% genes.intersec]

drop.sparse.gene.names <- colnames(GY)[colSums(GY) > 10]
Y <- Y[, colnames(Y) %in% drop.sparse.gene.names]
GY <- GY[rownames(GY) %in% drop.sparse.gene.names,
         colnames(GY) %in% drop.sparse.gene.names]

GY.names <- rownames(GY)
Y <- Y[, order(colnames(Y))]
GY <- GY[order(rownames(GY)), order(colnames(GY))]


assertthat::assert_that(all(colnames(Y) == colnames(GY)),
                        all(colnames(Y) == rownames(GY)),
                        all(rownames(GY) == colnames(GY)))


yeast <- list(X=X, Y=Y, GY=GY)
devtools::use_data(yeast, overwrite=TRUE)

write.table(X, "data/X.tsv", row.names=F, col.names=F, quote=F, sep="\t")
write.table(Y, "data/Y.tsv", row.names=F, col.names=T, quote=F, sep="\t")
write.table(GY, "data/GY.tsv", row.names=F, col.names=T, quote=F, sep="\t")
