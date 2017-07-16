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

GY.names <- rownames(GY)
Y <- Y[, order(colnames(Y))]
Y[1:5, 1:5]
GY <- GY[order(rownames(GY)), order(colnames(GY))]

assertthat::assert_that(all(colnames(Y) == colnames(GY)),
                        all(colnames(Y) == rownames(GY)),
                        all(rownames(GY) == colnames(GY)))

s <- edgenet(X=X, Y=Y, G.Y=GY, lambda=10, psigx=0, psigy=1, maxit=10000, thresh=10e-5)

