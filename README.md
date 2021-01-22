# Bfimpute
Bfimpute is a powerful imputation tool for scRNA-seq data that
recovers the dropout event by factorizing the count matrix into the product
of gene-specific and cell-specific feature matrices.

## Installation
You can also install the most recent updates of Bfimpute from github with:
```R
# install.packages("devtools")
devtools::install_github("Vivianstats/scImpute")
```

## Quick start
library(Bfimpute)
```R
# use the following commands to generate simulated data
if(!requireNamespace("splatter", quietly = TRUE)) stop('The package splatter was not installed')
if(!requireNamespace("scater", quietly = TRUE)) stop('The package scater was not installed')
sce <- scater::mockSCE()
params <- splatter::splatEstimate(sce)
params <- splatter::setParams(params, update = list(nGenes = 20000,
                                          group.prob = rep(0.125,8),
                                          de.prob = 0.08,
                                          de.facLoc = 0.3,
                                          de.facScale = 0.5,
                                          batchCells = 800,
                                          dropout.type = "experiment"))
sim <- splatter::splatSimulate(params, method = "groups")
counts <- sim@assays@data@listData[["counts"]]
# or you can use your own data and make its name counts
counts_imputed <- Bfimpute(counts, ccluster = "Seurat", label = NULL,
                           normalized = FALSE, S_G = NULL, S_C = NULL,
                           ncores = 5)
```
