# Bfimpute
Bfimpute is a powerful imputation tool for scRNA-seq data that
recovers the dropout event by factorizing the count matrix into the product
of gene-specific and cell-specific feature matrices.

## Installation
You can also install the most recent updates of Bfimpute from github with:
```R
# install.packages("devtools")
devtools::install_github("maiziezhoulab/Bfimpute")
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

## Parameters
### Required
- `counts` Expression count matrix with rows corresponding to genes and
columns corresponding to cells.
- `ccluster` The cluster approach: give `labeled` and corporate with
param `label` (see details in `label`) if the cells are
labeled, give the specific number `5` or `6` or ... if only the
cluster number is known, give `Seurat` and we will detect the clusters
on our own if lack of information, and of cause you can use your own cluster
method and give us a function with a matrix as input and a vector as output
(will only be used when `method` is set to `1`). Default is
`Seurat`.
- `S_G` Gene latent matrix with `D` rows and the columns
corresponding to genes. Default is `NULL` which means no Gene latent
matrix.
- `S_C` Cell latent matrix with `D` rows and the columns
corresponding to cells. Default is `NULL` which means no cell latent
matrix.

### Optional
- `method` Preprocessing choice: `1` for using Pre-QC, Pre-cluster
and dropout identification steps before imputation and `2` for just the
imputation. Default is `1`.
- `label` Cell cluster labels which can be a vector, data.frame, matrix
with one row or one column, and etc (will only be used when `ccluster`
is set to `labeled`). Default is `NULL`.
- `normalized` Whether the `counts` is raw or not. `FALSE` for
raw and `TRUE` for not. Default is `FALSE`.
- `D` Dimension of the latent factor. Default is `32`.
- `totalepoch` Total number of epochs. Default is `300`.
- `burnin` Number of burn-in epochs which are only be used as running
Markov chain. Default is `200`.
- `sn_max` Maximum adaptive precision. Default is `10`.
- `sn_init` Initial adaptive precision. Default is `1`.
- `threshold` The threshold on dropout probabilities. (will only be used
when `method` is set to `1`). Default is `0.5`.
- `ncores` Number of cores to use. Default is `5`.
- `out_dir` The path and folder to store the results. Default is
`"./Bfimpute/"`.
- `out_name` The file name which Bfimpute will save as. Default is
`"Bfimpute"`.
- `out_type` The file type which Bfimpute will save as: "csv", "txt",
"rds", and "all" for all the three types, or "none" for just returning
without saving. Default is `"all"`.

### Note
Larger `D`, `totalepoch`, and `burnin` will increase the accuracy of the
imputation, but it may take more efforts to run the process as a price.
On the contrary, smaller parameters will save your time.

