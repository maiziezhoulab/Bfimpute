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
- `ccluster` The cluster approach: give \code{labeled} and corporate with
param \code{label} (see details in \code{label}) if the cells are
labeled, give the specific number \code{5} or \code{6} or ... if only the
cluster number is known, give \code{Seurat} and we will detect the clusters
on our own if lack of information, and of cause you can use your own cluster
method and give us a function with a matrix as input and a vector as output
(will only be used when \code{method} is set to \code{1}). Default is
\code{Seurat}.
- `S_G` Gene latent matrix with \code{D} rows and the columns
corresponding to genes. Default is \code{NULL} which means no Gene latent
matrix.
- `S_C` Cell latent matrix with \code{D} rows and the columns
corresponding to cells. Default is \code{NULL} which means no cell latent
matrix.

### Optional
- `method` Preprocessing choice: \code{1} for using Pre-QC, Pre-cluster
and dropout identification steps before imputation and \code{2} for just the
imputation. Default is \code{1}.
- `label` Cell cluster labels which can be a vector, data.frame, matrix
with one row or one column, and etc (will only be used when \code{ccluster}
is set to \code{labeled}). Default is \code{NULL}.
- `normalized` Whether the \code{counts} is raw or not. \code{FALSE} for
raw and \code{TRUE} for not. Default is {FALSE}.
- `D` Dimension of the latent factor. Default is \code{32}.
- `totalepoch` Total number of epochs. Default is \code{300}.
- `burnin` Number of burn-in epochs which are only be used as running
Markov chain. Default is \code{200}.
- `sn_max` Maximum adaptive precision. Default is \code{10}.
- `sn_init` Initial adaptive precision. Default is \code{1}.
- `threshold` The threshold on dropout probabilities. (will only be used
when \code{method} is set to \code{1}). Default is \code{0.5}.
- `ncores` Number of cores to use. Default is \code{5}.
- `out_dir` The path and folder to store the results. Default is
\code{"./Bfimpute/"}.
- `out_name` The file name which Bfimpute will save as. Default is
\code{"Bfimpute"}.
- `out_type` The file type which Bfimpute will save as: "csv", "txt",
"rds", and "all" for all the three types, or "none" for just returning
without saving. Default is \code{"all"}.

### Note
Larger `D`, `totalepoch`, and `burnin` will increase the accuracy of the
imputation, but it may take more efforts to run the process as a price.
On the contrary, smaller parameters will save your time.

