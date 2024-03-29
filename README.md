# Bfimpute: A Bayesian factorization method to recover single- cell RNA sequencing data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5676122.svg)](https://doi.org/10.5281/zenodo.5676122)
[![repo size](https://img.shields.io/github/repo-size/maiziezhoulab/Bfimpute.svg)](https://github.com/maiziezhoulab/Bfimpute/archive/master.zip)

Bfimpute is a powerful imputation tool for scRNA-seq data that
recovers the dropout event by factorizing the count matrix into the product
of gene-specific and cell-specific feature matrices.

## Installation
You can also install the most recent updates of Bfimpute from github with:
```R
# install.packages("devtools")
devtools::install_github("maiziezhoulab/Bfimpute")
```

## Examples
In this section we will run the dataset from [Chu. et al (2016)](https://link.springer.com/article/10.1186/s13059-016-1033-x)
and the version with QC can be downloaded from
[here](https://drive.google.com/drive/folders/1C2rjTDy3Lvi4DE988FvGSOOCODVUyDI-?usp=sharing).

### with spectral clustering
```R
library(Bfimpute)
```
Set the folder direction of your data below:
```R
data_dir = "./"
```
For dataset we present, you can load the cell types, bulk and scRNAseq counts
matrices as followed:
```R
counts = read.csv(paste0(data_dir, "sc_qc.csv"), row.names = 1, header = T)
# counts[1:5,1:5]
bulk = read.csv(paste0(data_dir, "bulk_qc.csv"), row.names = 1, header = T)
# bulk[1:5,1:5]
cell_type = read.csv(paste0(data_dir, "cell_type.csv"), row.names = 1, header = T)
# unique(cell_type)
```
As we already know that there are `7` cell types in this dataset, we can run
`Bfimpute` using `Spectrum` or `specc` as shown:
```R
counts_imputed <- Bfimpute(counts, ccluster = c(7,"Spectrum"), ncores = 5)
# or
# counts_imputed <- Bfimpute(counts, ccluster = c(7,"specc"), ncores = 5)
```

### with other clustering methods
If other clustering methods implemented in the future is more powerful, we can
easily use them and replace our clustering step by building a function. The
input of function `ccluster` should be a matrix with points to cluster as
columns and rows as features, while the output is a vector of cluster result.

We present `kmeans` as an example here:
```R
N = 7 # number of cell types
ccluster = function(x){kmeans(t(x), centers = N)$cluster}
counts_imputed <- Bfimpute(counts = counts, ccluster = ccluster, ncores = 5)
```

### with cell labels instead of clustering
If the labeled cell type is given as it is now, we can use them to replace
clustering step and gain more accuracy:
```R
counts_imputed <- Bfimpute(counts, ccluster = "labeled", label = cell_type, ncores = 5)
```


### Gene or cell related information
If some gene side information is available for the raw counts, make it a
positive matrix and having columns in the same gene orders as in the row of the
count matrix. And we can use them to assist the imputation (`ccluster` here can
be changed to any options mentioned above):
```R
S_G = t(bulk)
counts_imputed <- Bfimpute(counts, ccluster = c(7,"Spectrum"), S_G = S_G, ncores = 5)
```
Similarly, if some cell side information is also available, make it a positive
matrix and having columns in the same cell orders as in the column of the count
matrix. And use the parameter `S_C` in the same way.

## Parameters
### Required
- `counts` Expression count matrix with rows corresponding to genes and
columns corresponding to cells.
- `ccluster` The cluster approach: give `labeled` and corporate with
param `label` (see details in `label`) if the cells are labeled;
give a specific number and a spectral clustering approach chosen from
`Spectrum, specc` otherwise; and of cause you can use your own cluster
method and give us a function with a matrix as input and a vector as output
Default is `c(7,"Spectrum")`.

### Optional
- `S_G` Gene latent matrix with `D` rows and the columns
corresponding to genes. Default is `NULL` which means no Gene latent
matrix.
- `S_C` Cell latent matrix with `D` rows and the columns
corresponding to cells. Default is `NULL` which means no cell latent
matrix.
- `label` Cell cluster labels which can be a vector, data.frame, matrix
with one row or one column, and etc (will only be used when `ccluster`
is set to `labeled`). Default is `NULL`.
- `normalized` Whether the `counts` is normalized. `TRUE` for
yes and `FALSE` for not. If `FALSE`, Bfimpute will perform CPM and
log10 transform with bias 1.01. If `TRUE`, Bfimpute will not perform CPM
or log10 transform. But we highly recommend you to use log-transform after
normalization on your own. Default is `FALSE`.
- `D` Dimension of the latent factor. Default is `32`.
- `totalepoch` Total number of epochs. Default is `300`.
- `burnin` Number of burn-in epochs which are only be used as running
Markov chain. Default is `200`.
- `sn_max` Maximum adaptive precision. Default is `10`.
- `sn_init` Initial adaptive precision. Default is `1`.
- `threshold` The threshold on dropout probabilities. Default is `0.5`.
- `ncores` Number of cores used in parallel computation. Default is `5`.
- `out_dir` The path and folder to store the results. Default is
`"./Bfimpute/"`.
- `out_name` The file name which Bfimpute will save as. Default is
`"Bfimpute"`.
- `out_type` The file type which Bfimpute will save as: "csv", "txt",
"rds", and "all" for all the three types, or "none" for just returning
without saving. Default is `"all"`.
- `returnGC` Whether return the G and C matrices of the final epoch. If
`TRUE`, `Bfimpute` will return a list which consists of the imputed
matrix `R_calculate`, `G`, and `C`. If `FALSE`, `Bfimpute` will return the
imputed matrix only. Default is `FALSE`.

### Note
Larger `D`, `totalepoch`, and `burnin` will increase the accuracy of the
imputation, but it may take more efforts to run the process as a price.
On the contrary, smaller parameters will save your time.

The tool is implemented with hyperparameters set to: `\mu_0=0`, `\beta_0=2`,
`\nu_0=D`, and `W_0=I` (the identity matrix).

## Reference
A Bayesian factorization method to recover single-cell RNA sequencing data. Z.-H. Wen, J. L. Langsam, L. Zhang, W. Shen, X. Zhou. Cell Reports Methods (2022) 2, 100133. https://doi.org/10.1016/j.crmeth.2021.100133

