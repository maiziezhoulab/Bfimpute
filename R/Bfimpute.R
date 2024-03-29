#' Run Bfimpute
#'
#' Uses Bfimpute to recover drop-out events
#'
#' Bfimpute is a powerful imputation tool for scRNA-seq data that
#' recovers the dropout event by factorizing the count matrix into the product
#' of gene-specific and cell-specific feature matrices. Bfimpute uses full
#' Bayesian inference to describe the latent information for genes and cells
#' and carries out a Markov chain Monte Carlo scheme which is able to easily
#' incorporates side information to train the model and perform the imputation.
#'
#' @param counts Expression count matrix with rows corresponding to genes and
#' columns corresponding to cells.
#' @param ccluster The cluster approach: give \code{labeled} and corporate with
#' param \code{label} (see details in \code{label}) if the cells are labeled;
#' give a specific number and a spectral clustering approach chosen from
#' \code{Spectrum, specc} otherwise; and of cause you can use your own cluster
#' method and give us a function with a matrix as input and a vector as output
#' Default is \code{c(7,"Spectrum")}.
#' @param label Cell cluster labels which can be a vector, data.frame, matrix
#' with one row or one column, and etc (will only be used when \code{ccluster}
#' is set to \code{labeled}). Default is \code{NULL}.
#' @param normalized Whether the \code{counts} is normalized. \code{TRUE} for
#' yes and \code{FALSE} for not. If \code{FALSE}, Bfimpute will perform CPM and
#' log10 transform with bias 1.01. If \code{TRUE}, Bfimpute will not perform CPM
#' or log10 transform. But we highly recommend you to use log-transform after
#' normalization on your own. Default is \code{FALSE}.
#' @param S_G Gene latent matrix with \code{D} rows and the columns
#' corresponding to genes. Default is \code{NULL} which means no Gene latent
#' matrix.
#' @param S_C Cell latent matrix with \code{D} rows and the columns
#' corresponding to cells. Default is \code{NULL} which means no cell latent
#' matrix.
#' @param D Dimension of the latent factor. Default is \code{32}.
#' @param totalepoch Total number of epochs. Default is \code{300}.
#' @param burnin Number of burn-in epochs which are only be used as running
#' Markov chain. Default is \code{200}.
#' @param sn_max Maximum adaptive precision. Default is \code{10}.
#' @param sn_init Initial adaptive precision. Default is \code{1}.
#' @param threshold The threshold on dropout probabilities. Default is \code{0.5}.
#' @param ncores Number of cores used in parallel computation. Default is \code{5}.
#' @param out_dir The path and folder to store the results. Default is
#' \code{"./Bfimpute/"}.
#' @param out_name The file name which Bfimpute will save as. Default is
#' \code{"Bfimpute"}.
#' @param out_type The file type which Bfimpute will save as: "csv", "txt",
#' "rds", and "all" for all the three types, or "none" for just returning
#' without saving. Default is \code{"all"}.
#' @param returnGC Whether return the G and C matrices of the final epoch. If
#' \code{TRUE}, \code{Bfimpute} will return a list which consists of the imputed
#' matrix \code{R_calculate}, \code{G}, and \code{C}. If \code{FALSE},
#' \code{Bfimpute} will return the imputed matrix only. Default is \code{FALSE}.
#'
#' @return Bfimpute returns the imputed matrix with the same dimension as
#' \code{counts}. And it also saves the imputed count matrix to the specific
#' direction with specific name and types (depending on \code{out_dir},
#' \code{out_name}, and \code{out_type}).
#' @export
#'
#' @import stats
#' @import utils
#' @import parallel
#' @import Spectrum
#' @import rsvd
#' @import mnormt
#' @importFrom kernlab specc
#'
#' @examples
#' \donttest{
#' library(Bfimpute)
#' # use the following commands to generate simulated data
#' if(!requireNamespace("splatter", quietly = TRUE)) stop('The package splatter was not installed')
#' if(!requireNamespace("scater", quietly = TRUE)) stop('The package scater was not installed')
#' sce <- scater::mockSCE()
#' params <- splatter::splatEstimate(sce)
#' params <- splatter::setParams(params, update = list(nGenes = 20000,
#'                                           group.prob = rep(0.2,5),
#'                                           de.prob = 0.08,
#'                                           de.facLoc = 0.3,
#'                                           de.facScale = 0.5,
#'                                           batchCells = 500,
#'                                           dropout.type = "experiment"))
#' sim <- splatter::splatSimulate(params, method = "groups")
#' counts <- sim@assays@data@listData[["counts"]]
#' # impute via Bfimpute
#' counts_imputed <- Bfimpute(counts, ccluster = c(5,"specc"), label = NULL,
#'                            normalized = FALSE, S_G = NULL, S_C = NULL,
#'                            ncores = 5)
#' }
#'
#' @author Zi-Hang Wen \email{wenzihang0506@gmail.com}
#'
#' @references Wen et al.
#' Bfimpute: A Bayesian factorization method to recover single-cell RNA
#' sequencing data. \emph{bioRxiv}, 2021. doi:
#' \url{https://doi.org/10.1101/2021.02.10.430649}.
#'
#' See more usage in \url{https://github.com/maiziezhoulab/Bfimpute}
#'
#'
Bfimpute <- function(counts, ccluster = c(7,"Spectrum"), label = NULL,
                     normalized = FALSE, S_G = NULL, S_C = NULL, D = 32,
                     totalepoch = 300, burnin = 200, sn_max = 10, sn_init = 1,
                     threshold = 0.5, ncores = 5, out_dir = "./Bfimpute/",
                     out_name = "Bfimpute", out_type = "all", returnGC = F){
  counts = as.matrix(counts)
  print("Running Bfimpute")
  #----------------------check---------------------#
  if(!is.null(S_G))
    if(ncol(S_G)!=nrow(counts))
      stop("The column number of S_G is different from the row number of counts.")
  if(!is.null(S_C))
    if(ncol(S_C)!=ncol(counts))
      stop("The column number of S_C is different from the column number of counts.")

  #------------------Normalization-----------------#
  if(normalized == FALSE){
    sum_for_cell = apply(counts,2,sum)
    sum_for_cell[sum_for_cell == 0] = 1
    counts = t(1e6 * t(counts) / sum_for_cell)
    counts = log10(counts + 1.01)                 ####
    infimum = log10(1.01)                 ####
  }else{
    counts = counts + 0.01
    infimum = 0.01
  }

  if(returnGC){
    global.G = list()
    global.C = list()
  }

  #--------------------Imputation-------------------#
  method = 1
  # if(method == 1){
    slist = scluster(counts, ccluster, label, infimum, threshold)
    R_calculate_list = mclapply(1:length(slist), function(cc){
      print(paste0("imputing cluster ", cc))
      R_cc = counts[slist[[cc]]$selectgs, slist[[cc]]$cells]
      counts_cc = matrix(0, nrow = nrow(R_cc), ncol = ncol(R_cc))
      counts_cc[slist[[cc]]$true_data] = R_cc[slist[[cc]]$true_data]

      S_G_cc = S_G[,slist[[cc]]$selectgs]
      if(length(S_G_cc) == length(slist[[cc]]$selectgs)){
        S_G_cc = t(S_G_cc)
      }
      S_C_cc = S_C[,slist[[cc]]$cells]
      if(length(S_C_cc) == length(slist[[cc]]$cells)){
        S_C_cc = t(S_C_cc)
      }

      return(Gibbs_sampling(counts_cc, S_G_cc, S_C_cc, D, totalepoch, burnin, sn_max, sn_init, method, returnGC))
    }, mc.cores = ncores)
    R_calculate = counts
    for(cc in 1:length(slist)){
      if(returnGC){
        R_calculate[slist[[cc]]$selectgs, slist[[cc]]$cells] = R_calculate_list[[cc]][[1]]

        G_temp = R_calculate_list[[cc]][[2]]
        C_temp = R_calculate_list[[cc]][[3]]
        if(is.null(rownames(counts))){
          colnames(G_temp) = slist[[cc]]$selectgs
        }else{
          colnames(G_temp) = rownames(counts)[slist[[cc]]$selectgs]
        }
        if(is.null(colnames(counts))){
          colnames(C_temp) = slist[[cc]]$cells
        }else{
          colnames(C_temp) = colnames(counts)[slist[[cc]]$cells]
        }
        global.G[[cc]] = G_temp
        names(global.G)[cc] = paste0("cluster_",cc)
        global.C[[cc]] = C_temp
        names(global.C)[cc] = paste0("cluster_",cc)
      }else{
        R_calculate[slist[[cc]]$selectgs, slist[[cc]]$cells] = R_calculate_list[[cc]]
      }
    }
  # }
  # if(method == 2){
  #   counts[counts == infimum] = 0
  #   R_calculate = Gibbs_sampling(counts, S_G, S_C, D, totalepoch, burnin, sn_max, sn_init, method, returnGC)
  # }
  R_calculate[R_calculate<infimum] = infimum

  #-----------------Denormalization----------------#
  if(normalized == FALSE){
    R_calculate = 10^R_calculate - 1.01                 ####
    R_calculate = t(1e-6 * t(R_calculate) * sum_for_cell)
  }
  else{
    R_calculate = R_calculate - 0.01
  }

  if(returnGC){
    Bfimpute.return = list(R_calculate = R_calculate, G = global.G, C = global.C)
  }else{
    Bfimpute.return = R_calculate
  }

  if(out_type != "none"
     & out_type != "all"
     & out_type != "csv"
     & out_type != "txt"
     & out_type != "rds"){
    warning("\"out_type\" can only be \"csv\", \"txt\", \"rds\", \"all\", or \"none\"!")
    return(Bfimpute.return)
  }
  if(out_type == "none"){
    return(Bfimpute.return)
  }
  dir.create(out_dir, recursive = TRUE)
  if(out_type == "all"){
    write.csv(R_calculate, file = paste0(out_dir, "/", out_name, ".csv"))
    write.table(R_calculate, file = paste0(out_dir, "/", out_name, ".txt"))
    saveRDS(R_calculate, file = paste0(out_dir, "/", out_name, ".rds"))
  }else if(out_type == "csv"){
    write.csv(R_calculate, file = paste0(out_dir, "/", out_name, ".csv"))
  }else if(out_type == "txt"){
    write.table(R_calculate, file = paste0(out_dir, "/", out_name, ".txt"))
  }else if(out_type == "rds"){
    saveRDS(R_calculate, file = paste0(out_dir, "/", out_name, ".rds"))
  }

  print("Finish imputing")
  return(Bfimpute.return)
}

