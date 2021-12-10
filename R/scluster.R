# This part of code is partly adopted from scImpute in doing the pca before
# clustering and in implementing Gamma-Normal mixture distribution model to
# detect real dropouts.

clustering <- function(counts, ccluster, label, infimum){
  I = nrow(counts)
  J = ncol(counts)

  selist = lapply(1:I, function(i) setdiff(counts[i, ], infimum))
  gmean = sapply(1:I,function(i){
    mean(selist[[i]])
  })
  gmean[is.na(gmean)] = 0
  gstd = sapply(1:I,function(i){
    sd(selist[[i]])
  })
  gstd[is.na(gstd)] = 0
  gcv = gstd/gmean
  gcv[is.na(gcv)] = 0
  selected = which(gmean >= 1 & gcv >= quantile(gcv, 0.25))
  if(length(selected) < 500){
    selected = 1:I}
  counts_selected = counts[selected, ]

  if(J < 5000){
    var_thre = 0.4
    pca = prcomp(t(counts_selected))
  }else{
    var_thre = 0.6
    pca = rpca(t(counts_selected), k = 1000, center = TRUE, scale = FALSE)
  }

  eigs = (pca$sdev)^2
  var_cum = cumsum(eigs)/sum(eigs)
  if(max(var_cum) <= var_thre){
    npc = length(var_cum)
  }else{
    npc = which.max(var_cum > var_thre)
    npc = max(npc, 15)                        ############
  }

  if(npc < 3){
    npc = 3
  }

  mat_pcs = t(pca$x[, 1:npc])

  dist_cells = matrix(0, nrow = J, ncol = J)
  for(id1 in 1:J){
    dist_cells[id1,] = c(sapply(1:id1, function(id2){
      sqrt(sum((mat_pcs[, id1] - mat_pcs[, id2])^2))
    }),rep(0,J-id1))
  }
  dist_cells = dist_cells + t(dist_cells)

  min_dist = sapply(1:J, function(i){
    min(dist_cells[i, -i])
  })
  iqr = quantile(min_dist, 0.75) - quantile(min_dist, 0.25)
  outliers = which(min_dist > 2.5 * quantile(min_dist, 0.75) -
                     1.5 * quantile(min_dist, 0.25))

  non_out = setdiff(1:J, outliers)

  clust = choose_method(mat_pcs, counts_selected, non_out, ccluster, label)

  return(clust)
}


update_pars <- function(xdata, wt){
  tp_s = sum(wt)
  tp_t = sum(wt * xdata)
  tp_u = sum(wt * log(xdata))
  tp_v = -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha = 20
  }else{
    alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){
      alpha = 20
    }else{
      alpha = uniroot(function(alpha, target){
        log(alpha) - digamma(alpha) - target
      },
      c(0.9, 1.1) * alpha0,
      target = tp_v,
      extendInt = "yes")$root
    }
  }
  beta = tp_s / tp_t * alpha
  return(c(alpha, beta))
}


mixing <- function(xdata, infimum){
  inits = rep(0, 5)
  inits[1] = sum(xdata == infimum)/length(xdata)
  if (inits[1] == 0) {inits[1] = 0.01}
  inits[2:3] = c(0.5, 1)
  xdata_rm = xdata[xdata > infimum]
  inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {inits[5] = 0}
  paramt = inits
  eps = 10
  iter = 0
  loglik_old = 0

  while(eps > 0.5) {
    pz1 = paramt[1] * dgamma(xdata, shape = paramt[2], rate = paramt[3])
    pz2 = (1 - paramt[1]) * dnorm(xdata, mean = paramt[4], sd = paramt[5])
    pz = pz1/(pz1 + pz2)
    pz[pz1 == 0] = 0
    wt = cbind(pz, 1 - pz)

    paramt[1] = sum(wt[, 1])/nrow(wt)
    paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])
    paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 2]))
    paramt[2:3] = update_pars(xdata, wt[,1])

    loglik = sum(log10(paramt[1] * dgamma(xdata, shape = paramt[2], rate = paramt[3])
                       + (1 - paramt[1]) * dnorm(xdata, mean = paramt[4], sd = paramt[5])))
    eps = (loglik - loglik_old)^2
    loglik_old = loglik
    iter = iter + 1
    if (iter > 100)
      break
  }
  return(paramt)
}


parameters <- function(counts, infimum){
    ngs = which(abs(rowSums(counts) - infimum * ncol(counts)) < 1e-10)
    parslist = lapply(1:nrow(counts), function(ii) {
      if (ii %in% ngs) {
        return(rep(NA, 5))
      }
      xdata = counts[ii, ]
      paramt = try(mixing(xdata, infimum), silent = TRUE)
      if (class(paramt) == "try-error"){
        paramt = rep(NA, 5)
      }
      return(paramt)
    })
    parslist = Reduce(rbind, parslist)
    colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
    return(parslist)
  }


select_genes <- function(parslist, subcount, infimum){
  selectgs = which( (rowSums(subcount) > infimum * ncol(subcount)) &
                         complete.cases(parslist) )
  if(length(selectgs) == 0) return(selectgs)
  # find out genes that violate assumption
  mu = parslist[, "mu"]
  sgene1 = which(mu <= ifelse(infimum == log10(1.01), log10(2.01), 1.01))

  dcheck1 = dgamma(mu+1, shape = parslist[, "alpha"], rate = parslist[, "beta"])
  dcheck2 = dnorm(mu+1, mean = parslist[, "mu"], sd = parslist[, "sigma"])
  sgene3 = which(dcheck1 >= dcheck2 & mu <= 1)

  sgene = union(sgene1, sgene3)
  selectgs = setdiff(selectgs, sgene)
  return(selectgs)
}


scluster <- function(counts, ccluster, label, infimum, threshold){
  clust = clustering(counts, ccluster, label, infimum)

  return(lapply(1:length(unique(clust[!is.na(clust)])), function(cc){
    cells = which(clust == cc)
    parslist = parameters(counts[, cells], infimum)
    selectgs = select_genes(parslist, counts[, cells], infimum)

    subc = counts[selectgs, cells, drop = FALSE]
    parslist = parslist[selectgs, , drop = FALSE]

    droprate = t(sapply(1:nrow(subc), function(i) {
      xdata = subc[i, ]
      paramt = parslist[i, ]
      pz1 = paramt[1] * dgamma(xdata, shape = paramt[2], rate = paramt[3])
      pz2 = (1 - paramt[1]) * dnorm(xdata, mean = paramt[4], sd = paramt[5])
      pz = pz1/(pz1 + pz2)
      pz[pz1 == 0] = 0
      wt = cbind(pz, 1 - pz)
      return(wt[, 1])
    }))

    meancheck = t(sapply(1:nrow(subc),function(nn){
      subc[nn,]>parslist[nn, "mu"]
    }))
    droprate[meancheck & droprate > threshold] = 0

    true_data = droprate <= threshold

    return(list(cells = cells, selectgs = selectgs, true_data = true_data))
  }))

}
