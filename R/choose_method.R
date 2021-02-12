choose_method <- function(mat_pcs, counts_selected, non_out, ccluster, label){
  # browser()
  okay = 0
  # function
  if(is.function(ccluster) & okay == 0){
    sres = ccluster(mat_pcs[,non_out])
    if(is.data.frame(sres) | is.matrix(sres) | is.array(sres) | is.vector(sres)){
      temp = sres
      temp = as.matrix(temp)
      temp = as.character(temp)
      temp1 = unique(temp)
      if(length(temp) == length(non_out)){
        temp2 = rep(0,length(temp))
        for(cc in 1:sum(!is.na(temp1))){
          temp2[temp == temp1[cc]] = cc
        }
        sres = as.numeric(temp2)
        okay = 1
        ccluster = "PRC"
      }
    }
    if(okay == 0)
      stop("There is something wrong with 'ccluster' function!")
  }
  # label
  if(sum(ccluster == "labeled") & okay == 0){
    if(is.null(label))
      stop("There is no label inputed!")
    if(!is.null(dim(label)) & sum(dim(label) != 1) > 1)
      stop("There are more than 1 dimensions which have more than 1 elements in 'label'!")
    label = as.matrix(label)
    label = as.character(label)
    unilabel = unique(label)
    if(length(label) != ncol(counts_selected))
      stop("The length of 'label' must be the same as the number of cells!")
    label_ = rep(0,length(label))
    for(cc in 1:sum(!is.na(unilabel))){
      label_[label == unilabel[cc]] = cc
    }
    clust = as.numeric(label_)
    okay = 2
  }
  # specific clusters
  if(length(ccluster) == 2 & okay == 0){
    k1 = as.numeric(ccluster[1])
    k2 = ccluster[2]
    if(k2 == "Spectrum"){
      spec_res = Spectrum(mat_pcs[,non_out],method=3,maxk=10,fixk=k1)
      sres = spec_res$assignments
      okay = 1
    }
    if(k2 == "specc"){
      set.seed(10)
      spec_res = specc(t(mat_pcs[, non_out]), centers = k1, kernel = "rbfdot")
      sres = spec_res
      okay = 1
    }
  }

  if(okay == 0)
    stop("Wrong input of ccluster!")
  if(okay == 1){
    clust = rep(NA, ncol(counts_selected))
    clust[non_out] = sres
  }

  num_clusters = sum(!is.na(unique(clust)))
  Knumbers = rep(0, num_clusters)
  for(ii in 1:num_clusters){
    Knumbers[ii] = sum(clust == ii, na.rm = TRUE)
  }
  print("cluster sizes:")
  print(Knumbers)

  return(clust)
}


seurat_clustering <- function(count_s){
  if(is.null(rownames(count_s))){
    rownames(count_s) = 1:nrow(count_s)
  }
  Sdata <- CreateSeuratObject(count_s)

  Sdata <- NormalizeData(Sdata)

  Sdata <- FindVariableFeatures(Sdata, selection.method = "vst", nfeatures = 2000)

  Sdata <- ScaleData(Sdata, features = rownames(Sdata))

  Sdata <- RunPCA(Sdata, features = VariableFeatures(object = Sdata))

  #Sdata <- JackStraw(Sdata, num.replicate = 100)
  #ElbowPlot(Sdata)

  #cell--1:15
  Sdata <- FindNeighbors(Sdata, dims = 1:15)
  # Sdata <- FindNeighbors(Sdata)
  Sdata <- FindClusters(Sdata, resolution = 0.5)

  new_clust <- as.numeric(Idents(Sdata))

  return(new_clust)
}
