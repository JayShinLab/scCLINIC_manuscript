library(Seurat)
##Following Seurat cell hashing tutorial https://satijalab.org/seurat/articles/hashing_vignette.html

# Read in UMI count matrix for RNA
hto12.umis <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/hto12_umi_mtx.rds")

# Read in HTO count matrix
hto12.htos <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/hto12_hto_mtx.rds")

# Select cell barcodes detected in both RNA and HTO
cells.use <- intersect(rownames(hto12.htos), colnames(hto12.umis))

# Create Seurat object and add HTO data
hto12 <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(hto12.umis[, cells.use]), sparse = T),
                            min.features = 300)
hto12[["HTO"]] <- CreateAssayObject(counts = t(x = hto12.htos[colnames(hto12), 1:12]))

# Normalize data
hto12 <- NormalizeData(hto12)
hto12 <- NormalizeData(hto12, assay = "HTO", normalization.method = "CLR")

hto12 <- HTODemux(hto12, assay = "HTO", positive.quantile = 0.99)
saveRDS(hto12, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/hto12.rds"))

#####Run different doublet detection - adapt from benchmarking code
#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }

#############################################################################################
##########################################################################
# Predict doublets based on doublet scores
##########################################################################
CallDoublets <- function(score, rate){
  
  #cat("Predict doublets...\n", file = stderr())
  num <- floor(length(score) * rate)
  threshold <- sort(score, decreasing = T)[num]
  pred <- score > threshold
  return(which(pred))
}

##########################################################################
# Call DoubletFinder to obtain doublet scores
##########################################################################
CallDoubletFinder_Update <- function(count, nfeatures=2000, PCs=10){
  
  rownames(count) <- as.character(1:dim(count)[1])
  colnames(count) <- as.character(1:dim(count)[2])
  seurat <- Seurat::CreateSeuratObject(count)
  seurat <- Seurat::NormalizeData(seurat, verbose = F)
  seurat <- Seurat::ScaleData(seurat, verbose = F)
  seurat <- Seurat::FindVariableFeatures(seurat, selection.method = "vst", nfeatures = nfeatures, verbose = F)
  seurat <- Seurat::RunPCA(seurat, verbose = F)
  
  sink('NUL')
  sink(stdout(), type = "message")
  
  tryCatch({
    sweep.vector <- DoubletFinder::paramSweep_v3(seurat, PCs = 1:PCs, sct = FALSE)
    sweep.table <- DoubletFinder::summarizeSweep(sweep.vector, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.table)
    
    pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
    pK <- as.numeric(levels(pK))[pK]
    seurat <- DoubletFinder::doubletFinder_v3(seurat, PCs = 1:PCs, pN = 0.25, pK = pK,
                                              nExp = 0.1, reuse.pANN = FALSE, sct = FALSE)
  },
  interrupt = function(e){
    sink(NULL, type="message")
    sink()
  }
  )
  sink(NULL, type="message")
  sink()
  metadata_names <- colnames(seurat@meta.data)
  matching_column <- grep("pANN", metadata_names, value = TRUE)
  score <- seurat@meta.data[, matching_column]
  return(score)
}

##########################################################################
# Call cxds, bcds, or hybrid to obtain doublet scores
##########################################################################
Callscds <- function(count, method, ntop.cxds=500, ntop.bcds=500){
  
  rownames(count) <- as.character(1:dim(count)[1])
  colnames(count) <- as.character(1:dim(count)[2])
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count))
  
  sink('NUL')
  
  tryCatch({
    if(method=='cxds'){
      sce <- scds::cxds(sce, ntop = ntop.cxds)
      CD <- SummarizedExperiment::colData(sce)
      score <- CD$cxds_score
      score <- as.numeric(score)  # remove vector names
      print(score)
    }
    
    if(method=='bcds'){
      sce <- scds::bcds(sce, ntop = ntop.bcds)
      CD <- SummarizedExperiment::colData(sce)
      score <- CD$bcds_score
      score <- as.numeric(score)  # remove vector names
    }
    
    if(method=='hybrid'){
      sce <- scds::cxds_bcds_hybrid(sce)
      CD <- SummarizedExperiment::colData(sce)
      score <- CD$hybrid_score
      score <- as.numeric(score)  # remove vector names
    }
  },
  interrupt = function(e){
    sink()
  }
  )
  sink()
  return(score)
}

#' Calculate doublet scores on single dataset
#'
#' Call different computational doublet-detection methods to calculate doublet scores on single dataset.
#' @param count A scRNA-seq count matrix.
#' @param method A name vector of doublet-detection methods.
#' @param n_neighbors The number of nearest neighbors in KNN classifier (Scrublet).
#' @param min_gene_variability_pctl The top percentile of highly variable genes (Scrublet).
#' @param n_prin_comps Number of principal components used to construct KNN classifer (Scrublet).
#' @param nfeatures Number of highly variable genes (DoubletFinder).
#' @param PCs Number of principal components used to construct KNN classifer (DoubletFinder).
#' @param nf Number of highly variable genes (scDblFinder).
#' @param includePCs The index of principal components to include in the predictors (scDblFinder).
#' @param max_depth Maximum depth of decision trees (scDblFinder).
#' @param k The number of nearest neighbors in KNN classifier (doubletCells).
#' @param d Number of principal components used to construct KNN classifer (doubletCells).
#' @param ntop.cxds Number of top variance genes to consider (cxds).
#' @param ntop.bcds Number of top variance genes to consider (bcds).
#' @param n_components Number of principal components used for clustering (DoubletDetection).
#' @param n_top_var_genes Number of highest variance genes to use (DoubletDetection).
#' @param n_iters Number of fit operations from which to collect p-values (DoubletDetection).
#'
#' @return A list of doublet scores calculated by each doublet-detection method.
#' @export
#'
#' @examples
#' score.list <- FindScores(count = count.list$`J293t-dm`, methods = c('cxds','bcds','hybrid'))
#'
FindScores <- function(count, methods,
                       # Scrublet
                       n_neighbors=round(0.5*sqrt(dim(count)[2])), min_gene_variability_pctl=85L, n_prin_comps=30L,
                       # DoubletFinder
                       nfeatures=2000, PCs=10,
                       # scDblFinder
                       nf=1000, includePCs=5, max_depth=5,
                       # doubletCells
                       k=50, d=50,
                       # cxds and bcds
                       ntop.cxds=500, ntop.bcds=500,
                       # DoubletDetection
                       n_components=30, n_top_var_genes=10000, n_iters=5
){
  score.list <- list()
  for(method in methods){
    tryCatch({
      if(method == 'DoubletFinder'){
        cat("Execute DoubletFinder...\n", file = stderr())
        score <- CallDoubletFinder_Update(count = count, nfeatures, PCs)
        score.list[[method]] <- score
      }
      
      if(method%in%c('cxds', 'bcds', 'hybrid')){
        cat(paste("Execute", method, "...\n"), file = stderr())
        score <- Callscds(count = count, method = method, ntop.cxds, ntop.bcds)
        score.list[[method]] <- score
      }
      
      if(method == 'Scrublet'){
        cat("Execute Scrublet...\n", file = stderr())
        score <- CallScrublet(count = count, n_neighbors, min_gene_variability_pctl, n_prin_comps)
        score.list[[method]] <- score
      }
      
      if(method == 'scDblFinder'){
        cat("Execute scDblFinder...\n", file = stderr())
        score <- CallscDblFinder(count = count, nf, includePCs, max_depth)
        score.list[[method]] <- score
      }
      
      if(method == 'DoubletDetection'){
        cat("Execute DoubletDetection...\n", file = stderr())
        score <- CallDoubletDetection(count = count, n_components, n_top_var_genes, n_iters)
        score.list[[method]] <- score
        cat('\n')
      }
      
      if(method == 'doubletCells'){
        cat("Execute doubletCells...\n", file = stderr())
        score <- CalldoubletCells(count = count, k=k, d=d)
        score.list[[method]] <- score
      }
    }, error = function(e){
      cat('\n', method, " ERROR :", conditionMessage(e), "\n")
      score.list[[method]] <<- NA
    }
    )
  }
  return(score.list)
}


#' Calculate doublet scores on multiple dataset
#'
#' Call different computational doublet-detection methods to calculate doublet scores on multiple datasets.
#' @param count.list A list of scRNA-seq count matrix.
#' @param method A name vector of doublet-detection methods.
#' @param n_neighbors The number of nearest neighbors in KNN classifier (Scrublet).
#' @param min_gene_variability_pctl The top percentile of highly variable genes (Scrublet).
#' @param n_prin_comps Number of principal components used to construct KNN classifer (Scrublet).
#' @param nfeatures Number of highly variable genes (DoubletFinder).
#' @param PCs Number of principal components used to construct KNN classifer (DoubletFinder).
#' @param nf Number of highly variable genes (scDblFinder).
#' @param includePCs The index of principal components to include in the predictors (scDblFinder).
#' @param max_depth Maximum depth of decision trees (scDblFinder).
#' @param k The number of nearest neighbors in KNN classifier (doubletCells).
#' @param d Number of principal components used to construct KNN classifer (doubletCells).
#' @param ntop.cxds Number of top variance genes to consider (cxds).
#' @param ntop.bcds Number of top variance genes to consider (bcds).
#' @param n_components Number of principal components used for clustering (DoubletDetection).
#' @param n_top_var_genes Number of highest variance genes to use (DoubletDetection).
#' @param n_iters Number of fit operations from which to collect p-values (DoubletDetection).
#'
#' @return A list of doublet scores calculated by each doublet-detection method on multiple datasets.
#' @export
#'
#' @examples
#' data.list <- ReadData(path = ".../real_datasets")
#' count.list <- data.list$count
#' methods <- c('Scrublet','doubletCells','cxds','bcds','hybrid','scDblFinder','DoubletDetection','DoubletFinder')
#' score.list.all <- FindScores.All(count.list, methods = methods)
#'
FindScores.All <- function(count.list, methods,
                           # Scrublet
                           n_neighbors=round(0.5*sqrt(dim(count)[2])), min_gene_variability_pctl=85L, n_prin_comps=30L,
                           # DoubletFinder
                           nfeatures=2000, PCs=10,
                           # scDblFinder
                           nf=1000, includePCs=5, max_depth=5,
                           # doubletCells
                           k=50, d=50,
                           # cxds and bcds
                           ntop.cxds=500, ntop.bcds=500,
                           # DoubletDetection
                           n_components=30, n_top_var_genes=10000, n_iters=5){
  score.list.all <- list()
  for(i in 1:length(count.list)){
    count <- count.list[[i]]; dim(count)
    data.name <- names(count.list)[i]; data.name
    cat('\nCalculating doublet scores on dataset: ', data.name, '......\n', file = stderr())
    score.list <- FindScores(count = count, methods = methods,
                             n_neighbors, min_gene_variability_pctl, n_prin_comps, nfeatures, PCs, nf, includePCs, max_depth,
                             k, d, ntop.cxds, ntop.bcds, n_components, n_top_var_genes, n_iters)
    score.list.all[[data.name]] <- score.list
  }
  return(score.list.all)
}


#' Call doublets on one dataset
#'
#' Call doublets based on doublet scores and a user-specified doublet rate.
#' @param score.list A list of doublet scores on one dataset.
#' @param rate A user-specified doublet rate (from 0 to 1).
#'
#' @return A list of doublet indices for different doublet-detection methods.
#' @export
#'
#' @examples
#' doublet.list <- FindDoublets(score.list = score.list, rate = .1)
#'
FindDoublets <- function(score.list, rate){
  
  options(warn=-1)
  doublet.list <- list()
  for(i in 1:length(score.list)){
    score <- score.list[[i]]
    method <- names(score.list)[i]
    if(is.na(score[1])){
      doublet.list[[method]] <- NA
    }else{
      index.doublet <- CallDoublets(score, rate)
      doublet.list[[method]] <- index.doublet
    }
  }
  options(warn=0)
  return(doublet.list)
}


#' Call doublets on multiple datasets
#'
#' Call doublets based on doublet scores and a user-specified doublet rate on multiple datasets.
#' @param score.list.all A list of doublet scores on multiple dataset.
#' @param rate A user-specified doublet rate (from 0 to 1).
#'
#' @return A list of doublet indices for different doublet-detection methods on multiple datasets.
#' @export
#'
#' @examples
#' doublet.list.all <- FindDoublets.All(score.list.all = score.list.all, rate = .1)
#'
FindDoublets.All <- function(score.list.all, rate){
  doublet.list.all <- list()
  for(i in 1:length(score.list.all)){
    score.list <- score.list.all[[i]]
    data.name <- names(score.list.all)[i]; data.name
    #cat('\nCall doublets on dataset:', data.name, '......\n')
    doublet.list <- FindDoublets(score.list = score.list, rate = rate)
    doublet.list.all[[data.name]] <- doublet.list
  }
  return(doublet.list.all)
}


#' Call doublets on multiple datasets and doublet rates
#'
#' Call doublets based on doublet scores and user-specified doublet rates on multiple datasets.
#' @param score.list.all A list of doublet scores on multiple dataset.
#' @param rates A user-specified vector of doublet rates (from 0 to 1).
#'
#' @return A list of doublet indices for different doublet-detection methods on multiple datasets and doublet rates.
#' @export
#'
#' @examples
#' doublet.list.all.rate <- FindDoublets.All.Rate(score.list.all = score.list.all, rates = seq(0.01, 0.25, 0.01))
#'
FindDoublets.All.Rate <- function(score.list.all, rates){
  doublet.list.all.rates <- list()
  cat('\nCall doublets ...\n',file = stderr())
  for(rate in rates){
    #print(rate)
    doublet.list.all <- FindDoublets.All(score.list.all, rate)
    doublet.list.all.rates[[as.character(rate)]] <- doublet.list.all
  }
  return(doublet.list.all.rates)
}


#' Calculate AUPRC or AUROC
#'
#' Calculate AUPRC and AUROC based on doublet score and annotation on single dataset.
#' @param score A vector of doublet scores.
#' @param label A vector of doublet annotations (0/1).
#' @param type A string of "AUPRC" or "AUROC".
#'
#' @return A number of AUPRC or AUROC.
#' @export
#'
FindAUC.Single <- function(score, label, type){
  
  options(warn=-1)
  fg <- score[label==1]
  bg <- score[label==0]
  if(type=='AUPRC'){
    pr <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    auc <- pr$auc.integral
  }
  if(type=='AUROC'){
    roc <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    auc <- roc$auc
  }
  options(warn=0)
  return(auc)
}



#' Calculate AUPRC or AUROC for one dataset
#'
#' Calculate AUPRC and AUROC based on doublet scores and annotations for different doublet-detection methods on one dataset.
#' @param score.list A list of doublet scores on one dataset.
#' @param label A vector of 0/1 doublet annotations.
#' @param type A character of "AUPRC" or "AUROC".
#'
#' @return A list of AUPRCs or AUROCs of different doublet-detection methods.
#' @export
#'
#' @examples
#' auprc.list <- FindAUC(score.list = score.list, label = label.list$`J293t-dm`, type = 'AUPRC')
#' auroc.list <- FindAUC(score.list = score.list, label = label.list$`J293t-dm`, type = 'AUROC')
#'
FindAUC <- function(score.list, label, type){
  
  options(warn=-1)
  auc.list <- list()
  for(i in 1:length(score.list)){
    score <- score.list[[i]]
    method <- names(score.list)[i]
    if(is.na(score[1])){
      auc.list[[method]] <- NA
    }else{
      fg <- score[label==1]
      bg <- score[label==0]
      if(type=='AUPRC'){
        pr <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
        auprc <- pr$auc.integral
        auc.list[[method]] <- auprc
      }
      if(type=='AUROC'){
        roc <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
        auroc <- roc$auc
        auc.list[[method]] <- auroc
      }
    }
  }
  options(warn=0)
  return(auc.list)
}


#' Calculate AUPRC or AUROC for mutiple datasets
#'
#' Calculate AUPRC and AUROC based on doublet scores and annotations for different doublet-detection methods on multiple datasets.
#' @param score.list.all A list of doublet scores on multiple datasets.
#' @param label.list A list of vectors of 0/1 doublet annotations.
#' @param type A character of "AUPRC" or "AUROC".
#'
#' @return A list of AUPRCs or AUROCs of different doublet-detection methods on multiple datasets.
#' @export
#'
#' @examples
#' auprc.list.all <- FindAUC.All(score.list.all, label.list, 'AUPRC')
#'
FindAUC.All <- function(score.list.all, label.list, type){
  auc.list.all <- list()
  for(i in 1:length(score.list.all)){
    score.list <- score.list.all[[i]]
    label <- label.list[[i]]; table(label)
    data.name <- names(score.list.all)[i]; data.name
    cat('\nCalculating', type, 'on dataset:', data.name, '......\n')
    auc.list <- FindAUC(score.list = score.list, label = label, type = type)
    auc.list.all[[data.name]] <- auc.list
  }
  return(auc.list.all)
}


#' Calculate precision, recall, or true negative rate (TNR) of identified doublets
#'
#' Calculate precision, recall, or TNR based on identified doublet indices and true doublet annotations.
#' @param doublet.list A list of doublet indices for different doublet-detection methods.
#' @param label A vector of 0/1 doublet annotations.
#' @param type A character of "precision", "recall", or "TNR".
#'
#' @return A list of precision, recall, or TNR for different doublet-detection methods.
#' @export
#'
#' @examples
#' acc.list <- FindACC(doublet.list = doublet.list, label = label.list$`J293t-dm`, type = 'precision')
#'
FindACC <- function(doublet.list, label, type){
  
  acc.list <- list()
  
  for(i in 1:length(doublet.list)){
    
    size <- length(label)
    pred.doublet <- doublet.list[[i]]
    pred.singlet <- setdiff(1:size, pred.doublet)
    truth.doublet <- which(label==1)
    truth.singlet <- which(label==0)
    method <- names(doublet.list)[i]
    tp <- intersect(pred.doublet, truth.doublet)
    fp <- setdiff(pred.doublet, truth.doublet)
    tn <- intersect(pred.singlet, truth.singlet)
    fn <- setdiff(pred.singlet, truth.singlet)
    
    if(type=='precision'){
      precision <- length(tp)/(length(tp) + length(fp))
      acc.list[[method]] <- precision
    }
    if(type=='recall'){
      recall <- length(tp)/(length(tp) + length(fn))
      acc.list[[method]] <- recall
    }
    if(type=='TNR'){
      TNR <- length(tn)/(length(tn) + length(fp))
      acc.list[[method]] <- TNR
    }
  }
  return(acc.list)
}


#' Calculate precision, recall, or true negative rate (TNR) of identified doublets on multiple datasets.
#'
#' Calculate precision, recall, or TNR based on identified doublet indices and true doublet annotations on multiple datasets.
#' @param doublet.list.all A list of doublet indices for different doublet-detection methods on multiple datasets.
#' @param label.list A list of vectors of 0/1 doublet annotations on multiple datasets.
#'
#' @return A list of precision, recall, or TNR for different doublet-detection methods on multiple datasets.
#' @export
#'
#' @examples
#' precision.list.all <- FindACC.All(doublet.list.all = doublet.list.all, label.list = label.list, type = 'precision')
#' recall.list.all <- FindACC.All(doublet.list.all = doublet.list.all, label.list = label.list, type = 'recall')
#' tnr.list.all <- FindACC.All(doublet.list.all = doublet.list.all, label.list = label.list, type = 'TNR')
#'
FindACC.All <- function(doublet.list.all, label.list, type){
  acc.list.all <- list()
  for(i in 1:length(doublet.list.all)){
    doublet.list <- doublet.list.all[[i]]
    label <- label.list[[i]]
    data.name <- names(doublet.list.all)[i]; data.name
    cat('\nCalculate', type, 'on dataset:', data.name, '......\n')
    acc.list <- FindACC(doublet.list = doublet.list, label = label, type = type)
    acc.list.all[[data.name]] <- acc.list
  }
  return(acc.list.all)
}



#' Remove doublets for single method on single dataset
#'
#' Remove doublets based on identified doublet indices for single method on single dataset.
#' @param count A scRNA-seq data matrix.
#' @param label A vector of cell type annotations.
#' @param doublets A vector of identified doublet indices.
#'
#' @return A list of scRNA-seq data matrix and cell type annotations after removing doublets.
#' @export
#'
#' @examples
#'
RemoveDoublets <- function(count, label, doublets){
  count.removal <- count[,-doublets]; dim(count.removal); dim(count)
  label.removal <- label[-doublets]
  data.removal <- list(count=count.removal, label=label.removal)
  return(data.removal)
}


#' Remove doublets for multiple methods on single dataset
#'
#' Remove doublets based on identified doublet indices for multiple methods on single dataset.
#' @param count A scRNA-seq data matrix.
#' @param label A vector of cell type annotations.
#' @param doublet.list A list of identified doublet indices of different methods.
#'
#' @return A list of scRNA-seq data matrix and cell type annotations after removing doublets by different methods.
#' @export
#'
#' @examples
#' data.removal.list <- RemoveDoublets.Method(count = data.de$count, label = data.de$label.cluster, doublet.list = doublet.list)
#'
RemoveDoublets.Method <- function(count, label, doublet.list){
  data.removal.list <- list()
  for(method in names(doublet.list)){
    #method <- names(doublet.list)[1]; method
    data.removal <- RemoveDoublets(count = count, label = label, doublets = doublet.list[[method]])
    data.removal.list[[method]] <- data.removal
  }
  return(data.removal.list)
}


#' Remove doublets for multiple methods on multiple datasets
#'
#' Remove doublets based on identified doublet indices for different doublet-detection methods on multiple datasets.
#' @param count.list A list of scRNA-seq count matrices.
#' @param label.list A list doublet annoations (0/1).
#' @param doublet.list.all A list of identified doublet indices of different methods on multiple datasets.
#'
#' @return A list of scRNA-seq count matrices and doublet annotations after removing doublets.
#' @export
#'
#' @examples
#' data.removal.all <- RemoveDoublets.All(count.list=count.list, label.list=label.list, doublet.list.all=doublet.list.all)
#'
RemoveDoublets.all <- function(count.list, label.list, doublet.list.all){
  count.removal.list <- list()
  label.removal.list <- list()
  for(name in names(count.list)){
    #print(name)
    count <- count.list[[name]]; dim(count)
    label <- label.list[[name]]; table(label)
    for(method in names(doublet.list.all[[1]])){
      #method <- names(doublet.list.all[[1]])[1]
      #print(method)
      doublets <- doublet.list.all[[name]][[method]]
      data.removal <- RemoveDoublets(count, label, doublets)
      count.removal.list[[name]][[method]] <- data.removal$count
      label.removal.list[[name]][[method]] <- data.removal$label
    }
  }
  data.removal.all <- list(count=count.removal.list, label=label.removal.list)
  return(data.removal.all)
}



#' Remove doublets for multiple methods on multiple datasets and doublet rates
#'
#' Remove doublets based on identified doublet indices for different doublet-detection methods on multiple datasets and doublet rates.
#' @param count.list A list of scRNA-seq count matrices.
#' @param label.list A list doublet annoations (0/1).
#' @param doublet.list.all.rate A list of identified doublet indices of different methods on multiple datasets and doublet rates.
#'
#' @return A list of scRNA-seq count matrices and doublet annotations after removing doublets based on multiple doublet rates.
#' @export
#'
#' @examples
#' doublet.list.all.rate <- FindDoublets.All.Rate(score.list.all = score.list.all, rates = seq(0.01, 0.25, 0.01))
#'
RemoveDoublets.All.Rate <- function(count.list, label.list, doublet.list.all.rate){
  data.removal.all.rate <- list()
  cat("Remove identified doublets ...\n",  file = stderr())
  for(rate in names(doublet.list.all.rate)){
    doublet.list.all <- doublet.list.all.rate[[rate]]
    data.removal.all <- RemoveDoublets.all(count.list, label.list, doublet.list.all)
    data.removal.all.rate[[rate]] <- data.removal.all
  }
  return(data.removal.all.rate)
}
#########

#### Running Different Doublet Detection Algorithms and DEAlgo
#Calculate all
synwithcontam.recluster <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/hto12.rds")
count_lst <- list()
count_lst[["Syn.1"]] <- synwithcontam.recluster@assays$RNA@counts
methods <- c('cxds','bcds','hybrid','DoubletFinder')
score.list.all <- FindScores.All(count_lst, methods)

saveRDS(score.list.all,"~/DEAlgoManuscript/Manuscript_Figures/Fig2B/Fig2B_DFsScore.rds")

library(dealgolorg)
library(dplyr)
library(pracma)
library(Seurat)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(viridisLite)
library(reshape2)
library(gridExtra)

Name <-  "DEAlgo_Fig2B_hto12"
Input <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2B/hto12.rds"
Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2B/"
filteredmatrix=NA
rawmatrix=NA
resol=0.8
overlapRatioList=c(0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
OverlapRatio=0.25
ISThreshold=0
Cutoff=0
gene_n=150

obj <- STEP1A_GlobalMarkers(Input,Output,Name,resol)

obj <- STEP1B_MergingCluster(obj,Output,Name,resol,overlapRatioList,gene_n)

obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n)

obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold)
saveRDS(obj,"~/DEAlgoManuscript/Manuscript_Figures/Fig2B/step1d.rds")
STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n)

obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/step1d.rds")
obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,Cutoff,filteredmatrix,rawmatrix)

PlotContaminationPattern(obj,Output,Name,OverlapRatio)

PlotContaminationPattern <- function(obj,Outdir,Name,OverlapRatio=0.5,CELLANNOTATION = FALSE, verbose = TRUE){
  ###No Return
  message("Plot Contamination Pattern.")
  
  if (CELLANNOTATION){
    message("Using user-annotated clusters.")
    obj@meta.data$annotation_index <- paste0("M",as.numeric(factor(obj@meta.data[,OverlapRatio])))#change the manual annotation to index, eg. CellType0 CellType1 CellType2 CellType3 to 1 2 3 4
    OverlapRatio <- "annotation_index" #User manual cellannotation
  }else{
    OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
    message(paste0("Using ",OverlapRatio))
  }
  
  folder_path_Step2_Output <- paste0(Outdir,Name,"_Step2/Output_",OverlapRatio,"/")
  
  #Load artifacts information and DEAlgo result rds
  contamgeneinfo <- read.csv(paste0(folder_path_Step2_Output,"overlaplst_filtered_",OverlapRatio,"_ContaminationInfo.csv"),na.strings = c("", "NA"))
  
  #DEAlgo Subcluster ID
  contamgeneinfo$MajorSub <- paste0(contamgeneinfo$globalcluster,contamgeneinfo$cluster_local)
  #Summarize ES score and their source of major cluster for each subclusters
  result <- contamgeneinfo %>%
    group_by(MajorSub) %>%
    summarize(
      cluser_reflst = list(cluster_ref),
      dp1lst = list(dp1)
    ) %>%
    ungroup()
  
  # Function to calculate the average ES score (dp1lst) for each source of artifacts (cluster_reflst)
  tabulate_Cluster_Contribution_Score <- function(cluster_reflst, dp1lst) {
    components <- unlist(cluster_reflst)
    dp1_values <- unlist(dp1lst)
    CCS_Matrix <- tapply(dp1_values, components, mean, na.rm = TRUE)
    return(CCS_Matrix)
  }
  
  # For each subclusters (each row in result), calculate the average ES score for each source of artifacts
  tabulate_CCS <- mapply(tabulate_Cluster_Contribution_Score, result$cluser_reflst, result$dp1lst, SIMPLIFY = FALSE)
  
  # List of all source of artifacts which contaminated major cluster X
  lst_of_source_of_artifacts <- unique(unlist(lapply(tabulate_CCS, names)))
  
  # Create a matrix, each row represent one subclusters and each column present each source of artifacts, to store the CCS for each major clusters
  CCS_Matrix <- matrix(NA, nrow = length(tabulate_CCS), ncol = length(lst_of_source_of_artifacts), dimnames = list(NULL, lst_of_source_of_artifacts))
  
  # store the CCS for each major clusters in the matrix
  for (i in seq_along(tabulate_CCS)) {
    CCS_Matrix[i, names(tabulate_CCS[[i]])] <- tabulate_CCS[[i]]
  }
  
  # Replace NA with 0
  CCS_Matrix[is.na(CCS_Matrix)] <- 0
  
  # Convert to dataframe for plotting
  CCS_Matrix <- as.data.frame(CCS_Matrix)
  CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
  data <- CCS_Matrix %>%
    pivot_longer(cols = -MajorSub, names_to = "Component", values_to = "value") #dependency tidyr
  
  for (clusterx in unique(obj@meta.data[,OverlapRatio])){
    # Filter out other clusters, only keep cluster X which wish to display
    filtered_data <- data %>%
      dplyr::filter(startsWith(MajorSub, clusterx))
    
    if (nrow(filtered_data) != 0){
      message(paste0("Plot Contamination Pattern ",clusterx))
      # Stacked Plot of CCS vs Subclusters
      p1 <- ggplot(filtered_data, aes(fill=Component, y=value, x=MajorSub)) +
        geom_bar(position="stack", stat="identity")+
        labs(title = paste0(clusterx),
             x = "Subclusters",
             y = "Contamination Score") +
        scale_fill_manual(values = viridis(length(lst_of_source_of_artifacts)),
                          name = "Source of artifacts") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(size = 8),  # Adjust x-axis title size
              axis.title.y = element_text(size = 8))
      
      ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScoreStackPlot",clusterx,".png"), p1, height = 8, width = 8, dpi = 300)
      
      # Heatmap
      # filtered_data$value[startsWith(as.character(filtered_data$MajorSub), substr(as.character(filtered_data$Component), 1, 2))] <- 0
      # COnvert dataframe
      heatmap_data <- dcast(filtered_data, MajorSub ~ Component, value.var = "value")
      heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts
      
      rownames(heatmap_data) <- heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
      heatmap_data$MajorSub <- NULL #remove MajorSub columns
      
      #heatmap_data$Total <- rowSums(heatmap_data,na.rm = TRUE)
      
      # Summarize DEAlgo Contamination Score (CS) of each subclusters
      CS_table <- obj@meta.data %>%
        group_by(DEAlgo_ClusterID) %>%
        summarize(
          DiffDP1_sum = mean(DiffDP1_sum), #DiffDP1_sum of each cells within each subclusters are same value...
        ) %>%
        ungroup()
      
      #Convert to plotting dataframe
      heatmap_data <- as.data.frame(heatmap_data)  # Convert heatmap_data to a dataframe
      heatmap_data <- heatmap_data %>%
        mutate(Score = CS_table$DiffDP1_sum[match(rownames(heatmap_data), CS_table$DEAlgo_ClusterID)]) #Score
      
      # Convert data to matrix
      heatmap_matrix <- as.matrix(heatmap_data)
      
      # Calculate the color scale for heatmap based on quantile breaks
      quantile_breaks <- function(xs, n = 10) {
        breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
        breaks[!duplicated(breaks)]
      }
      mat_breaks <- quantile_breaks(heatmap_matrix, n = 100)
      
      # Create heatmap
      p2 <- pheatmap(heatmap_matrix,
                     cluster_rows = FALSE,  # Do not cluster rows
                     cluster_cols = FALSE,  # Do not cluster columns
                     main = paste0(clusterx),
                     na_col = "grey",  # Fill missing values with grey
                     color             = viridis(length(mat_breaks) - 1),
                     breaks            = mat_breaks,
                     labels_row = rownames(heatmap_matrix),
                     labels_col = colnames(heatmap_matrix),
                     show_rownames = TRUE,  # Show row names
                     show_colnames = TRUE,
                     angle_col = 0)
      ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScoreHeatmap",clusterx,".png"), p2, height = 2.8, width = 8, dpi = 300)
      
      # Plot multicolumn bar plot
      multicolplot <- data.frame((heatmap_matrix))
      
      rownames_sorted <- rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))),decreasing = T)]
      
      # Reorder the heatmap matrix rows according to the sorted row names
      multicolplot <- multicolplot[rownames_sorted, ]
      
      colnames_sorted <- colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))),decreasing = F)]
      
      # Reorder the heatmap matrix rows according to the sorted row names
      multicolplot <- multicolplot[,colnames_sorted]
      
      multicolplot$Category <- as.character(rownames(multicolplot))
      
      multicolplot$Category <- factor(multicolplot$Category, levels = rownames_sorted)
      
      x_limits <- c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
      custom_breaks <- seq( x_limits[1], x_limits[2], length.out = 3)
      text_size = 5
      
      plot_list <- list()
      for (i in seq(colnames_sorted)) {
        clus <- colnames_sorted[i]
        p <- ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
          geom_bar(stat = "identity", fill = c("#d55046","#8d8d0f","#b488bb","#FDE725")[i]) + # viridis(length(colnames_sorted))[i]) +
          labs(title = paste(clus), y ="", x ="") +  # Set y-axis label for the first plot only
          theme_minimal()+
          theme(panel.grid = element_blank(),axis.line = element_line(color = "black"),
                axis.text = element_text(size = text_size),
                axis.title = element_text(size = text_size),
                plot.title = element_text(size = text_size)) +
          coord_flip()+
          #ylim(x_limits)+
          scale_y_continuous(breaks = custom_breaks, limits = x_limits)+
          if (i != 1){
            theme(axis.text.y = element_blank())
          }
        # Add the plot to the list
        plot_list[[i]] <- p
      }
      
      widths <- c(6.2, rep(5, length(plot_list) - 1))
      
      g1 <- grid.arrange(grobs = plot_list, ncol = length(plot_list),right = "",widths = widths)
      ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScore",clusterx,".png"),g1,  height = 2.5, width = 8, dpi = 300)
      
    }
  }
  
  ###Plot Source of Artifacts Contamination Patterns
  #Plotting and visualize the DEAlgo score of each source of artifacts on the major cluster's UMAP,
  p1 <- FeaturePlot(obj,features = colnames(CCS_Matrix)[colnames(CCS_Matrix) != "MajorSub"])
  ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationPattern_SourceOfArtifacts",".png"),p1,  height = 20, width = 20, dpi = 300)
  
  folder_path_Step2_L1R <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/")
  folder_path_Step2_L1R_Marker <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/Marker/")
  
  qcsclst <- sort(unique(obj@meta.data[,OverlapRatio]))#list of major clusters seurat ID
  cluster_to_consider <- list()
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".csv")
    if (file_name %in% list.files(folder_path_Step2_L1R_Marker)){
      cluster_to_consider <- unlist(c(cluster_to_consider,i))
    }
  }
  
  for (i in cluster_to_consider){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))
      recluster@meta.data <- obj@meta.data[rownames(recluster@meta.data),]#copy-paste the metadata of updated seurat object to subcluster's metadata
      #Plotting and visualize the DEAlgo score of each source of artifacts on the subcluster's UMAP,
      p1 <- FeaturePlot(recluster, features = colnames(CCS_Matrix)[!colnames(CCS_Matrix) %in% c("MajorSub", i)])
      ggsave(filename = paste0(folder_path_Step2_Output,OverlapRatio,"_cluster_",i,"_ContaminationPattern_SourceofArtifacts.png"), p1, height = 10, width = 10, dpi = 300)
    }
  }
}
obj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/DEAlgo_Fig2B_hto12_Step2/Output_Overlap_Ratio_0.25/","DEAlgoResult.rds"))
PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = FALSE)

  
#DoubletFinder
library(DoubletFinder)
pbmcLog <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","hto12.rds"))
pbmcLog <- NormalizeData(pbmcLog)
pbmcLog <- FindVariableFeatures(pbmcLog, selection.method = "vst", nfeatures = 2000)
pbmcLog <- ScaleData(pbmcLog)
pbmcLog <- RunPCA(pbmcLog)
pbmcLog <- RunUMAP(pbmcLog, dims = 1:10)
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pbmc <- paramSweep_v3(pbmcLog, PCs = 1:30, sct = F)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

p1 <- ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

saveRDS(bcmvn_pbmc, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","hto12_pk.rds"))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
pbmcLog <- FindNeighbors(pbmcLog, dims = 1:10, reduction = "pca")
pbmcLog <- FindClusters(pbmcLog)
annotations <- pbmcLog@meta.data$seurat_clusters
DimPlot(pbmcLog, reduction = 'umap', group.by = "seurat_clusters")
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1788*nrow(pbmcLog@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# nExp_poi.adj
# [1] 1329
# > nExp_poi
# [1] 1465
# run doubletFinder 
pbmcLog.Doublet <- doubletFinder_v3(pbmcLog, 
                                   PCs = 1:30, 
                                   pN = 0.25, 
                                   pK = pK, 
                                   nExp = nExp_poi.adj,
                                   reuse.pANN = FALSE, sct = F)

colnam <- colnames(pbmcLog.Doublet@meta.data)
# visualize doublets
p2 <- DimPlot(pbmcLog.Doublet, reduction = 'umap', group.by = grep("DF.",colnam,value = T))#paste0("DF.classifications_0.25_",pK,"_",)
p3 <- DimPlot(pbmcLog.Doublet, reduction = 'umap', group.by = "seurat_clusters")
p4 <- FeaturePlot(pbmcLog.Doublet,features = grep("pANN",colnam,value = T),reduction = 'umap')

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","DoubletFinder.png"), p1+p2+p3+p4, height = 20, width = 20, dpi = 300)
saveRDS(pbmcLog.Doublet, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","DoubletFinderScore.rds"))

# Define a function to check for repeated instances in a string
singletlevel = 3
homoheterochecker <- function(string) {
  if(string == "Negative"){
    return("Negative")
  }
  else{
    words <- unlist(strsplit(string, "_"))
    cell_type <- gsub("-.*", "", words)
    if(length(cell_type) == 1) {
      return("Singlet")
    } 
    else{
      if(length(cell_type) > 1 & length(unique(cell_type)) == 1) {
        return("Homotypic")
      } 
      else {
        return("Heterotypic")
      }
    }
  }
}
###Calculate AUC

#read doublet
DFScore <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/Fig2B_DFsScore.rds")##different doublet detection algorithms

doublet_b4qc <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","DoubletFinderScore.rds"))##DoubletFinder result

dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/DEAlgo_Fig2B_hto12_Step2/Output_Overlap_Ratio_0.25/","DEAlgoResult.rds"))##DEAlgo result

hto12obj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/hto12.rds"))
dealgoseuobj$HTO_classification <- hto12obj$HTO_classification

# hto12obj$HomoHetero <- unlist(lapply(hto12obj$HTO_classification,homoheterochecker))
# hto12obj$Status <- unlist(lapply(hto12obj$HTO_classification, function(x) gsub("-.", "", x)))
# 
# hto12obj$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
# 
# max_value <- max(hto12obj$DiffDP1_sum, na.rm = TRUE)
# 
# # Replace NA values (low IS cluster remove in Step1) with the maximum contamination score in Step2
# hto12obj$DiffDP1_sum[is.na(hto12obj$DiffDP1_sum)] <- max_value *2
# 
# hto12obj$AUC <- NA
# hto12obj$AUC <- ifelse(hto12obj$HomoHetero == "Singlet", 0,hto12obj$AUC)
# hto12obj$AUC <- ifelse(hto12obj$HomoHetero == "Heterotypic", 1,hto12obj$AUC)
# hto12obj$AUC <- ifelse(hto12obj$HomoHetero == "Homotypic", 2,hto12obj$AUC)
# hto12obj$AUC <- ifelse(hto12obj$HomoHetero == "Negative", 3,hto12obj$AUC)
# 
# doublet_b4qc$cxds <- DFScore$Syn.1$cxds
# doublet_b4qc$bcds <- DFScore$Syn.1$bcds
# doublet_b4qc$hybrid <- DFScore$Syn.1$hybrid
# doublet_b4qc$DoubletFinder <- DFScore$Syn.1$DoubletFinder
# doublet_b4qc$DFbfqc_pANN <- doublet_b4qc@meta.data[, grep("pANN", colnames(doublet_b4qc@meta.data), value = TRUE)]
# doublet_b4qc$DFbfqc <- doublet_b4qc@meta.data[, grep("DF.classifications", colnames(doublet_b4qc@meta.data), value = TRUE)]
# 
# hto12obj$cxds <- doublet_b4qc$cxds
# hto12obj$bcds <- doublet_b4qc$bcds
# hto12obj$hybrid <- doublet_b4qc$hybrid
# hto12obj$DoubletFinder <- doublet_b4qc$DoubletFinder
# hto12obj$DFbfqc_pANN <- doublet_b4qc$DFbfqc_pANN
# hto12obj$DFbfqc <- doublet_b4qc$DFbfqc

dealgoseuobj$HomoHetero <- unlist(lapply(dealgoseuobj$HTO_classification,homoheterochecker))
dealgoseuobj$Status <- unlist(lapply(dealgoseuobj$HTO_classification, function(x) gsub("-.", "", x)))

dealgoseuobj$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum

max_value <- max(dealgoseuobj$DiffDP1_sum, na.rm = TRUE)

# Replace NA values (low IS cluster remove in Step1) with the maximum contamination score in Step2
dealgoseuobj$DiffDP1_sum[is.na(dealgoseuobj$DiffDP1_sum)] <- max_value

dealgoseuobj$AUC <- NA
dealgoseuobj$AUC <- ifelse(dealgoseuobj$HomoHetero == "Singlet", 0,dealgoseuobj$AUC)
dealgoseuobj$AUC <- ifelse(dealgoseuobj$HomoHetero == "Heterotypic", 1,dealgoseuobj$AUC)
dealgoseuobj$AUC <- ifelse(dealgoseuobj$HomoHetero == "Homotypic", 2,dealgoseuobj$AUC)
dealgoseuobj$AUC <- ifelse(dealgoseuobj$HomoHetero == "Negative", 3,dealgoseuobj$AUC)

doublet_b4qc$cxds <- DFScore$Syn.1$cxds
doublet_b4qc$bcds <- DFScore$Syn.1$bcds
doublet_b4qc$hybrid <- DFScore$Syn.1$hybrid
doublet_b4qc$DoubletFinder <- DFScore$Syn.1$DoubletFinder
doublet_b4qc$DFbfqc_pANN <- doublet_b4qc@meta.data[, grep("pANN", colnames(doublet_b4qc@meta.data), value = TRUE)]
doublet_b4qc$DFbfqc <- doublet_b4qc@meta.data[, grep("DF.classifications", colnames(doublet_b4qc@meta.data), value = TRUE)]

dealgoseuobj$cxds <- doublet_b4qc$cxds
dealgoseuobj$bcds <- doublet_b4qc$bcds
dealgoseuobj$hybrid <- doublet_b4qc$hybrid
dealgoseuobj$DoubletFinder <- doublet_b4qc$DoubletFinder
dealgoseuobj$DFbfqc_pANN <- doublet_b4qc$DFbfqc_pANN
dealgoseuobj$DFbfqc <- doublet_b4qc$DFbfqc

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$dealgocontamclus <- ifelse(dealgoseuobj$DEAlgo_Contaminated == "Artifact", dealgoseuobj$DEAlgo_ClusterID, NA)

for (i in c("M0","M1","M2","M4")){
  file_name <- paste0("Overlap_Ratio_0.25","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/DEAlgo_Fig2B_hto12_Step2/Overlap_Ratio_0.25_recluster/",file_name))
  
  recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
  recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
  recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
  recluster$DEAlgo_Contaminated <- dealgoseuobj$DEAlgo_Contaminated
  recluster$dealgocontamclus <- dealgoseuobj$dealgocontamclus 
  
  p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
  p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
  p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE,  sizes.highlight = 0.1)
  p5 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE,  sizes.highlight = 0.1)
  
  p3 <- DimPlot(recluster, group.by = "dealgocontamclus")
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","Fig2B_cluster_",i,"_ID.png"), p7+p5+p8+p4+p3, height = 10, width = 15, dpi = 300)
  
}

###Calculate AUROC
auroclst <- list()
aucinput <- dealgoseuobj@meta.data[dealgoseuobj$HomoHetero %in% c("Singlet","Heterotypic"),]

score <- aucinput$DiffDP1_sum
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr1 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df1 <- data.frame(Recall = pr1$curve[,1], Precision = pr1$curve[,2], Dataset = "scCLINIC")

auroc <- pr1$auc#pr$auc.integral
print(paste0("scCLINIC AUROC:",auroc))
auroclst <- c(auroclst,scCLINIC = auroc)

score <- aucinput$DFbfqc_pANN
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr2 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df2 <- data.frame(Recall = pr2$curve[,1], Precision = pr2$curve[,2], Dataset = "DoubletFinder")

auroc <- pr2$auc#pr$auc.integral
print(paste0("DoubletFinder AUROC:",auroc))
auroclst <- c(auroclst,DoubletFinder = auroc)
#CXDS
score <- aucinput$cxds
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr3 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df3 <- data.frame(Recall = pr3$curve[,1], Precision = pr3$curve[,2], Dataset = "cxds")

auroc <- pr3$auc#pr$auc.integral
print(paste0("cxds AUROC:",auroc))
auroclst <- c(auroclst,cxds = auroc)
#BCDS
score <- aucinput$bcds
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr4 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df4 <- data.frame(Recall = pr4$curve[,1], Precision = pr4$curve[,2], Dataset = "bcds")

auroc <- pr4$auc#pr$auc.integral
print(paste0("bcds AUROC:",auroc))
auroclst <- c(auroclst,bcds = auroc)
#Hybrid
score <- aucinput$hybrid
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr5 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df5 <- data.frame(Recall = pr5$curve[,1], Precision = pr5$curve[,2], Dataset = "hybrid")

auroc <- pr5$auc#pr$auc.integral
print(paste0("hybrid AUROC:",auroc))
auroclst <- c(auroclst,hybrid = auroc)
#DF
score <- aucinput$DoubletFinder
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr6 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df6 <- data.frame(Recall = pr6$curve[,1], Precision = pr6$curve[,2], Dataset = "DoubletFinder2")

auroc <- pr6$auc#pr$auc.integral
print(paste0("DoubletFinder.BenchmarkingPaper AUROC:",auroc))
auroclst <- c(auroclst,DoubletFinder.BenchmarkingPaper = auroc)
# Combine the dataframes
pr_combined <- rbind(pr_df1, pr_df2,pr_df3,pr_df4,pr_df5)

# Plot using ggplot2
p1 <- ggplot(pr_combined, aes(x = Recall, y = Precision, color = Dataset)) +
  geom_line() +
  labs(title = "AUROC", x = "FPR", y = "TPR") +
  theme_minimal() +
  scale_color_manual(values = c("#53a0b9", "#9671c3", "#d55046","#be883d","#69a75f"))+  # Customize colors if needed
  theme(panel.grid = element_blank())

#Save AUROC Score
auroclst_df <- data.frame(
  Method = names(auroclst),
  AUROC = unlist(auroclst)
)
rownames(auroclst_df) <- 1:nrow(auroclst_df)
write.table(auroclst_df,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","auroc.csv"), sep = ",")

#Calculate AUPRC
auprclst <- list()
score <- aucinput$DiffDP1_sum
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr1 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df1 <- data.frame(Recall = pr1$curve[,1], Precision = pr1$curve[,2], Dataset = "scCLINIC")

auprc <- pr1$auc.integral#pr$auc.integral
print(paste0("scCLINIC AUPRC:",auprc))
auprclst <- c(auprclst,scCLINIC = auprc)

score <- aucinput$DFbfqc_pANN
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr2 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df2 <- data.frame(Recall = pr2$curve[,1], Precision = pr2$curve[,2], Dataset = "DoubletFinder")

auprc <- pr2$auc.integral#pr$auc.integral
print(paste0("DoubletFinder AUPRC:",auprc))
auprclst <- c(auprclst,DoubletFinder = auprc)
#CXDS
score <- aucinput$cxds
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr3 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df3 <- data.frame(Recall = pr3$curve[,1], Precision = pr3$curve[,2], Dataset = "cxds")

auprc <- pr3$auc.integral#pr$auc.integral
print(paste0("cxds AUPRC:",auprc))
auprclst <- c(auprclst,cxds = auprc)
#BCDS
score <- aucinput$bcds
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr4 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df4 <- data.frame(Recall = pr4$curve[,1], Precision = pr4$curve[,2], Dataset = "bcds")

auprc <- pr4$auc.integral#pr$auc.integral
print(paste0("bcds AUPRC:",bcds = auprc))
auprclst <- c(auprclst,bcds = auprc)
#Hybrid
score <- aucinput$hybrid
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr5 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df5 <- data.frame(Recall = pr5$curve[,1], Precision = pr5$curve[,2], Dataset = "hybrid")

auprc <- pr5$auc.integral#pr$auc.integral
print(paste0("hybrid AUPRC:",auprc))
auprclst <- c(auprclst,hybrid = auprc)
#DF
score <- aucinput$DoubletFinder
label <- aucinput$AUC
fg <- score[label==1]
bg <- score[label==0]

pr6 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr_df6 <- data.frame(Recall = pr6$curve[,1], Precision = pr6$curve[,2], Dataset = "DoubletFinder2")

auprc <- pr6$auc.integral#pr$auc.integral
print(paste0("DoubletFinder.BenchmarkingPaper AUPRC:",auprc))
auprclst <- c(auprclst,DoubletFinder.BenchmarkingPaper = auprc)
# Combine the dataframes
pr_combined <- rbind(pr_df1, pr_df2,pr_df3,pr_df4,pr_df5)

# Plot using ggplot2
p2 <- ggplot(pr_combined, aes(x = Recall, y = Precision, color = Dataset)) +
  geom_line() +
  labs(title = "AUPRC", x = "Recall", y = "Precision") +
  theme_minimal() +
  scale_color_manual(values =  c("#53a0b9", "#9671c3", "#d55046","#be883d","#69a75f"))+  # Customize colors if needed
  theme(panel.grid = element_blank())

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","Fig2B.png"), plot = p1+p2, height = 5, width = 12, dpi = 300)

#Save AUPRC Score
auprclst_df <- data.frame(
  Method = names(auprclst),
  AUPRC = unlist(auprclst)
)
rownames(auprclst_df) <- 1:nrow(auprclst_df)
write.table(auprclst_df,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","auprc.csv"), sep = ",")

#UPSET PLOT
metadata <- dealgoseuobj@meta.data
GroundHetero <- rownames(metadata[metadata$HomoHetero == "Heterotypic",])
# length(GroundHetero)
cxdscall <- rownames(metadata[order(metadata$cxds, decreasing = TRUE), ][1:length(GroundHetero), ])
bcdscall <- rownames(metadata[order(metadata$bcds, decreasing = TRUE), ][1:length(GroundHetero), ])
Hybridcall <- rownames(metadata[order(metadata$hybrid, decreasing = TRUE), ][1:length(GroundHetero), ])
DFcall <- rownames(metadata[order(metadata$DFbfqc_pANN, decreasing = TRUE), ][1:length(GroundHetero), ])
DEAlgocall <- rownames(metadata[order(metadata$DiffDP1_sum, decreasing = TRUE), ][1:length(GroundHetero), ])

#install.packages("UpSetR")
library(UpSetR)

# example of list input (list of named vectors)
listInput <- list(cxds = cxdscall, bcds = bcdscall, Hybrid = Hybridcall, DoubletFinder = DFcall, scCLINIC = DEAlgocall, GroundTruth = GroundHetero)
png(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/upsetplot_ngroundtruth.png"), width = 2000, height = 2000,res=300)
print(upset(fromList(listInput), order.by = "freq", nsets = 6))
dev.off()

# DEAlgoPredict <- rownames(metadata[metadata$DEAlgo_Contaminated == TRUE,])
# 
# cxdscall <- rownames(metadata[order(metadata$cxds, decreasing = TRUE), ][1:533, ])
# bcdscall <- rownames(metadata[order(metadata$bcds, decreasing = TRUE), ][1:533, ])
# Hybridcall <- rownames(metadata[order(metadata$hybrid, decreasing = TRUE), ][1:533, ])
# DFcall <- rownames(metadata[order(metadata$DFbfqc_pANN, decreasing = TRUE), ][1:533, ])
# 
# listInput1 <- list(cxds = cxdscall, bcds = bcdscall, Hybrid = Hybridcall, DoubletFinder = DFcall, DEAlgo = DEAlgoPredict, GroundTruth = GroundHetero)
# png(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/upsetplot_.png"), width = 2000, height = 2000,res=300)
# print(upset(fromList(listInput1), order.by = "freq", nsets = 6))
# dev.off()

##Performance based on default parameters
#DF Hybrid bcds cxds called doublet based on percentage expected provided by 10x genomics
#8193 (total cell before qc) × 0.065 (doublet rate from 10x) -> 533 (expected doublet)
cxdscall <- rownames(metadata[order(metadata$cxds, decreasing = TRUE), ][1:length(GroundHetero), ])
bcdscall <- rownames(metadata[order(metadata$bcds, decreasing = TRUE), ][1:533, ])
Hybridcall <- rownames(metadata[order(metadata$hybrid, decreasing = TRUE), ][1:533, ])
DFcall <- rownames(metadata[order(metadata$DFbfqc_pANN, decreasing = TRUE), ][1:533, ])
DEAlgoPredict <- rownames(metadata[as.numeric(as.character(metadata$DEAlgocluster_Contam)) < singletlevel,]) #contaminated cells = level 1 and 2
listInput1 <- list(cxds = cxdscall, bcds = bcdscall, Hybrid = Hybridcall, DoubletFinder = DFcall, scCLINIC = DEAlgoPredict, GroundTruth = GroundHetero)
png(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/upsetplot_defaultperformance.png"), width = 2000, height = 1500,res=300)
print(upset(fromList(listInput1), order.by = "freq", nsets = 6))
dev.off()


####Visualize Plot
folder_path_Step2_Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2B/DEAlgo_Fig2B_hto12_Step2/Overlap_Ratio_0.25_recluster/"
folder_path_Step2_L1R_Marker <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2B/DEAlgo_Fig2B_hto12_Step2/Overlap_Ratio_0.25_recluster/Marker"
res <- "Overlap_Ratio_0.25"
qcsclst <- sort(unique(dealgoseuobj@meta.data[,res]))
cluster_to_consider <- list()
for (i in sort(qcsclst)){
  file_name <- paste0(res,"_cluster_",i,".csv")
  if (file_name %in% list.files(folder_path_Step2_L1R_Marker)){
    cluster_to_consider <- unlist(c(cluster_to_consider,i))
  }
}
library(gridExtra)
for (i in cluster_to_consider){
  file_name <- paste0(res,"_cluster_",i,".rds")
  if (file_name %in% list.files(folder_path_Step2_Output)){
    recluster <- readRDS(paste0(folder_path_Step2_Output,file_name))
    
    #recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
    recluster$DFbfqc_pANN <- dealgoseuobj$DFbfqc_pANN
    #recluster$DFafqc <- dealgoseuobj$DFafqc
    recluster$DFbfqc <- dealgoseuobj$DFbfqc
    recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
    recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
    recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
    recluster$DEAlgo_Contaminated <- ifelse(as.numeric(levels(recluster$DEAlgocluster_Contam))[recluster$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
    #recluster$DEAlgoref <- dealgoseuobj$DEAlgoref
    recluster$HTO_classification <- hto12obj$HTO_classification
    recluster$HomoHetero <- unlist(lapply(recluster$HTO_classification,homoheterochecker))
    recluster$Status <- unlist(lapply(recluster$HTO_classification, function(x) gsub("-.", "", x)))
    
    
    p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE, sizes.highlight = 0.05)
    p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
    p6 <- DimPlot(recluster, reduction = "umap",group.by = "HomoHetero",raster=FALSE, sizes.highlight = 0.05,pt.size = 0.05)
    p3 <- DimPlot(recluster, reduction = "umap",group.by = "Status",raster=FALSE, sizes.highlight = 0.05,pt.size = 0.05)
    p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE, sizes.highlight = 0.05)
    #p5 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgoref",raster=FALSE, sizes.highlight = 0.05)
    
    #p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
    p2 <- FeaturePlot(recluster, reduction = "umap",features = "DFbfqc_pANN")
    p0 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE, label = F, sizes.highlight = 0.1,pt.size = 0.05)
    p10 <- DimPlot(recluster, reduction = "umap",group.by = "DFbfqc",raster=FALSE, sizes.highlight = 0.05,cols= c("red","grey"),pt.size = 0.05)
    
    # Combine plots into a 3x3 grid
    pa <- grid.arrange(p7,p0,p8,p4,p6,p3,p10,p2, ncol = 3)
    
    # Save the combined plot
    ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/", res, "_cluster_", i, "_rawscore_Validate.png"), pa,
           height = 15, width = 15, dpi = 300)
    
    #ggsave(filename = paste0(folder_path_Step2_Output,res,"_cluster_",i,"_rawscore_Validate.png"), p2+p10+p3+p4+p5+p6+p7+p8+p9, height = 45, width = 5, dpi = 300)
    
  }
}
##M2 Published
i <- "M2"
file_name <- paste0(res,"_cluster_",i,".rds")
recluster <- readRDS(paste0(folder_path_Step2_Output,file_name))

recluster$DFbfqc <- dealgoseuobj$DFbfqc
recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
recluster$HTO_classification <- dealgoseuobj$HTO_classification
recluster$DEAlgo_Contaminated <- ifelse(as.numeric(levels(recluster$DEAlgocluster_Contam))[recluster$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
recluster$dealgocontamclus <- ifelse(recluster$DEAlgo_Contaminated == "Artifact", recluster$DEAlgo_ClusterID, NA)
recluster$Status <- unlist(lapply(recluster$HTO_classification, function(x) gsub("-.", "", x)))

p3 <- DimPlot(recluster, reduction = "umap",group.by = "Status",raster=FALSE, sizes.highlight = 0.05,pt.size = 0.05, cols = 
                c(  
                  "#c4944a",
                  "#d55046",
                  "#11C3C7",
                  "#1e88e5",
                  "#b488bb",
                  "#8d8d0f",  # Unique color 3 (changed)
                  "grey"
                  
                ))
p4 <- DimPlot(recluster, reduction = "umap",group.by = "dealgocontamclus",raster=FALSE, sizes.highlight = 0.05,pt.size = 0.05, cols = c("#8d8d0f","#b488bb","#d55046"))
p10 <- DimPlot(recluster, reduction = "umap",group.by = "DFbfqc",raster=FALSE, sizes.highlight = 0.05,cols= c("red","grey"),pt.size = 0.05)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/", "Published", "_cluster_", i, "_rawscore_Validate.png"), p3+p4+p10,
       height = 5, width = 15, dpi = 300)


originalseu <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/DEAlgo_Fig2B_hto12_Step1/1b_resol_0.8.rds")
originalseu$HTO_classification <- hto12obj$HTO_classification
originalseu$HomoHetero <- unlist(lapply(originalseu$HTO_classification,homoheterochecker))
dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")

p2 <- DimPlot(originalseu,group.by = "HomoHetero", reduction = "umap",raster=FALSE)

p3 <- FeaturePlot(dealgoseuobj,features = "DiffDP1_sum", reduction = "umap",raster=FALSE, label = F) #& scale_color_gradient(limits = c(0, 1))
p6 <- FeaturePlot(dealgoseuobj, reduction = "umap",features = "DFbfqc_pANN")
p8 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DFbfqc",raster=FALSE, label = F, sizes.highlight = 0.1)


p9 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "HomoHetero",raster=FALSE, label = F, sizes.highlight = 0.1)
p10 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "Status",raster=FALSE, label = F, sizes.highlight = 0.1)
p11 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE, label = F, sizes.highlight = 0.1)
p12 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE, label = F, sizes.highlight = 0.1)


p13 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE, label = F, sizes.highlight = 0.1)
p14 <- DimPlot(originalseu, reduction = "umap",group.by = "RNA_snn_res.0.8",raster=FALSE, label = T, sizes.highlight = 0.1, label.size = 7)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","Overlap_Ratio_0.25","_Fig3a.png"), p14+p12+p3+p11+p9+p10+p8+p6+p2, height = 20, width = 22, dpi = 300)


originalseu$DEAlgo_Contaminated <- dealgoseuobj$DEAlgo_Contaminated
originalseu$DEAlgo_Contaminated <- ifelse(is.na(originalseu$DEAlgo_Contaminated),"Low Quality Cell",originalseu$DEAlgo_Contaminated)
originalseu$DFbfqc <- doublet_b4qc$DFbfqc#dealgoseuobj$DFbfqc
originalseu$HomoHetero <- dealgoseuobj$HomoHetero

p2 <- DimPlot(originalseu,group.by = "DEAlgo_Contaminated", reduction = "umap",raster=FALSE, cols = c("red","black","grey"))
p8 <- DimPlot(originalseu, reduction = "umap",group.by = "DFbfqc",raster=FALSE, label = F, sizes.highlight = 0.1, cols = c("red","grey"))
p9 <- DimPlot(originalseu, reduction = "umap",group.by = "HomoHetero",raster=FALSE, label = F, sizes.highlight = 0.1, cols = c(
  "red",
  "#71b84e",
  "black",
  "grey"))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2B/","Overlap_Ratio_0.25","_Fig2BFullUmap.png"), p9+p2+p8, height = 8, width = 24, dpi = 300)


