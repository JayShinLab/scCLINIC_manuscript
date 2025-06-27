#!/usr/bin/env Rscript
data_path = commandArgs(trailingOnly = TRUE)

if (length(data_path) == 0) {
  stop("Please specify path to manuscript data.", call. = FALSE)
}

library(scCLINIC)
library(dplyr)
library(pracma)
library(Seurat)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(viridisLite)
library(reshape2)
library(gridExtra)
library(cowplot)
library(R.utils)
library(fs)
library(patchwork)
library(purrr)
library(DoubletCollection)

text_size <- 12
font_type <- "Arial"
line_width <- 1.2
bar_width <- 2.5

manuscript_colors <-
  c(
    '#53A85F',    '#F1BB72',    '#57C3F3',    '#D6E7A3',    '#3A6963',
    '#E63863',    '#E95C59',    '#E59CC4',    '#AB3282',    '#23452F',
    '#BD956A',    '#8C549C',    '#585658',    '#9FA3A8',    '#E0D4CA',
    '#5F3D69',    '#C5DEBA',    '#58A4C3',    '#E4C755',    '#F7F398',
    '#AA9A59',    '#E39A35',    '#C1E6F3',    '#6778AE',    '#91D0BE',
    '#B53E2B',    '#476D87',    '#712820',    '#DCC1DD',    '#CCE0F5',
    '#CCC9E6',    '#625D9E',    '#68A180',    '#968175',    '#E5D2DD'
  )

outdir <- paste0(data_path, "/scCLINIC_Figures/")

create_folder_if_not_exists(outdir)

create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
}

create_folder_if_not_exists(paste0(data_path,"/Fig2A/"))
create_folder_if_not_exists(paste0(data_path,"/Fig2B/"))
##### Figure 2
### Figure 2A

##################################################################
### Add Rscript for generating scDesign3 dataset to this folder###
library(Seurat)
library(scCustomize)
data.list <- readRDS(paste0(data_path,"/sim_type.rds"))
#Please download sim_type.rds from the synthetic_datasets.zip archive available at https://zenodo.org/records/4562782.

celltype5 <- data.list$count$"5"
celltype5_label <- data.list$label$"5"
celltype5 <- CreateSeuratObject(celltype5)
celltype5@meta.data$label <- celltype5_label
celltype5_singlet <- subset(celltype5,subset = label == "singlet")

create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
}

celltype5_singlet <- standard_seurat_clustering(celltype5_singlet,res = 0.8)
celltype5_singlet$CellType <- paste0("CellType_",celltype5_singlet$RNA_snn_res.0.8)
saveRDS(celltype5_singlet, paste0(data_path,"/Fig2A/Syn_5_Cell_Type_original_1666Singlets.rds"))

#Scdesign3
# remotes::install_github("csoneson/DuoClustering2018")

library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(scran)
library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)

theme_set(theme_bw())

# Load the PBMC dataset
obj <- readRDS(paste0(data_path,"/Fig2A/Syn_5_Cell_Type_original_1666Singlets.rds"))

objconvert <- as.SingleCellExperiment(obj)
colData(objconvert)$cell_type = as.factor(colData(objconvert)$CellType)
colData(objconvert)$library = colSums(counts(objconvert))
rowData(objconvert)$id = rownames(obj)

set.seed(123)
example_simu <- scdesign3(
  sce = objconvert,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = "library",
  mu_formula = "cell_type + offset(log(library))",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "pbmcmapply",
  important_feature = "auto",
  ncell = 10000
)

saveRDS(example_simu,paste0(data_path,"/Fig2A/OUTPUT_Syn_5_Cell_Type_10000.rds"))

example_simu <- readRDS(paste0(data_path,"/Fig2A/OUTPUT_Syn_5_Cell_Type_10000.rds"))

#####Porpotion Dataset
library(scCustomize)
library(DoubletCollection)
library(Matrix)
###Contaminated cell generator
Samplingcell <- function(cluster1,samplingrate){
  cells.to.sample <- ncol(cluster1)
  sampled.cells <- sample(x = cells.to.sample, size = samplingrate, replace = F)
  cluster.sample <- cluster1[,sampled.cells]
  return (cluster.sample)
}

SetProportion <- function(expression_matrix, proportiontocreate) {
  for (i in 1:ncol(expression_matrix)) {
    expression_matrix[, i] <- expression_matrix[, i]/100*proportiontocreate
  }
  return(expression_matrix)
}


set.seed(123) #Set seed for reproducibility
random_integers <- sample(1:100, 20, replace = FALSE)
proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
DoubletsPerCluster = 100

#######(OPTIONAL) can create batch script for different seedidx to run simultaneously
for (seedidx in seq(5)){
  Syn.seurat.sub <- readRDS(paste0(data_path,"/Fig2A/OUTPUT_Syn_5_Cell_Type_10000.rds"))
  Syn.seurat.sub <- CreateSeuratObject(
    counts = Syn.seurat.sub$new_count,
    meta.data = Syn.seurat.sub$new_covariate
  )###This is added for seurat4
  
  Syn.seurat.sub$CellType<- Syn.seurat.sub$cell_type     #changed this for Seurat4
  
  seed1 <- random_integers[seedidx]
  seed2 <- random_integers[seedidx+10]
  outdir_1 <- paste0(paste0(data_path,"/Fig2A/Syn_5_Cell_Type_",seedidx,"/"))
  create_folder_if_not_exists(outdir_1)
  
  for (contamper in proportion_lst){
    falsecell_lst <- list()
    ncol_lst <- list()
    rmcell_lst <- list()
    ctlst <- unique(Syn.seurat.sub$CellType)          #changed this
    for (celltype in sort(ctlst)){
      cluster1 <- subset(Syn.seurat.sub, subset = CellType == celltype)   #"celltype" needs to be changed
      set.seed(seed1)
      sampled <- Samplingcell(cluster1, DoubletsPerCluster)
      rmcell <- colnames(sampled)
      rmcell_lst <- c(rmcell_lst,rmcell)
      sampled.count <- sampled@assays$RNA@counts
      ncol_lst <- c(ncol_lst,ncol(sampled.count))
      sampled.count <- SetProportion(sampled.count,contamper)
      falsecell_lst <- c(falsecell_lst, Matrix(sampled.count, sparse = T))
    }
    
    falsecell_lst2 <- list()
    ncol_lst2 <- list()
    for (celltype in sort(ctlst)){
      cluster1 <- subset(Syn.seurat.sub, subset = CellType == celltype)
      set.seed(seed2)
      sampled <- Samplingcell(cluster1, DoubletsPerCluster)#5cell type got 4 type of contamination
      rmcell <- colnames(sampled)
      rmcell_lst <- c(rmcell_lst,rmcell)
      sampled.count <- sampled@assays$RNA@counts
      ncol_lst2 <- c(ncol_lst2,ncol(sampled.count))
      falsecell_lst2 <- c(falsecell_lst2, sampled.count)
    }
    
    
    no_cell <- min(unlist(ncol_lst))/4
    cocontam <- list()
    combinations <- combn(seq(length(ncol_lst)), 2)
    indextrack <- c(1,1,1,1,1)
    indextrack1 <- c(1,1,1,1,1)
    for (idx in seq(ncol(combinations))){
      x <- combinations[1,idx]
      y <- combinations[2,idx]
      empty_matrix <- Matrix(0, nrow = nrow(falsecell_lst[[x]]), ncol = no_cell, sparse = TRUE)
      cellname <- paste("Contam",x-1,".",y-1,".",seq(no_cell),".",contamper)
      colnames(empty_matrix) <- cellname
      rownames(empty_matrix) <- rownames(falsecell_lst2[[x]])
      for (i in seq(no_cell)) {
        empty_matrix[,i] <- falsecell_lst2[[y]][,indextrack[y]]+falsecell_lst[[x]][,indextrack1[x]]
        indextrack[y] <- indextrack[y] + 1
        indextrack1[x] <- indextrack1[x] + 1
      }
      cocontam <- c(cocontam,empty_matrix)
    }
    
    for (idx in seq(ncol(combinations))){
      x <- combinations[2,idx]
      y <- combinations[1,idx]
      empty_matrix <- Matrix(0, nrow = nrow(falsecell_lst[[x]]), ncol = no_cell, sparse = TRUE)
      cellname <- paste("Contam",x-1,".",y-1,".",seq(no_cell),".",contamper)
      colnames(empty_matrix) <- cellname
      rownames(empty_matrix) <- rownames(falsecell_lst2[[x]])
      for (i in seq(no_cell)) {
        empty_matrix[,i] <- falsecell_lst2[[y]][,indextrack[y]]+falsecell_lst[[x]][,indextrack1[x]]
        indextrack[y] <- indextrack[y] + 1
        indextrack1[x] <- indextrack1[x] + 1
      }
      cocontam <- c(cocontam,empty_matrix)
    }
    
    seurat_object_list <- lapply(unlist(cocontam), CreateSeuratObject)
    merged_object <- Merge_Seurat_List(list_seurat = seurat_object_list)
    
    synwithcontam <- merge(Syn.seurat.sub,merged_object)
    synwithcontam <- synwithcontam[,!colnames(synwithcontam)%in% unlist(rmcell_lst)]
    synwithcontam.recluster <- standard_seurat_clustering(synwithcontam,res = 0.8)
  
    synwithcontam.recluster$Ground <- "Singlets"
    synwithcontam.recluster$Ground <- ifelse(grepl("Contam", colnames(synwithcontam.recluster)),
                                             colnames(synwithcontam.recluster),
                                             synwithcontam.recluster$Ground)
    
    groundlist <- list()
    for (j in synwithcontam.recluster$Ground){
      if (j != "Singlets"){
        groundlist <- c(groundlist,paste0(strsplit(j, " ")[[1]][2]))
      } else{
        groundlist <- c(groundlist,j)
      }
    }
    
    synwithcontam.recluster$x <- unlist(groundlist)
    
    groundlist <- list()
    for (j in synwithcontam.recluster$Ground){
      if (j != "Singlets"){
        groundlist <- c(groundlist,paste0(strsplit(j, " ")[[1]][4]))
      } else{
        groundlist <- c(groundlist,j)
      }
    }
    
    synwithcontam.recluster$y <- unlist(groundlist)
    
    synwithcontam.recluster$CellType <- ifelse(is.na(synwithcontam.recluster$CellType),"Artifact",synwithcontam.recluster$CellType)
    
    saveRDS(synwithcontam.recluster,paste0(outdir_1,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))
    
  }} 
#######(OPTIONAL) batch script end here

#########################################################################
####Adapt from Xi, N., & Li, J. (2021). Benchmarking Computational Doublet-Detection Methods for Single-Cell RNA Sequencing Data.
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

proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
DoubletsPerCluster = 100

#######(OPTIONAL) can create batch script for different seedidx to run simultaneously
for (seedidx in seq(5)){
  outdir_1 <- paste0(data_path,"/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
  for (contamper in proportion_lst){
    synwithcontam.recluster <- readRDS(paste0(outdir_1,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))   
    count_lst <- list()
    count_lst[["Syn.1"]] <- synwithcontam.recluster@assays$RNA@counts
    methods <- c('cxds','bcds','hybrid','DoubletFinder')
    score.list.all <- FindScores.All(count_lst, methods)
    
    saveRDS(score.list.all,paste0(outdir_1,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,"_DFScore.rds"))
    
  }}
#######(OPTIONAL) batch script end here

####scCLINIC

#######(OPTIONAL) can create batch script for different seedidx to run simultaneously

proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
DoubletsPerCluster = 100

for (seedidx in seq(5)){
  outdir_1 <- paste0(data_path,"/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
for (contamper in proportion_lst){
  
  Name <- paste0("Syn_5_Cell_Type_scCLINIC_",seedidx,"_",contamper)
  Input <- paste0(outdir_1,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds")
  Output <- outdir_1

  resol="Manual"
  OverlapRatio="source"
  ISThreshold=0
  gene_n=150
  CELLANNOTATION = TRUE
  
  obj <- readRDS(Input)
  
  obj$source <- ifelse(is.na(obj$cell_type), paste("CellType",obj$y),obj$cell_type)
  
  obj$source <- gsub("CellType ","CellType_",obj$source)
  
  obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)
  
  obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold,CELLANNOTATION = TRUE)
  
  saveRDS(obj,paste0(Output,Name,"_",contamper,"_step1d.rds"))
  
  STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)
  
  obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)
  
  PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)
}
}

#######(OPTIONAL) batch script end here

##################################################################
indir <- paste0(data_path, "/Fig2A/Syn_5_Cell_Type_1/")

Percent100 <-
  readRDS(paste0(indir, "Syn_5_Cell_Type_CP_", 100, "_DR_100.rds"))
Percent50 <-
  readRDS(paste0(indir, "Syn_5_Cell_Type_CP_", 50, "_DR_100.rds"))
Percent10 <-
  readRDS(paste0(indir, "Syn_5_Cell_Type_CP_", 10, "_DR_100.rds"))
Percent1 <-
  readRDS(paste0(indir, "Syn_5_Cell_Type_CP_", 1, "_DR_100.rds"))

Percent100$CellType <- Percent100$cell_type
Percent50$CellType <- Percent50$cell_type
Percent10$CellType <- Percent10$cell_type
Percent1$CellType <- Percent1$cell_type

Percent100$CellType[is.na(Percent100$CellType)] <- "Artifact"
Percent50$CellType[is.na(Percent50$CellType)] <- "Artifact"
Percent10$CellType[is.na(Percent10$CellType)] <- "Artifact"
Percent1$CellType[is.na(Percent1$CellType)] <- "Artifact"

#Extract CellType vector
cell_types <- Percent100@meta.data$CellType

#Reorder 'Artifact' last
new_levels <- c(setdiff(unique(cell_types), "Artifact"), "Artifact")

#Apply factor releveling
Percent100@meta.data$CellType <- factor(cell_types, levels = new_levels)
Percent50@meta.data$CellType <- factor(cell_types, levels = new_levels)
Percent10@meta.data$CellType <- factor(cell_types, levels = new_levels)
Percent1@meta.data$CellType <- factor(cell_types, levels = new_levels)

A1 <-
  DimPlot(Percent100,
          group.by = "CellType",
          cols = manuscript_colors,
          pt.size = 0.1) &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.3, -0.3),
    legend.direction = "horizontal"
  ) &
  ggtitle("Artifact Proportion 100%") &
  theme(
    text = element_text(size = text_size, family = font_type),
    legend.title = element_text(size = text_size, family = font_type)
  ) &
  guides(
    override.aes = list(size = 5),
    nrow = 1,
    color = guide_legend(nrow = 1)
  )
A2 <-
  DimPlot(Percent50,
          group.by = "CellType",
          cols = manuscript_colors,
          pt.size = 0.1) &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  ) &
  ggtitle("50%") &
  theme(
    text = element_text(size = text_size, family = font_type),
    legend.title = element_text(size = text_size, family = font_type)
  )
A3 <-
  DimPlot(Percent10,
          group.by = "CellType",
          cols = manuscript_colors,
          pt.size = 0.1) &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  ) &
  ggtitle("10%") &
  theme(
    text = element_text(size = text_size, family = font_type),
    legend.title = element_text(size = text_size, family = font_type)
  )
A4 <-
  DimPlot(Percent1,
          group.by = "CellType",
          cols = manuscript_colors,
          pt.size = 0.1) &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  ) &
  ggtitle("1%") &
  theme(
    text = element_text(size = text_size, family = font_type),
    legend.title = element_text(size = text_size, family = font_type)
  )

####################################################
### Calculating AUROC and AUPRC#####################
final_result_lst <- list()
final_result_lst2 <- list()
for (seedidx in seq(5)){
  outdirDFs <- paste0(data_path,"/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
  score_Methods <- list()
  label_lst <- list()
  for (contamper in proportion_lst){
    outdir_scCLINIC <- paste0(data_path,"/Fig2A/Syn_5_Cell_Type_",seedidx,"/","Syn_5_Cell_Type_scCLINIC_",seedidx, "_",contamper,"_Step2/Output_annotation_index/")
    
    syn_obj <- readRDS(paste0(outdir_scCLINIC,"scCLINICResult.rds"))#!#scCLINICResult.rds
    
    originalrds <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))
    
    max_scCLINICScore <- max(syn_obj@meta.data$scCLINICScore)#!#scCLINICScore
    originalrds$scCLINICScore <- syn_obj$scCLINICScore#!#scCLINICScore
    originalrds$scCLINICScore <- ifelse(is.na(originalrds$scCLINICScore) ,max_scCLINICScore,originalrds$scCLINICScore)
    
    cname <- colnames(originalrds)
    contamlabel <- grepl("Contam", cname)
    contamlabel <- as.integer(contamlabel)
    
    label_lst[[paste0("SynDataset",contamper)]] <- contamlabel
    
    scoreDFs <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,"_DFScore.rds"))
    
    score_scCLINIC <- originalrds@meta.data$scCLINICScore
    
    score_Methods[[paste0("SynDataset",contamper)]]["scCLINIC"] <- list(score_scCLINIC)
    score_Methods[[paste0("SynDataset",contamper)]]["cxds"] <- list(scoreDFs$Syn.1$cxds)
    score_Methods[[paste0("SynDataset",contamper)]]["bcds"] <- list(scoreDFs$Syn.1$bcds)
    score_Methods[[paste0("SynDataset",contamper)]]["hybrid"] <- list(scoreDFs$Syn.1$hybrid)
    score_Methods[[paste0("SynDataset",contamper)]]["DoubletFinder"] <- list(scoreDFs$Syn.1$DoubletFinder)
    
  }
  auprc.list.all <- FindAUC.All(score_Methods, label_lst, 'AUPRC')
  
  auroc.list.all <- FindAUC.All(score_Methods, label_lst, 'AUROC')
  
  # transform the output of FindAUC.All to a data frame for visualization
  
  result.auprc <- DoubletCollection::ListToDataframe(auprc.list.all, 'boxplot')
  
  result.auroc <- DoubletCollection::ListToDataframe(auroc.list.all, 'boxplot')
  
  # write.table(result.auprc,".csv",sep = ",")
  # write.table(result.auroc,".csv",sep = ",")
  
  result.auprc$dataset <- as.integer(gsub("SynDataset", "", result.auprc$dataset))
  result.auroc$dataset <- as.integer(gsub("SynDataset", "", result.auroc$dataset))
  
  if (length(final_result_lst)==0){
    final_result_lst <- result.auroc
    final_result_lst2 <- result.auprc
  }
  final_result_lst <- cbind(final_result_lst,result.auroc)
  final_result_lst2 <- cbind(final_result_lst2,result.auprc)
}

final_result_lst <- final_result_lst[, !duplicated(t(final_result_lst))] #?
final_result_lst2 <- final_result_lst2[, !duplicated(t(final_result_lst2))] #?

value_columns <- final_result_lst[grep("^value", names(final_result_lst))]

# Take the average of each column
average_values <- rowMeans(value_columns, na.rm = TRUE)
std_values <- apply(value_columns, 1, sd, na.rm = TRUE)

final_result_lst$average <- average_values
final_result_lst$std <- std_values

#PRC
value_columns <- final_result_lst2[grep("^value", names(final_result_lst2))]

# Take the average of each column
average_values <- rowMeans(value_columns, na.rm = TRUE)
std_values <- apply(value_columns, 1, sd, na.rm = TRUE)

final_result_lst2$average <- average_values
final_result_lst2$std <- std_values

write.table(final_result_lst,file = paste0(data_path,"/Fig2A/AUROC_Fig2A.csv"), sep = ",")
write.table(final_result_lst2,file = paste0(data_path,"/Fig2A/AUPRC_Fig2A.csv"), sep = ",")

####################################################
final_result_lst <-
  read.csv(
    paste0(data_path,
      "/Fig2A/AUROC_Fig2A.csv"
    ),
    na.strings = c("", "NA")
  )
final_result_lst2 <-
  read.csv(
    paste0(data_path,
      "/Fig2A/AUPRC_Fig2A.csv"
    ),
    na.strings = c("", "NA")
  )

Figure_2B <-
  ggplot(final_result_lst, aes(x = dataset, y = average, color = method)) +
  #geom_boxplot() +
  geom_line(size = line_width) +
  theme_minimal() +
  #geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_errorbar(aes(ymin = average - std, ymax = average + std), width = bar_width) +  # Add error bars
  scale_color_manual(values = manuscript_colors) + # Specify fill colors
  labs(
    title = "",
    x = "Proportion of Artifacts",
    y = "AUROC",
    color = "Algorithm"
  ) +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme(
    legend.position = c(0.7, 0.3),
    text = element_text(size = text_size, family = font_type),
    legend.text =  element_text(size = text_size, family = font_type),
    legend.key.size = unit(1, 'lines')
  ) &
  theme(
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.line = element_line(color = "black", size = 0.5)
  )

Figure_2C <-
  ggplot(final_result_lst2, aes(x = dataset, y = average, color = method)) +
  #geom_boxplot() +
  geom_line(size = line_width) +
  theme_minimal() +
  #geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_errorbar(aes(ymin = average - std, ymax = average + std), width = bar_width) +  # Add error bars
  scale_color_manual(values = manuscript_colors) + # Specify fill colors
  labs(
    title = "",
    x = "Proportion of Artifacts",
    y = "AUPRC",
    color = "Algorithm"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    legend.position = c(0.8, 0.3),
    text = element_text(size = text_size, family = font_type),
    legend.text =  element_text(size = text_size, family = font_type),
    legend.key.size = unit(1, 'lines')
  ) &
  theme(
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.line = element_line(color = "black", size = 0.5)
  )
################################################
### B. CITE-Seq dataset #############
##Following Seurat cell hashing tutorial https://satijalab.org/seurat/articles/hashing_vignette.html

# Read in UMI count matrix for RNA
hto12.umis <- readRDS(paste0(data_path,"/hto12_umi_mtx.rds"))
#Please download hto12_umi_mtx.rds from https://www.dropbox.com/scl/fo/tygiouyv6spn8x0coyyau/AK3HDK42JJbkNMjQa53mtqA?rlkey=3urtt7msejtbwnhflkxopm6zr&e=1&st=yrf1txrq&dl=0
#archive available at Seurat cell hashing tutorial

# Read in HTO count matrix
hto12.htos <- readRDS(paste0(data_path,"/hto12_hto_mtx.rds"))
#Please download hto12_hto_mtx.rds from https://www.dropbox.com/scl/fo/tygiouyv6spn8x0coyyau/AK3HDK42JJbkNMjQa53mtqA?rlkey=3urtt7msejtbwnhflkxopm6zr&e=1&st=yrf1txrq&dl=0
#archive available at Seurat cell hashing tutorial

# Select cell barcodes detected in both RNA and HTO
cells.use <- intersect(rownames(hto12.htos), colnames(hto12.umis))

# Create Seurat object and add HTO data
hto12 <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(hto12.umis[, cells.use]), sparse = T),min.features = 300)

hto12[["HTO"]] <- CreateAssayObject(counts = t(x = hto12.htos[colnames(hto12), 1:12]))

# Normalize data
hto12 <- NormalizeData(hto12)

hto12 <- NormalizeData(hto12, assay = "HTO", normalization.method = "CLR")

hto12 <- HTODemux(hto12, assay = "HTO", positive.quantile = 0.99)

saveRDS(hto12, paste0(data_path,"/Fig2B/hto12.rds"))

###### Apply scCLINIC and DoubletFinder ########
Name <-  "Fig2B_hto12"
Input <- paste0(data_path,"/Fig2B/hto12.rds")
Output <- paste0(data_path,"/Fig2B/")

resol=0.8
overlapRatioList=c(0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
OverlapRatio=0.25
ISThreshold=0
gene_n=150

obj <- STEP1A_GlobalMarkers(Input,Output,Name,resol)

obj <- STEP1B_MergingCluster(obj,Output,Name,resol,overlapRatioList,gene_n)

obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n)

obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold)

saveRDS(obj,paste0(data_path,"/Fig2B/scCLINIC_step1d.rds"))

STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n)

obj <- readRDS(paste0(data_path,"/Fig2B/scCLINIC_step1d.rds")) #change typo

obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n)

PlotContaminationPattern(obj,Output,Name,OverlapRatio)

#DoubletFinder
library(DoubletFinder)

hto_obj <- readRDS(paste0(data_path,"/Fig2B/","hto12.rds"))
hto_obj <- NormalizeData(hto_obj)
hto_obj <- FindVariableFeatures(hto_obj, selection.method = "vst", nfeatures = 2000)
hto_obj <- ScaleData(hto_obj)
hto_obj <- RunPCA(hto_obj)
hto_obj <- RunUMAP(hto_obj, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pbmc <- paramSweep_v3(hto_obj, PCs = 1:30, sct = F)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

p1 <- ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
hto_obj <- FindNeighbors(hto_obj, dims = 1:10, reduction = "pca")
hto_obj <- FindClusters(hto_obj)
annotations <- hto_obj@meta.data$seurat_clusters
DimPlot(hto_obj, reduction = 'umap', group.by = "seurat_clusters")
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.06*nrow(hto_obj@meta.data))  ##NOT 0.1788 Assuming 6% doublet formation rate - for 8193 cells recovered
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# nExp_poi.adj
# [1] 1329
# > nExp_poi
# [1] 1465
# run doubletFinder 
hto_obj.Doublet <- doubletFinder_v3(hto_obj, 
                                    PCs = 1:30, 
                                    pN = 0.25, 
                                    pK = pK, 
                                    nExp = nExp_poi.adj,
                                    reuse.pANN = FALSE, sct = F)

saveRDS(hto_obj.Doublet, paste0(data_path,"/Fig2B/DoubletFinderScore.rds"))

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
################################################
#Check the stats...
#read doublet
doublet_b4qc <- readRDS(paste0(data_path,"/Fig2B/DoubletFinderScore.rds"))##DoubletFinder result

scclinicseuobj <- readRDS(paste0(data_path,"/Fig2B/Fig2B_hto12_Step2/Output_Overlap_Ratio_0.25/scCLINICResult.rds"))

hto12obj <- readRDS(paste0(data_path,"/Fig2B/hto12.rds"))

scclinicseuobj$HTO_classification <- hto12obj$HTO_classification

scclinicseuobj$HomoHetero <- unlist(lapply(scclinicseuobj$HTO_classification,homoheterochecker))

scclinicseuobj$Status <- unlist(lapply(scclinicseuobj$HTO_classification, function(x) gsub("-.", "", x)))

originalseu <- readRDS(paste0(data_path,"/Fig2B/Fig2B_hto12_Step1/1b_resol_0.8.rds"))

originalseu$HTO_classification <- hto12obj$HTO_classification

originalseu$HomoHetero <- unlist(lapply(originalseu$HTO_classification,homoheterochecker))

scclinicseuobj$scCLINIC_prediction <- ifelse(as.numeric(levels(scclinicseuobj$scCLINIC_Level))[scclinicseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")

originalseu$scCLINIC_prediction <- scclinicseuobj$scCLINIC_prediction

originalseu$scCLINIC_prediction <- ifelse(is.na(originalseu$scCLINIC_prediction),"Low Quality Cell",originalseu$scCLINIC_prediction)

originalseu$DFbfqc <- doublet_b4qc@meta.data[, grep("DF.classifications", colnames(doublet_b4qc@meta.data), value = TRUE)]

originalseu$HomoHetero <- scclinicseuobj$HomoHetero

saveRDS(originalseu,paste0(data_path,"/Fig2B/originalseu_plot.rds"))

originalseu <-
  readRDS(paste0(data_path,"/Fig2B/originalseu_plot.rds"))
originalseu$scCLINIC_prediction <- 				###SCCLINIC
  ifelse(
    originalseu$scCLINIC_prediction == "Low Quality Cell",	###SCCLINIC
    "Low Identity Score",
    originalseu$scCLINIC_prediction				###SCCLINIC
  )
  
D2 <-
  DimPlot(
    originalseu,
    group.by = "scCLINIC_prediction",				###SCCLINIC
    reduction = "umap",
    raster = FALSE,
    cols = c("#E63863", "#53A85F", "#9FA3A8")
  ) & labs(color = "scCLINIC") &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.1, 0.2)
  ) & ggtitle(NULL) #axis.text = element_blank(),

SA1 <-
  DimPlot(
    originalseu,
    reduction = "umap",
    group.by = "DFbfqc",
    raster = FALSE,
    label = F,
    sizes.highlight = 0.1,
    cols = c("#E63863", "#9FA3A8")
  )
D1 <-
  DimPlot(
    originalseu,
    reduction = "umap",
    group.by = "HomoHetero",
    raster = FALSE,
    label = F,
    sizes.highlight = 0.1,
    cols = c("#E63863",
             "#F1BB72",
             "#53A85F",
             "#9FA3A8")
  ) & labs(color = "CITE-seq") &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.1, 0.2)
  ) & ggtitle(NULL)
  
step1d <-
  readRDS(paste0(data_path,"/Fig2B/scCLINIC_step1d.rds")) #correct typo

# Replace elements in step1d dataset
step1d$Overlap_Ratio_0.25 <-
  dplyr::recode(
    step1d$Overlap_Ratio_0.25,
    "M1" = "HEK_M1",
    "M2" = "THP1_M2",
    "M3" = "K562_M3",
    "M4" = "KG1_M4"
  )

E1 <-
  DimPlot(
    step1d,
    reduction = 'umap',
    group.by = "Overlap_Ratio_0.25",
    label = TRUE,
    repel = TRUE,
    cols = manuscript_colors
  ) &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  ) & ggtitle(NULL)
LabelClusters(E1, id = "Overlap_Ratio_0.25") + theme(legend.text = element_text(size = text_size, family = font_type))

##############################################
### Plot M3 Subcluster###########
##############################################
recluster <- readRDS(paste0(data_path,"/Fig2B/Fig2B_hto12_Step2/Overlap_Ratio_0.25_recluster/Overlap_Ratio_0.25_cluster_M3.rds"))

recluster$DFbfqc <- doublet_b4qc@meta.data[, grep("DF.classifications", colnames(doublet_b4qc@meta.data), value = TRUE)]
recluster$scCLINIC_ClusterID <- scclinicseuobj$scCLINIC_ClusterID
recluster$scCLINIC_Level <- scclinicseuobj$scCLINIC_Level
recluster$HTO_classification <- scclinicseuobj$HTO_classification
recluster$scCLINIC_prediction <- ifelse(as.numeric(levels(recluster$scCLINIC_Level))[recluster$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
recluster$scCLINIC_Artifact_Cluster <- ifelse(recluster$scCLINIC_prediction == "Artifact", recluster$scCLINIC_ClusterID, NA)
recluster$Status <- unlist(lapply(recluster$HTO_classification, function(x) gsub("-.", "", x)))

saveRDS(recluster,paste0(data_path,"/Fig2B/recluster_plot.rds"))

recluster <-
  readRDS(paste0(data_path,"/Fig2B/recluster_plot.rds"))

E2 <-
  DimPlot(
    recluster,
    reduction = "umap",
    group.by = "Status",
    raster = FALSE,
    sizes.highlight = 0.05,
    pt.size = 0.05,
    cols =
      c(
        '#F3B1A0',
        '#53A85F',
        '#F1BB72',
        '#3A6963',
        '#57C3F3',
        '#D6E7A3',
        '#E63863'
        
      )
  ) & labs(color = "CITE-Seq") &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.65, 0.25)
  ) & ggtitle("K562_M3 Subcluster")


# Replace elements in step1d dataset
recluster$scCLINIC_Artifact_Cluster <-				###SCCLINIC
  ifelse(
    is.na(recluster$scCLINIC_Artifact_Cluster),			###SCCLINIC
    "Singlet",
    sub("^M3_", "", recluster$scCLINIC_Artifact_Cluster)		###SCCLINIC
  )

E3 <-
  DimPlot(
    recluster,
    reduction = "umap",
    group.by = "scCLINIC_Artifact_Cluster",			###SCCLINIC
    raster = FALSE,
    sizes.highlight = 0.05,
    pt.size = 0.05,
    cols = c('#D6E7A3', '#57C3F3', '#53A85F', '#F1BB72'),
    label = TRUE,
    repel = TRUE
  ) &
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  ) & ggtitle("K562_M3 scCLINIC")
LabelClusters(E3, id = "scCLINIC_Artifact_Cluster") + 		###SCCLINIC
	theme(legend.text = element_text(size = text_size, family = font_type))

SA2 <-
  DimPlot(
    recluster,
    reduction = "umap",
    group.by = "DFbfqc",
    raster = FALSE,
    sizes.highlight = 0.05,
    cols = c("red", "grey"),
    pt.size = 0.05
  )
  
scCLINICseuobj <-
  readRDS(
    paste0(
      data_path,"/Fig2B/Fig2B_hto12_Step2/Output_Overlap_Ratio_0.25/scCLINICResult.rds"	###SCCLINIC
    )
  )

#THP1
ELANE_S <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    features = c("ELANE"),
    pt.size = 0.05
  ) &
  scale_color_gradientn(
    colors = c("grey", '#E63863'),
    limits = c(0, 5),
    oob = scales::squish
  )
ELANE_G <-
  FeaturePlot(scCLINICseuobj,
              reduction = "umap",
              features = c("ELANE")) &
  scale_color_gradientn(
    colors = c("grey", '#E63863'),
    limits = c(0, 5),
    oob = scales::squish
  )

#HEK
CDKN2A_S <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    features = c("CDKN2A"),
    pt.size = 0.05
  ) &
  scale_color_gradientn(
    colors = c("grey", '#E63863'),
    limits = c(0, 4),
    oob = scales::squish
  )
CDKN2A_G <-
  FeaturePlot(scCLINICseuobj,
              reduction = "umap",
              features = c("CDKN2A")) &
  scale_color_gradientn(
    colors = c("grey", '#E63863'),
    limits = c(0, 4),
    oob = scales::squish
  )

#KG1
HLADRA_S <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    features = c("HLA-DRA"),
    pt.size = 0.05
  ) &
  scale_color_gradientn(
    colors = c("grey", '#E63863'),
    limits = c(0, 5),
    oob = scales::squish
  )
HLADRA_G <-
  FeaturePlot(scCLINICseuobj,
              reduction = "umap",
              features = c("HLA-DRA")) &
  scale_color_gradientn(
    colors = c("grey", '#E63863'),
    limits = c(0, 5),
    oob = scales::squish
  )

Figure_2G <-
  wrap_elements(((HLADRA_G |
                    ELANE_G |
                    CDKN2A_G) &
                   theme(plot.title = element_text(
                     face = "italic",
                     size = text_size,
                     family = font_type
                   ))
  ) / ((HLADRA_S |
          ELANE_S |
          CDKN2A_S) &
         ggtitle(NULL)) &
    theme(
      text = element_text(size = text_size, family = font_type),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank()
    )
  )

#Load artifacts information and scCLINIC result rds
contamgeneinfo <-
  read.csv(
    paste0(
      data_path,"/Fig2B/Fig2B_hto12_Step2/Output_Overlap_Ratio_0.25/ArtifactsInfo.csv"	###SCCLINIC
    ),
    na.strings = c("", "NA")
  )

#scCLINIC Subcluster ID
contamgeneinfo$MajorSub <-
  paste0(contamgeneinfo$Major_Cluster,
         "_",
         contamgeneinfo$Subcluster)
#Summarize ES score and their source of major cluster for each subclusters
result <- contamgeneinfo %>%
  group_by(MajorSub) %>%
  summarize(cluser_reflst = list(Source_of_Artifacts_SoA),
            dp1lst = list(Enrichment_Score)) %>%
  ungroup()

# Function to calculate the average ES score (dp1lst) for each source of artifacts (cluster_reflst)
tabulate_Cluster_Contribution_Score <-
  function(cluster_reflst, dp1lst) {
    components <- unlist(cluster_reflst)
    dp1_values <- unlist(dp1lst)
    CCS_Matrix <- tapply(dp1_values, components, mean, na.rm = TRUE)
    return(CCS_Matrix)
  }

# For each subclusters (each row in result), calculate the average ES score for each source of artifacts
tabulate_CCS <-
  mapply(
    tabulate_Cluster_Contribution_Score,
    result$cluser_reflst,
    result$dp1lst,
    SIMPLIFY = FALSE
  )

# List of all source of artifacts which contaminated major cluster X
lst_of_source_of_artifacts <-
  unique(unlist(lapply(tabulate_CCS, names)))

# Create a matrix, each row represent one subclusters and each column present each source of artifacts, to store the CCS for each major clusters
CCS_Matrix <-
  matrix(
    NA,
    nrow = length(tabulate_CCS),
    ncol = length(lst_of_source_of_artifacts),
    dimnames = list(NULL, lst_of_source_of_artifacts)
  )

# store the CCS for each major clusters in the matrix
for (i in seq_along(tabulate_CCS)) {
  CCS_Matrix[i, names(tabulate_CCS[[i]])] <- tabulate_CCS[[i]]
}

# Replace NA with 0
CCS_Matrix[is.na(CCS_Matrix)] <- 0

# Convert to dataframe for plotting
CCS_Matrix <- as.data.frame(CCS_Matrix)
CCS_Matrix$MajorSub <-
  result$MajorSub #Named each rows with their Subclusters ID
data <- CCS_Matrix %>%
  pivot_longer(cols = -MajorSub,
               names_to = "Component",
               values_to = "value") #dependency tidyr

#for (clusterx in unique(obj@meta.data[,OverlapRatio])){
# Filter out other clusters, only keep cluster X which wish to display
clusterx <- "M3"
filtered_data <- data %>%
  dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

heatmap_data <-
  dcast(filtered_data, MajorSub ~ Component, value.var = "value")
heatmap_data <-
  heatmap_data[,!colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

rownames(heatmap_data) <-
  heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL #remove MajorSub columns

#heatmap_data$Total <- rowSums(heatmap_data,na.rm = TRUE)

# Summarize scCLINIC Score (CS) of each subclusters
CS_table <- scCLINICseuobj@meta.data %>%
  group_by(scCLINIC_ClusterID) %>% 							###SCCLINIC
  summarize(scCLINICScore = mean(scCLINICScore), #scCLINICScore of each cells within each subclusters are same value... #
            scCLINICScore)  %>% ungroup()
            
#Convert to plotting dataframe
heatmap_data <-
  as.data.frame(heatmap_data)  # Convert heatmap_data to a dataframe
heatmap_data <- heatmap_data %>%
  mutate(scCLINIC_Score = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)])  ###SCCLINIC

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <- data.frame((heatmap_matrix))

colnames(multicolplot) <- dplyr::recode(
  colnames(multicolplot),
  "M1" = "HEK_M1",
  "M2" = "THP1_M2",
  "M3" = "K562_M3",
  "M4" = "KG1_M4"
)
rownames(multicolplot) <-
  gsub(".*_", "", rownames(multicolplot))

rownames_sorted <-
  rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
multicolplot <- multicolplot[rownames_sorted,]

colnames_sorted <-
  colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
multicolplot <- multicolplot[, colnames_sorted]


multicolplot$Category <- as.character(rownames(multicolplot))
multicolplot$Category <-
  factor(multicolplot$Category, levels = rownames_sorted)

x_limits <-
  c(0, ceiling(max(multicolplot[,-ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <- seq(x_limits[1], x_limits[2], length.out = 3)

# Create individual plots without y-axis text
plot_list <- list()
for (i in seq_along(colnames_sorted)) {
  clus <- colnames_sorted[i]
  p <-
    ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(stat = "identity",
             fill = c('#53A85F', '#D6E7A3', '#57C3F3', '#E63863')[i]) + 

    labs(title = paste(clus),
         y = "",
         x = "Cluster") +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "grey", size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = text_size, family = font_type),
      axis.title = element_text(size = text_size, family = font_type),
      plot.title = element_text(size = text_size, family = font_type),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    coord_flip() +
    scale_y_continuous(breaks = custom_breaks, limits = x_limits)
  
  plot_list[[i]] <- p
}

# Remove y-axis labels, title, and ticks from all but the left-most plot
plot_list[2:length(plot_list)] = plot_list[2:length(plot_list)] %>%
  map(
    ~ .x + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  )


# Combine the y-axis labels plot with the other plots
Figure_2F <-
  wrap_elements(
    plot_grid(
      plotlist = plot_list,
      ncol = length(plot_list),
      align = "v"
    ) + draw_label(
      "K562_M3 Artifact Enrichment Score",
      fontface = 'bold',
      y = 1.1,
      fontfamily = font_type,
      size = text_size
    ) + draw_label(
      "Enrichment Score",
      fontface = NULL,
      x = 0.15,
      y = 0,
      fontfamily = font_type,
      size = text_size
    )
  )


Figure_2A <-
  wrap_elements((A1 | A2 | A3 | A4) &
                  theme(
                    legend.spacing.x = unit(0.2, "lines"),
                    text = element_text(size = text_size, family = font_type)
                  ))
Figure_2D <-
  wrap_elements((D1 | D2) &
                  theme(
                    legend.key.size = unit(0, 'lines'),
                    legend.spacing.x = unit(0.2, "lines"),
                    text = element_text(size = text_size, family = font_type)
                  ))
Figure_2E <-
  wrap_elements((E1 | E2 | E3) &
                  theme(
                    legend.key.size = unit(0, 'lines'),
                    legend.spacing.x = unit(0.2, "lines"),
                    text = element_text(size = text_size, family = font_type),
                    legend.text = element_text(size = text_size, family = font_type)
                  )
  )

thm <-
  theme(text = element_text(
    size = text_size,
    family = font_type,
    face = "bold"
  ))
figure_2 <-
  (((Figure_2A | Figure_2B | Figure_2C) + plot_layout(widths = c(0.6, 0.2, 0.2), ncol = 3)) /
     (Figure_2D | Figure_2E) /
     (Figure_2F | Figure_2G) + plot_layout(heights = c(0.3, 0.4, 0.4), nrow = 3) +
      plot_annotation(tag_levels = "A", theme = thm)
  )
            
ggsave(
  filename = paste0(outdir, "Figure_2.png"),
  plot = figure_2,
  height = 12,
  width = 18,
  dpi = 300
)

            
