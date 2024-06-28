library(Seurat)
library(scCustomize)
data.list <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/sim_type.rds") # Do not directly# copy this line (“directory” is where “real_datasets.zip” is downloaded# and uncompressed)
#data from https://zenodo.org/records/4562782
celltype5 <- data.list$count$"5"
celltype5_label <- data.list$label$"5"
celltype5 <- CreateSeuratObject(celltype5)
celltype5@meta.data$label <- celltype5_label
celltype5_singlet <- subset(celltype5,subset = label == "singlet")

standard_seurat_clustering <- function(obj, res = 0.1, npc = 30,VF = TRUE, Verbose = FALSE){
  message("Clustering started.")
  if (ncol(obj) < npc){
    npc = ncol(obj)
  }
  obj <- NormalizeData(obj, verbose = Verbose)
  if (VF){
    obj <- FindVariableFeatures(obj, verbose = Verbose)
  }
  obj <- ScaleData(obj, verbose = Verbose)
  obj <- RunPCA(obj,npcs = npc, verbose = Verbose)
  obj <- FindNeighbors(obj, dims = 1:npc, reduction = "pca", verbose = Verbose)
  obj <- FindClusters(obj, resolution = as.numeric(res), verbose = Verbose)
  obj <- RunUMAP(obj, dims = 1:npc, verbose = Verbose)
  message("Clustering completed.")
  return(obj)
}

create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
}

celltype5_singlet <- standard_seurat_clustering(celltype5_singlet,res = 0.8)
celltype5_singlet$CellType <- paste0("CellType_",celltype5_singlet$RNA_snn_res.0.8)
saveRDS(celltype5_singlet, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_original_1666Singlets.rds"))

#Scdesign3
# remotes::install_github("csoneson/DuoClustering2018"

library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(scran)
library(tidyverse)
theme_set(theme_bw())

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_original_1666Singlets.rds")

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

saveRDS(example_simu,"~/DEAlgoManuscript/Manuscript_Figures/Fig2A/OUTPUT_Syn_5_Cell_Type_10000.rds")

example_simu <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/OUTPUT_Syn_5_Cell_Type_10000.rds")
#Validate Plot
sce <- objconvert
logcounts(sce) <- log1p(counts(sce))

simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))
colData(simu_sce)$library <- colSums(counts(simu_sce))

df1 = colData(sce) %>% as_tibble() %>% select(library) %>% mutate(Method = "Reference")
df2 = colData(simu_sce) %>% as_tibble()  %>% select(library) %>% mutate(Method = "scDesign3")
df = rbind(df1,df2)
p1 <- ggplot(df, aes(x = Method, y = library, fill = Method)) +
  geom_violin(color = "black") +
  geom_point(aes(color = Method), size = 0.1, alpha = 0.2, position = position_jitter(width = 0.3)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal() +
  scale_fill_manual(values = c("#69a75f", "#9671c3")) + # Specify fill colors
  scale_color_manual(values = c("black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "Average Expression Level") +
  theme(legend.position = "none") # Remove legend

simu_sce_seu <- CreateSeuratObject(counts = example_simu$new_count,meta.data = example_simu$new_covariate)
simu_sce_seu <- standard_seurat_clustering(simu_sce_seu,0.8)

p2 <- DimPlot(simu_sce_seu,group.by = "cell_type", cols = c(
  "#7c3b36",
  "#71b84e",
  "#7c4ccb",
  "#c4944a",
  "#cf4ba6",
  "#597b4a",
  "#d55046",
  "#53a0b9",
  "#4e336c",
  "#b488bb"
))

saveRDS(simu_sce_seu,"~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_10000.rds")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/2OUTPUT_Syn_5_Cell_Type_10000.png"), p1+p2, height = 5, width = 10, dpi = 300)

obj$Dataset <- "Reference"
simu_sce_seu$Dataset <- "scdesign3"
obj$cell_type <- obj$RNA_snn_res.0.8
mergeseuobj <- Merge_Seurat_List(c(simu_sce_seu,obj))

mergeseuobj <- standard_seurat_clustering(mergeseuobj,res = 0.8)

p1 <- DimPlot(mergeseuobj,group.by = "Dataset", cols = c(
  "#7c3b36",
  "#71b84e",
  "#7c4ccb",
  "#c4944a",
  "#cf4ba6",
  "#597b4a",
  "#d55046",
  "#53a0b9",
  "#4e336c",
  "#b488bb"
),pt.size = 0.01,alpha = 0.5)
p2 <- DimPlot(mergeseuobj,group.by = "cell_type", cols = c(
  "#7c3b36",
  "#71b84e",
  "#7c4ccb",
  "#c4944a",
  "#cf4ba6",
  "#597b4a",
  "#d55046",
  "#53a0b9",
  "#4e336c",
  "#b488bb"
))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/3OUTPUT_Syn_5_Cell_Type_10000.png"), p1+p2, height = 5, width = 10, dpi = 300)


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

# Set seed for reproducibility <- create batch script named Fig2A_published_batch1.R
set.seed(123)
random_integers <- sample(1:100, 20, replace = FALSE)
proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
DoubletsPerCluster = 100
###
for (seedidx in seq(5)){
  Syn.seurat.sub <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_10000.rds")
  Syn.seurat.sub$CellType<- Syn.seurat.sub$cell_type

  seed1 <- random_integers[seedidx]
  seed2 <- random_integers[seedidx+10]
  outdir <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
  create_folder_if_not_exists(outdir)
  
  for (contamper in proportion_lst){
    falsecell_lst <- list()
    ncol_lst <- list()
    rmcell_lst <- list()
    ctlst <- unique(Syn.seurat.sub@meta.data$CellType)
    for (celltype in sort(ctlst)){
      cluster1 <- subset(Syn.seurat.sub, subset = CellType == celltype)
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
        groundlist <- c(groundlist,paste0(strsplit(j, " ")[[1]][2]))#,"_",strsplit(j, " ")[[1]][4]))
      } else{
        groundlist <- c(groundlist,j)
      }
    }
    
    synwithcontam.recluster$x <- unlist(groundlist)
    
    groundlist <- list()
    for (j in synwithcontam.recluster$Ground){
      if (j != "Singlets"){
        groundlist <- c(groundlist,paste0(strsplit(j, " ")[[1]][4]))#,"_",strsplit(j, " ")[[1]][4]))
      } else{
        groundlist <- c(groundlist,j)
      }
    }
    
    synwithcontam.recluster$y <- unlist(groundlist)
    
    synwithcontam.recluster$CellType <- ifelse(is.na(synwithcontam.recluster$CellType),"Doublet",synwithcontam.recluster$CellType)

    p1 <- DimPlot(synwithcontam.recluster,group.by = "CellType",cols = c("#727cce",
                                                                         "#b4943e",
                                                                         "#c15ca5",
                                                                         "#60a862",
                                                                         "#cb5a4c",
                                                                         "#cecece"),pt.size = 0.1)
    p2 <- DimPlot(synwithcontam.recluster,group.by = "x",cols = c("#727cce",
                                                                  "#b4943e",
                                                                  "#c15ca5",
                                                                  "#60a862",
                                                                  "#cb5a4c",
                                                                  "#cecece"),pt.size = 0.1)
    p3 <- DimPlot(synwithcontam.recluster,group.by = "y",cols = c("#727cce",
                                                                  "#b4943e",
                                                                  "#c15ca5",
                                                                  "#60a862",
                                                                  "#cb5a4c",
                                                                  "#cecece"),pt.size = 0.1)
    
    ggsave(filename = paste0(outdir,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".png"), p1+p2+p3, height = 5, width = 15, dpi = 300)
    
    saveRDS(synwithcontam.recluster,paste0(outdir,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))
    
  }} 
# <- create batch script named Fig2A_published_batch1.R end here

# <- create batch script named Fig2A_published_batch2.R

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

proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
DoubletsPerCluster = 100
for (seedidx in seq(5)){
  outdir <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
  for (contamper in proportion_lst){
    synwithcontam.recluster <- readRDS(paste0(outdir,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))   
    count_lst <- list()
    count_lst[["Syn.1"]] <- synwithcontam.recluster@assays$RNA@counts
    methods <- c('cxds','bcds','hybrid','DoubletFinder')
    score.list.all <- FindScores.All(count_lst, methods)
    
    saveRDS(score.list.all,paste0(outdir,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,"_DFScore.rds"))
    
  }}
# <- create batch script named Fig2A_published_batch2.R end here

# <- create batch script named Fig2A_published_batch3.R

####scCLINIC

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


proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
DoubletsPerCluster = 100

outdir <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
for (contamper in proportion_lst){
  
  Name <- paste0("Syn_5_Cell_Type_scCLINIC_",seedidx,"_",contamper)
  Input <- paste0(outdir,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds")
  Output <- outdir
  filteredmatrix=NA
  rawmatrix=NA
  resol="Manual"
  OverlapRatio="source"
  Cutoff=0
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
  
  obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,Cutoff,filteredmatrix,rawmatrix,CELLANNOTATION = TRUE)
  
  PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)
}


# <- create batch script named Fig2A_published_batch3.R end here
####

# #scCLINIC vs others doublet detection algorithms, plot AUPRC and AUROC (each runs)
# for (seedidx in seq(5)){
#   outdirDFs <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
#   score_DEAlgo <- list()
#   label_lst <- list()
#   for (contamper in proportion_lst){
#     outdirDEAlgo <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/","Syn_5_Cell_Type_scCLINIC_",seedidx, "_",contamper,"_Step2/Output_annotation_index/")
#     
#     pbmc_demuUMAP_Harmony <- readRDS(paste0(outdirDEAlgo,"DEAlgoResult.rds"))
#     
#     originalrds <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))
#     
#     #Assign filter cells scclinic score as the maximum sccinic score in the dataset
#     max_DiffDP1_sum <- max(pbmc_demuUMAP_Harmony@meta.data$DiffDP1_sum)
#     originalrds$DiffDP1_sum <- pbmc_demuUMAP_Harmony$DiffDP1_sum
#     originalrds$DiffDP1_sum <- ifelse(is.na(originalrds$DiffDP1_sum) ,max_DiffDP1_sum,originalrds$DiffDP1_sum)
#     
#     #Ground truth, cells with name "Contam" are artifacts
#     cname <- colnames(originalrds)
#     contamlabel <- grepl("Contam", cname)
#     contamlabel <- as.integer(contamlabel)
#     
#     label_lst[[paste0("SynDataset",contamper)]] <- contamlabel
# 
#     scoreDFs <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,"_DFScore.rds"))
#     
#     score_dp1_sum <- originalrds@meta.data$DiffDP1_sum
#     
#     score_DEAlgo[[paste0("SynDataset",contamper)]]["scCLINIC"] <- list(score_dp1_sum)
#     score_DEAlgo[[paste0("SynDataset",contamper)]]["cxds"] <- list(scoreDFs$Syn.1$cxds)
#     score_DEAlgo[[paste0("SynDataset",contamper)]]["bcds"] <- list(scoreDFs$Syn.1$bcds)
#     score_DEAlgo[[paste0("SynDataset",contamper)]]["hybrid"] <- list(scoreDFs$Syn.1$hybrid)
#     score_DEAlgo[[paste0("SynDataset",contamper)]]["DoubletFinder"] <- list(scoreDFs$Syn.1$DoubletFinder)
#     
#   }
#   auprc.list.all <- FindAUC.All(score_DEAlgo, label_lst, 'AUPRC')
#   
#   auroc.list.all <- FindAUC.All(score_DEAlgo, label_lst, 'AUROC')
#   
#   # transform the output of FindAUC.All to a data frame for visualization
#   
#   result.auprc <- DoubletCollection::ListToDataframe(auprc.list.all, 'boxplot')
#   
#   result.auroc <- DoubletCollection::ListToDataframe(auroc.list.all, 'boxplot')
#   
#   # write.table(result.auprc,".csv",sep = ",")
#   # write.table(result.auroc,".csv",sep = ",")
# 
#   result.auprc$dataset <- as.integer(gsub("SynDataset", "", result.auprc$dataset))
#   result.auroc$dataset <- as.integer(gsub("SynDataset", "", result.auroc$dataset))
#   
#   p1 <- ggplot(result.auroc, aes(x = dataset, y = value, color = method)) +
#     geom_line() +
#     theme_minimal() +
#     geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
#     scale_color_manual(values = c("#69a75f", "#9671c3", "#be883d","#2d5736", "#442d67", "#6e4d1f")) + # Specify fill colors
#     labs(title = "",
#          x = "",
#          y = "AUROC")+
#     scale_y_continuous(limits = c(0, 1))
#   
#   p2 <- ggplot(result.auprc, aes(x = dataset, y = value, color = method)) +
#     geom_line() +
#     theme_minimal() +
#     geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
#     scale_color_manual(values = c("#69a75f", "#9671c3", "#be883d","#2d5736", "#442d67", "#6e4d1f")) + # Specify fill colors
#     labs(title = "",
#          x = "",
#          y = "AUPRC")+
#     scale_y_continuous(limits = c(0, 1))
#   
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_10k_Result_",seedidx,".png"), p1+p2, height = 5, width = 10, dpi = 300)
# }

#scCLINIC vs others doublet detection algorithms, plot AUPRC and AUROC (5 runs)
final_result_lst <- list()
final_result_lst2 <- list()
for (seedidx in seq(5)){
  outdirDFs <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
  score_DEAlgo <- list()
  label_lst <- list()
  for (contamper in proportion_lst){
    outdirDEAlgo <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/","Syn_5_Cell_Type_scCLINIC_",seedidx, "_",contamper,"_Step2/Output_annotation_index/")
    
    pbmc_demuUMAP_Harmony <- readRDS(paste0(outdirDEAlgo,"DEAlgoResult.rds"))
    
    originalrds <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))
    
    max_DiffDP1_sum <- max(pbmc_demuUMAP_Harmony@meta.data$DiffDP1_sum)
    originalrds$DiffDP1_sum <- pbmc_demuUMAP_Harmony$DiffDP1_sum
    originalrds$DiffDP1_sum <- ifelse(is.na(originalrds$DiffDP1_sum) ,max_DiffDP1_sum,originalrds$DiffDP1_sum)
    
    cname <- colnames(originalrds)
    contamlabel <- grepl("Contam", cname)
    contamlabel <- as.integer(contamlabel)
    
    label_lst[[paste0("SynDataset",contamper)]] <- contamlabel
    
    scoreDFs <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,"_DFScore.rds"))
    
    score_dp1_sum <- originalrds@meta.data$DiffDP1_sum
    
    score_DEAlgo[[paste0("SynDataset",contamper)]]["scCLINIC"] <- list(score_dp1_sum)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["cxds"] <- list(scoreDFs$Syn.1$cxds)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["bcds"] <- list(scoreDFs$Syn.1$bcds)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["hybrid"] <- list(scoreDFs$Syn.1$hybrid)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["DoubletFinder"] <- list(scoreDFs$Syn.1$DoubletFinder)
    
  }
  auprc.list.all <- FindAUC.All(score_DEAlgo, label_lst, 'AUPRC')
  
  auroc.list.all <- FindAUC.All(score_DEAlgo, label_lst, 'AUROC')
  
  # transform the output of FindAUC.All to a data frame for visualization
  
  result.auprc <- DoubletCollection::ListToDataframe(auprc.list.all, 'boxplot')
  
  result.auroc <- DoubletCollection::ListToDataframe(auroc.list.all, 'boxplot')
  
  # write.table(result.auprc,".csv",sep = ",")
  # write.table(result.auroc,".csv",sep = ",")
  
  result.auprc$dataset <- as.integer(gsub("SynDataset", "", result.auprc$dataset))
  result.auroc$dataset <- as.integer(gsub("SynDataset", "", result.auroc$dataset))
  
  p1 <- ggplot(result.auroc, aes(x = dataset, y = value, color = method)) +
    geom_line() +
    theme_minimal() +
    geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
    scale_color_manual(values = c("#69a75f", "#9671c3", "#be883d","#2d5736", "#442d67", "#6e4d1f")) + # Specify fill colors
    labs(title = "",
         x = "",
         y = "AUROC")+
    scale_y_continuous(limits = c(0, 1))
  
  p2 <- ggplot(result.auprc, aes(x = dataset, y = value, color = method)) +
    geom_line() +
    theme_minimal() +
    geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
    scale_color_manual(values = c("#69a75f", "#9671c3", "#be883d","#2d5736", "#442d67", "#6e4d1f")) + # Specify fill colors
    labs(title = "",
         x = "",
         y = "AUPRC")+
    scale_y_continuous(limits = c(0, 1))
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_10k_Result_",seedidx,".png"), p1+p2, height = 5, width = 10, dpi = 300)
  
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

p1 <- ggplot(final_result_lst, aes(x = dataset, y = average, color = method)) +
  #geom_boxplot() +
  geom_line() +
  theme_minimal() +
  #geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_errorbar(aes(ymin = average - std, ymax = average + std), width = 2) +  # Add error bars
  scale_color_manual(values = c("#53a0b9", "#9671c3", "#d55046","#be883d","#69a75f")) + # Specify fill colors
  labs(title = "",
       x = "Proportion of Artifacts",
       y = "AUROC")+
  scale_y_continuous(limits = c(0.5, 1))

p2 <- ggplot(final_result_lst2, aes(x = dataset, y = average, color = method)) +
  #geom_boxplot() +
  geom_line() +
  theme_minimal() +
  #geom_point(aes(color = method), size = 2, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_errorbar(aes(ymin = average - std, ymax = average + std), width = 2) +  # Add error bars
  scale_color_manual(values = c("#53a0b9", "#9671c3", "#d55046","#be883d","#69a75f")) + # Specify fill colors
  labs(title = "",
       x = "Proportion of Artifacts",
       y = "AUPRC")+
  scale_y_continuous(limits = c(0, 1))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_10k_Result_",".png"), p1+p2, height = 2.5, width = 10, dpi = 300)

###Calculate ROC PRC Curve for each proportion in each run
DoubletsPerCluster <- 100
proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
for (seedidx in seq(5)){
  outdirDFs <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/")
  for (contamper in proportion_lst){
    score_DEAlgo <- list()
    label_lst <- list()
    
    outdirDEAlgo <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_",seedidx,"/","Syn_5_Cell_Type_scCLINIC_",seedidx, "_",contamper,"_Step2/Output_annotation_index/")
    
    pbmc_demuUMAP_Harmony <- readRDS(paste0(outdirDEAlgo,"DEAlgoResult.rds"))
    
    originalrds <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))
    
    #
    max_DiffDP1_sum <- max(pbmc_demuUMAP_Harmony@meta.data$DiffDP1_sum)
    originalrds$DiffDP1_sum <- pbmc_demuUMAP_Harmony$DiffDP1_sum
    originalrds$DiffDP1_sum <- ifelse(is.na(originalrds$DiffDP1_sum) ,max_DiffDP1_sum,originalrds$DiffDP1_sum)
    
    #Ground truth
    cname <- colnames(originalrds)
    contamlabel <- grepl("Contam", cname)
    contamlabel <- as.integer(contamlabel)
    
    label_lst[[paste0("SynDataset",contamper)]] <- contamlabel
    
    scoreDFs <- readRDS(paste0(outdirDFs,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,"_DFScore.rds"))
    
    score_dp1_sum <- originalrds@meta.data$DiffDP1_sum
    
    score_DEAlgo[[paste0("SynDataset",contamper)]]["scCLINIC"] <- list(score_dp1_sum)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["cxds"] <- list(scoreDFs$Syn.1$cxds)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["bcds"] <- list(scoreDFs$Syn.1$bcds)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["hybrid"] <- list(scoreDFs$Syn.1$hybrid)
    score_DEAlgo[[paste0("SynDataset",contamper)]]["DoubletFinder"] <- list(scoreDFs$Syn.1$DoubletFinder)
    
    ###AUC START HERE

    auroclst <- list()
    
    score <- score_DEAlgo[[1]]$scCLINIC
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr1 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df1 <- data.frame(Recall = pr1$curve[,1], Precision = pr1$curve[,2], Dataset = "scCLINIC")
    
    auroc <- pr1$auc#pr$auc.integral
    print(paste0("scCLINIC AUROC:",auroc))
    auroclst <- c(auroclst,scCLINIC = auroc)
    
    score <- score_DEAlgo[[1]]$DoubletFinder
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr2 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df2 <- data.frame(Recall = pr2$curve[,1], Precision = pr2$curve[,2], Dataset = "DoubletFinder")
    
    auroc <- pr2$auc#pr$auc.integral
    print(paste0("DoubletFinder AUROC:",auroc))
    auroclst <- c(auroclst,DoubletFinder = auroc)
    #CXDS
    score <- score_DEAlgo[[1]]$cxds
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr3 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df3 <- data.frame(Recall = pr3$curve[,1], Precision = pr3$curve[,2], Dataset = "cxds")
    
    auroc <- pr3$auc#pr$auc.integral
    print(paste0("cxds AUROC:",auroc))
    auroclst <- c(auroclst,cxds = auroc)
    #BCDS
    score <- score_DEAlgo[[1]]$bcds
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr4 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df4 <- data.frame(Recall = pr4$curve[,1], Precision = pr4$curve[,2], Dataset = "bcds")
    
    auroc <- pr4$auc#pr$auc.integral
    print(paste0("bcds AUROC:",auroc))
    auroclst <- c(auroclst,bcds = auroc)
    #Hybrid
    score <- score_DEAlgo[[1]]$hybrid
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr5 <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df5 <- data.frame(Recall = pr5$curve[,1], Precision = pr5$curve[,2], Dataset = "hybrid")
    
    auroc <- pr5$auc#pr$auc.integral
    print(paste0("hybrid AUROC:",auroc))
    auroclst <- c(auroclst,hybrid = auroc)

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
    write.table(auroclst_df,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/","ROCPRC_",seedidx, "_",contamper,"_auroc.csv"), sep = ",")
    
    #Calculate AUPRC
    auprclst <- list()
    score <- score_DEAlgo[[1]]$scCLINIC
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr1 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df1 <- data.frame(Recall = pr1$curve[,1], Precision = pr1$curve[,2], Dataset = "scCLINIC")
    
    auprc <- pr1$auc.integral#pr$auc.integral
    print(paste0("scCLINIC AUPRC:",auprc))
    auprclst <- c(auprclst,scCLINIC = auprc)
    
    score <- score_DEAlgo[[1]]$DoubletFinder
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr2 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df2 <- data.frame(Recall = pr2$curve[,1], Precision = pr2$curve[,2], Dataset = "DoubletFinder")
    
    auprc <- pr2$auc.integral#pr$auc.integral
    print(paste0("DoubletFinder AUPRC:",auprc))
    auprclst <- c(auprclst,DoubletFinder = auprc)
    #CXDS
    score <- score_DEAlgo[[1]]$cxds
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr3 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df3 <- data.frame(Recall = pr3$curve[,1], Precision = pr3$curve[,2], Dataset = "cxds")
    
    auprc <- pr3$auc.integral#pr$auc.integral
    print(paste0("cxds AUPRC:",auprc))
    auprclst <- c(auprclst,cxds = auprc)
    #BCDS
    score <- score_DEAlgo[[1]]$bcds
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr4 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df4 <- data.frame(Recall = pr4$curve[,1], Precision = pr4$curve[,2], Dataset = "bcds")
    
    auprc <- pr4$auc.integral#pr$auc.integral
    print(paste0("bcds AUPRC:",bcds = auprc))
    auprclst <- c(auprclst,auprc)
    #Hybrid
    score <- score_DEAlgo[[1]]$hybrid
    label <- label_lst[[1]]
    fg <- score[label==1]
    bg <- score[label==0]
    
    pr5 <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_df5 <- data.frame(Recall = pr5$curve[,1], Precision = pr5$curve[,2], Dataset = "hybrid")
    
    auprc <- pr5$auc.integral#pr$auc.integral
    print(paste0("hybrid AUPRC:",auprc))
    auprclst <- c(auprclst,hybrid = auprc)

    # Combine the dataframes
    pr_combined <- rbind(pr_df1, pr_df2,pr_df3,pr_df4,pr_df5)
    
    # Plot using ggplot2
    p2 <- ggplot(pr_combined, aes(x = Recall, y = Precision, color = Dataset)) +
      geom_line() +
      labs(title = "AUPRC", x = "Recall", y = "Precision") +
      theme_minimal() +
      scale_color_manual(values =  c("#53a0b9", "#9671c3", "#d55046","#be883d","#69a75f"))+  # Customize colors if needed
      theme(panel.grid = element_blank())
    
    ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/","ROCPRC_",seedidx, "_",contamper,"Fig2A.png"), plot = p1+p2, height = 5, width = 12, dpi = 300)
    
    #Save AUPRC Score
    auprclst_df <- data.frame(
      Method = names(auprclst),
      AUROC = unlist(auprclst)
    )
    rownames(auprclst_df) <- 1:nrow(auprclst_df)
    write.table(auprclst_df,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/","ROCPRC_",seedidx, "_",contamper,"_auprc.csv"), sep = ",")
    
  }}

####Visualize Plot
folder_path_Step2_Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_scCLINIC_1_100_Step2/Output_annotation_index/"
folder_path_Step2_L1R_Marker <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_scCLINIC_1_100_Step2/annotation_index_recluster/Marker/"
folder_path_Step2 <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_scCLINIC_1_100_Step2/annotation_index_recluster/"
df_score <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_CP_100_DR_100_DFScore.rds"))
dealgoseuobj <- readRDS(paste0(folder_path_Step2_Output,"DEAlgoResult.rds"))
oriseuobj <- dealgoseuobj #readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_CP_100_DR_100.rds"))

cname <- colnames(oriseuobj)
contamlabel <- grepl("Contam", cname)
contamlabel <- as.integer(contamlabel)
contamlabel <- ifelse(contamlabel == 0, "Singlet", "Artifact")
oriseuobj@meta.data$GroundTruth <- contamlabel


oriseuobj$DFbfqc_pANN <- df_score$Syn.1$DoubletFinder

sorted_cells <- sort(oriseuobj$DFbfqc_pANN, decreasing = TRUE)
top_500_cells <- names(sorted_cells)[1:500]

oriseuobj$DFbfqc <- "Singlet"
oriseuobj@meta.data$DFbfqc[rownames(oriseuobj@meta.data) %in% top_500_cells]<- "Doublet"

res <- "annotation_index"
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
  if (file_name %in% list.files(folder_path_Step2)){
    recluster <- readRDS(paste0(folder_path_Step2,file_name))
    
    #recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
    recluster$DFbfqc_pANN <- oriseuobj$DFbfqc_pANN
    recluster$GroundTruth <- oriseuobj$GroundTruth
    #recluster$DFafqc <- dealgoseuobj$DFafqc
    recluster$DFbfqc <- oriseuobj$DFbfqc
    recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
    recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
    recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
    recluster$DEAlgoref <- dealgoseuobj$DEAlgoref
    
    p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE, sizes.highlight = 0.05)
    p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
    p6 <- DimPlot(recluster, reduction = "umap",group.by = "GroundTruth",raster=FALSE, sizes.highlight = 0.05)
    #p3 <- DimPlot(recluster, reduction = "umap",group.by = "Status",raster=FALSE, sizes.highlight = 0.05)
    p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE, sizes.highlight = 0.05)
    #p5 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgoref",raster=FALSE, sizes.highlight = 0.05)
    
    #p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
    p2 <- FeaturePlot(recluster, reduction = "umap",features = "DFbfqc_pANN")
    #p0 <- DimPlot(recluster, reduction = "umap",group.by = "DFafqc",raster=FALSE, label = T, sizes.highlight = 0.1)
    p10 <- DimPlot(recluster, reduction = "umap",group.by = "DFbfqc",raster=FALSE, sizes.highlight = 0.05)
    
    # Combine plots into a 3x3 grid
    pa <- grid.arrange(p2, p10, p4, p6, p7, p8, ncol = 3)
    
    # Save the combined plot
    ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/", res, "_cluster_", i, "_scCLINIC_1_100.png"), pa,
           height = 10, width = 17, dpi = 300)
    
    #ggsave(filename = paste0(folder_path_Step2_Output,res,"_cluster_",i,"_rawscore_Validate.png"), p2+p10+p3+p4+p5+p6+p7+p8+p9, height = 45, width = 5, dpi = 300)
    
  }
}

p2 <- DimPlot(oriseuobj,group.by = "GroundTruth", reduction = "umap",raster=FALSE)

p3 <- FeaturePlot(dealgoseuobj,features = "DiffDP1_sum", reduction = "umap",raster=FALSE, label = F) #& scale_color_gradient(limits = c(0, 1))
p6 <- FeaturePlot(oriseuobj, reduction = "umap",features = "DFbfqc_pANN")
p8 <- DimPlot(oriseuobj, reduction = "umap",group.by = "DFbfqc",raster=FALSE, label = F, sizes.highlight = 0.1)

p11 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE, label = F, sizes.highlight = 0.1)
#p12 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgoref",raster=FALSE, label = F, sizes.highlight = 0.1)


p13 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE, label = F, sizes.highlight = 0.1)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/",res,"_Fig2A_scCLINIC_1_100.png"), p2+p3+p6+p8+p11+p13, height = 10, width = 17, dpi = 300)


###1
####Visualize Plot
library(ggplot2)
library(Seurat)
folder_path_Step2_Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_scCLINIC_1_1_Step2/Output_annotation_index/"
folder_path_Step2_L1R_Marker <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_scCLINIC_1_1_Step2/annotation_index_recluster/Marker/"
folder_path_Step2 <- "~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_scCLINIC_1_1_Step2/annotation_index_recluster/"
df_score <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_CP_1_DR_100_DFScore.rds"))
dealgoseuobj <- readRDS(paste0(folder_path_Step2_Output,"DEAlgoResult.rds"))
oriseuobj <- dealgoseuobj #readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/Syn_5_Cell_Type_1/Syn_5_Cell_Type_CP_1_DR_100.rds"))

cname <- colnames(oriseuobj)
contamlabel <- grepl("Contam", cname)
contamlabel <- as.integer(contamlabel)
contamlabel <- ifelse(contamlabel == 0, "Singlet", "Artifact")
oriseuobj@meta.data$GroundTruth <- contamlabel


oriseuobj$DFbfqc_pANN <- df_score$Syn.1$DoubletFinder

sorted_cells <- sort(oriseuobj$DFbfqc_pANN, decreasing = TRUE)
top_500_cells <- names(sorted_cells)[1:500]

oriseuobj$DFbfqc <- "Singlet"
oriseuobj@meta.data$DFbfqc[rownames(oriseuobj@meta.data) %in% top_500_cells]<- "Doublet"

res <- "annotation_index"
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
  if (file_name %in% list.files(folder_path_Step2)){
    recluster <- readRDS(paste0(folder_path_Step2,file_name))
    
    #recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
    recluster$DFbfqc_pANN <- oriseuobj$DFbfqc_pANN
    recluster$GroundTruth <- oriseuobj$GroundTruth
    #recluster$DFafqc <- dealgoseuobj$DFafqc
    recluster$DFbfqc <- oriseuobj$DFbfqc
    recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
    recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
    recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
    recluster$DEAlgoref <- dealgoseuobj$DEAlgoref
    
    p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE, sizes.highlight = 0.05)
    p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
    p6 <- DimPlot(recluster, reduction = "umap",group.by = "GroundTruth",raster=FALSE, sizes.highlight = 0.05)
    #p3 <- DimPlot(recluster, reduction = "umap",group.by = "Status",raster=FALSE, sizes.highlight = 0.05)
    p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE, sizes.highlight = 0.05)
    #p5 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgoref",raster=FALSE, sizes.highlight = 0.05)
    
    #p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
    p2 <- FeaturePlot(recluster, reduction = "umap",features = "DFbfqc_pANN")
    #p0 <- DimPlot(recluster, reduction = "umap",group.by = "DFafqc",raster=FALSE, label = T, sizes.highlight = 0.1)
    p10 <- DimPlot(recluster, reduction = "umap",group.by = "DFbfqc",raster=FALSE, sizes.highlight = 0.05)
    
    # Combine plots into a 3x3 grid
    pa <- grid.arrange(p2, p10, p4, p6, p7, p8, ncol = 3)
    
    # Save the combined plot
    ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/", res, "_cluster_", i, "_scCLINIC_1_1.png"), pa,
           height = 10, width = 17, dpi = 300)
    
    #ggsave(filename = paste0(folder_path_Step2_Output,res,"_cluster_",i,"_rawscore_Validate.png"), p2+p10+p3+p4+p5+p6+p7+p8+p9, height = 45, width = 5, dpi = 300)
    
  }
}

p2 <- DimPlot(oriseuobj,group.by = "GroundTruth", reduction = "umap",raster=FALSE)

p3 <- FeaturePlot(dealgoseuobj,features = "DiffDP1_sum", reduction = "umap",raster=FALSE, label = F) #& scale_color_gradient(limits = c(0, 1))
p6 <- FeaturePlot(oriseuobj, reduction = "umap",features = "DFbfqc_pANN")
p8 <- DimPlot(oriseuobj, reduction = "umap",group.by = "DFbfqc",raster=FALSE, label = F, sizes.highlight = 0.1)

p11 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE, label = F, sizes.highlight = 0.1)
#p12 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgoref",raster=FALSE, label = F, sizes.highlight = 0.1)


p13 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE, label = F, sizes.highlight = 0.1)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig2A/",res,"_Fig2A_scCLINIC_1_1.png"), p2+p3+p6+p8+p11+p13, height = 10, width = 17, dpi = 300)
