seedidx = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(seedidx)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
seedidx <- as.integer(seedidx)

#####Porpotion Dataset
library(scCustomize)
library(DoubletCollection)
library(Matrix)
library(Seurat)
library(ggplot2)

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

# Set seed for reproducibility
set.seed(123)
random_integers <- sample(1:100, 20, replace = FALSE)
proportion_lst <- c(100,90,80,70,60,50,40,30,20,10,5,1)
DoubletsPerCluster = 100
###
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
  
  ggsave(filename = paste0(outdir,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".jpeg"), p1+p2+p3, height = 5, width = 15, dpi = 300)
  
  saveRDS(synwithcontam.recluster,paste0(outdir,"Syn_5_Cell_Type_CP_",contamper,"_DR_",DoubletsPerCluster,".rds"))
  
}
