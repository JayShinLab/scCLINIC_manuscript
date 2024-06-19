seedidx = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(seedidx)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
seedidx <- as.integer(seedidx)

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
  
  # obj <- readRDS(Input)
  # 
  # obj$source <- ifelse(is.na(obj$cell_type), paste("CellType",obj$y),obj$cell_type)
  # 
  # obj$source <- gsub("CellType ","CellType_",obj$source)
  # 
  # obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)
  # 
  # obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold,CELLANNOTATION = TRUE)
  # saveRDS(obj,paste0(Output,Name,"_",contamper,"_step1d.rds"))
  obj <- readRDS(paste0(Output,Name,"_",contamper,"_step1d.rds"))
  
  STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)
  
  obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,Cutoff,filteredmatrix,rawmatrix,CELLANNOTATION = TRUE)
  
  PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)
}
