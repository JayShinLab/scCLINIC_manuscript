###Fig4 is adipose, Fig5 is pbmc covid
#Data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129363, stored in ~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatureMetabolism
#Reference https://www.nature.com/articles/s42255-019-0152-6 Nature Metabolism
###Fig4B_LoadDataCluster18.R
library(ggrepel)
library(dplyr)
library(scCustomize)
library(enrichR)
library(ggplot2)
library(pheatmap)
pacman::p_load(fgsea,tidyverse,pheatmap,openxlsx)
Hgmt<-gmtPathways("~/DEAlgoManuscript/Manuscript_Figures/h.all.v2023.2.Hs.symbols.gmt")#url("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.v2023.2.Hs.symbols.gmt"))
GOgmt<-gmtPathways("~/DEAlgoManuscript/Manuscript_Figures/c5.go.v2023.2.Hs.symbols.gmt")#(#url("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"))

saveRDS(GOgmt,"~/DEAlgoManuscript/Manuscript_Figures/GOgmt.rds")
saveRDS(Hgmt,"~/DEAlgoManuscript/Manuscript_Figures/Hgmt.rds")

library(Seurat)
adinat <- Read10X("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatureMetabolism") 
adinatobj <- CreateSeuratObject(adinat)

adinatobj$SampleIDX <- sapply(strsplit(rownames(adinatobj@meta.data), "-"), "[[", 2)
adinatobj$nCount_RNA <- colSums(x = adinatobj, slot = "counts")  # nCount_RNA
adinatobj$nFeatureRNA <- colSums(x = GetAssayData(object = adinatobj, slot = "counts") > 0)  # nFeatureRNA
adinatobj[["percent.mt"]] <- PercentageFeatureSet(adinatobj, pattern = "^MT-")


adinatobjQC <- subset(adinatobj, 
                      nCount_RNA >= 200 & nFeatureRNA >= 200 & percent.mt <= 20)

#DF
library(DoubletFinder)

subobj <- adinatobjQC

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu_kidney <- subobj
seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

sweep.vector <- DoubletFinder::paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
sweep.table <- DoubletFinder::summarizeSweep(sweep.vector, GT = FALSE)
bcmvn <- DoubletFinder::find.pK(sweep.table)

saveRDS(bcmvn, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/adinatobjQC_bcmvn.rds"))

pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
pK <- as.numeric(levels(pK))[pK]
seu_kidney <- DoubletFinder::doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = pK,
                                              nExp = 0.1, reuse.pANN = FALSE, sct = FALSE)

metadata_names <- colnames(seu_kidney@meta.data)
matching_column <- grep("pANN", metadata_names, value = TRUE)
score <- seu_kidney@meta.data[, matching_column]

saveRDS(score, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/adinatobjQC_DFscore.rds"))


adinatobjQC$DFScore <- score
saveRDS(adinatobjQC, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/adinatobjQC.rds"))


#DEAlgo

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

Name <- "AdiposeNatMet"
Input <- "~/DEAlgoManuscript/Manuscript_Figures/Fig4/adinatobjQC.rds"
Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig4/"
resol=0.8
overlapRatioList=c(0.1,0.25,0.5,0.75,0.9)
OverlapRatio=0.5
ISThreshold=0
gene_n=150

obj <- STEP1A_GlobalMarkers(Input,Output,Name,resol)
obj <- STEP1B_MergingCluster(obj,Output,Name,resol,overlapRatioList,gene_n)
obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold)
STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n)
PlotContaminationPattern(obj,Output,Name,OverlapRatio)

###Fig4B_ConditionGSEA.R
singletlevel = 6
dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Output_Overlap_Ratio_0.5/scCLINICResult.rds"))
doublet_afqc <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/adinatobjQC.rds"))

metatable_unique<- read.table("~/DEAlgoManuscript/Manuscript_Figures/Fig4/ClinicalMetaData.csv",sep = ",")

dealgoseuobj$SampleIDX <- sapply(strsplit(rownames(dealgoseuobj@meta.data), "-"), "[[", 2)
indices <- match(dealgoseuobj@meta.data$SampleIDX, metatable_unique$x.SampleIDX)

# Assign the values from metatable_unique to the corresponding rows in dealgoseuobj@meta.data
dealgoseuobj@meta.data$ncells <- metatable_unique$ncells[indices]
dealgoseuobj@meta.data$Condition <- metatable_unique$y.Condition[indices]
dealgoseuobj@meta.data$Tissue <- metatable_unique$y.Tissue[indices]
dealgoseuobj@meta.data$SampleName <- metatable_unique$y.SampleName[indices]


#7.5% of 26540 is 1990
top_indices <- order(doublet_afqc@meta.data$DFScore, decreasing = TRUE)[1190]
cutoff <- doublet_afqc@meta.data$DFScore[top_indices]

dealgoseuobj$DFafqc_pANN <- doublet_afqc$DFScore
dealgoseuobj$DFafqc <- NA
dealgoseuobj$DFafqc <- ifelse(doublet_afqc$DFScore < cutoff, "Singlet", "Doublet")

qcsclst <- sort(unique(dealgoseuobj@meta.data[,"Overlap_Ratio_0.5"]))

OriginalALL <- dealgoseuobj

dealgoseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(dealgoseuobj$scCLINIC_Level))[dealgoseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$scCLINIC_contamID <- ifelse(dealgoseuobj$scCLINIC_artifact == "Artifact", dealgoseuobj$scCLINIC_ClusterID, NA)

for (i in c("M3","M4","M5","M6","M7","M8")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
  recluster$scCLINIC_ClusterID <- dealgoseuobj$scCLINIC_ClusterID
  recluster$scCLINICScore <- dealgoseuobj$scCLINICScore
  recluster$scCLINIC_Level <- dealgoseuobj$scCLINIC_Level
  recluster$scCLINIC_artifact <- dealgoseuobj$scCLINIC_artifact
  recluster$scCLINIC_contamID <- dealgoseuobj$scCLINIC_contamID 
  
  p7 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
  p8 <- FeaturePlot(recluster, reduction = "umap",features = "scCLINICScore")
  p4 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,  sizes.highlight = 0.1)
  p5 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_artifact",raster=FALSE,  sizes.highlight = 0.1)
  
  p3 <- DimPlot(recluster, group.by = "scCLINIC_contamID")
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4_cluster_",i,"_ID.png"), p7+p5+p8+p4+p3, height = 10, width = 17, dpi = 300)
  
}


DEAlgoSinglet <- subset(dealgoseuobj,subset = scCLINIC_artifact == "Singlet")

DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")

OriginalALL$ConditionTissue <- paste0(OriginalALL$Condition, OriginalALL$Tissue)
p1 <- DimPlot(OriginalALL,group.by = "ConditionTissue",cols = c(
  "#7c3b36",
  "#71b84e",
  "#cf4ba6",
  "#53a0b9",
  "#4e336c",
  "#b488bb"
))


ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/ConditionTissue.png"), p1, height = 8, width = 10, dpi = 300)

p1 <- DimPlot(OriginalALL,group.by = "Overlap_Ratio_0.5",cols = c(
  "#597b4a",
  "#71b84e",
  "#7c4ccb",
  "#c4944a",
  "#cf4ba6",
  "#53a0b9",
  "#d55046",
  "#4e336c",
  "#b488bb"
))


ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/UMAPFig4B.png"), p1, height = 8, width = 10, dpi = 300)


#Diabetic + SAT
DiaSATOri <- subset(OriginalALL, subset = Condition == "Diabetic" & Tissue == "SAT")
DiaSATDE <- subset(DEAlgoSinglet, subset = Condition == "Diabetic" & Tissue == "SAT")
DiaSATDF <- subset(DFSinglet, subset = Condition == "Diabetic" & Tissue == "SAT")
#NonDiabetic + SAT
NonSATOri <- subset(OriginalALL, subset = Condition == "NonDiabetic" & Tissue == "SAT")
NonSATDE <- subset(DEAlgoSinglet, subset = Condition == "NonDiabetic" & Tissue == "SAT")
NonSATDF <- subset(DFSinglet, subset = Condition == "NonDiabetic" & Tissue == "SAT")
#Diabetic + VAT
DiaVATOri <- subset(OriginalALL, subset = Condition == "Diabetic" & Tissue == "VAT")
DiaVATDE <- subset(DEAlgoSinglet, subset = Condition == "Diabetic" & Tissue == "VAT")
DiaVATDF <- subset(DFSinglet, subset = Condition == "Diabetic" & Tissue == "VAT")
#NonDiabetic + VAT
NonVATOri <- subset(OriginalALL, subset = Condition == "NonDiabetic" & Tissue == "VAT")
NonVATDE <- subset(DEAlgoSinglet, subset = Condition == "NonDiabetic" & Tissue == "VAT")
NonVATDF <- subset(DFSinglet, subset = Condition == "NonDiabetic" & Tissue == "VAT")

p1 <- DimPlot(DiaSATOri,group.by = "Overlap_Ratio_0.5")
p2 <- DimPlot(DiaSATDE,group.by = "Overlap_Ratio_0.5")
p3 <- DimPlot(DiaSATDF,group.by = "Overlap_Ratio_0.5")
p4 <- DimPlot(NonSATOri,group.by = "Overlap_Ratio_0.5")
p5 <- DimPlot(NonSATDE,group.by = "Overlap_Ratio_0.5")
p6 <- DimPlot(NonSATDF,group.by = "Overlap_Ratio_0.5")
p7 <- DimPlot(DiaVATOri,group.by = "Overlap_Ratio_0.5")
p8 <- DimPlot(DiaVATDE,group.by = "Overlap_Ratio_0.5")
p9 <- DimPlot(DiaVATDF,group.by = "Overlap_Ratio_0.5")
p10 <- DimPlot(NonVATOri,group.by = "Overlap_Ratio_0.5")
p11 <- DimPlot(NonVATDE,group.by = "Overlap_Ratio_0.5")
p12 <- DimPlot(NonVATDF,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaSAT.png"), p1+p2+p3, height = 5, width = 15, dpi = 300)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/NonSAT.png"), p4+p5+p6, height = 5, width = 15, dpi = 300)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaVAT.png"), p7+p8+p9, height = 5, width = 15, dpi = 300)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/NonVAT.png"), p10+p11+p12, height = 5, width = 15, dpi = 300)


#Azimuth
obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig4/step1d.rds")
obj_azimuth <- Azimuth::RunAzimuth(obj,reference = "adiposeref")

for (i in c("M1","M2","M3","M4","M5","M6","M7","M8")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
  recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
  recluster$DFafqc <- dealgoseuobj$DFafqc
  recluster$scCLINIC_ClusterID <- dealgoseuobj$scCLINIC_ClusterID
  recluster$scCLINICScore <- dealgoseuobj$scCLINICScore
  recluster$scCLINIC_Level <- dealgoseuobj$scCLINIC_Level
  recluster$scCLINIC_artifact <- dealgoseuobj$scCLINIC_artifact
  
  recluster$azimuthl1 <- obj_azimuth$predicted.celltype.l1
  recluster$azimuthl2 <- obj_azimuth$predicted.celltype.l2
  
  p7 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
  p8 <- FeaturePlot(recluster, reduction = "umap",features = "scCLINICScore")
  p4 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,  sizes.highlight = 0.1)
  p5 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_artifact",raster=FALSE,  sizes.highlight = 0.1)
  p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
  p0 <- DimPlot(recluster, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1, cols = c("red","grey"))
  
  p2 <- DimPlot(recluster, group.by = "azimuthl1", cols = c(
    "#7c3b36",
    "#71b84e",
    "#7c4ccb",
    "#c4944a",
    "#cf4ba6",
    "#597b4a",
    "#d55046",
    "#53a0b9",
    "#4e336c",
    "#b488bb",
    "#8d8d0f",  # Unique color 1
    "#3b8d89",  # Unique color 2
    "#e51c23",  # Unique color 3 (changed)
    "#1e88e5",  # Unique color 4 (changed)
    "#ff9800"   # Unique color 5 (changed)  # Unique color 5
  ))
  
  p3 <- DimPlot(recluster, group.by = "azimuthl2")
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4_cluster_",i,"_rawscore_Validate.png"), p7+p5+p8+p4+p0+p1+p2+p3, height = 15, width = 17, dpi = 300)
  
}

dealgoseuobj$azimuthl1 <- obj_azimuth$predicted.celltype.l1
dealgoseuobj$azimuthl2 <- obj_azimuth$predicted.celltype.l2

p0 <- DimPlot(dealgoseuobj,group.by = "azimuthl1", reduction = "umap",raster=FALSE, cols = c(
  "#7c3b36",
  "#71b84e",
  "#7c4ccb",
  "#c4944a",
  "#cf4ba6",
  "#597b4a",
  "#d55046",
  "#53a0b9",
  "#4e336c",
  "#b488bb",
  "#8d8d0f",  # Unique color 1
  "#3b8d89",  # Unique color 2
  "#e51c23",  # Unique color 3 (changed)
  "#1e88e5",  # Unique color 4 (changed)
  "#ff9800"   # Unique color 5 (changed)  # Unique color 5
))
p1 <- DimPlot(dealgoseuobj,group.by = "azimuthl2", reduction = "umap",raster=FALSE)
p3 <- FeaturePlot(dealgoseuobj,features = "scCLINICScore", reduction = "umap",raster=FALSE, label = F) #& scale_color_gradient(limits = c(0, 1))
p5 <- FeaturePlot(dealgoseuobj, reduction = "umap",features = "DFafqc_pANN")
p7 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1)
p11 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,  sizes.highlight = 0.1)
p12 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "scCLINIC_artifact",raster=FALSE,  sizes.highlight = 0.1)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4_Azimuth.png"), p12 + p3 + p11 + p7 + p5 +p0 + p1, height = 25, width = 23, dpi = 300)

saveRDS(obj_azimuth,"~/DEAlgoManuscript/Manuscript_Figures/Fig4/Fig4_Azimuth.rds")

#Compare DiaSAT vs NonSAT
#Batch Script start Rscript ~/Script/script/8_Fig4_published_batch1.R, "M0","M2","M6","M10"