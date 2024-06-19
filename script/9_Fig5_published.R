#Fig 5
#data downloaded from https://www.covid19cellatlas.org/index.patient.html Blish Lab Peripheral Blood Mononuclear Cells (PBMCs) Nature Medicine, stored in ~/DEAlgoManuscript/Manuscript_Figures/Fig5
#Reference https://www.nature.com/articles/s41591-020-0944-y


covid19natmed <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/blish_covid19.rds")
covid19natmed.updated = UpdateSeuratObject(object = covid19natmed)
saveRDS(covid19natmed.updated,"~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemed.rds")

####Fig4C_covid19.R
library(Seurat)
adinatobj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemed.rds")
adinat <- adinatobj@assays$RNA$counts

#DF
library(DoubletFinder)

subobj <- CreateSeuratObject(adinat)

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

saveRDS(bcmvn, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC_bcmvn.rds"))

pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
pK <- as.numeric(levels(pK))[pK]
seu_kidney <- DoubletFinder::doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = pK,
                                              nExp = 0.1, reuse.pANN = FALSE, sct = FALSE)

metadata_names <- colnames(seu_kidney@meta.data)
matching_column <- grep("pANN", metadata_names, value = TRUE)
score <- seu_kidney@meta.data[, matching_column]

saveRDS(score, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC_DFscore.rds"))


adinatobj$DFScore <- score
saveRDS(adinatobj, paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC.rds"))


#DEAlgo
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

Name <- "PBMCCOVID_NatureMed"
Input <- "~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC.rds"
Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig5/"
filteredmatrix=NA
rawmatrix=NA
resol=0.8
overlapRatioList=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
OverlapRatio=0.5
ISThreshold=0
Cutoff=0
gene_n=150

obj <- STEP1A_GlobalMarkers(Input,Output,Name,resol)

obj <- STEP1B_MergingCluster(obj,Output,Name,resol,overlapRatioList,gene_n)

obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n)

obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold)

saveRDS(obj,"~/DEAlgoManuscript/Manuscript_Figures/Fig5/step1d.rds")
STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/step1d.rds")


obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,Cutoff,filteredmatrix,rawmatrix)

PlotContaminationPattern(obj,Output,Name,OverlapRatio)



###Start Analysis
library(ggrepel)
library(dplyr)
library(scCustomize)
library(enrichR)
library(ggplot2)
library(pheatmap)
pacman::p_load(fgsea,tidyverse,pheatmap,openxlsx)

#read doublet
GOgmt <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/GOgmt.rds")
Hgmt <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Hgmt.rds")
singletlevel = 4


dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Output_Overlap_Ratio_0.5/","DEAlgoResult.rds"))
doublet_afqc <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC.rds"))

metatable_unique<- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemed.rds"))

dealgoseuobj@meta.data <- cbind(metatable_unique@meta.data,dealgoseuobj@meta.data)

# # Add prefix to column names of metatable_unique
# colnames(metatable_unique@meta.data) <- paste0("Ori_", colnames(metatable_unique@meta.data))
# 
# # Add prefix to column names of dealgoseuobj
# colnames(dealgoseuobj@meta.data) <- paste0("DEAlgo_", colnames(dealgoseuobj@meta.data))
# 
# # Merge the two data frames
# dealgoseuobj@meta.data <- cbind(metatable_unique@meta.data, dealgoseuobj@meta.data)

#7.5% of 44721 is 3354
top_indices <- order(doublet_afqc@meta.data$DFScore, decreasing = TRUE)[3354]
cutoff <- doublet_afqc@meta.data$DFScore[top_indices]

dealgoseuobj$DFafqc_pANN <- doublet_afqc$DFScore
dealgoseuobj$DFafqc <- NA
dealgoseuobj$DFafqc <- ifelse(doublet_afqc$DFScore < cutoff, "Singlet", "Doublet")

qcsclst <- sort(unique(dealgoseuobj@meta.data[,"Overlap_Ratio_0.5"]))


#Healthy vs Control
OriginalALL <- dealgoseuobj

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
DEAlgoSinglet <- subset(dealgoseuobj,subset = DEAlgo_Contaminated == "Singlet")

DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")



#Diabetic + SAT
DiaSATOri <- subset(OriginalALL, subset = Status == "COVID")
DiaSATDE <- subset(DEAlgoSinglet, subset = Status == "COVID")
DiaSATDF <- subset(DFSinglet, subset = Status == "COVID")
#NonDiabetic + SAT
NonSATOri <- subset(OriginalALL, subset = Status == "Healthy")
NonSATDE <- subset(DEAlgoSinglet, subset = Status == "Healthy")
NonSATDF <- subset(DFSinglet, subset = Status == "Healthy")

p1 <- DimPlot(DiaSATOri,group.by = "Overlap_Ratio_0.5",cols = c(
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
p2 <- DimPlot(DiaSATDE,group.by = "Overlap_Ratio_0.5",cols = c(
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
p3 <- DimPlot(DiaSATDF,group.by = "Overlap_Ratio_0.5",cols = c(
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
p4 <- DimPlot(NonSATOri,group.by = "Overlap_Ratio_0.5",cols = c(
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
p5 <- DimPlot(NonSATDE,group.by = "Overlap_Ratio_0.5",cols = c(
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
p6 <- DimPlot(NonSATDF,group.by = "Overlap_Ratio_0.5",cols = c(
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

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID.jpeg"), p1+p2+p3, height = 5, width = 17, dpi = 300)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/Healthy.jpeg"), p4+p5+p6, height = 5, width = 17, dpi = 300)

#Azimuth
obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/step1d.rds")
obj_azimuth <- Azimuth::RunAzimuth(obj,reference = "pbmcref")

for (i in c("M6","M16","M20")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
  recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
  recluster$DFafqc <- dealgoseuobj$DFafqc
  recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
  recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
  recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
  recluster$DEAlgo_Contaminated <- dealgoseuobj$DEAlgo_Contaminated
  
  recluster$azimuthl1 <- obj_azimuth$predicted.celltype.l1
  recluster$azimuthl2 <- obj_azimuth$predicted.celltype.l2
  
  p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
  p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
  p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE,  sizes.highlight = 0.1)
  p5 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE,  sizes.highlight = 0.1)
  p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
  p0 <- DimPlot(recluster, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1)
  
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
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig5_cluster_",i,"_rawscore_Validate.png"), p7+p5+p8+p4+p0+p1+p2+p3, height = 10, width = 17, dpi = 300)
  
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
p3 <- FeaturePlot(dealgoseuobj,features = "DiffDP1_sum", reduction = "umap",raster=FALSE, label = F) #& scale_color_gradient(limits = c(0, 1))
p5 <- FeaturePlot(dealgoseuobj, reduction = "umap",features = "DFafqc_pANN")
p7 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1)
p11 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE,  sizes.highlight = 0.1)
p12 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE,  sizes.highlight = 0.1)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig5_Azimuth.png"), p12 + p3 + p11 + p7 + p5 +p0 + p1, height = 20, width = 23, dpi = 300)

saveRDS(obj_azimuth,"~/DEAlgoManuscript/Manuscript_Figures/Fig5/Fig5_Azimuth.rds")

#Marker COVID VS HEALTHY
#Fig5_published_batch1.R, "M6","M16","M20"