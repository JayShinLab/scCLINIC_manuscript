cell_type = commandArgs(trailingOnly=TRUE) #"M0","M2","M6","M10"
# test if there is at least one argument: if not, return an error
if (length(cell_type)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
###Fig4B_ConditionGSEA.R
library(ggrepel)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(pheatmap)
pacman::p_load(fgsea,tidyverse,pheatmap,openxlsx)
library(Seurat)

GOgmt <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/GOgmt.rds")
Hgmt <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Hgmt.rds")

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

DEAlgoSinglet <- subset(dealgoseuobj,subset = scCLINIC_artifact == "Singlet")


DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")

OriginalALL$ConditionTissue <- paste0(OriginalALL$Condition, OriginalALL$Tissue)

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


#Compare DiaSAT vs NonSAT
#Batch Script start Fig4_published_batch1.R
subDiaSATOri <- subset(DiaSATOri,subset = Overlap_Ratio_0.5 == cell_type)
subDiaSATDE <- subset(DiaSATDE,subset = Overlap_Ratio_0.5 == cell_type)
subDiaSATDF <- subset(DiaSATDF,subset = Overlap_Ratio_0.5 == cell_type)
subNonSATOri <- subset(NonSATOri,subset = Overlap_Ratio_0.5 == cell_type)
subNonSATDE <- subset(NonSATDE,subset = Overlap_Ratio_0.5 == cell_type)
subNonSATDF <- subset(NonSATDF,subset = Overlap_Ratio_0.5 == cell_type)
#Ori
mergeseurat <- Merge_Seurat_List(c(subDiaSATOri,subNonSATOri),c("DiaSATOri","NonSATOri"))
mergeseurat <- DietSeurat(mergeseurat, counts = TRUE, data = TRUE, scale.data = FALSE)
mergeseurat$Group <- sub("_.*", "", rownames(mergeseurat@meta.data))
genes.percent.expression <- rowMeans(mergeseurat[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(object = mergeseurat) <- "Group"
markers<-FindMarkers(mergeseurat,ident.1 ="DiaSATOri", ident.2 = "NonSATOri", 
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)  ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATOrivsNonSATOri_",cell_type,".rds"))

#DE
mergeseurat <- Merge_Seurat_List(c(subDiaSATDE,subNonSATDE),c("DiaSATDE","NonSATDE"))
mergeseurat <- DietSeurat(mergeseurat, counts = TRUE, data = TRUE, scale.data = FALSE)
mergeseurat$Group <- sub("_.*", "", rownames(mergeseurat@meta.data))
genes.percent.expression <- rowMeans(mergeseurat[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(object = mergeseurat) <- "Group"
markers<-FindMarkers(mergeseurat,ident.1 ="DiaSATDE", ident.2 = "NonSATDE", 
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)   ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATDEvsNonSATDE_",cell_type,".rds"))
#DF
mergeseurat <- Merge_Seurat_List(c(subDiaSATDF,subNonSATDF),c("DiaSATDF","NonSATDF"))
mergeseurat <- DietSeurat(mergeseurat, counts = TRUE, data = TRUE, scale.data = FALSE)
mergeseurat$Group <- sub("_.*", "", rownames(mergeseurat@meta.data))
genes.percent.expression <- rowMeans(mergeseurat[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(object = mergeseurat) <- "Group"
markers<-FindMarkers(mergeseurat,ident.1 ="DiaSATDF", ident.2 = "NonSATDF", 
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)   ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATDFvsNonSATDF_",cell_type,".rds"))
 



dbused <- "Hgmt"
#CELL TYPE SPECIFIC GSEA
markersDEAlgo <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATDEvsNonSATDE_",cell_type,".rds"))
markersDF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATDFvsNonSATDF_",cell_type,".rds"))
markersOri <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATOrivsNonSATOri_",cell_type,".rds"))

#markers$gene <- rownames(markers)
#markers <- markers[ markers$avg_log2FC >0,]
markersDEAlgo$stat=markersDEAlgo$avg_log2FC
markersDF$stat=markersDF$avg_log2FC
markersOri$stat=markersOri$avg_log2FC
#replace Inf with max finite value
markersDEAlgo$stat<-replace(markersDEAlgo$stat,markersDEAlgo$stat>max(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)]),max(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)])*10)
markersDEAlgo$stat<-replace(markersDEAlgo$stat,markersDEAlgo$stat<min(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)]),min(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)])*10)

markersDF$stat<-replace(markersDF$stat,markersDF$stat>max(markersDF$stat[is.finite(markersDF$stat)]),max(markersDF$stat[is.finite(markersDF$stat)])*10)
markersDF$stat<-replace(markersDF$stat,markersDF$stat<min(markersDF$stat[is.finite(markersDF$stat)]),min(markersDF$stat[is.finite(markersDF$stat)])*10)

markersOri$stat<-replace(markersOri$stat,markersOri$stat>max(markersOri$stat[is.finite(markersOri$stat)]),max(markersOri$stat[is.finite(markersOri$stat)])*10)
markersOri$stat<-replace(markersOri$stat,markersOri$stat<min(markersOri$stat[is.finite(markersOri$stat)]),min(markersOri$stat[is.finite(markersOri$stat)])*10)


Hlst <- list()
markerssub <- markersOri
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(Hgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaSATvsNonSAT_Ori_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(OriResult=H,Hlst)

markerssub <- markersDF
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(Hgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaSATvsNonSAT_DFQC_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DFResult=H,Hlst)

markerssub <- markersDEAlgo
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(Hgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaSATvsNonSAT_DEAlgo_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DEAlgoResult=H,Hlst)

saveRDS(Hlst,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","DiaSATvsNonSAT_Result_",dbused,"_",cell_type,".rds"))

dbused <- "GOgmt"
  
#CELL TYPE SPECIFIC GSEA
markersDEAlgo <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATDEvsNonSATDE_",cell_type,".rds"))
markersDF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATDFvsNonSATDF_",cell_type,".rds"))
markersOri <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/GSEA_Marker_DiaSATOrivsNonSATOri_",cell_type,".rds"))

#markers$gene <- rownames(markers)
#markers <- markers[ markers$avg_log2FC >0,]
markersDEAlgo$stat=markersDEAlgo$avg_log2FC
markersDF$stat=markersDF$avg_log2FC
markersOri$stat=markersOri$avg_log2FC
#replace Inf with max finite value
markersDEAlgo$stat<-replace(markersDEAlgo$stat,markersDEAlgo$stat>max(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)]),max(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)])*10)
markersDEAlgo$stat<-replace(markersDEAlgo$stat,markersDEAlgo$stat<min(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)]),min(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)])*10)

markersDF$stat<-replace(markersDF$stat,markersDF$stat>max(markersDF$stat[is.finite(markersDF$stat)]),max(markersDF$stat[is.finite(markersDF$stat)])*10)
markersDF$stat<-replace(markersDF$stat,markersDF$stat<min(markersDF$stat[is.finite(markersDF$stat)]),min(markersDF$stat[is.finite(markersDF$stat)])*10)

markersOri$stat<-replace(markersOri$stat,markersOri$stat>max(markersOri$stat[is.finite(markersOri$stat)]),max(markersOri$stat[is.finite(markersOri$stat)])*10)
markersOri$stat<-replace(markersOri$stat,markersOri$stat<min(markersOri$stat[is.finite(markersOri$stat)]),min(markersOri$stat[is.finite(markersOri$stat)])*10)


Hlst <- list()
markerssub <- markersOri
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaSATvsNonSAT_Ori_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(OriResult=H,Hlst)

markerssub <- markersDF
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaSATvsNonSAT_DFQC_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DFResult=H,Hlst)

markerssub <- markersDEAlgo
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/DiaSATvsNonSAT_DEAlgo_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DEAlgoResult=H,Hlst)

saveRDS(Hlst,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","DiaSATvsNonSAT_Result_",dbused,"_",cell_type,".rds"))

