library(ggrepel)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(pheatmap)
pacman::p_load(fgsea,tidyverse,pheatmap,openxlsx)
library(Seurat)
#read doublet
GOgmt <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/GOgmt.rds")
Hgmt <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Hgmt.rds")
singletlevel = 4
dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","DEAlgoResult.rds"))
doublet_afqc <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/seu_kidney_codeocean_updated_withDFScore.rds"))

dealgoseuobj$DFafqc_pANN <- doublet_afqc@meta.data[,  grep("pANN", colnames(doublet_afqc@meta.data), value = TRUE)]

dealgoseuobj$DFafqc <- doublet_afqc@meta.data[,  grep("DF.classifications", colnames(doublet_afqc@meta.data), value = TRUE)]#doublet_afqc$DF.classifications_0.25_0.09_913

#Healthy vs Control
OriginalALL <- dealgoseuobj

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
DEAlgoSinglet <- subset(dealgoseuobj,subset = DEAlgo_Contaminated == "Singlet")

DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")

DPT <- subset(dealgoseuobj,subset = CT.Park == "Stroma")
DFSinglet <- subset(DPT,subset = DFafqc == "Singlet")
DEAlgoSinglet <- subset(DPT,subset = DEAlgo_Contaminated == "Singlet")
All <- DPT

mergeseurat <- Merge_Seurat_List(c(All,DEAlgoSinglet,DFSinglet),c("Ori.","DEAlgo.","DF."))

mergeseurat <- DietSeurat(mergeseurat, counts = TRUE, data = TRUE, scale.data = FALSE)

mergeseurat$Group <- sub("\\..*", "", rownames(mergeseurat@meta.data))

genes.percent.expression <- rowMeans(mergeseurat[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])

obj <- mergeseurat
Idents(object = obj) <- "Group"
markers<-FindMarkers(obj,ident.1 ="DEAlgo", ident.2 = "Ori",
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)  ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/GSEA_Marker_DEAlgo_Ori_Stroma.rds"))

markers<-FindMarkers(obj,ident.1 ="DF", ident.2 = "Ori",
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)  ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/GSEA_Marker_DF_Ori_Stroma.rds"))
  
#CELL TYPE SPECIFIC GSEA
markersDEAlgo <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/GSEA_Marker_DEAlgo_Ori_Stroma.rds"))
markersDF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/GSEA_Marker_DF_Ori_Stroma.rds"))
#markers$gene <- rownames(markers)
#markers <- markers[ markers$avg_log2FC >0,]
markersDEAlgo$stat=markersDEAlgo$avg_log2FC
markersDF$stat=markersDF$avg_log2FC
#replace Inf with max finite value
markersDEAlgo$stat<-replace(markersDEAlgo$stat,markersDEAlgo$stat>max(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)]),max(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)])*10)
markersDEAlgo$stat<-replace(markersDEAlgo$stat,markersDEAlgo$stat<min(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)]),min(markersDEAlgo$stat[is.finite(markersDEAlgo$stat)])*10)

markersDF$stat<-replace(markersDF$stat,markersDF$stat>max(markersDF$stat[is.finite(markersDF$stat)]),max(markersDF$stat[is.finite(markersDF$stat)])*10)
markersDF$stat<-replace(markersDF$stat,markersDF$stat<min(markersDF$stat[is.finite(markersDF$stat)]),min(markersDF$stat[is.finite(markersDF$stat)])*10)



dbused <- "Hgmt"
Hlst <- list()

markerssub <- markersDF
rank<-setNames(markerssub$stat,toupper(rownames(markerssub)))
rank <- rank[names(rank) != "PISD"]
H<-fgsea(Hgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
Hlst <- c(DFResult=H,Hlst)

markerssub <- markersDEAlgo
rank<-setNames(markerssub$stat,toupper(rownames(markerssub)))
rank <- rank[names(rank) != "PISD"]
H<-fgsea(Hgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
Hlst <- c(DEAlgoResult=H,Hlst)

saveRDS(Hlst,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_","Stroma_X_Ori",".rds"))
