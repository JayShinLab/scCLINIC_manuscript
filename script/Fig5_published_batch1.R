cell_type = commandArgs(trailingOnly=TRUE) #"M2","M16","M20" Monocyte, Platelet, Neutrophil
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
singletlevel = 4


dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Output_Overlap_Ratio_0.5/","DEAlgoResult.rds"))
doublet_afqc <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC.rds"))

metatable_unique<- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemed.rds"))

dealgoseuobj@meta.data <- cbind(metatable_unique@meta.data,dealgoseuobj@meta.data)

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

#All COVID
obj <- DiaSATOri
genes.percent.expression <- rowMeans(obj[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(obj)="Overlap_Ratio_0.5"
markers<-FindMarkers(obj,ident.1 =cell_type,
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)  ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_Ori_COVID_",cell_type,".rds"))

#DE COVID
obj <- DiaSATDE
genes.percent.expression <- rowMeans(obj[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(obj)="Overlap_Ratio_0.5"
markers<-FindMarkers(obj,ident.1 =cell_type,
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)  ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DE_COVID_",cell_type,".rds"))

#DF COVID
obj <- DiaSATDF
genes.percent.expression <- rowMeans(obj[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(obj)="Overlap_Ratio_0.5"
markers<-FindMarkers(obj,ident.1 =cell_type,
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)  ##Output and rank all filtered genes
saveRDS(markers,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DF_COVID_",cell_type,".rds"))


#
dbused <- "Hgmt"
#CELL TYPE SPECIFIC GSEA
markersDEAlgo <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DE_COVID_",cell_type,".rds"))
markersDF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DF_COVID_",cell_type,".rds"))
markersOri <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_Ori_COVID_",cell_type,".rds"))

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
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_Ori_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(OriResult=H,Hlst)

markerssub <- markersDF
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(Hgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_DFQC_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DFResult=H,Hlst)

markerssub <- markersDEAlgo
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(Hgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_DEAlgo_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DEAlgoResult=H,Hlst)

saveRDS(Hlst,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_Result_",dbused,"_",cell_type,".rds"))

dbused <- "GOgmt"

#CELL TYPE SPECIFIC GSEA
markersDEAlgo <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DE_COVID_",cell_type,".rds"))
markersDF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DF_COVID_",cell_type,".rds"))
markersOri <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_Ori_COVID_",cell_type,".rds"))

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
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_Ori_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(OriResult=H,Hlst)

markerssub <- markersDF
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_DFQC_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DFResult=H,Hlst)

markerssub <- markersDEAlgo
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)
write.table(as.matrix(H),file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_DEAlgo_FGSEA_",dbused,"_",cell_type,".csv"),sep = ",")

Hlst <- c(DEAlgoResult=H,Hlst)

saveRDS(Hlst,paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID_Result_",dbused,"_",cell_type,".rds"))

