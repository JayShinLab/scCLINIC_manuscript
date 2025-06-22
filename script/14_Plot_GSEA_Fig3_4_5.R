#PBMC Dataset
#Antigen mispresentation in B cells and neutrophils pathway changes
library(ggrepel)
library(dplyr)
library(scCustomize)
library(enrichR)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(readr)
library(forcats)
#read doublet
dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Output_Overlap_Ratio_0.5/","scCLINICResult.rds"))
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

dbused <- "GOgmt"
cell_type <- "M10" #Neutrophil
overlap <- c("GOBP_REGULATION_OF_RESPONSE_TO_STRESS","GOBP_REGULATION_OF_INTRACELLULAR_SIGNAL_TRANSDUCTION","GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS","GOBP_REGULATION_OF_IMMUNE_RESPONSE","GOBP_POSITIVE_REGULATION_OF_SIGNALING","GOBP_POSITIVE_REGULATION_OF_MULTICELLULAR_ORGANISMAL_PROCESS","GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS")

Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]

DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]

Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]

pathwaydisplay2 <- unique(c(Ori_p,DF_p,DEAlgo_p))

# DF_n <- data[data$Group == "DF" & data$NES < 0 & data$padj < 0.05, ]
# DF_n <- DF_n[head(order(DF_n$padj), n=10), "Gene"]
# 
# DEAlgo_n <- data[data$Group == "DEAlgo" & data$NES < 0 & data$padj < 0.05, ]
# DEAlgo_n <- DEAlgo_n[head(order(DEAlgo_n$padj), n=10), "Gene"]
# 
# Ori_n <- data[data$Group == "Ori" & data$NES < 0 & data$padj < 0.05, ]
# Ori_n <- Ori_n[head(order(Ori_n$padj), n=10), "Gene"]
# 
# pathwaydisplay3 <- unique(c(Ori_n,DF_n,DEAlgo_n))

data <- data[data$Gene %in% unique(pathwaydisplay2),]

#data <- data[data$Gene %in% names(tubuleGO),]

# plot: dot plot
data$pstat <- 1  # Initialize the new column

# Assign values based on conditions
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
padj_break <- levels(data$pstat)
data$Gene_short <- substr(data$Gene, 1, 50)
p1 <- ggplot(data = data, aes(x = Group, y = Gene_short, color = NES, size = pstat)) + 
  geom_point() +
  scale_color_viridis_c(option = "viridis", breaks = NES_break, limits = c(min(NES_break), max(NES_break))) +
  scale_size_manual(values = c(1, 2, 3, 4), breaks = padj_break) +  # Use padj_break as breaks argument
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle(paste0(cell_type))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Neutrophil_GSEA_Result_",dbused,"_",cell_type,"GOBPTOP.png"), p1, height = 10, width = 10, dpi = 300)


data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

p1 <- ggplot(data = data, aes(x = -log10(padj), y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0("Neutrophil"))

p2 <- ggplot(data = data, aes(x = NES, y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  ggtitle(paste0("Neutrophil"))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVIDPBMC_Neutrophil_",dbused,"_",cell_type,"_1.png"), p1, height = 4, width = 7, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVIDPBMC_Neutrophil_",dbused,"_",cell_type,"_2.png"), p2, height = 4, width = 5, dpi = 300)


clusterx <- "M10"
Positive_regulation_of_immune_system_process_plot <- c("IL7R","CCL5","MS4A1","HLA-DPB1","IGHM","IGLC3","RPS19")
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))
cluster2$status <- dealgoseuobj$Status

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot) & scale_color_gradientn(colors = c("grey","red"),limit = c(0,8))
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p5 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters", split.by = "status")+ theme(legend.position = "right") 
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limit = c(0,8))
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_5.png"), p5, height = 10, width = 10, dpi = 300)

clusterx <- "M10"
Positive_regulation_of_immune_system_process_plot <- c("HBB", "HBA1","HBA2")
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))
cluster2$status <- dealgoseuobj$Status

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot) & scale_color_gradientn(colors = c("grey","red"),limit = c(0,8), oob = scales::squish)
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p5 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters", split.by = "status")+ theme(legend.position = "right") 
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limit = c(0,8), oob = scales::squish)
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_5.png"), p5, height = 10, width = 10, dpi = 300)

###B Cells
dbused <- "GOgmt"
cell_type <- "M3"
overlap <- c("GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY","GOBP_B_CELL_ACTIVATION","GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY","GOBP_B_CELL_PROLIFERATION")

Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Healthy_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]

DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]

Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]

pathwaydisplay2 <- unique(c(Ori_p,DF_p,DEAlgo_p))

data <- data[data$Gene %in% unique(pathwaydisplay2),]

# plot: dot plot
data$pstat <- 1  # Initialize the new column

# Assign values based on conditions
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
padj_break <- levels(data$pstat)
data$Gene_short <- substr(data$Gene, 1, 50)
p1 <- ggplot(data = data, aes(x = Group, y = Gene_short, color = NES, size = pstat)) + 
  geom_point() +
  scale_color_viridis_c(option = "viridis", breaks = NES_break, limits = c(min(NES_break), max(NES_break))) +
  scale_size_manual(values = c(1, 2, 3, 4), breaks = padj_break) +  # Use padj_break as breaks argument
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle(paste0(cell_type))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Healthy_BCells_GSEA_Result_",dbused,"_",cell_type,"GOBPTOP.png"), p1, height = 10, width = 10, dpi = 300)

data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

p1 <- ggplot(data = data, aes(x = -log10(padj), y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0("B Cells"))

p2 <- ggplot(data = data, aes(x = NES, y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  ggtitle(paste0("B Cells"))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","HealthyPBMC_B_CELLS_",dbused,"_",cell_type,"_1.png"), p1, height = 2.5, width = 7, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","HealthyPBMC_B_CELLS_",dbused,"_",cell_type,"_2.png"), p2, height = 2.5, width = 5, dpi = 300)


# clusterx <- 5
# tables <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DEAlgoQC_05H_GOgmt_",clusterx,".rds"))
# lstcalc <- list()
# for (pathplot in overlap){
#   lstcalc <- c(setNames(list(unlist(tables[tables$pathway == pathplot,]$leadingEdge)[1:20]), as.character(pathplot)), lstcalc)
# }
# 
# cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Output/DEAlgo_PBMCNatMed_Fig4C_Step2/DEAlgo_G0OR_0.5_recluster/DEAlgo_G0OR_0.5_cluster_",clusterx,".rds"))
# cluster2$status <- dealgoseuobj$Ori_Status
# #cluster2 <- subset(cluster2, subset = status == "Healthy")
# for (sigpath in names(lstcalc)){
#   p1 <- FeaturePlot(cluster2, reduction = "umap",features =lstcalc[[sigpath]])
#   p2 <- VlnPlot(cluster2,features =lstcalc[[sigpath]],group.by = "seurat_clusters")
#   p5 <- VlnPlot(cluster2,features =lstcalc[[sigpath]],group.by = "seurat_clusters", split.by = "status")+ theme(legend.position = "right") 
#   p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =lstcalc[[sigpath]])
#   p4 <- VlnPlot(dealgoseuobj,features =lstcalc[[sigpath]],group.by = "DEAlgo_DEAlgo_G0OR_0.5")
#   
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",sigpath,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",sigpath,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",sigpath,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",sigpath,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",sigpath,"_",clusterx,"_5.png"), p5, height = 10, width = 10, dpi = 300)
# }

clusterx <- "M3"
Positive_regulation_of_immune_system_process_plot <- c("CD247","CD8A","TRAC","CD3D","MS4A1","CD79A","CD22")
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))
cluster2$status <- dealgoseuobj$Status

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limits = c(0, 5))
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p5 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters", split.by = "status")+ theme(legend.position = "right") 
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_5.png"), p5, height = 10, width = 10, dpi = 300)

###Platelet-monocyte aggregation
cell_type <- "M2"
dbused <- "GOgmt"
Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]

DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]

Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]

pathwaydisplay2 <- unique(c(Ori_p,DF_p,DEAlgo_p))

data <- data[data$Gene %in% unique(pathwaydisplay2),]

# plot: dot plot
data$pstat <- 1  # Initialize the new column

# Assign values based on conditions
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
padj_break <- levels(data$pstat)
data$Gene_short <- substr(data$Gene, 1, 50)
p1 <- ggplot(data = data, aes(x = Group, y = Gene_short, color = NES, size = pstat)) + 
  geom_point() +
  scale_color_viridis_c(option = "viridis", breaks = NES_break, limits = c(min(NES_break), max(NES_break))) +
  scale_size_manual(values = c(1, 2, 3, 4), breaks = padj_break) +  # Use padj_break as breaks argument
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle(paste0(cell_type))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Monocyte_GSEA_Result_",dbused,"_",cell_type,"GOBPTOP.png"), p1, height = 10, width = 10, dpi = 300)

cell_type <- "M8"
dbused <- "GOgmt"
Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]

DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]

DEAlgo50_p <- data[data$Group == "DEAlgo+" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo50_p <- DEAlgo50_p[head(order(DEAlgo50_p$padj), n=50), "Gene"]

Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]

pathwaydisplay3 <- unique(c(Ori_p,DF_p,DEAlgo_p,DEAlgo50_p))

overlap <- intersect(pathwaydisplay2, pathwaydisplay3)
print(overlap)


Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)

data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
data <- data[data$Group != "DEAlgo+",]
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

p1 <- ggplot(data = data, aes(x = -log10(padj), y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0("Platelet Cells"))

p2 <- ggplot(data = data, aes(x = NES, y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  ggtitle(paste0("Platelet Cells"))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVIDPBMC_PlateletMonocyte_",dbused,"_",cell_type,"_1.png"), p1, height = 5, width = 7, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVIDPBMC_PlateletMonocyte_",dbused,"_",cell_type,"_2.png"), p2, height = 5, width = 5, dpi = 300)

# clusterx <- "M16"
# tables <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/GSEA_Marker_DEAlgoQC_05C_GOgmt_",clusterx,".rds"))
# lstcalc <- list()
# for (pathplot in overlap){
#   lstcalc <- c(setNames(list(unlist(tables[tables$pathway == pathplot,]$leadingEdge)[1:20]), as.character(pathplot)), lstcalc)
# }
# 
# cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Output/DEAlgo_PBMCNatMed_Fig4C_Step2/DEAlgo_G0OR_0.5_recluster/DEAlgo_G0OR_0.5_cluster_",clusterx,".rds"))
# cluster2$status <- dealgoseuobj$Ori_Status
# #cluster2 <- subset(cluster2, subset = status == "Healthy")
# for (sigpath in names(lstcalc)){
#   p1 <- FeaturePlot(cluster2, reduction = "umap",features =lstcalc[[sigpath]])
#   p2 <- VlnPlot(cluster2,features =lstcalc[[sigpath]],group.by = "seurat_clusters")
#   p5 <- VlnPlot(cluster2,features =lstcalc[[sigpath]],group.by = "seurat_clusters", split.by = "status")+ theme(legend.position = "right") 
#   p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =lstcalc[[sigpath]])
#   p4 <- VlnPlot(dealgoseuobj,features =lstcalc[[sigpath]],group.by = "DEAlgo_DEAlgo_G0OR_0.5")
#   
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","CFig4B_",sigpath,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","CFig4B_",sigpath,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","CFig4B_",sigpath,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","CFig4B_",sigpath,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
#   ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","CFig4B_",sigpath,"_",clusterx,"_5.png"), p5, height = 10, width = 10, dpi = 300)
# }

clusterx <- "M8"
Positive_regulation_of_immune_system_process_plot <- c("FCN1","CD14","AIF1","CYBB","MS4A1","TRAC","CD8A","CD247")
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))
cluster2$status <- dealgoseuobj$Status

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limits = c(0, 5))
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p5 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters", split.by = "status")+ theme(legend.position = "right") 
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limits = c(0, 5))
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_5.png"), p5, height = 10, width = 10, dpi = 300)


for (clusterx in qcsclst){
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))
cluster2$status <- dealgoseuobj$Status

p1 <- DimPlot(cluster2,group.by = "status")
p4 <- DimPlot(dealgoseuobj,group.by = "Status")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_","Status","_",clusterx,"_1.png"), p1, height = 5, width = 5, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig4B_","Status","_",clusterx,"_4.png"), p4, height = 5, width = 5, dpi = 300)

}

###Platelet-monocyte aggregation
cell_type <- "M2"
dbused <- "GOgmt"
Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]

DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]

Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]

pathwaydisplay2 <- unique(c(Ori_p,DF_p,DEAlgo_p))

data <- data[data$Gene %in% unique(pathwaydisplay2),]

# plot: dot plot
data$pstat <- 1  # Initialize the new column

# Assign values based on conditions
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
padj_break <- levels(data$pstat)
data$Gene_short <- substr(data$Gene, 1, 50)
p1 <- ggplot(data = data, aes(x = Group, y = Gene_short, color = NES, size = pstat)) + 
  geom_point() +
  scale_color_viridis_c(option = "viridis", breaks = NES_break, limits = c(min(NES_break), max(NES_break))) +
  scale_size_manual(values = c(1, 2, 3, 4), breaks = padj_break) +  # Use padj_break as breaks argument
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle(paste0(cell_type))

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Monocyte_GSEA_Result_",dbused,"_",cell_type,"GOBPTOP.png"), p1, height = 10, width = 10, dpi = 300)

cell_type <- "M10"
dbused <- "GOgmt"
Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]

DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]

DEAlgo50_p <- data[data$Group == "DEAlgo+" & data$NES > 0 & data$padj < 0.05, ]
DEAlgo50_p <- DEAlgo50_p[head(order(DEAlgo50_p$padj), n=50), "Gene"]

Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]

pathwaydisplay3 <- unique(c(Ori_p,DF_p,DEAlgo_p,DEAlgo50_p))

overlap <- intersect(pathwaydisplay2, pathwaydisplay3)
print(overlap)


Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVID_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)

data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
data <- data[data$Group != "DEAlgo+",]
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

p1 <- ggplot(data = data, aes(x = -log10(padj), y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0("Platelet Cells"))

p2 <- ggplot(data = data, aes(x = NES, y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  ggtitle(paste0("Platelet Cells"))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVIDPBMC_NeutrophilsMonocyte_",dbused,"_",cell_type,"_1.png"), p1, height = 20, width = 7, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","COVIDPBMC_NeutrophilsMonocyte_",dbused,"_",cell_type,"_2.png"), p2, height = 20, width = 5, dpi = 300)

# #Healthy
# ###Platelet-monocyte aggregation
# cell_type1 <- "2"
# Ori <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_Ori_05H_",dbused,"_",cell_type1,".rds"))
# DF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DFQC_05H_",dbused,"_",cell_type1,".rds"))
# DE <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DEAlgoQC_05H_",dbused,"_",cell_type1,".rds"))
# DE50 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DEAlgoQC50_05H_",dbused,"_",cell_type1,".rds"))
# 
# data <- data.frame(
#   Group = c(rep("Ori", length(Ori$padj)),rep("DEAlgo", length(DE$padj)),rep("DEAlgo+", length(DE50$padj)),rep("DF", length(DF$padj))),
#   padj = c(Ori$padj,DE$padj,DE50$padj,DF$padj),
#   #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
#   Gene = c(Ori$pathway,DE$pathway,DE50$pathway,DF$pathway),
#   NES = c(Ori$NES,DE$NES,DE50$NES,DF$NES)
#   
#   #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
#   #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
# )
# data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]
# 
# DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
# DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]
# 
# DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
# DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]
# 
# DEAlgo50_p <- data[data$Group == "DEAlgo+" & data$NES > 0 & data$padj < 0.05, ]
# DEAlgo50_p <- DEAlgo50_p[head(order(DEAlgo50_p$padj), n=50), "Gene"]
# 
# Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
# Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]
# 
# pathwaydisplay2 <- unique(c(Ori_p,DF_p,DEAlgo_p,DEAlgo50_p))
# 
# cell_type <- 17
# Ori <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_Ori_05H_",dbused,"_",cell_type,".rds"))
# DF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DFQC_05H_",dbused,"_",cell_type,".rds"))
# DE <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DEAlgoQC_05H_",dbused,"_",cell_type,".rds"))
# DE50 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DEAlgoQC50_05H_",dbused,"_",cell_type,".rds"))
# 
# data <- data.frame(
#   Group = c(rep("Ori", length(Ori$padj)),rep("DEAlgo", length(DE$padj)),rep("DEAlgo+", length(DE50$padj)),rep("DF", length(DF$padj))),
#   padj = c(Ori$padj,DE$padj,DE50$padj,DF$padj),
#   #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
#   Gene = c(Ori$pathway,DE$pathway,DE50$pathway,DF$pathway),
#   NES = c(Ori$NES,DE$NES,DE50$NES,DF$NES)
#   
#   #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
#   #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
# )
# data <- data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]
# 
# DF_p <- data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
# DF_p <- DF_p[head(order(DF_p$padj), n=50), "Gene"]
# 
# DEAlgo_p <- data[data$Group == "DEAlgo" & data$NES > 0 & data$padj < 0.05, ]
# DEAlgo_p <- DEAlgo_p[head(order(DEAlgo_p$padj), n=50), "Gene"]
# 
# DEAlgo50_p <- data[data$Group == "DEAlgo+" & data$NES > 0 & data$padj < 0.05, ]
# DEAlgo50_p <- DEAlgo50_p[head(order(DEAlgo50_p$padj), n=50), "Gene"]
# 
# Ori_p <- data[data$Group == "Ori" & data$NES > 0 & data$padj < 0.05, ]
# Ori_p <- Ori_p[head(order(Ori_p$padj), n=50), "Gene"]
# 
# pathwaydisplay3 <- unique(c(Ori_p,DF_p,DEAlgo_p,DEAlgo50_p))
# 
# overlap <- intersect(pathwaydisplay2, pathwaydisplay3)
# print(overlap)
# 
# Ori <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_Ori_05H_",dbused,"_",cell_type,".rds"))
# DF <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DFQC_05H_",dbused,"_",cell_type,".rds"))
# DE <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DEAlgoQC_05H_",dbused,"_",cell_type,".rds"))
# DE50 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","GSEA_Marker_DEAlgoQC50_05H_",dbused,"_",cell_type,".rds"))
# 
# data <- data.frame(
#   Group = c(rep("Ori", length(Ori$padj)),rep("DEAlgo", length(DE$padj)),rep("DEAlgo+", length(DE50$padj)),rep("DF", length(DF$padj))),
#   padj = c(Ori$padj,DE$padj,DE50$padj,DF$padj),
#   Gene = c(Ori$pathway,DE$pathway,DE50$pathway,DF$pathway),
#   NES = c(Ori$NES,DE$NES,DE50$NES,DF$NES)
# )
# data <- data[data$Gene %in% unique(overlap),]
# data$pstat <- 1  # Initialize the new column
# data$pstat[data$padj <= 0.05] <- 2
# data$pstat[data$padj <= 0.01] <- 3
# data$pstat[data$padj <= 0.001] <- 4
# data$pstat[is.na(data$padj)] <- NA
# 
# NES_break <- c(-4,-2,0,2,4)
# # Convert pstat to a factor
# data$pstat <- factor(data$pstat)
# 
# # Modify the size breaks to match the levels of pstat
# data$Gene_short <- data$Gene
# data <- data[data$Group != "DEAlgo+",]
# substitute_underscores <- function(string) {
#   gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
# }
# 
# # Apply the function to Gene_short
# data$Gene_short <- substitute_underscores(data$Gene_short)
# 
# p1 <- ggplot(data = data, aes(x = -log10(padj), y = Gene_short, fill = Group)) + 
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
#   geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
#   theme_bw() + 
#   ylab("") + 
#   xlab("-log10(padj)") + 
#   ggtitle(paste0("B Cells"))
# 
# p2 <- ggplot(data = data, aes(x = NES, y = Gene_short, fill = Group)) + 
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
#   theme_bw() + 
#   ylab("") + 
#   xlab("Normalized Enrichment Score (NES)") + 
#   ggtitle(paste0("B Cells"))
# ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","HealthyPBMC_PlateletMonocyte_",dbused,"_",cell_type,".png"), p1+p2, height = 7, width = 15, dpi = 300)


####Adipose tissue dataset
####Fig 4B
dbused <- "GOgmt"
cell_type <- "M4" #NK T B Cell
overlap <- c("GOBP_LIPID_LOCALIZATION","GOBP_FATTY_ACID_TRANSPORT","GOBP_ERK1_AND_ERK2_CASCADE","GOBP_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND","GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM","GOBP_VASCULATURE_DEVELOPMENT","GOBP_TUBE_MORPHOGENESIS","GOBP_RESPONSE_TO_ZINC_ION","GOBP_RESPONSE_TO_COPPER_ION","GOBP_INTRACELLULAR_ZINC_ION_HOMEOSTASIS")

Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","DiaSATvsNonSAT_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <- data[data$Gene %in% unique(overlap),]

data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
data$Gene_short <- factor(data$Gene_short, levels = overlap)
data <- data[order(data$Group, data$Gene_short), ]
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)
data$Gene_short <- factor(data$Gene_short, levels = unique(data$Gene_short))
p1 <- ggplot(data = data, aes(x = -log10(padj), y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0(""))

p2 <- ggplot(data = data, aes(x = NES, y = Gene_short, fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  scale_x_continuous(breaks = seq(0, max(data$NES), by = 1)) + # Set x-axis breaks to 0, 1, 2, 3, ...
  ggtitle(paste0(""))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","P7_",dbused,"_",cell_type,"_1.png"), p1, height = 5, width = 7, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","P7_",dbused,"_",cell_type,"_2.png"), p2, height = 5, width = 5, dpi = 300)

clusterx <- cell_type
dbused <- "GOgmt"
tables <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/", "DiaSATvsNonSAT_Result_",dbused,"_",cell_type,".rds"))
lstcalc <- list()
for (pathplot in overlap){
  leading_edge <- unlist(tables$OriResult.leadingEdge[which(tables$OriResult.pathway == pathplot)])[1:20]
  leading_edge <- leading_edge[!is.na(leading_edge)]
  lstcalc <- c(setNames(list(leading_edge), as.character(pathplot)), lstcalc)
}
dealgoseuobj<- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Output_Overlap_Ratio_0.5/scCLINICResult.rds"))
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))
#cluster2 <- subset(cluster2, subset = status == "Healthy")
for (sigpath in unique(names(lstcalc))){
  p1 <- FeaturePlot(cluster2, reduction = "umap",features =lstcalc[[sigpath]])
  p2 <- VlnPlot(cluster2,features =lstcalc[[sigpath]],group.by = "seurat_clusters")
  p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =lstcalc[[sigpath]])
  p4 <- VlnPlot(dealgoseuobj,features =lstcalc[[sigpath]],group.by = "Overlap_Ratio_0.5")
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",sigpath,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",sigpath,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",sigpath,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",sigpath,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
}

clusterx <- cell_type
Positive_regulation_of_immune_system_process_plot <- c("APOD","APOE","PLA2G2A","DCN","CD36","GNLY","NKG7","CCL5")
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limits = c(0, 6))
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)

clusterx <- cell_type
Positive_regulation_of_immune_system_process_plot <- c("MT2A" ,  "MT1E"   ,"MT1M"  , "MT1A"  ,"MT1F"  ,"MT1G"  ,   "MT1X" )
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limits = c(0, 6), oob = scales::squish)
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)

clusterx <- cell_type
Positive_regulation_of_immune_system_process_plot <- c("CD3E", "CD3D", "CD3G","TRAC","CD79A","CD79B","NKG7")
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limits = c(0, 6))
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)

#Estrogen Response Late Pathway in all celltype
####Fig 4B
dbused <- "Hgmt"
overlap <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_ESTROGEN_RESPONSE_LATE")

# Function to read RDS files and extract necessary information
read_and_extract <- function(dbused, step) {
  file_path <- paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/",
                      "DiaSATvsNonSAT_Result_", dbused, "_", step, ".rds")
  data <- readRDS(file_path)
  data.frame(
    Group = c(rep("Ori", length(data$OriResult.padj)),
              rep("DEAlgo", length(data$DEAlgoResult.padj)),
              rep("DF", length(data$DFResult.padj))),
    padj = c(data$OriResult.padj, data$DEAlgoResult.padj, data$DFResult.padj),
    Gene = c(data$OriResult.pathway, data$DEAlgoResult.pathway, data$DFResult.pathway),
    NES = c(data$OriResult.NES, data$DEAlgoResult.NES, data$DFResult.NES)
  )
}

# Read and combine data for all 0, 2, 5, and 8
Hgmt_0 <- read_and_extract("Hgmt", "M1")
Hgmt_2 <- read_and_extract("Hgmt", "M3")
Hgmt_6 <- read_and_extract("Hgmt", "M4")
Hgmt_10 <- read_and_extract("Hgmt", "M5")

Hgmt_0$Celltype <- "M1"
Hgmt_2$Celltype <- "M3"
Hgmt_6$Celltype <- "M4"
Hgmt_10$Celltype <- "M5"

# Combine data into a single dataframe
combined_data <-rbind(Hgmt_0, Hgmt_6) #rbind(Hgmt_0, Hgmt_6) #rbind(Hgmt_0, Hgmt_2, Hgmt_6, Hgmt_10)
data <- combined_data
data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- paste0(data$Gene," (",data$Celltype,")")
data <- data[data$Group != "DEAlgo+",]
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

data <- data %>%
  mutate(Gene = factor(Gene, levels = overlap)) %>%
  arrange(Group, Gene)

data <- data[! data$Gene_short %in% c("HALLMARK_TNFA_SIGNALING_\nVIA_NFKB (M3)","HALLMARK_TNFA_SIGNALING_\nVIA_NFKB (M4)","HALLMARK_TNFA_SIGNALING_\nVIA_NFKB (M5)"),]
p1 <- ggplot(data = data, aes(x = -log10(padj), y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0(""))

p2 <- ggplot(data = data, aes(x = NES, y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f", "#be883d"),labels = c("scCLINIC", "DF", "All cells")) +
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  scale_x_continuous(breaks = seq(0, max(data$NES), by = 1)) + # Set x-axis breaks to 0, 1, 2, 3, ...
  ggtitle(paste0(""))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Hallmark_Fig4B1p1b.png"), p1, height = 2, width = 7, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Hallmark_Fig4B1p2b.png"), p2, height = 2, width = 5, dpi = 300)

HALLMARK_ESTROGEN_RESPONSE_LATE <- c("CXCL14", "PTGER3", "MEST", "DNAJC1", "BTG3", "PDLIM3", "XBP1", "ELOVL5", "CPE", "BAG1", "ASS1", "CHPT1", "DCXR", "NXT1", "CYP26B1", "FABP5", "JAK1", "EGR3", "TPBG", "KLF4", "IL6ST", "FDFT1", "UGDH", "PERP", "CKB", "RABEP1", "TOB1")
HALLMARK_ESTROGEN_RESPONSE_LATE <- c("CAV1","CXCL14", "CXCL12", "KLF4", "IGFBP4", "CD9", "DUSP2", "JAK1", "IL6ST", "ADD3", "BLVRB", "SGK1", "ATP2B4", "DYNLT3", "FOS")#c("CXCL14", "CAV1", "CXCL12", "KLF4", "IGFBP4", "CD9", "DUSP2", "JAK1", "IL6ST", "ADD3", "BLVRB", "SGK1", "ATP2B4", "DYNLT3", "FOS")
HALLMARK_ESTROGEN_RESPONSE_LATE <- c("KLF4","CAV1","CXCL14", "CXCL12", "IGFBP4", "CD9","FOS" )

for (clusterx in c("M1","M3","M4","M5")){
  cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",clusterx,".rds"))
  
p1 <- FeaturePlot(cluster2, reduction = "umap",features =HALLMARK_ESTROGEN_RESPONSE_LATE) & scale_color_gradientn(colors = c("grey","red"),limits = c(0, 6))
p2 <- VlnPlot(cluster2,features =HALLMARK_ESTROGEN_RESPONSE_LATE,group.by = "seurat_clusters")
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =HALLMARK_ESTROGEN_RESPONSE_LATE)
p4 <- VlnPlot(dealgoseuobj,features =HALLMARK_ESTROGEN_RESPONSE_LATE,group.by = "Overlap_Ratio_0.5")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",HALLMARK_ESTROGEN_RESPONSE_LATE,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",HALLMARK_ESTROGEN_RESPONSE_LATE,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",HALLMARK_ESTROGEN_RESPONSE_LATE,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4B_",HALLMARK_ESTROGEN_RESPONSE_LATE,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)
}



#####Kidney
####Fig 4B
dbused <- "Hgmt"
cell_type <- "ProxTubule"
overlap <- c("HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_BILE_ACID_METABOLISM")

Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("Ori", length(Hlst$OriResult.padj)),rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$OriResult.padj,Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$OriResult.pathway,Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$OriResult.NES,Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)

data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

data <- data %>%
  mutate(Gene = factor(Gene, levels = overlap)) %>%
  arrange(Group, Gene)

p1 <- ggplot(data = data, aes(x = -log10(padj), y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f","#be883d"),labels = c("scCLINIC", "DF", "All cells")) + #,labels = c("scCLINIC", "DF", "All cells")
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0("Prox.Tubule"))

p2 <- ggplot(data = data, aes(x = NES, y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f","#be883d"),labels = c("scCLINIC", "DF", "All cells")) + #,labels = c("scCLINIC", "DF", "All cells")
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  ggtitle(paste0("Prox.Tubule"))+
  scale_x_continuous(breaks = seq(0, ceiling(max(data$NES)), by = 1))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_",cell_type,"_1.png"), p1, height = 2.5, width = 7, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_",cell_type,"_2.png"), p2, height = 2.5, width = 5, dpi = 300)

dbused <- "Hgmt"
cell_type <- "Stroma"
overlap <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_BILE_ACID_METABOLISM")

Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_",cell_type,".rds"))

data <- data.frame(
  Group = c(rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)

data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

data <- data %>%
  mutate(Gene = factor(Gene, levels = overlap)) %>%
  arrange(Group, Gene)

p1 <- ggplot(data = data, aes(x = -log10(padj), y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f"),labels = c("All Cells vs scCLINIC", "All Cells vs DF")) + #,labels = c("scCLINIC", "DF", "All cells")
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0("Stroma"))

p2 <- ggplot(data = data, aes(x = NES, y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f"),labels = c("All Cells vs scCLINIC", "All Cells vs DF")) + #,labels = c("scCLINIC", "DF", "All cells")
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  ggtitle(paste0("Stroma"))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_",cell_type,".png"), p1+p2, height = 2.5, width = 15, dpi = 300)

###X_Ori
dbused <- "Hgmt"
cell_type <- "Stroma"
overlap <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_BILE_ACID_METABOLISM")

Hlst <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_",cell_type,"_X_Ori.rds"))

data <- data.frame(
  Group = c(rep("DEAlgo", length(Hlst$DEAlgoResult.padj)),rep("DF", length(Hlst$DFResult.padj))),
  padj = c(Hlst$DEAlgoResult.padj,Hlst$DFResult.padj),
  #Value = c(OriH$NES,DEAlgoH$NES,DFH$NES),
  Gene = c(Hlst$DEAlgoResult.pathway,Hlst$DFResult.pathway),
  NES = c(Hlst$DEAlgoResult.NES,Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$DEAlgoResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$DEAlgoResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)

data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2,0,2,4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
substitute_underscores <- function(string) {
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <- substitute_underscores(data$Gene_short)

data <- data %>%
  mutate(Gene = factor(Gene, levels = overlap)) %>%
  arrange(Group, Gene)

p1 <- ggplot(data = data, aes(x = -log10(padj), y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f"),labels = c("scCLINIC vs All Cells", "DF vs All Cells")) + #,labels = c("scCLINIC", "DF", "All cells")
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", size = 1) +  # Add dotted vertical line
  theme_bw() + 
  ylab("") + 
  xlab("-log10(padj)") + 
  ggtitle(paste0("Stroma"))

p2 <- ggplot(data = data, aes(x = NES, y = fct_inorder(Gene_short), fill = Group)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#9671c3","#69a75f"),labels = c("scCLINIC vs All Cells", "DF vs All Cells")) + #,labels = c("scCLINIC", "DF", "All cells")
  theme_bw() + 
  ylab("") + 
  xlab("Normalized Enrichment Score (NES)") + 
  ggtitle(paste0("Stroma"))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","GSEA_Result_",dbused,"_",cell_type,"_X_Ori.png"), p1+p2, height = 2.5, width = 15, dpi = 300)


Positive_regulation_of_immune_system_process_plot <- c("Slc23a1","Slc27a2","Abcd3","Fmo1","Acox1","Rdh16","Acox3")

clusterx <- "Stroma"
cluster2 <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/annotation_index_recluster/annotation_index_cluster_M8.rds"))
dealgoseuobj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/DEAlgoResult.rds")
p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","red"),limits = c(0, 4), oob = scales::squish)
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "DEAlgo_ClusterID")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)

Positive_regulation_of_immune_system_process_plot <- c("Lyn","Cybb","Plaur","Ifitm2","Ifitm3","Csf2rb","Csf2ra","Cd44")

p1 <- FeaturePlot(cluster2, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)& scale_color_gradientn(colors = c("grey","blue"),limits = c(0, 4), oob = scales::squish)
p2 <- VlnPlot(cluster2,features =Positive_regulation_of_immune_system_process_plot,group.by = "seurat_clusters")
p3 <- FeaturePlot(dealgoseuobj, reduction = "umap",features =Positive_regulation_of_immune_system_process_plot)
p4 <- VlnPlot(dealgoseuobj,features =Positive_regulation_of_immune_system_process_plot,group.by = "DEAlgo_ClusterID")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_1.png"), p1, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_2.png"), p2, height = 10, width = 10, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_3.png"), p3, height = 15, width = 15, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig4B_",Positive_regulation_of_immune_system_process_plot,"_",clusterx,"_4.png"), p4, height = 10, width = 10, dpi = 300)

#Azimuth
obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig3/step1d.rds")
obj_azimuth <- Azimuth::RunAzimuth(obj,reference = "kidneyref")


file_name <- paste0("annotation_index","_cluster_","M8",".rds")
recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/annotation_index_recluster/",file_name))

recluster$azimuthl1 <- obj_azimuth$predicted.annotation.l1
recluster$azimuthl2 <- obj_azimuth$predicted.annotation.l2

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

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig3_cluster_","M8","_azimuth.png"), p2+p3, height = 5, width = 12, dpi = 300)

cell1 <- WhichCells(recluster, expression = azimuthl1 == "Proximal Tubule")
cell2 <- WhichCells(recluster, expression = azimuthl2 == "Non-classical monocyte")

p1 <- DimPlot(recluster, cells.highlight = list("Proximal Tubule" = cell1, "Non-classical monocyte" = cell2), 
            cols.highlight = c("#d55046", "#1e88e5")) +
            scale_color_manual(values = c("Proximal Tubule" = "#d55046", "Non-classical monocyte" = "#1e88e5"),
            labels = c("Proximal Tubule\n(scCLINIC +ve)", "Non-classical monocyte\n(DF +ve)")) +
            guides(color = guide_legend(override.aes = list(size = 5))) +
            theme(legend.title = element_blank())
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig3_cluster_","M8","_StromaContam.png"), p1, height = 5, width = 7, dpi = 300)

saveRDS(obj_azimuth,"~/DEAlgoManuscript/Manuscript_Figures/Fig3/Fig3_Azimuth.rds")

