#DEAlgo Published Fig3b_PCT

#Load and save Robj
# load("/data/seu_kidney_codeocean.Robj") <- from https://codeocean.com/capsule/5650599/tree/v1
# saveRDS(seu_kidney,"~/DEAlgoManuscript/Manuscript_Figures/Fig3/seu_kidney_codeocean.rds")
#Convert Old seurat (s4 object) to new seurat
# pbmc_demuUMAP_Harmony <- readRDS(Indir)
# seurat_object <- CreateSeuratObject(counts = pbmc_demuUMAP_Harmony@raw.data)
# seurat_object@meta.data$CT.Park <- pbmc_demuUMAP_Harmony@meta.data$CT.Park
# saveRDS(seurat_object,"~/DEAlgoManuscript/Manuscript_Figures/Fig3/seu_kidney_codeocean_updated.rds")

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

Name <- "kidneymouse_science"
Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig3/"
filteredmatrix=NA
rawmatrix=NA
resol="Manual"
OverlapRatio="CT.Park"
ISThreshold = 0
Cutoff=0
gene_n=150
CELLANNOTATION = TRUE

obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig3/seu_kidney_codeocean_updated.rds")

obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)

obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold,CELLANNOTATION = TRUE)

saveRDS(obj,"~/DEAlgoManuscript/Manuscript_Figures/Fig3/step1d.rds")

STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)

obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig3/step1d.rds")

obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,Cutoff,filteredmatrix,rawmatrix,CELLANNOTATION = TRUE)

PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)


#DF
###################################################
## Step 2: Perform and summarize parameter sweep ##
###################################################
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
library(DoubletFinder)
seu_kidney <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig3/seu_kidney_codeocean_updated.rds")
seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

print("Performing pN-pK parameter sweep...")
sweep.res.list_kidney <- paramSweep_v3(seu_kidney)

print("Finding optimal pK using BCmvn maximization...")
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
pK <- bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]
pK <- as.numeric(levels(pK))[pK]
saveRDS(bcmvn_kidney,"~/DEAlgoManuscript/Manuscript_Figures/Fig3/bcmv_kidney.rds")
pK_kidney <- 0.09 ## Visually-discerned from find.pK graphical output

######################################################################
## Step 3: Model homotypic doublet proportions, set pANN thresholds ##
######################################################################
print("Modeling homotypic doublet proportions...")
annotations <- seu_kidney@meta.data$CT.Park
homotypic.prop <- modelHomotypic(annotations)

print("Setting pANN thresholds...")
poi.dfr <- 0.075
nExp_poi <- round(poi.dfr*ncol(seu_kidney))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##################################################################
## Step 4: Run DoubletFinder using (un)adjusted pANN thresholds ##
##################################################################
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = pK_kidney, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

saveRDS(seu_kidney,"~/DEAlgoManuscript/Manuscript_Figures/Fig3/seu_kidney_codeocean_updated_withDFScore.rds")



##Figure3BFinalNEWVERSION.R
##Fig3B Compare DF vs DEAlgo
#read doublet
res <- "annotation_index"
dealgoseuobj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/DEAlgoResult.rds")
doublet_afqc <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig3/seu_kidney_codeocean_updated_withDFScore.rds")

dealgoseuobj@meta.data$DFafqc_pANN <- doublet_afqc@meta.data[, grep("pANN",  colnames(doublet_afqc@meta.data), value = TRUE)]

dealgoseuobj@meta.data$DFafqc <- doublet_afqc@meta.data[, grep("DF.",  colnames(doublet_afqc@meta.data), value = TRUE)]

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < 4, "Artifact", "Singlet")


library(ggrepel)
library(dplyr)

qcsclst <- sort(unique(dealgoseuobj@meta.data[,res]))
cluster_to_consider <- list()
for (i in sort(qcsclst)){
  file_name <- paste0(res,"_cluster_",i,".csv")
  if (file_name %in% list.files("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/annotation_index_recluster/Marker/")){
    cluster_to_consider <- unlist(c(cluster_to_consider,i))
  }
}
for (i in cluster_to_consider){
  file_name <- paste0(res,"_cluster_",i,".rds")
  if (file_name %in% list.files("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/annotation_index_recluster/")){
    recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/annotation_index_recluster/",file_name))
    
    recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
    recluster$DFafqc <- dealgoseuobj$DFafqc
    recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
    recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
    recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
    recluster$DEAlgo_Contaminated <- dealgoseuobj$DEAlgo_Contaminated
    
    if (all(is.na(unique(recluster$DEAlgocluster_Contam)))){
      p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
      p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
      p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
      p0 <- DimPlot(recluster, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1)
      
      ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig3_cluster_",i,"_rawscore_Validate.png"), p0+p1+p7+p8, height = 10, width = 12, dpi = 300)
    }
    else{
      p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
      p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
      p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE,  sizes.highlight = 0.1)
      p5 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE,  sizes.highlight = 0.1)
      p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
      p0 <- DimPlot(recluster, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1)
      
      ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig3_cluster_",i,"_rawscore_Validate.png"), p7+p5+p8+p4+p0+p1, height = 10, width = 17, dpi = 300)
      
    }
  }
}

p1 <- DimPlot(dealgoseuobj,group.by = "CT.Park", reduction = "umap",raster=FALSE)
p3 <- FeaturePlot(dealgoseuobj,features = "DiffDP1_sum", reduction = "umap",raster=FALSE, label = F) #& scale_color_gradient(limits = c(0, 1))
p5 <- FeaturePlot(dealgoseuobj, reduction = "umap",features = "DFafqc_pANN")
p7 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1)
p11 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE,  sizes.highlight = 0.1)
p12 <- DimPlot(dealgoseuobj, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE,  sizes.highlight = 0.1)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/","Fig3b.png"), p1 + p12 + p3 + p11 + p7 + p5, height = 15, width = 30, dpi = 300)


celltype <- dealgoseuobj@meta.data$CT.Park
state_counts <- table(celltype)
print(state_counts)
# celltype
# Coll.Duct.IC     Coll.Duct.PC     Coll.Duct.TC Dist.Prox.Tubule             Endo    Loop.of.Henle      Prox.Tubule           Stroma 
# 461              238               39             1979              273              296             8158              733
state_percentages <- prop.table(state_counts) * 100
# Print the percentages
print(state_percentages)
# celltype
# Coll.Duct.IC     Coll.Duct.PC     Coll.Duct.TC Dist.Prox.Tubule             Endo    Loop.of.Henle      Prox.Tubule           Stroma 
# 3.7858257        1.9545044        0.3202759       16.2519504        2.2419315        2.4308122       66.9951548        6.0195450

table(dealgoseuobj$DEAlgo_Contaminated)
# Artifact  Singlet 
# 853    11324  
table(dealgoseuobj$DFafqc)
# Doublet Singlet 
# 913   11264 
#

#NOTE: Coll.Duct.IC (M1)    Coll.Duct.PC  (M2)  Coll.Duct.TC (M3) Dist.Prox.Tubule    (M4)     Endo (M5)  Loop.of.Henle (M6)   Prox.Tubule  (M7)    Stroma (M8)
DEAlgoContam1 <- c("M1S4","M2S6","M4S0","M5S0","M6S5","M8S2")#Stroma (M8), Endo (M5), Dist.Prox.Tubule (M4), Coll.Duct.IC (M1), Loop of Henle (M6)

Contam1Cell <- subset(dealgoseuobj,subset = DEAlgo_ClusterID %in% DEAlgoContam1)
noncontam1 <- subset(dealgoseuobj,subset = L0 %in% c("M1","M2","M4","M5","M6","M8"))

negative <- noncontam1[,!(noncontam1$DEAlgo_ClusterID %in% DEAlgoContam1)]
Contam1CellDF <- subset(Contam1Cell,subset = DFafqc == "Doublet")
Contam1CellNDF <- subset(Contam1Cell,subset = DFafqc == "Singlet")


noncontam1$grouping <- NA
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(negative), "Singlet",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellDF), "DF",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellNDF), "NDF",noncontam1$grouping)


# expall <- AverageExpression(noncontam1,group.by = "grouping",slot = "data")
#pseudo_ifnb <- AggregateExpression(noncontam1, assays = "RNA", group.by = "grouping")#IF return.seurat = T mean normalized again? so should I
#pseudodata <- pseudo_ifnb@assays$RNA$data


negativeexp <- AverageExpression(negative,group.by = "all",slot = "data")
Contam1CellDFexp <- AverageExpression(Contam1CellDF,group.by = "all",slot = "data")
Contam1CellNDFexp <- AverageExpression(Contam1CellNDF,group.by = "all",slot = "data")
Contam1CellNDFexp_split <- AverageExpression(Contam1CellNDF,group.by = "DEAlgo_ClusterID",slot = "data")

contam1gene <- read.csv("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/overlaplst_filtered_annotation_index_ContaminationInfo.csv",na.strings = c("", "NA"))
PTMarkergene <- read.csv("~/DEAlgoManuscript/PTMarkerGene.csv",na.strings = c("", "NA"))
contam1genelist <- unique(contam1gene[contam1gene$cluster_ref=="M7",]$gene)
PTMarkergenelist <- unique(PTMarkergene$genes)

intersection <- intersect(contam1genelist, PTMarkergenelist)

# filtered_PTMarkergenelist <- PTMarkergenelist[!PTMarkergenelist %in% "Fth1"]
# negative_res <- negativeexp$RNA[filtered_PTMarkergenelist,]
# Contam1CellDFexp_res <- Contam1CellDFexp$RNA[filtered_PTMarkergenelist,]
# Contam1CellNDFexp_res <- Contam1CellNDFexp$RNA[filtered_PTMarkergenelist,]
# 
# Contam1CellNDFexp_split_res <- Contam1CellNDFexp_split$RNA[filtered_PTMarkergenelist,]

negative_res <- negativeexp$RNA[intersection,]
Contam1CellDFexp_res <- Contam1CellDFexp$RNA[intersection,]
Contam1CellNDFexp_res <- Contam1CellNDFexp$RNA[intersection,]

# expallres <- expall$RNA[filtered_PTMarkergenelist,]
# # Get gene names from row names
# Convert sparse matrix to data frame
# expallres_df <- as.data.frame(as(Contam1CellNDFexp_split_res, "matrix"))
# expallres_df$Gene <- rownames(expallres_df)
# 
# # Reshape data frame
# reshaped_df <- reshape2::melt(expallres_df, id.vars = "Gene", variable.name = "Condition", value.name = "Value")

# Load necessary libraries
#69a75f
#9671c3
#be883d
#cc5366

library(ggplot2)
# Create data frames for each condition
data_negative <- data.frame(Gene = names(negative_res), Value = negative_res)
data2 <- data.frame(Gene = names(Contam1CellDFexp_res), Value = Contam1CellDFexp_res)
data3 <- data.frame(Gene = names(Contam1CellNDFexp_res), Value = Contam1CellNDFexp_res)

# Combine all data frames
all_data <- rbind(data_negative, data3, data2)
all_data$Condition <- factor(rep(c("Singlet", "DEAlgo", "DF+DEAlgo"), each = length(names(negative_res))),levels = c("DF+DEAlgo", "DEAlgo","Singlet"))
# all_data <- rbind(all_data,reshaped_df)
# Create the violin plot with scatter plot
partial_contamination_data <- subset(all_data, Condition == "DEAlgo" & Value > 6)

p1 <-ggplot(all_data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(color = alpha("black",0.25)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  geom_point(aes(color = Condition), size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_text_repel(data = partial_contamination_data, aes(label = Gene), 
                  hjust = 0.5, vjust = 0, segment.color = "transparent", size = 2,position = position_jitter(width = 0.2))+
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linetype = "dashed", color = alpha("black",0.5), size = 0.5, position = position_dodge(width = 0.75)) + # Add average line
  
  theme_minimal() +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "Proximal Tubule Marker Average Expression Level")+
  
  theme(legend.position = "none") # Remove legend



ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_A.png"), p1, height = 5, width = 5, dpi = 300)

# Perform pairwise t-tests
pairwise_result <- pairwise.t.test(all_data$Value, all_data$Condition)

write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_intersect.csv"),sep = ",")


partial_contamination_data <- subset(all_data, Condition == "Partial Contamination" & Value > 5)
p1 <- ggplot(all_data, aes(x = Condition, y = Value, fill = Condition, group = Gene)) +
  geom_line(alpha = 0.2)+
  geom_point(size = 2, alpha = 0.4)+
  geom_text_repel(data = partial_contamination_data, aes(label = Gene),
                  hjust = 0.5, vjust = 0.5, segment.color = "transparent")+
  labs(title = "", x = "", y = "Average Expression Level") +
  scale_color_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  theme_minimal()+
  theme(legend.position = "none")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_geomline.png"), p1, height = 5, width = 5, dpi = 300)

#Proved is not Doublet
p <- VlnPlot(noncontam1,features = "nCount_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal() +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nCount_RNA") +
  theme(legend.position = "none")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nCount.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nCount_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nCount.csv"),sep = ",")

p <- VlnPlot(noncontam1,features = "nFeature_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal() +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nFeature_RNA") +
  theme(legend.position = "none") # Remove legend
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nFeature_RNA.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nFeature_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nFeature_RNA.csv"),sep = ",")

#Partial contamination by multiple cell type in Endo (M5)
DEAlgoContam1 <- c("M5S5","M5S4")#Pure DBSCAN G7L2 G6L2 G3L2 with other G2L3 (V) G4L2 (X)
Contam1Cell <- subset(dealgoseuobj,subset = DEAlgo_ClusterID %in% DEAlgoContam1)
noncontam1 <- subset(dealgoseuobj,subset = L0 %in% c("M5"))

negative <- noncontam1[,as.numeric(levels(noncontam1$DEAlgocluster_Contam))[noncontam1$DEAlgocluster_Contam] >= 4] 
positive <- noncontam1[,as.numeric(levels(noncontam1$DEAlgocluster_Contam))[noncontam1$DEAlgocluster_Contam] < 4]

Contam1CellDF <- subset(Contam1Cell,subset = DFafqc == "Doublet")
Contam1CellNDF <- subset(Contam1Cell,subset = DFafqc == "Singlet")

noncontam1$grouping <- "NA"
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(negative), "Singlet",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellDF), "DF",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellNDF), "NDF",noncontam1$grouping)

negativeexp <- AverageExpression(negative,group.by = "all",slot = "data")
Contam1CellDFexp <- AverageExpression(Contam1CellDF,group.by = "all",slot = "data")
Contam1CellNDFexp <- AverageExpression(Contam1CellNDF,group.by = "all",slot = "data")

Contam1Cellexp <- AverageExpression(Contam1Cell,group.by = "all",slot = "data")
positiveexp <- AverageExpression(positive,group.by = "all",slot = "data")
contamgeneinfo <- read.csv("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/overlaplst_filtered_annotation_index_ContaminationInfo.csv",na.strings = c("", "NA"))

genelist1 <- list()
for (ref in c("M1","M2","M4","M6")){
  contamXgene <- contamgeneinfo[contamgeneinfo$cluster_ref==ref & contamgeneinfo$globalduplicate > 1,]
  contamXgenelist <- unique(contamXgene$gene)
  genelist1 <- c(contamXgenelist,genelist1)
}
genelist <- list()
reflist <- list()
for (ref in c("CollDuctIC","CollDuctPC","DistProxTubule","LoopofHenle")){
  markergene <- read.csv(paste0("~/DEAlgoManuscript/KidneyMarkerCTPark/",ref,".csv"),na.strings = c("", "NA"))
  markergenelist <- unique(markergene$genes)
  genelist <- c(markergenelist,genelist)
  
  reflist <- c(rep(length(markergenelist),x = ref),reflist)
}

intersection <- intersect(unlist(genelist1), unlist(genelist))

markergenelist_filtered <- intersection[!intersection %in% "mt-Co1"]

# Find indices of elements in genelist that are in the intersection
indices <- match(markergenelist_filtered, unlist(genelist))

# Filter reflist corresponding to the intersection
filtered_reflist <- reflist[indices]

negative_res <- negativeexp$RNA[markergenelist_filtered,]
Contam1Cellexpres <- Contam1Cellexp$RNA[markergenelist_filtered,]
#positiveexpres <- positiveexp$RNA[markergenelist_filtered,]

data_negative <- data.frame(Gene = names(negative_res), Value = negative_res, Ref = unlist(filtered_reflist))
data2 <- data.frame(Gene = names(Contam1Cellexpres), Value = Contam1Cellexpres, Ref = unlist(filtered_reflist))
#data2 <- data.frame(Gene = names(positiveexpres), Value = positiveexpres, Ref = unlist(filtered_reflist))

negative_res <- negativeexp$RNA[markergenelist_filtered,]
Contam1CellDFexp_res <- Contam1CellDFexp$RNA[markergenelist_filtered,]
Contam1CellNDFexp_res <- Contam1CellNDFexp$RNA[markergenelist_filtered,]

data_negative <- data.frame(Gene = names(negative_res), Value = negative_res, Ref = unlist(filtered_reflist))
data2 <- data.frame(Gene = names(Contam1CellDFexp_res), Value = Contam1CellDFexp_res, Ref = unlist(filtered_reflist))
data3 <- data.frame(Gene = names(Contam1CellNDFexp_res), Value = Contam1CellNDFexp_res, Ref = unlist(filtered_reflist))

# Combine all data frames
all_data <- rbind(data_negative, data2, data3)
all_data$Condition <- factor(rep(c("Singlet", "DF+DEAlgo","DEAlgo"), each = length(names(negative_res))),levels = c("DF+DEAlgo","DEAlgo","Singlet"))
partial_contamination_data <- subset(all_data, Condition == "DEAlgo" & Value > 6)
#all_data <- rbind(all_data,reshaped_df)
# Create the violin plot with scatter plot
p1 <-ggplot(all_data, aes(x = Condition, y = Value, fill = Condition, color = as.factor(Ref), group = Condition)) +
  geom_violin(color = alpha("black", 0.25)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white", outlier.shape = NA) +
  geom_text_repel(data = partial_contamination_data, aes(label = Gene), segment.color = "transparent", size = 2,position = position_jitter(width = 0.2))+
  geom_point(size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linetype = "dashed", color = alpha("black",0.5), size = 0.5, position = position_dodge(width = 0.75)) + # Add average line
  
  scale_color_manual(values = c("#648ace", "#cc5143", "#58a865","#9671c3")) +  # Custom colors based on Ref
  scale_fill_manual(values = c("white","white", "white")) + # Specify fill colors
  theme_minimal() +
  labs(title = "", x = "", y = "CD-PC, LOH, CD-IC Marker Average Expression Level") +
  theme(legend.position = "none") # Remove legend


pairwise_result <- pairwise.t.test(all_data$Value, all_data$Condition)

write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_marker_multipleDFDEAlgo.csv"),sep = ",")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_B.png"), p1, height = 5, width = 5, dpi = 300)

partial_contamination_data <- subset(all_data, Condition == "Partial Contamination" & Value > 5)
p1 <- ggplot(all_data, aes(x = Condition, y = Value, fill = Condition, color = as.factor(Ref), group = Gene)) +
  geom_line(alpha = 0.4)+
  geom_point(size = 1, alpha = 0.4)+
  geom_text_repel(data = partial_contamination_data, aes(label = Gene),
                  hjust = -0.5, segment.color = "transparent")+
  labs(title = "", x = "", y = "Average Expression Level") +
  scale_color_manual(values = c("#648ace", "#cc5143", "#58a865","#9671c3")) + # Specify fill colors
  theme_minimal()+
  theme(legend.position = "none")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_marker_multipleDFDEAlgogeomline.png"), p1, height = 5, width = 5, dpi = 300)

noncontam1 <- subset(noncontam1,subset = grouping %in% c("DF","NDF","Singlet"))
p1 <- VlnPlot(noncontam1,features = "nCount_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal() +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nCount_RNA") +
  theme(legend.position = "none") # Remove legend
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nCount_RNAMultiple.png"), p1, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nCount_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nCount_RNAMultiple.csv"),sep = ",")

p1 <- VlnPlot(noncontam1,features = "nFeature_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal() +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nFeature_RNA") +
  theme(legend.position = "none") # Remove legend
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nFeature_RNAMultiple.png"), p1, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nFeature_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3B_DFNDEAlgoP_nFeature_RNAMultiple.csv"),sep = ",")



#Triplet Contam
color4marker <- c("#e63c29", "#e63c29","#9146ec","#ed9307","#e63c29", "#008000","#9146ec")
markertoseelist <- c("Apela","Npnt","Tspan13","Slc12a3","Scnn1b","Slc5a3","Itpr2")
limitlst <- c(5,4,4,6,4,4,4)
plot_list <- list()
plot_list1 <- list()
for (markertosee in seq(length(markertoseelist))){
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/annotation_index_recluster/","annotation_index_cluster_","M5",".rds"))
  p2<- FeaturePlot(dealgoseuobj, reduction = "umap",features = markertoseelist[markertosee], pt.size = 0.05)+scale_color_gradientn( colours = c('lightgrey', color4marker[markertosee]),  limits = c(0, limitlst[markertosee])) 
  p1 <- FeaturePlot(recluster, reduction = "umap",features = markertoseelist[markertosee], pt.size = 0.1)+scale_color_gradientn( colours = c('lightgrey', color4marker[markertosee]),  limits = c(0, limitlst[markertosee]))
  plot_list[[markertosee]] <- p1
  plot_list1[[markertosee]] <- p2
}

all_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = ceiling(length(markertoseelist)/2))
all_plots1 <- grid.arrange(grobs = plot_list1, nrow = 2, ncol = ceiling(length(markertoseelist)/2))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig33W.png"), plot = all_plots, height = 6, width = 13, dpi = 300)
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/","Fig3BW.png"), plot = all_plots1, height = 12, width = 23, dpi = 300)


library(gridExtra)
celltypelist <- unique(dealgoseuobj$L0)[1:7]
# c("Slc34a1","Acsm2","Slc27a2","Dnase1","Miox","Pck1","Ass1","Akr1c21","Errfi1","Lrp2","Sord","Slc4a4")
# c("Umod","Egf","Wfdc15b","Mt1","Ppp1r1a","Ly6a","Slc12a3","Calb1","Wnk1","Oxct1","Apela")
# markergenelist[markergenelist$celltype=="5",]$genes
for (markertosee in c("Slc34a1","Acsm2","Slc27a2","Dnase1","Miox","Pck1","Ass1","Akr1c21","Errfi1","Lrp2","Sord","Slc4a4")){
  #markertosee <- "Prox1" #"Pdpn" #Lyve not found
  plot_list <- list()
  for (curcelltype in celltypelist[1:7]){
    
    recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/annotation_index_recluster/","annotation_index_cluster_",curcelltype,".rds"))
    
    p1 <- FeaturePlot(recluster, reduction = "umap",features = markertosee, pt.size = 0.1)#DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE, label = T, sizes.highlight = 0.1)
    #ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/",curcelltype,"_marker_",markertosee,".png"), p1, height = 5, width = 5, dpi = 300)
    plot_list[[curcelltype]] <- p1
  }
  
  all_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = ceiling(length(celltypelist)/2))
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/",markertosee,".png"), plot = all_plots, height = 5, width = 10, dpi = 300)
}




obj <- dealgoseuobj
Outdir <- "~/DEAlgoManuscript/Manuscript_Figures/Fig3/"
Name <- "kidneymouse_science"
OverlapRatio="CT.Park"
CELLANNOTATION = TRUE
###No Return
message("Plot Contamination Pattern.")

if (CELLANNOTATION){
  message("Using user-annotated clusters.")
  obj@meta.data$annotation_index <- paste0("M",as.numeric(factor(obj@meta.data[,OverlapRatio])))#change the manual annotation to index, eg. CellType0 CellType1 CellType2 CellType3 to 1 2 3 4
  OverlapRatio <- "annotation_index" #User manual cellannotation
}else{
  OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
  message(paste0("Using ",OverlapRatio))
}

folder_path_Step2_Output <- paste0(Outdir,Name,"_Step2/Output_",OverlapRatio,"/")

#Load artifacts information and DEAlgo result rds
contamgeneinfo <- read.csv(paste0(folder_path_Step2_Output,"overlaplst_filtered_",OverlapRatio,"_ContaminationInfo.csv"),na.strings = c("", "NA"))

#DEAlgo Subcluster ID
contamgeneinfo$MajorSub <- paste0(contamgeneinfo$globalcluster,contamgeneinfo$cluster_local)
#Summarize ES score and their source of major cluster for each subclusters
result <- contamgeneinfo %>%
  group_by(MajorSub) %>%
  summarize(
    cluser_reflst = list(cluster_ref),
    dp1lst = list(dp1)
  ) %>%
  ungroup()

# Function to calculate the average ES score (dp1lst) for each source of artifacts (cluster_reflst)
tabulate_Cluster_Contribution_Score <- function(cluster_reflst, dp1lst) {
  components <- unlist(cluster_reflst)
  dp1_values <- unlist(dp1lst)
  CCS_Matrix <- tapply(dp1_values, components, mean, na.rm = TRUE)
  return(CCS_Matrix)
}

# For each subclusters (each row in result), calculate the average ES score for each source of artifacts
tabulate_CCS <- mapply(tabulate_Cluster_Contribution_Score, result$cluser_reflst, result$dp1lst, SIMPLIFY = FALSE)

# List of all source of artifacts which contaminated major cluster X
lst_of_source_of_artifacts <- unique(unlist(lapply(tabulate_CCS, names)))

# Create a matrix, each row represent one subclusters and each column present each source of artifacts, to store the CCS for each major clusters
CCS_Matrix <- matrix(NA, nrow = length(tabulate_CCS), ncol = length(lst_of_source_of_artifacts), dimnames = list(NULL, lst_of_source_of_artifacts))

# store the CCS for each major clusters in the matrix
for (i in seq_along(tabulate_CCS)) {
  CCS_Matrix[i, names(tabulate_CCS[[i]])] <- tabulate_CCS[[i]]
}

# Replace NA with 0
CCS_Matrix[is.na(CCS_Matrix)] <- 0

# Convert to dataframe for plotting
CCS_Matrix <- as.data.frame(CCS_Matrix)
CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
data <- CCS_Matrix %>%
  pivot_longer(cols = -MajorSub, names_to = "Component", values_to = "value") #dependency tidyr

filtered_data <- data #keep all clusters

# Heatmap
# filtered_data$value[startsWith(as.character(filtered_data$MajorSub), substr(as.character(filtered_data$Component), 1, 2))] <- 0
# COnvert dataframe
heatmap_data <- dcast(filtered_data, MajorSub ~ Component, value.var = "value")
rownames(heatmap_data) <- heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL #remove MajorSub columns

#heatmap_data$Total <- rowSums(heatmap_data,na.rm = TRUE)

# Summarize DEAlgo Contamination Score (CS) of each subclusters
CS_table <- obj@meta.data %>%
  group_by(DEAlgo_ClusterID) %>%
  summarize(
    DiffDP1_sum = mean(DiffDP1_sum), #DiffDP1_sum of each cells within each subclusters are same value...
  ) %>%
  ungroup()

#Convert to plotting dataframe
heatmap_data <- as.data.frame(heatmap_data)  # Convert heatmap_data to a dataframe
heatmap_data <- heatmap_data %>%
  mutate(Score = CS_table$DiffDP1_sum[match(rownames(heatmap_data), CS_table$DEAlgo_ClusterID)]) #Score

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Calculate the color scale for heatmap based on quantile breaks
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(heatmap_matrix, n = 100)

# Create heatmap
heatmap_matrix_t <- t(heatmap_matrix)
p2 <- pheatmap(heatmap_matrix_t,
               cluster_rows = FALSE,  # Do not cluster rows
               cluster_cols = FALSE,  # Do not cluster columns
               na_col = "grey",  # Fill missing values with grey
               color             = viridis(length(mat_breaks) - 1),
               breaks            = mat_breaks,
               labels_row = rownames(heatmap_matrix_t),
               labels_col = colnames(heatmap_matrix_t),
               show_rownames = TRUE,  # Show row names
               show_colnames = TRUE,
               angle_col = 90,
               fontsize_row = 16,  # Adjust row font size
               fontsize_col = 16)
ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScoreHeatmap_Fig3",".png"), p2, height = 4, width = 30, dpi = 300)

# Plot multicolumn bar plot
multicolplot <- data.frame((heatmap_matrix))
multicolplot$Category <- as.character(rownames(multicolplot))
plot_list <- list()

x_limits <- c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <- seq( x_limits[1], x_limits[2], length.out = 3)
text_size = 5
for (i in seq_along(colnames(heatmap_matrix))) {
  clus <- colnames(heatmap_matrix)[i]
  p <- ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(stat = "identity", fill = viridis(length(colnames(heatmap_matrix)))[i]) +
    labs(title = paste(clus), y ="", x ="") +  # Set y-axis label for the first plot only
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.line = element_line(color = "black"),
          axis.text = element_text(size = text_size),
          axis.title = element_text(size = text_size),
          plot.title = element_text(size = text_size))+
    coord_flip()+
    #ylim(x_limits)+
    scale_y_continuous(breaks = custom_breaks, limits = x_limits) +
    if (i != 1){
      theme(axis.text.y = element_blank())
    }
  
  # Add the plot to the list
  plot_list[[i]] <- p
}
g1 <- grid.arrange(grobs = plot_list, ncol = length(plot_list),right = "")
ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScore_MultiStackFig3",".png"),g1,  height = 5, width = 8, dpi = 300)

