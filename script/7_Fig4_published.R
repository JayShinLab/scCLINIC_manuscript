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

Name <- "AdiposeNatureMetabolism"
Input <- "~/DEAlgoManuscript/Manuscript_Figures/Fig4/adinatobjQC.rds"
Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig4/"
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
saveRDS(obj,"~/DEAlgoManuscript/Manuscript_Figures/Fig4/step1d.rds")
STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig4/step1d.rds")
obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,Cutoff,filteredmatrix,rawmatrix)

PlotContaminationPattern(obj,Output,Name,OverlapRatio)

PlotContaminationPattern <- function(obj,Outdir,Name,OverlapRatio=0.5,CELLANNOTATION = FALSE, verbose = TRUE){
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
  
  for (clusterx in unique(obj@meta.data[,OverlapRatio])){
    # Filter out other clusters, only keep cluster X which wish to display
    filtered_data <- data %>%
      dplyr::filter(startsWith(MajorSub, clusterx))
    
    if (nrow(filtered_data) != 0){
      message(paste0("Plot Contamination Pattern ",clusterx))
      # Stacked Plot of CCS vs Subclusters
      p1 <- ggplot(filtered_data, aes(fill=Component, y=value, x=MajorSub)) +
        geom_bar(position="stack", stat="identity")+
        labs(title = paste0(clusterx),
             x = "Subclusters",
             y = "Contamination Score") +
        scale_fill_manual(values = viridis(length(lst_of_source_of_artifacts)),
                          name = "Source of artifacts") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(size = 8),  # Adjust x-axis title size
              axis.title.y = element_text(size = 8))
      
      ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScoreStackPlot",clusterx,".png"), p1, height = 8, width = 8, dpi = 300)
      
      # Heatmap
      # filtered_data$value[startsWith(as.character(filtered_data$MajorSub), substr(as.character(filtered_data$Component), 1, 2))] <- 0
      # COnvert dataframe
      heatmap_data <- dcast(filtered_data, MajorSub ~ Component, value.var = "value")
      heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts
      
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
      p2 <- pheatmap(heatmap_matrix,
                     cluster_rows = FALSE,  # Do not cluster rows
                     cluster_cols = FALSE,  # Do not cluster columns
                     main = paste0(clusterx),
                     na_col = "grey",  # Fill missing values with grey
                     color             = viridis(length(mat_breaks) - 1),
                     breaks            = mat_breaks,
                     labels_row = rownames(heatmap_matrix),
                     labels_col = colnames(heatmap_matrix),
                     show_rownames = TRUE,  # Show row names
                     show_colnames = TRUE,
                     angle_col = 0)
      ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScoreHeatmap",clusterx,".png"), p2, height = 2.8, width = 8, dpi = 300)
      
      # Plot multicolumn bar plot
      multicolplot <- data.frame((heatmap_matrix))
      
      rownames_sorted <- rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))),decreasing = T)]
      
      # Reorder the heatmap matrix rows according to the sorted row names
      multicolplot <- multicolplot[rownames_sorted, ]
      
      colnames_sorted <- colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))),decreasing = F)]
      
      # Reorder the heatmap matrix rows according to the sorted row names
      multicolplot <- multicolplot[,colnames_sorted]
      
      multicolplot$Category <- as.character(rownames(multicolplot))
      
      multicolplot$Category <- factor(multicolplot$Category, levels = rownames_sorted)
      
      x_limits <- c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
      custom_breaks <- seq( x_limits[1], x_limits[2], length.out = 3)
      text_size = 5
      
      plot_list <- list()
      for (i in seq(colnames_sorted)) {
        clus <- colnames_sorted[i]
        p <- ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
          geom_bar(stat = "identity", fill = viridis(length(colnames_sorted))[i]) +
          labs(title = paste(clus), y ="", x ="") +  # Set y-axis label for the first plot only
          theme_minimal()+
          theme(panel.grid = element_blank(),axis.line = element_line(color = "black"),
                axis.text = element_text(size = text_size),
                axis.title = element_text(size = text_size),
                plot.title = element_text(size = text_size)) +
          coord_flip()+
          #ylim(x_limits)+
          scale_y_continuous(breaks = custom_breaks, limits = x_limits)+
          if (i != 1){
            theme(axis.text.y = element_blank())
          }
        # Add the plot to the list
        plot_list[[i]] <- p
      }
      
      widths <- c(6.2, rep(5, length(plot_list) - 1))
      
      g1 <- grid.arrange(grobs = plot_list, ncol = length(plot_list),right = "",widths = widths)
      ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationScore",clusterx,".png"),g1,  height = 2.5, width = 8, dpi = 300)
      
    }
  }
  
  ###Plot Source of Artifacts Contamination Patterns
  #Plotting and visualize the DEAlgo score of each source of artifacts on the major cluster's UMAP,
  p1 <- FeaturePlot(obj,features = colnames(CCS_Matrix)[colnames(CCS_Matrix) != "MajorSub"])
  ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationPattern_SourceOfArtifacts",".png"),p1,  height = 20, width = 20, dpi = 300)
  
  folder_path_Step2_L1R <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/")
  folder_path_Step2_L1R_Marker <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/Marker/")
  
  qcsclst <- sort(unique(obj@meta.data[,OverlapRatio]))#list of major clusters seurat ID
  cluster_to_consider <- list()
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".csv")
    if (file_name %in% list.files(folder_path_Step2_L1R_Marker)){
      cluster_to_consider <- unlist(c(cluster_to_consider,i))
    }
  }
  
  for (i in cluster_to_consider){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))
      recluster@meta.data <- obj@meta.data[rownames(recluster@meta.data),]#copy-paste the metadata of updated seurat object to subcluster's metadata
      #Plotting and visualize the DEAlgo score of each source of artifacts on the subcluster's UMAP,
      p1 <- FeaturePlot(recluster, features = colnames(CCS_Matrix)[!colnames(CCS_Matrix) %in% c("MajorSub", i)])
      ggsave(filename = paste0(folder_path_Step2_Output,OverlapRatio,"_cluster_",i,"_ContaminationPattern_SourceofArtifacts.png"), p1, height = 10, width = 10, dpi = 300)
    }
  }
}
obj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatureMetabolism_Step2/Output_Overlap_Ratio_0.5/DEAlgoResult.rds"))
PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = FALSE)

###Fig4B_ConditionGSEA.R
singletlevel = 6
dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatureMetabolism_Step2/Output_Overlap_Ratio_0.5/DEAlgoResult.rds"))
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

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$dealgocontamclus <- ifelse(dealgoseuobj$DEAlgo_Contaminated == "Artifact", dealgoseuobj$DEAlgo_ClusterID, NA)

for (i in c("M0","M2","M6","M10")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatureMetabolism_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
  recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
  recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
  recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
  recluster$DEAlgo_Contaminated <- dealgoseuobj$DEAlgo_Contaminated
  recluster$dealgocontamclus <- dealgoseuobj$dealgocontamclus 
  
  p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
  p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
  p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE,  sizes.highlight = 0.1)
  p5 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE,  sizes.highlight = 0.1)
  
  p3 <- DimPlot(recluster, group.by = "dealgocontamclus")
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4_cluster_",i,"_ID.png"), p7+p5+p8+p4+p3, height = 10/3*2, width = 17, dpi = 300)
  
}


DEAlgoSinglet <- subset(dealgoseuobj,subset = DEAlgo_Contaminated == "Singlet")


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

for (i in c("M0","M2","M6","M10")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/AdiposeNatureMetabolism_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
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
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4_cluster_",i,"_rawscore_Validate.png"), p7+p5+p8+p4+p0+p1+p2+p3, height = 10, width = 17, dpi = 300)
  
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

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig4/","Fig4_Azimuth.png"), p12 + p3 + p11 + p7 + p5 +p0 + p1, height = 20, width = 23, dpi = 300)

saveRDS(obj_azimuth,"~/DEAlgoManuscript/Manuscript_Figures/Fig4/Fig4_Azimuth.rds")

#Compare DiaSAT vs NonSAT
#Batch Script start Rscript ~/Script/script/8_Fig4_published_batch1.R, "M0","M2","M6","M10"