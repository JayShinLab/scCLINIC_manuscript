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
      dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))
    
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
      
      widths <- c(7, rep(5, length(plot_list) - 1))
      
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
obj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Output_Overlap_Ratio_0.5/","DEAlgoResult.rds"))
PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = FALSE)



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

p1 <- DimPlot(dealgoseuobj,group.by = "Status",cols = c(
  "#e51c23",  # Unique color 3 (changed)
  "#1e88e5"
))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/Condition.png"), p1, height = 5, width = 7, dpi = 300)

covid19natmed.updated <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemed.rds")
dealgoseuobj$cell.type <- covid19natmed.updated$cell.type 
dealgoseuobj$cell.type.fine <- covid19natmed.updated$cell.type.fine 
dealgoseuobj$cell.type.coarse <- covid19natmed.updated$cell.type.coarse 

p1 <- DimPlot(dealgoseuobj,group.by = "cell.type.coarse",cols = c(
  "#8B0000",
  "#71b84e",
  "#7c4ccb",
  "#c4944a",
  "black",
  "#597b4a",
  "#d55046",
  "grey",
  "#4e336c",
  "#b488bb",
  "#8d8d0f",  # Unique color 1
  "#3b8d89",  # Unique color 2
  "#e51c23",  # Unique color 3 (changed)
  "#1e88e5",  # Unique color 4 (changed)
  "#ff9800"   # Unique color 5 (changed)  # Unique color 5
))
p1 <- DimPlot(dealgoseuobj,group.by = "cell.type",cols = c(
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
  "#e51c23",
  "#ffa07a",
  "#8d8d0f",  # Unique color 1
  "#3b8d89",  # Unique color 2
  "grey",  # Unique color 3 (changed)
  "#1e88e5",  # Unique color 4 (changed)
  "#ff9800" ,  # Unique color 5 (changed)  # Unique color 5
  "black" ,# Light Salmon
  "#20b2aa" ,# Light Sea Green
  "#dda0dd" # Plum
  
))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/cell.type.png"), p1, height = 5, width = 7, dpi = 300)


p1 <- DimPlot(dealgoseuobj,group.by = "Overlap_Ratio_0.5",cols = c(
  "#1e88e5",
  "#71b84e",
  "#7c4ccb",
  "#c4944a",
  "black",
  "#597b4a",
  "#d55046",
  "grey",
  "#4e336c",
  "#b488bb",
  "#8d8d0f",  # Unique color 1
  "#3b8d89",  # Unique color 2
  "#e51c23",  # Unique color 3 (changed)
  "#dda0dd",  # Unique color 4 (changed)
  "#ff9800"   # Unique color 5 (changed)  # Unique color 5
))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/globalumap.png"), p1, height = 5, width = 7, dpi = 300)


#Healthy vs Control
OriginalALL <- dealgoseuobj

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
DEAlgoSinglet <- subset(dealgoseuobj,subset = DEAlgo_Contaminated == "Singlet")

DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$dealgocontamclus <- ifelse(dealgoseuobj$DEAlgo_Contaminated == "Artifact", dealgoseuobj$DEAlgo_ClusterID, NA)

for (i in c("M20","M6","M16","M13","M22")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
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
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig5_cluster_",i,"_ID.png"), p7+p5+p8+p4+p3, height = 10/3*2, width = 17, dpi = 300)
  
}

p1 <- DimPlot(dealgoseuobj,group.by = "DEAlgo_Contaminated")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/scCLINICresult.png"), p1, height = 5, width = 7, dpi = 300)


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

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/COVID.png"), p1+p2+p3, height = 5, width = 17, dpi = 300)

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/Healthy.png"), p4+p5+p6, height = 5, width = 17, dpi = 300)

#Azimuth
obj <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/step1d.rds")
obj_azimuth <- Azimuth::RunAzimuth(obj,reference = "pbmcref")

for (i in c("M6","M16","M20","M13","M22")){
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






###For Poster
#read doublet
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

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
DEAlgoSinglet <- subset(dealgoseuobj,subset = DEAlgo_Contaminated == "Singlet")

DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$dealgocontamclus <- ifelse(dealgoseuobj$DEAlgo_Contaminated == "Artifact", dealgoseuobj$DEAlgo_ClusterID, NA)

#T NK cells
DEAlgoContam1 <- c("M17S5","M13S3","M20S0","M21S3","M22S2","M23S4")#Stroma (M8), Endo (M5), Dist.Prox.Tubule (M4), Coll.Duct.IC (M1), Loop of Henle (M6)

Contam1Cell <- subset(dealgoseuobj,subset = DEAlgo_ClusterID %in% DEAlgoContam1)
noncontam1 <- subset(dealgoseuobj,subset = L0 %in% c("M17","M13","M20","M21","M22","M23"))

negative <- noncontam1[,!(noncontam1$DEAlgo_ClusterID %in% DEAlgoContam1)]
Contam1CellDF <- subset(Contam1Cell,subset = DFafqc == "Doublet")
Contam1CellNDF <- subset(Contam1Cell,subset = DFafqc == "Singlet")


noncontam1$grouping <- NA
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(negative), "Singlet",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellDF), "DF+scCLINIC",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellNDF), "scCLINIC",noncontam1$grouping)


# expall <- AverageExpression(noncontam1,group.by = "grouping",slot = "data")
#pseudo_ifnb <- AggregateExpression(noncontam1, assays = "RNA", group.by = "grouping")#IF return.seurat = T mean normalized again? so should I
#pseudodata <- pseudo_ifnb@assays$RNA$data

library(Seurat)
negativeexp <- AverageExpression(negative,group.by = "all",slot = "data")
Contam1CellDFexp <- AverageExpression(Contam1CellDF,group.by = "all",slot = "data")
Contam1CellNDFexp <- AverageExpression(Contam1CellNDF,group.by = "all",slot = "data")
Contam1CellNDFexp_split <- AverageExpression(Contam1CellNDF,group.by = "DEAlgo_ClusterID",slot = "data")

contam1gene <- read.csv("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Output_Overlap_Ratio_0.5/overlaplst_filtered_Overlap_Ratio_0.5_ContaminationInfo.csv",na.strings = c("", "NA"))
contam1gene <- contam1gene[contam1gene$ES_local >0 ,]
contam1genelist <- unique(contam1gene[contam1gene$cluster_ref=="M0",]$gene)

intersection <- contam1genelist

negative_res <- negativeexp$RNA[intersection,]
Contam1CellDFexp_res <- Contam1CellDFexp$RNA[intersection,]
Contam1CellNDFexp_res <- Contam1CellNDFexp$RNA[intersection,]

library(ggplot2)
library(ggrepel)
# Create data frames for each condition
data_negative <- data.frame(Gene = names(negative_res), Value = negative_res)
data2 <- data.frame(Gene = names(Contam1CellDFexp_res), Value = Contam1CellDFexp_res)
data3 <- data.frame(Gene = names(Contam1CellNDFexp_res), Value = Contam1CellNDFexp_res)

# Combine all data frames
all_data <- rbind(data_negative, data3, data2)
all_data$Condition <- factor(rep(c("Singlet", "scCLINIC", "DF+scCLINIC"), each = length(names(negative_res))),levels = c("DF+scCLINIC", "scCLINIC","Singlet"))
# all_data <- rbind(all_data,reshaped_df)
# Create the violin plot with scatter plot
partial_contamination_data <- subset(all_data, Condition == "scCLINIC" & Value > 7.5)

p1 <-ggplot(all_data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(color = alpha("black",0.25)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  geom_point(aes(color = Condition), size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_text_repel(data = partial_contamination_data, aes(label = Gene), 
                  hjust = 0.5, vjust = 1, segment.color = "transparent", size = 4,position = position_jitter(width = 0.2))+
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linetype = "dashed", color = alpha("black",0.5), size = 0.5, position = position_dodge(width = 0.75)) + # Add average line
  
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "T NK Cell Artifact Average Expression Level")+
  
  theme(legend.position = "none") # Remove legend



ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5A.png"), p1, height = 5, width = 5, dpi = 300)

# Perform pairwise t-tests
pairwise_result <- pairwise.t.test(all_data$Value, all_data$Condition)

write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5A.csv"),sep = ",")

#Proved is not Doublet
p <- VlnPlot(noncontam1,features = "nCount_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nCount_RNA") +
  theme(legend.position = "none")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Ancount.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nCount_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Ancount.csv"),sep = ",")

p <- VlnPlot(noncontam1,features = "nFeature_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nFeature_RNA") +
  theme(legend.position = "none") # Remove legend
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Anfeature.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nFeature_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Anfeature.csv"),sep = ",")


#Monocyte M2 cells
DEAlgoContam1 <- c("M16S2","M16S4","M13S1","M22S2")#Stroma (M8), Endo (M5), Dist.Prox.Tubule (M4), Coll.Duct.IC (M1), Loop of Henle (M6)

Contam1Cell <- subset(dealgoseuobj,subset = DEAlgo_ClusterID %in% DEAlgoContam1)
noncontam1 <- subset(dealgoseuobj,subset = L0 %in% c("M16","M13","M22"))

negative <- noncontam1[,!(noncontam1$DEAlgo_ClusterID %in% DEAlgoContam1)]
Contam1CellDF <- subset(Contam1Cell,subset = DFafqc == "Doublet")
Contam1CellNDF <- subset(Contam1Cell,subset = DFafqc == "Singlet")


noncontam1$grouping <- NA
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(negative), "Singlet",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellDF), "DF+scCLINIC",noncontam1$grouping)
noncontam1$grouping <- ifelse(colnames(noncontam1) %in% colnames(Contam1CellNDF), "scCLINIC",noncontam1$grouping)


# expall <- AverageExpression(noncontam1,group.by = "grouping",slot = "data")
#pseudo_ifnb <- AggregateExpression(noncontam1, assays = "RNA", group.by = "grouping")#IF return.seurat = T mean normalized again? so should I
#pseudodata <- pseudo_ifnb@assays$RNA$data

library(Seurat)
negativeexp <- AverageExpression(negative,group.by = "all",slot = "data")
Contam1CellDFexp <- AverageExpression(Contam1CellDF,group.by = "all",slot = "data")
Contam1CellNDFexp <- AverageExpression(Contam1CellNDF,group.by = "all",slot = "data")
Contam1CellNDFexp_split <- AverageExpression(Contam1CellNDF,group.by = "DEAlgo_ClusterID",slot = "data")

contam1gene <- read.csv("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Output_Overlap_Ratio_0.5/overlaplst_filtered_Overlap_Ratio_0.5_ContaminationInfo.csv",na.strings = c("", "NA"))
contam1gene <- contam1gene[contam1gene$ES_local >0,]
contam1genelist <- unique(contam1gene[contam1gene$cluster_ref=="M2",]$gene)

intersection <- contam1genelist

negative_res <- negativeexp$RNA[intersection,]
Contam1CellDFexp_res <- Contam1CellDFexp$RNA[intersection,]
Contam1CellNDFexp_res <- Contam1CellNDFexp$RNA[intersection,]

library(ggplot2)
library(ggrepel)
# Create data frames for each condition
data_negative <- data.frame(Gene = names(negative_res), Value = negative_res)
data2 <- data.frame(Gene = names(Contam1CellDFexp_res), Value = Contam1CellDFexp_res)
data3 <- data.frame(Gene = names(Contam1CellNDFexp_res), Value = Contam1CellNDFexp_res)

# Combine all data frames
all_data <- rbind(data_negative, data3, data2)
all_data$Condition <- factor(rep(c("Singlet", "scCLINIC", "DF+scCLINIC"), each = length(names(negative_res))),levels = c("DF+scCLINIC", "scCLINIC","Singlet"))
# all_data <- rbind(all_data,reshaped_df)
# Create the violin plot with scatter plot
partial_contamination_data <- subset(all_data, Condition == "scCLINIC" & Value > 12)

library(ggbreak)
library(ggplot2)
library(ggrepel)
library(ggbreak)
# 
# p1 <- ggplot(all_data, aes(x = Condition, y = Value, fill = Condition)) +
#   geom_violin(color = alpha("black", 0.25)) +
#   geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white", outlier.shape = NA) + # Add white boxplot within each violin
#   geom_point(aes(color = Condition), size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
#   geom_text_repel(data = partial_contamination_data, aes(label = Gene), 
#                   hjust = 0.5, vjust = 1, segment.color = "transparent", size = 4, position = position_jitter(width = 0.2)) +
#   stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
#   stat_summary(fun = mean, geom = "line", aes(group = 1), linetype = "dashed", color = alpha("black", 0.5), size = 0.5, position = position_dodge(width = 0.75)) + # Add average line
#   theme_minimal(base_size = 14) +
#   scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
#   scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
#   labs(title = "",
#        x = "",
#        y = "CD14 Monocyte Artifact Average Expression Level") +
#   theme(legend.position = "none") + # Remove legend
#   scale_y_break(c(60, 90), scales = 0.5) # Set y-axis break between 60 and 90
# 
# ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/", "SupplementaryFig5B.png"), p1, height = 5, width = 5, dpi = 300)


p1 <-ggplot(all_data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(color = alpha("black",0.25)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  geom_point(aes(color = Condition), size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_text_repel(data = partial_contamination_data, aes(label = Gene), 
                  hjust = 0.5, vjust = 1, segment.color = "transparent", size = 4,position = position_jitter(width = 0.2))+
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linetype = "dashed", color = alpha("black",0.5), size = 0.5, position = position_dodge(width = 0.75)) + # Add average line
  
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "CD14 Monocyte Artifact Average Expression Level")+
  
  theme(legend.position = "none") 


ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5B.png"), p1, height = 5, width = 5, dpi = 300)

p1 <-ggplot(all_data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(color = alpha("black",0.25)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  geom_point(aes(color = Condition), size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_text_repel(data = partial_contamination_data, aes(label = Gene), 
                  hjust = 0.5, vjust = 1, segment.color = "transparent", size = 4,position = position_jitter(width = 0.2))+
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linetype = "dashed", color = alpha("black",0.5), size = 0.5, position = position_dodge(width = 0.75)) + # Add average line
  
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "CD14 Monocyte Artifact Average Expression Level")+
  
  theme(legend.position = "none") 

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5B_nolimit.png"), p1, height = 5, width = 5, dpi = 300)

# Perform pairwise t-tests
pairwise_result <- pairwise.t.test(all_data$Value, all_data$Condition)

write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5B.csv"),sep = ",")

#Proved is not Doublet
p <- VlnPlot(noncontam1,features = "nCount_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nCount_RNA") +
  theme(legend.position = "none")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Bncount.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nCount_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Bncount.csv"),sep = ",")

p <- VlnPlot(noncontam1,features = "nFeature_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#69a75f", "#9671c3", "#be883d")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nFeature_RNA") +
  theme(legend.position = "none") # Remove legend
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Bnfeature.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(noncontam1$nFeature_RNA, noncontam1$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Bnfeature.csv"),sep = ",")

#Monocyte gene
library(gridExtra)
celltypelist <- c("M2","M16","M13","M22")
for (markertosee in c("S100A9","S100A8","FCN1","LYZ","VCAN")){
  plot_list <- list()
  for (curcelltype in celltypelist){
    
    recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Overlap_Ratio_0.5_recluster/","Overlap_Ratio_0.5_cluster_",curcelltype,".rds"))
    
    p1 <- FeaturePlot(recluster, reduction = "umap",features = markertosee, pt.size = 0.1) & scale_color_gradientn(colors = c("grey","red"),limits = c(0, 6), oob = scales::squish)#,limits = c(0, 5)DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_Contaminated",raster=FALSE, label = T, sizes.highlight = 0.1)
    #ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/",curcelltype,"_marker_",markertosee,".png"), p1, height = 5, width = 5, dpi = 300)
    plot_list[[curcelltype]] <- p1
  }
  
  all_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = ceiling(length(celltypelist)/2))
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/",markertosee,".png"), plot = all_plots, height = 5, width = 5, dpi = 300)
}



#B cells supplementary Fig5
obj_azimuth <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/Fig5_Azimuth.rds")

singletlevel = 7

dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Output_Overlap_Ratio_0.5/","DEAlgoResult.rds"))
doublet_afqc <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC.rds"))

metatable_unique<- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemed.rds"))

dealgoseuobj@meta.data <- cbind(metatable_unique@meta.data,dealgoseuobj@meta.data)

#16.014% of 44721 is 7162
top_indices <- order(doublet_afqc@meta.data$DFScore, decreasing = TRUE)[7162]
cutoff <- doublet_afqc@meta.data$DFScore[top_indices]

dealgoseuobj$DFafqc_pANN <- doublet_afqc$DFScore
dealgoseuobj$DFafqc <- NA
dealgoseuobj$DFafqc <- ifelse(doublet_afqc$DFScore < cutoff, "Singlet", "Doublet")

dealgoseuobj$DEAlgo_Contaminated <- ifelse(as.numeric(levels(dealgoseuobj$DEAlgocluster_Contam))[dealgoseuobj$DEAlgocluster_Contam] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$dealgocontamclus <- ifelse(dealgoseuobj$DEAlgo_Contaminated == "Artifact", dealgoseuobj$DEAlgo_ClusterID, NA)



i <- "M6"
file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Overlap_Ratio_0.5_recluster/",file_name))

recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
recluster$DFafqc <- dealgoseuobj$DFafqc
recluster$DEAlgo_ClusterID <- dealgoseuobj$DEAlgo_ClusterID
recluster$DiffDP1_sum <- dealgoseuobj$DiffDP1_sum
recluster$DEAlgocluster_Contam <- dealgoseuobj$DEAlgocluster_Contam
recluster$DEAlgo_Contaminated <- dealgoseuobj$DEAlgo_Contaminated

recluster$dealgocontamclus <- dealgoseuobj$dealgocontamclus

recluster$azimuthl1 <- obj_azimuth$predicted.celltype.l1
recluster$azimuthl2 <- obj_azimuth$predicted.celltype.l2

p7 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgo_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
p8 <- FeaturePlot(recluster, reduction = "umap",features = "DiffDP1_sum")
p4 <- DimPlot(recluster, reduction = "umap",group.by = "DEAlgocluster_Contam",raster=FALSE,  sizes.highlight = 0.1)
p5 <- DimPlot(recluster, reduction = "umap",group.by = "dealgocontamclus",raster=FALSE,  sizes.highlight = 0.1)
p1 <- FeaturePlot(recluster, reduction = "umap",features = "DFafqc_pANN")
p0 <- DimPlot(recluster, reduction = "umap",group.by = "DFafqc",raster=FALSE,  sizes.highlight = 0.1, cols = c("grey"))

p2 <- DimPlot(recluster, group.by = "azimuthl1", cols = c(
  "#3b8d89",
  "#ff9800",
  "#e51c23",
  "#c4944a",
  "#cf4ba6",
  "#d55046",
  "#4e336c"))

p3 <- DimPlot(recluster, group.by = "azimuthl2")

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig5_cluster_",i,"_rawscore_Validate.png"), p7+p5+p8+p4+p0+p1+p2+p3, height = 10, width = 12, dpi = 300)
  

