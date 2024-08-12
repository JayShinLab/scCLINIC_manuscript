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

Name <- "PBMCCOVIDNatMed"
Input <- "~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC.rds"
Output <- "~/DEAlgoManuscript/Manuscript_Figures/Fig5/"
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

Outdir <- Output
STEP2B_ContaminationScore<-function(obj,Outdir,Name,resol=0.8,OverlapRatio=0.5,gene_n=150,CELLANNOTATION = FALSE){
  message("Step2B started.")
  folder_path_Step1 <- paste0(Outdir,Name,"_Step1/")
  ##Step 2###################################################################Creating Output Folders
  if (CELLANNOTATION){
    message("Using user-annotated clusters.")
    obj@meta.data$annotation_index <- as.numeric(factor(obj@meta.data[,OverlapRatio]))#change the manual annotation to index, eg. CellType0 CellType1 CellType2 CellType3 to 1 2 3 4
    OverlapRatio <- "annotation_index" #User manual cellannotation
  }else{
    OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
    message(paste0("Using ",OverlapRatio))
  }
  clus.names <- obj@meta.data[[OverlapRatio]]
  if (any(!(grepl("^[a-zA-Z]|^\\.[^0-9]", clus.names)))) {
    clus.names <- ifelse(
      !(grepl("^[a-zA-Z]|^\\.[^0-9]", clus.names)),
      paste0("M", clus.names),
      clus.names
    )
    obj@meta.data[[OverlapRatio]] <- clus.names
    message = paste0("Appending `M` to cluster names to ensure valid variable names")
  }
  
  folder_path_Step2 <- paste0(Outdir,Name,"_Step2/")
  folder_path_Step2_L1R <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/")
  folder_path_Step2_L1R_Marker <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/Marker/")
  folder_path_Step2_Output <- paste0(Outdir,Name,"_Step2/Output_",OverlapRatio,"/")
  
  ###scCLINIC
  totalcells <- ncol(obj)
  avegloballst <- AverageExpression(obj,group.by = OverlapRatio,slot = "counts")#avg expression of each major clusters
  gcell_table <- table(obj@meta.data[,OverlapRatio])#number of cells of each major clusters
  qcsclst <- sort(unique(obj@meta.data[,OverlapRatio]))#list of major clusters seurat ID
  cluster_to_consider <- list()
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".csv")
    if (file_name %in% list.files(folder_path_Step2_L1R_Marker)){
      cluster_to_consider <- unlist(c(cluster_to_consider,i))
    }
  }
  
  #Load the top150 global markers
  n_topgene <- gene_n#################################ONLY DEPENDENT
  Global.markers <- read.table(paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,"_GlobalMarker.csv"), sep= ",", header = T, row.names = 1)
  #Global.markers <- Global.markers[Global.markers$p_val_adj <= 0.05,]
  Global.markers <- Global.markers[Global.markers$avg_log2FC > 0,]
  Global.markers <- Global.markers[Global.markers$pct.1 > Global.markers$pct.2,]
  Global.markers$ES <- Global.markers$avg_log2FC*(Global.markers$pct.1-Global.markers$pct.2)
  Global.markers <- Global.markers[order(Global.markers$ES, decreasing = TRUE), ]
  Global.markers <- Global.markers[Global.markers$cluster %in% cluster_to_consider,]
  
  Global.markers %>%
    group_by(cluster) %>%
    top_n(n = n_topgene, wt = ES) %>%
    ungroup()-> Global_Top20
  
  #Identify and calculate the Identity Score of the Artifacts in each Subclusters and calculate the scCLINIC score for each subclusters
  overlaplst <- list()
  
  for (i in cluster_to_consider){#For each major cluster (i)
    message(paste0("Calculating Enrichment Score for ",i))
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (paste0(sub("\\.rds$", "", file_name),".csv") %in% list.files(folder_path_Step2_L1R_Marker)){ #Load each subcluster markers
      sc.markers <- read.table(paste0(folder_path_Step2_L1R_Marker,sub("\\.rds$", "", file_name),".csv"), sep = ",", header = T, row.names = 1)#readRDS(paste0(folder_path_Step2_L1R_Marker,file_name))
      subc <- readRDS(paste0(folder_path_Step2_L1R,OverlapRatio,"_cluster_",i,".rds"))
      avegenelst <- AverageExpression(object = subc, group.by = "scCLINIC_subcluster",slot = "counts")#Avg marker expression of each subclusters
      ncell_table <- table(subc@meta.data$"scCLINIC_subcluster")#Number of cells of each subclusters
      
      #Prefilter the findallmarker result, only keep positive markers with pct1 > pct2, ensure no negative ES
      if (length(sc.markers)!=0){
        sc.markers$ES <- sc.markers$avg_log2FC*(sc.markers$pct.1-sc.markers$pct.2) #Calculate ES score for each markers
        sc.markers$ES <- ifelse(sc.markers$avg_log2FC < 0, 0, sc.markers$ES) #Set marker ES to 0, if marker avglog2FC < 0
        sc.markers$ES <- ifelse(sc.markers$pct.1 < sc.markers$pct.2, 0, sc.markers$ES) #Set marker ES to 0, if marker pct.1 - pct.2 < 0
        
        sc.markers <- sc.markers[order(sc.markers$ES, decreasing = TRUE), ]
        
        localsclst <- unique(sc.markers$cluster)
        for (j in sort(localsclst)){#For each subcluster (j) in major cluster (i)
          OwnGene <- Global_Top20[Global_Top20$cluster %in% i,]#Global marker which are present in major cluster i
          Global_sub_Top20 <- anti_join(Global_Top20, OwnGene, by = "gene")# Global marker genes which are not present in major cluster i
          
          sc_sub_Top20 <- sc.markers[sc.markers$cluster %in% j,]# Subcluster marker (j)
          
          overlap <- sc_sub_Top20[sc_sub_Top20$gene %in% Global_sub_Top20$gene,]#Subcluster marker (j) which also present in other major clusters
          
          if (nrow(overlap)>0){
            overlap$globalcluster <- i #global cluster = major cluster i
            contam_global <- Global_sub_Top20[Global_sub_Top20$gene %in% sc_sub_Top20$gene,] #extract the global marker information, including avglog2FC, pct1, pct2
            merged_df <- merge(overlap, contam_global, by = "gene", suffixes = c("_local", "_ref")) #merge both information, _local = subclusters marker information, _ref = global marker information
            merged_df$ncells <- ncell_table[[j]] #number of cells in subcluster j
            merged_df$refcells <- gcell_table[[i]] #number of cells in that major cluster where the global marker present
            if (length(merged_df$gene) == 1){
              merged_df$aveexp_local <- avegenelst$RNA[merged_df$gene,as.character(merged_df$cluster_local)] #avg expression of that global marker in that major cluster where the global marker present
              merged_df$aveexp_ref <- avegloballst$RNA[merged_df$gene,as.character(merged_df$cluster_ref)]  #avg expression of that subcluster marker in subcluster j
            }
            else{
              merged_df$aveexp_local <- avegenelst$RNA[merged_df$gene,as.character(j)]
              merged_df$aveexp_ref <- diag(as.matrix(avegloballst$RNA[merged_df$gene,as.character(merged_df$cluster_ref)]))#as.character(merged_df$cluster_ref)
            }
            overlaplst <- rbind(overlaplst,merged_df)
          }
        }
      }}}
  message("Calculating Contamination Scores.")
  #Notes: Subcluster (j) markers which not present in major cluster (i) aka artifacts
  overlaplst.filtered <- overlaplst
  overlaplst.filtered$cluster_ref <- overlaplst.filtered$cluster_ref
  overlaplst.filtered$cluster_local <- overlaplst.filtered$cluster_local
  overlaplst.filtered$globalcluster <- overlaplst.filtered$globalcluster
  overlaplst.filtered$dp1 <- overlaplst.filtered$ES_local  #Identity Score of the artifact aka ES_local
  
  write.table(overlaplst.filtered,file = paste0(folder_path_Step2_Output,"ArtifactsInfo.csv"),sep = ",")
  
  referenceoverlap <- overlaplst.filtered[, c("dp1","cluster_local" , "globalcluster", "cluster_ref")]
  overlaplst.filtered <-  overlaplst.filtered[, c("gene","dp1","p_val_local", "avg_log2FC_local", "pct.1_local", "pct.2_local", "p_val_adj_local", "cluster_local" ,   "ES_local", "globalcluster", "ncells", "aveexp_local")]
  
  #Remove duplicated artifacts contributed by different major clusters, eg. Artifact A present across major cluster 1 2 3, in overlaplst.filtered appeared three times...
  dup_rows <- duplicated(overlaplst.filtered)
  overlaplst.filtered <- overlaplst.filtered[!dup_rows, ]
  
  #Summarize the artifacts per subclusters, calculate the scCLINIC score for each subclusters
  contamdfall <- overlaplst.filtered %>%
    group_by(cluster_local, globalcluster) %>%
    summarise(count = n(),mean(dp1),mean(ncells))#sum the Identity Score for each artifacts in the subcluster aka Subcluster scCLINIC Score, count the number of unique artifacts in the subcluster
  
  contamdfall$cluster_local <- as.matrix(contamdfall)[,"cluster_local"]
  #contamdfall$cluster_local <- as.integer(contamdfall$cluster_local)
  #contamdfall$globalcluster <- as.integer(contamdfall$globalcluster)#Summary including the Subclusters ID (globalcluster & cluster_local) and the respective scCLINIC Score (mean(dp1))
  
  #Classified subclusters to either contaminated or non-contaminated subclusters based on elbow point
  contamdfall$scCLINICScore <- contamdfall$`mean(dp1)`
  contamdfall_backup <- contamdfall
  contamdfall <- contamdfall[order(contamdfall$"scCLINICScore", decreasing = TRUE), ]
  
  #write.table(contamdfall,file = paste0(folder_path_Step2_Output,"ArtifactsInfo.csv"),sep = ",")
  
  x1 <- as.numeric(rownames(contamdfall))#Subclusters instances
  y1 <- contamdfall$"scCLINICScore"#scCLINIC Score of respective subclusters
  z1 <- paste0(contamdfall$globalcluster,"_",contamdfall$cluster_local)#Respective subclusters ID, eg. M1S2 = Subcluster 2 in Major Cluster 1
  ncell1 <- contamdfall$`mean(ncells)`#no. of Cells of respective subclusters
  
  ycopy <- y1
  # Step 3: Calculate the first and second derivatives
  # y1 = y1#sgolayfilt(y1, p = 3, n = 5) #Optional Smoothing
  # first_derivative <- diff(y1) / diff(x1)
  # second_derivative <- diff(first_derivative) / diff(x1[-length(x1)])
  # # Step 5: Identify local maxima of the second derivative
  # local_maxima <- findpeaks(second_derivative, minpeakheight = 0.01, minpeakdistance = 5)
  # maxima_x <- sort(x1[local_maxima[,2]],decreasing = F)
  # maxima_y <- sort(y1[local_maxima[,2]] ,decreasing = T)
  
  #Maxima
  local_maxima <- findpeaks(diff(diff(y1)), npeaks = 7,minpeakdistance = 5, sortstr = TRUE)
  maxima_x <- sort(local_maxima[,2], decreasing = F)
  maxima_y <- sort(y1[local_maxima[,2]] ,decreasing = T)
  
  #write.table(local_maxima,file = paste0(folder_path_Step2_Output,"second_derivative.csv"),sep = ",")
  
  #Optimal maxima
  optimal_maxima <- findpeaks(diff(diff(diff(y1))), npeaks = 7,minpeakdistance = 5, sortstr = TRUE)
  opmaxima_x <- optimal_maxima[,2]+1
  #opmaxima_y <- sort(y1[optimal_maxima[,2]+1] ,decreasing = T)
  
  # Function to find the closest point in maxima_x for a given value
  find_closest <- function(value, maxima_x) {
    # Compute the absolute differences
    differences <- abs(maxima_x - value)
    # Find the index of the minimum difference
    closest_index <- which.min(differences)
    # Return the closest value
    return(maxima_x[closest_index])
  }
  
  # Apply the function to each element in opmaxima_x
  closest_points <- sapply(opmaxima_x, find_closest, maxima_x = maxima_x)
  #write.table(optimal_maxima,file = paste0(folder_path_Step2_Output,"third_derivative.csv"),sep = ",")
  
  #Plot scCLINIC Score vs Percentage Droplets (%)
  y1 <- ycopy #Plot using raw score not smoothing (if smoothing is applied)
  data <- data.frame(x = x1, y = y1,z = z1,ncelllst = cumsum(ncell1)/totalcells*100)#ncelllst = accumulative percentage droplets of each subclusters
  # Calculate the scaling factor
  range_y1 <- range(data$y)
  range_ncell1_accumulated <- range(data$ncelllst)
  scaling_factor <- diff(range_y1) / diff(range_ncell1_accumulated)
  
  # Base plot
  p <- ggplot(data, aes(x = x1, y = y1, label=z1)) +
    geom_point() +
    geom_line(aes(y = y), col = "#53a0b9") +
    geom_line(aes(y = ncelllst * scaling_factor), col = "#d55046") +
    geom_text(aes(label = z), angle = 90, hjust = 0.5, vjust = -0.5) +
    labs(title = paste0("Contamination Score (no. of cells = ", totalcells, ")"),
         x = "Subclusters (ID)", y = "Contamination Score") +
    scale_y_continuous(
      name = "Contamination Score",
      sec.axis = sec_axis(~./scaling_factor, name = "Percentage Droplets (%)")
    ) +
    theme_minimal() +
    theme(
      axis.line.y = element_line(color = "#53a0b9"),
      axis.text.y = element_text(color = "#53a0b9"),
      axis.line.y.right = element_line(color = "#d55046"),
      axis.text.y.right = element_text(color = "#d55046"),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(color = "black")
    )+
    annotate("text", x =  x1[length(x1)], y = y1[1], label = length(maxima_x)+1, hjust = 1, col = "black")
  
  
  # Create list to store additional layers
  additional_layers <- list()
  
  nelbowidx <- 0
  # Loop to create vertical and horizontal lines and annotations
  for (nelbow in maxima_x) {
    nelbowidx <- nelbowidx + 1
    nelbow_y <- y1[nelbow]
    additional_layers <- c(
      additional_layers,
      geom_vline(xintercept = nelbow, linetype = "dashed", color = "black"),
      #geom_hline(yintercept = nelbow_y, linetype = "dashed", color = "#53a0b9"),
      #geom_hline(yintercept = data$ncelllst[nelbow] * scaling_factor, linetype = "dashed", color = "#d55046"),
      # geom_segment(x = nelbow, y = data$ncelllst[nelbow] * scaling_factor,
      #                  xend = x1[length(x1)], yend = data$ncelllst[nelbow] * scaling_factor,
      #              linetype = "dashed", color = "#d55046"),
      #annotate("text", x = 0, y = nelbow_y, label = round(nelbow_y, 3), vjust = -1, col = "#53a0b9"),
      annotate("text", x = nelbow, y = data$ncelllst[nelbow] * scaling_factor, label = floor(data$ncelllst[nelbow]*1000)/1000, vjust = 1,hjust=-0.1, col = "#d55046"),
      annotate("text", x = nelbow, y = y1[1], label = nelbowidx, hjust = 1.5, col = "black")
    )
  }
  
  # Add the additional layers to the plot
  p <- p + additional_layers +
    geom_vline(xintercept = closest_points[1], linetype = "dashed", color = "#d55046")+
    geom_vline(xintercept = closest_points[2], linetype = "dashed", color = "#7c3b36")
  
  ggsave(filename = paste0(folder_path_Step2_Output,"scCLINICScoreKneePlot",".png"), p, height = 10, width = 10, dpi = 300)
  
  #Summarize the artifacts per subclusters per source of artifacts, calculate the scCLINIC score for each subclusters and each source of artifacts
  summary_clusterref <- referenceoverlap %>% #sum the Identity Score for each artifacts in the subcluster per source of artifacts aka Subcluster scCLINIC Score per source of artifacts, count the number of unique artifacts in the subcluster per source of artifacts
    group_by(cluster_local, globalcluster,cluster_ref) %>%
    summarise(count = n(),mean(dp1))
  
  summary_clusterref$cluster_local <- as.matrix(summary_clusterref)[,"cluster_local"]
  #summary_clusterref$cluster_local <- as.integer(summary_clusterref$cluster_local)
  summary_clusterref$cluster_ref <- as.matrix(summary_clusterref)[,"cluster_ref"]
  #summary_clusterref$cluster_ref <- as.integer(summary_clusterref$cluster_ref)
  #summary_clusterref$globalcluster <- as.integer(summary_clusterref$globalcluster)
  summary_clusterref$DiffDP1 <- summary_clusterref$`mean(dp1)`
  
  #Extend the source of artifacts information to each contaminated subclusters
  GlobalPlusCellSpecific <- inner_join(contamdfall, summary_clusterref,
                                       by = c("cluster_local", "globalcluster"))
  GlobalPlusCellSpecific <- GlobalPlusCellSpecific[,c("cluster_local", "globalcluster","scCLINICScore", "cluster_ref","DiffDP1")]
  
  
  ################################
  #Store scCLINIC Result in Metadata of the original seurat object
  #scCLINIC Score (scCLINICScore), Major Cluster ID (L0), Subcluster ID (L1),
  #scCLINIC level (scCLINIC_Level)
  message("Saving results in Seurat Object.")
  obj@meta.data$"scCLINICScore" <- 0
  obj@meta.data$"L0" <- NA
  obj@meta.data$"L1" <- NA
  obj@meta.data$"scCLINIC_Level" <- length(maxima_y)+1
  
  contaminatedinfo <- GlobalPlusCellSpecific
  contaminatedinfo <- contaminatedinfo[order(contaminatedinfo$"DiffDP1", decreasing = TRUE), ]
  
  refinfo <- contaminatedinfo %>% #Summarized the source of artifacts for each subclusters
    group_by(cluster_local, globalcluster) %>%
    summarise(cluster_ref = list(cluster_ref))
  
  contamV <- contamdfall_backup#Including the scCLINIC score for each subclusters (for both contaminated and non-contaminated subclusters)
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))
      
      obj$L0 <- ifelse(colnames(obj) %in% colnames(recluster), i ,obj@meta.data$L0)#Store L0
      
      for (clusterid in unique(recluster$"scCLINIC_subcluster")){
        obj$L1 <- ifelse(colnames(obj) %in% rownames(recluster@meta.data[recluster@meta.data$"scCLINIC_subcluster" == clusterid,]), clusterid ,obj@meta.data$L1)#Store L1
      }
      
      subclusterlst <- contamV[contamV$globalcluster == i,]$cluster_local
      if (length(subclusterlst)>0){
        for (subclusterid in subclusterlst){
          recluster.filterout <- rownames(recluster@meta.data[recluster@meta.data$scCLINIC_subcluster == subclusterid,]) #cell name of each contaminated subclusters
          curDP1_sum <- contamV[contamV$globalcluster == i & contamV$cluster_local == subclusterid,]$scCLINICScore #curDP1_sum is the scCLINIC Score for the current subcluster (subclusterid)
          obj@meta.data[recluster.filterout,]$scCLINICScore <- rep(curDP1_sum, length(recluster.filterout))#Store scCLINIC Score for each subclusters
        }
      }
    }
  }
  
  #Store the scCLINIC_Level for each subclusters
  contam_level = length(maxima_y)
  for (maxima_y_cur in sort(maxima_y, decreasing = FALSE)){
    obj$scCLINIC_Level <- ifelse(obj$scCLINICScore >= maxima_y_cur,contam_level, obj$scCLINIC_Level)
    contam_level = contam_level - 1
  }
  obj@meta.data$"scCLINIC_ClusterID" <- paste0(obj@meta.data$"L0","_",obj@meta.data$"L1")#Store scCLINIC_ClusterID
  
  #Plotting and visualize the scCLINIC results on the UMAP,
  palette15 <- c(
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
    "#8d8d0f",
    "#3b8d89",
    "#e51c23",
    "#1e88e5",
    "#ff9800"
  )
  obj$scCLINIC_Level <- factor(obj$scCLINIC_Level, levels = sort(unique(obj$scCLINIC_Level),decreasing=F))
  col <- setNames(palette15, levels(obj$scCLINIC_Level))
  
  p4 <- DimPlot(obj, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,pt.size = 0.1,cols = col)#rev(viridis(length(unique(obj$scCLINIC_Level)))))
  p6 <- DimPlot(obj, reduction = "umap",group.by = OverlapRatio,raster=FALSE,pt.size = 0.1)
  p8 <- FeaturePlot(obj, reduction = "umap",features = "scCLINICScore",pt.size = 0.1)
  ggsave(filename = paste0(folder_path_Step2_Output,"scCLINICResult.png"), p6+p4+p8, height = 5, width = 15, dpi = 300)
  
  ##Store score in Subcluster (sync with updated seurat object), OPTIONAL
  for (i in cluster_to_consider){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))
      recluster@meta.data <- obj@meta.data[rownames(recluster@meta.data),]#copy-paste the metadata of updated seurat object to subcluster's metadata
      #Plotting and visualize the scCLINIC results on the subcluster's UMAP,
      p4 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,  sizes.highlight = 0.1,cols=col)
      p7 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
      p8 <- FeaturePlot(recluster, reduction = "umap",features = "scCLINICScore")
      ggsave(filename = paste0(folder_path_Step2_Output,i,"_scCLINICResult.png"), p7+p4+p8, height = 5, width = 15, dpi = 300)
    }
  }
  
  ###Step4 Interpretation
  #Load the contamination information
  file_name <- paste0("ArtifactsInfo.csv")
  curcluster <- read.table(paste0(folder_path_Step2_Output,file_name), sep= ",", header = T, row.names = 1)
  source_of_artifact <- unique(curcluster$cluster_ref)
  
  #
  artifact_info_lst <- t(data.frame(row.names = c(colnames(curcluster),"multiplehit","globalduplicate","nduplicate")))
  for (source_of_artifact_i in source_of_artifact){#For each source of artifacts (source_of_artifact_i)
    artifact_info <- curcluster[curcluster$cluster_ref== source_of_artifact_i,]#Subset each artifacts from the same source (source_of_artifact_i)
    artifact_info$multiplehit <- duplicated(artifact_info$gene) | duplicated(artifact_info$gene, fromLast = T)# Identify artifacts from the same source which present multiple times (multiplehit = TRUE)
    
    #Investigate the artifacts from the same source which present multiple times
    dupligenetable <- artifact_info[artifact_info$multiplehit==T,]
    artifact_info$globalduplicate <- 1
    artifact_info$nduplicate <- 1
    for (gene_dupli in unique(dupligenetable$gene)){
      ndupli <- length(dupligenetable[dupligenetable$gene == gene_dupli,]$gene)#ndupli: Number of times the artifacts from the same source are present.
      globaldupli <- length(unique(dupligenetable[dupligenetable$gene == gene_dupli,]$globalcluster))#globaldupli: The major clusters of contaminated subclusters where the artifacts are present.
      #Store both information in artifact_info
      artifact_info$nduplicate <- ifelse(artifact_info$gene == gene_dupli, ndupli ,artifact_info$nduplicate)
      artifact_info$globalduplicate <- ifelse(artifact_info$gene == gene_dupli, globaldupli ,artifact_info$globalduplicate)
    }
    artifact_info_lst <- rbind(artifact_info_lst, artifact_info)#rbind each source of artifacts' artifact_info
  }
  
  artifact_info <- artifact_info_lst
  
  #Change the name of the columns
  colnames(artifact_info) <- c(
    
    "Artifact_Gene" ,   "p_val"  ,    "avg_log2FC" ,"pct.1"  ,    "pct.2"  ,
    "p_val_adj",  "Subcluster" ,   "Enrichment_Score"    ,     "Major_Cluster" ,   "SoA_p_val"  ,
    "SoA_avg_log2FC" ,  "SoA_pct.1"  ,      "SoA_pct.2"   ,     "SoA_p_val_adj" ,   "Source_of_Artifacts_SoA" ,
    "SoA_Enrichment_Score"         ,  "No. cell"      ,     "SoA_No. cell"      ,   "Avg_Exp"   ,  "SoA_Avg_Exp"   ,
    "dp1"          ,    "multiplehit"   ,   "globalduplicate" , "nduplicate"
  )
  
  #Remove certain columns
  artifact_info <- artifact_info[, !colnames(artifact_info) %in% c("dp1", "multiplehit", "globalduplicate", "nduplicate")]
  
  # Reorder the columns in artifact_info
  new_order <- c(
    
    "Artifact_Gene",
    
    "Major_Cluster", "Subcluster", "No. cell", "Enrichment_Score",  "avg_log2FC", "pct.1", "pct.2", "p_val","p_val_adj", "Avg_Exp",
    
    "Source_of_Artifacts_SoA","SoA_No. cell", "SoA_Enrichment_Score","SoA_avg_log2FC", "SoA_pct.1", "SoA_pct.2", "SoA_p_val","SoA_p_val_adj","SoA_Avg_Exp"
    
  )
  artifact_info <- artifact_info[, new_order, drop = FALSE]
  
  file_name <- paste0("ArtifactsInfo.csv")
  write.table(artifact_info,file = paste0(folder_path_Step2_Output,file_name),sep = ",")
  
  ###Add the scCLINIC Score for each source of artifacts into the seurat object metadata
  #scCLINIC Subcluster ID
  artifact_info$MajorSub <- paste0(artifact_info$Major_Cluster,"_",artifact_info$Subcluster)
  #Summarize ES score and their source of major cluster for each subclusters
  result <- artifact_info %>%
    group_by(MajorSub) %>%
    summarize(
      cluser_reflst = list(Source_of_Artifacts_SoA),
      dp1lst = list(Enrichment_Score)
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
  
  CCS_Matrix <- as.data.frame(CCS_Matrix)
  CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
  #Rowname = MajorSub and remove the column
  rownames(CCS_Matrix) <- CCS_Matrix$MajorSub
  CCS_Matrix <- CCS_Matrix[, !colnames(CCS_Matrix) %in% "MajorSub"]
  #Add each source of artifacts into metadata as seperate column
  for (colname in colnames(CCS_Matrix)){
    obj@meta.data[[colname]]<-0
  }
  #Assign each source of artifacts' scCLINIC score into metadata for each subclusters
  for (rowccs in rownames(CCS_Matrix)){
    for (colname in colnames(CCS_Matrix)){
      obj@meta.data[[colname]] <- ifelse(obj$scCLINIC_ClusterID == rowccs, CCS_Matrix[rowccs,colname] ,obj@meta.data[[colname]])#Store L0
    }
  }
  #Tabulated subclusters information, including their scCLINIC_ClusterID, scCLINIC_Level, PercentageDroplets, ElbowPoint, n_artifacts, n_cells, scCLINIC score, CCS
  
  scoreinfo <- obj@meta.data[,c(c("scCLINIC_ClusterID","scCLINIC_Level"),colnames(CCS_Matrix))]
  scoreinfo <- scoreinfo[!duplicated(scoreinfo), ]
  scoreinfo <- scoreinfo[complete.cases(scoreinfo), ]
  rownames(scoreinfo) <- NULL
  
  contamdfall$scCLINIC_ClusterID <- paste0(contamdfall$globalcluster,"_",contamdfall$cluster_local)
  contamdfall$PercentageDroplets <- contamdfall$`mean(ncells)`/totalcells*100
  contamdfall <- merge(contamdfall, scoreinfo, by = "scCLINIC_ClusterID", all = TRUE)
  contamdfall$ElbowPoint <- 0
  contam_level <- length(maxima_y)
  for (maxima_y_cur in sort(maxima_y,decreasing = F)){
    contamdfall$ElbowPoint <- ifelse(contamdfall$scCLINIC_Level == contam_level,maxima_y_cur,contamdfall$ElbowPoint)
    contam_level <- contam_level -1
  }
  contamdfall <- transform(contamdfall, n_artifacts = count, n_cells = `mean(ncells)`)
  contamdfall <- subset(contamdfall, select = -c(count  , mean.dp1., mean.ncells.))
  
  #Remove certain columns
  contamdfall <- contamdfall[, !colnames(contamdfall) %in% c("cluster_local", "globalcluster","n_artifacts")]
  
  #Change the name of the columns
  colnames(contamdfall)[colnames(contamdfall) == "n_cells"] <- "No. cell"
  colnames(contamdfall)[colnames(contamdfall) == "PercentageDroplets"] <- "No. cell (%)"
  colnames(contamdfall)[colnames(contamdfall) == "ElbowPoint"] <- "scCLINIC_Elbow_Point"
  # Identify columns with "^M" in their names
  cols_with_M <- grep("\\M", colnames(contamdfall), value = TRUE)
  # Assign new column names to the data frame
  colnames(contamdfall)[colnames(contamdfall) %in% cols_with_M] <- paste0("CCS","_", cols_with_M)
  
  # Reorder the columns in artifact_info
  new_order <- c(
    "scCLINIC_ClusterID" ,  "scCLINICScore" , "No. cell (%)" ,"No. cell",  "scCLINIC_Level", "scCLINIC_Elbow_Point",grep("\\CCS", colnames(contamdfall), value = TRUE)
  )
  
  contamdfall <- contamdfall[, new_order, drop = FALSE]
  
  write.table(contamdfall,file = paste0(folder_path_Step2_Output,"scCLINICScore.csv"),sep = ",",row.names = FALSE)
  
  #Tidy up Metadata (Remove unnecessary columns)
  obj_metadata <- obj@meta.data
  metadata_saved <- c(OverlapRatio,"scCLINICScore","scCLINIC_Level","scCLINIC_ClusterID",grep("^M",colnames(obj_metadata),value = T))
  obj_metadata <- obj_metadata[,metadata_saved]
  # Get the column names
  colnames <- colnames(obj_metadata)
  # Add "CCS_" in front of column names that start with "M"
  colnames(obj_metadata) <- ifelse(grepl("^M", colnames), paste0("CCS_", colnames), colnames)
  obj@meta.data <- obj_metadata
  
  #Save the updated seurat object (including scCLINIC results)
  saveRDS(obj,file = paste0(folder_path_Step2_Output,"scCLINICResult.rds"))
  
  message("Step2B completed.")
  
  return(obj)
}


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


dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Output_Overlap_Ratio_0.5/","scCLINICResult.rds"))
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
  "grey",
  "#7c4ccb",
  "#c4944a",
  "#3b8d89",
  "#597b4a",
  "#d55046",
  "#71b84e",
  "#4e336c",
  "#b488bb",
  "#8d8d0f",  # Unique color 1
  "black",  # Unique color 2
  "#e51c23",  # Unique color 3 (changed)
  "#dda0dd",  # Unique color 4 (changed)
  "#ff9800"   # Unique color 5 (changed)  # Unique color 5
))
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/globalumap.png"), p1, height = 5, width = 7, dpi = 300)


#Healthy vs Control
OriginalALL <- dealgoseuobj

dealgoseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(dealgoseuobj$scCLINIC_Level))[dealgoseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
DEAlgoSinglet <- subset(dealgoseuobj,subset = scCLINIC_artifact == "Singlet")

DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")

dealgoseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(dealgoseuobj$scCLINIC_Level))[dealgoseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$scCLINIC_contamID <- ifelse(dealgoseuobj$scCLINIC_artifact == "Artifact", dealgoseuobj$scCLINIC_ClusterID, NA)

for (i in c("M8","M10")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
  recluster$scCLINIC_ClusterID <- dealgoseuobj$scCLINIC_ClusterID
  recluster$scCLINICScore <- dealgoseuobj$scCLINICScore
  recluster$scCLINIC_Level <- dealgoseuobj$scCLINIC_Level
  recluster$scCLINIC_artifact <- dealgoseuobj$scCLINIC_artifact
  recluster$scCLINIC_contamID <- dealgoseuobj$scCLINIC_contamID 
  recluster$Status <- dealgoseuobj$Status 
  
  p1 <- DimPlot(recluster,group.by = "Status",cols = c(
    "#e51c23",  # Unique color 3 (changed)
    "#1e88e5"
  ))
  
  p7 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
  p8 <- FeaturePlot(recluster, reduction = "umap",features = "scCLINICScore")
  p4 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,  sizes.highlight = 0.1)
  p5 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_artifact",raster=FALSE,  sizes.highlight = 0.1)
  
  p3 <- DimPlot(recluster, group.by = "scCLINIC_contamID")
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig5_cluster_",i,"_ID.png"), p7+p5+p8+p4+p3+p1, height = 10, width = 17, dpi = 300)
  
}

p1 <- DimPlot(dealgoseuobj,group.by = "scCLINIC_artifact")
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

for (i in c("M8","M10")){
  file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
  recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/",file_name))
  
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
  
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig5_cluster_",i,"_rawscore_Validate.png"), p7+p5+p8+p4+p0+p1+p2+p3, height = 15, width = 17, dpi = 300)
  
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

dealgoseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(dealgoseuobj$scCLINIC_Level))[dealgoseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
DEAlgoSinglet <- subset(dealgoseuobj,subset = scCLINIC_artifact == "Singlet")

DFSinglet <- subset(dealgoseuobj,subset = DFafqc == "Singlet")

dealgoseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(dealgoseuobj$scCLINIC_Level))[dealgoseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$scCLINIC_contamID <- ifelse(dealgoseuobj$scCLINIC_artifact == "Artifact", dealgoseuobj$scCLINIC_ClusterID, NA)

#CD14 Monocyte cells
DEAlgoContam1 <- c("M16S2","M16S4","M13S1","M22S2")
Contam1Cell <- subset(dealgoseuobj,subset = scCLINIC_ClusterID %in% DEAlgoContam1)
Contam1CellDF <- subset(Contam1Cell,subset = DFafqc == "Doublet") # scCLINIC + DF
Contam1CellNDF <- subset(Contam1Cell,subset = DFafqc == "Singlet") # scCLINIC

negative_others <- dealgoseuobj[,dealgoseuobj$scCLINIC_artifact == "Singlet" & dealgoseuobj$DFafqc == "Singlet" & !(dealgoseuobj$L0 == "M2")]
negative_PT <- dealgoseuobj[,dealgoseuobj$scCLINIC_artifact == "Singlet" & dealgoseuobj$DFafqc == "Singlet" & (dealgoseuobj$L0 == "M2")]


dealgoseuobj$grouping <- NA
dealgoseuobj$grouping <- ifelse(colnames(dealgoseuobj) %in% colnames(negative_PT), "CD14Monocyte\n-Singlet",dealgoseuobj$grouping)
dealgoseuobj$grouping <- ifelse(colnames(dealgoseuobj) %in% colnames(Contam1CellDF), "DF+scCLINIC",dealgoseuobj$grouping)
dealgoseuobj$grouping <- ifelse(colnames(dealgoseuobj) %in% colnames(Contam1CellNDF), "scCLINIC",dealgoseuobj$grouping)
dealgoseuobj$grouping <- ifelse(colnames(dealgoseuobj) %in% colnames(negative_others), "Others\n-Singlet",dealgoseuobj$grouping)
ncountobjplot <- dealgoseuobj[,! is.na(dealgoseuobj$grouping)]

negativeexp <- AverageExpression(negative_others,group.by = "all",slot = "counts")
negativePTexp <- AverageExpression(negative_PT,group.by = "all",slot = "counts")
Contam1CellDFexp <- AverageExpression(Contam1CellDF,group.by = "all",slot = "counts")
Contam1CellNDFexp <- AverageExpression(Contam1CellNDF,group.by = "all",slot = "counts")
Contam1CellNDFexp_split <- AverageExpression(Contam1CellNDF,group.by = "scCLINIC_ClusterID",slot = "counts")

contam1gene <- read.csv("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVID_NatureMed_Step2/Output_Overlap_Ratio_0.5/overlaplst_filtered_Overlap_Ratio_0.5_ContaminationInfo.csv",na.strings = c("", "NA"))
contam1gene <- contam1gene[contam1gene$ES_local >0 ,]
contam1genelist <- unique(contam1gene[contam1gene$cluster_ref=="M2",]$gene)

intersection <- contam1genelist

negative_res <- negativeexp$RNA[intersection,]
negativePT_res <- negativePTexp$RNA[intersection,]
Contam1CellDFexp_res <- Contam1CellDFexp$RNA[intersection,]
Contam1CellNDFexp_res <- Contam1CellNDFexp$RNA[intersection,]

library(ggplot2)
library(ggrepel)
# Create data frames for each condition
data_negative <- data.frame(Gene = names(negative_res), Value = negative_res)
data_negative_PT <- data.frame(Gene = names(negativePT_res), Value = negativePT_res)
data2 <- data.frame(Gene = names(Contam1CellDFexp_res), Value = Contam1CellDFexp_res)
data3 <- data.frame(Gene = names(Contam1CellNDFexp_res), Value = Contam1CellNDFexp_res)

# Combine all data frames
all_data <- rbind(data_negative, data_negative_PT, data3, data2)
all_data$Log2Value <- log2(all_data$Value + 1) # + 1 to avoid 0 , log(0) = -infinity
all_data$Condition <- factor(rep(c("Others\n-Singlet", "CD14Monocyte\n-Singlet","scCLINIC", "DF+scCLINIC"), each = length(names(negative_res))),levels = c("CD14Monocyte\n-Singlet","Others\n-Singlet","scCLINIC","DF+scCLINIC"))
# all_data <- rbind(all_data,reshaped_df)
# Create the violin plot with scatter plot
partial_contamination_data <- subset(all_data, Condition == "scCLINIC" & Log2Value > 2)

p1 <-ggplot(all_data, aes(x = Condition, y = Log2Value, fill = Condition)) +
  geom_violin(color = alpha("black",0.25)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  geom_point(aes(color = Condition), size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
  geom_text_repel(data = partial_contamination_data, aes(label = Gene), 
                  hjust = 0.5, vjust = 1, segment.color = "transparent", size = 4,position = position_jitter(width = 0.2))+
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linetype = "dashed", color = alpha("black",0.5), size = 0.5, position = position_dodge(width = 0.75)) + # Add average line
  
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#d55046", "#be883d","#9671c3","#69a75f")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black","black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "CD14 Monocyte Marker\nlog2(Average Expression Level)")+
  
  theme(legend.position = "none") # Remove legend



ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5A.png"), p1, height = 5, width = 5, dpi = 300)

# Perform pairwise t-tests
pairwise_result <- pairwise.t.test(all_data$Log2Value, all_data$Condition)

write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5A.csv"),sep = ",")

ncountobjplot$grouping <- factor(ncountobjplot$grouping, levels = c("CD14Monocyte\n-Singlet", "Others\n-Singlet", "scCLINIC", "DF+scCLINIC"))

#Proved is not Doublet
p <- VlnPlot(ncountobjplot,features = "nCount_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#d55046", "#be883d","#9671c3","#69a75f")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nCount_RNA") +
  theme(legend.position = "none")
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Ancount.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(ncountobjplot$nCount_RNA, ncountobjplot$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Ancount.csv"),sep = ",")

p <- VlnPlot(ncountobjplot,features = "nFeature_RNA",group.by = "grouping")+
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.1), alpha = 1, fill = "white",outlier.shape = NA) + # Add white boxplot within each violin
  #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#d55046", "#be883d","#9671c3","#69a75f")) + # Specify fill colors
  scale_color_manual(values = c("black", "black", "black", "black")) + # Specify scatter plot point colors
  labs(title = "",
       x = "",
       y = "nFeature_RNA") +
  theme(legend.position = "none") # Remove legend
ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Anfeature.png"), p, height = 5, width = 5, dpi = 300)

pairwise_result <- pairwise.t.test(ncountobjplot$nFeature_RNA, ncountobjplot$grouping)
write.table(pairwise_result$p.value,file = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","SupplementaryFig5Anfeature.csv"),sep = ",")

###Arh
c("M17S5","M13S3","M20S0","M21S3","M22S2","M23S4")#Stroma (M8), Endo (M5), Dist.Prox.Tubule (M4), Coll.Duct.IC (M1), Loop of Henle (M6)



Contam1Cell <- subset(dealgoseuobj,subset = scCLINIC_ClusterID %in% DEAlgoContam1)
noncontam1 <- subset(dealgoseuobj,subset = L0 %in% c("M17","M13","M20","M21","M22","M23"))

negative <- noncontam1[,!(noncontam1$scCLINIC_ClusterID %in% DEAlgoContam1)]
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
Contam1CellNDFexp_split <- AverageExpression(Contam1CellNDF,group.by = "scCLINIC_ClusterID",slot = "data")

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

Contam1Cell <- subset(dealgoseuobj,subset = scCLINIC_ClusterID %in% DEAlgoContam1)
noncontam1 <- subset(dealgoseuobj,subset = L0 %in% c("M16","M13","M22"))

negative <- noncontam1[,!(noncontam1$scCLINIC_ClusterID %in% DEAlgoContam1)]
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
negativeexp <- AverageExpression(negative,group.by = "all",slot = "counts")
Contam1CellDFexp <- AverageExpression(Contam1CellDF,group.by = "all",slot = "counts")
Contam1CellNDFexp <- AverageExpression(Contam1CellNDF,group.by = "all",slot = "counts")
Contam1CellNDFexp_split <- AverageExpression(Contam1CellNDF,group.by = "scCLINIC_ClusterID",slot = "counts")

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
    
    p1 <- FeaturePlot(recluster, reduction = "umap",features = markertosee, pt.size = 0.1) & scale_color_gradientn(colors = c("grey","red"),limits = c(0, 6), oob = scales::squish)#,limits = c(0, 5)DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_artifact",raster=FALSE, label = T, sizes.highlight = 0.1)
    #ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig3/kidneymouse_science_Step2/Output_annotation_index/",curcelltype,"_marker_",markertosee,".png"), p1, height = 5, width = 5, dpi = 300)
    plot_list[[curcelltype]] <- p1
  }
  
  all_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = ceiling(length(celltypelist)/2))
  ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/",markertosee,".png"), plot = all_plots, height = 5, width = 5, dpi = 300)
}



#B cells supplementary Fig5
obj_azimuth <- readRDS("~/DEAlgoManuscript/Manuscript_Figures/Fig5/Fig5_Azimuth.rds")

singletlevel = 7

dealgoseuobj <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Output_Overlap_Ratio_0.5/","scCLINICResult.rds"))
doublet_afqc <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemedQC.rds"))

metatable_unique<- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/covid19naturemed.rds"))

dealgoseuobj@meta.data <- cbind(metatable_unique@meta.data,dealgoseuobj@meta.data)

#16.014% of 44721 is 7162
top_indices <- order(doublet_afqc@meta.data$DFScore, decreasing = TRUE)[7162]
cutoff <- doublet_afqc@meta.data$DFScore[top_indices]

dealgoseuobj$DFafqc_pANN <- doublet_afqc$DFScore
dealgoseuobj$DFafqc <- NA
dealgoseuobj$DFafqc <- ifelse(doublet_afqc$DFScore < cutoff, "Singlet", "Doublet")

dealgoseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(dealgoseuobj$scCLINIC_Level))[dealgoseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
dealgoseuobj$scCLINIC_contamID <- ifelse(dealgoseuobj$scCLINIC_artifact == "Artifact", dealgoseuobj$scCLINIC_ClusterID, NA)



i <- "M3"
file_name <- paste0("Overlap_Ratio_0.5","_cluster_",i,".rds")
recluster <- readRDS(paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/",file_name))

recluster$DFafqc_pANN <- dealgoseuobj$DFafqc_pANN
recluster$DFafqc <- dealgoseuobj$DFafqc
recluster$scCLINIC_ClusterID <- dealgoseuobj$scCLINIC_ClusterID
recluster$scCLINICScore <- dealgoseuobj$scCLINICScore
recluster$scCLINIC_Level <- dealgoseuobj$scCLINIC_Level
recluster$scCLINIC_artifact <- dealgoseuobj$scCLINIC_artifact

recluster$scCLINIC_contamID <- dealgoseuobj$scCLINIC_contamID

recluster$azimuthl1 <- obj_azimuth$predicted.celltype.l1
recluster$azimuthl2 <- obj_azimuth$predicted.celltype.l2

recluster$Status <- dealgoseuobj$Status 

p10 <- DimPlot(recluster,group.by = "Status",cols = c(
  "#e51c23",  # Unique color 3 (changed)
  "#1e88e5"
))

p7 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
p8 <- FeaturePlot(recluster, reduction = "umap",features = "scCLINICScore")
p4 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,  sizes.highlight = 0.1)
p5 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_contamID",raster=FALSE,  sizes.highlight = 0.1)
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

ggsave(filename = paste0("~/DEAlgoManuscript/Manuscript_Figures/Fig5/","Fig5_cluster_",i,"_rawscore_Validate.png"), p7+p5+p8+p4+p0+p1+p2+p3+p10, height = 15, width = 17, dpi = 300)


