#!/usr/bin/env Rscript
library(patchwork)
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
library(ggpubr)
library(viridis)
library(purrr)
library(ggforce)

text_size <- 12
font_type <- "Arial"#"Times New Roman"
line_width <- 1.2
bar_width <- 2.5

manuscript_colors <-
  c(
    '#53A85F',    '#F1BB72',    '#57C3F3',    '#D6E7A3',    '#3A6963',
    '#E63863',    '#E95C59',    '#E59CC4',    '#AB3282',    '#23452F',
    '#BD956A',    '#8C549C',    '#585658',    '#9FA3A8',    '#E0D4CA',
    '#5F3D69',    '#C5DEBA',    '#58A4C3',    '#E4C755',    '#F7F398',
    '#AA9A59',    '#E39A35',    '#C1E6F3',    '#6778AE',    '#91D0BE',
    '#B53E2B',    '#476D87',    '#712820',    '#DCC1DD',    '#CCE0F5',
    '#CCC9E6',    '#625D9E',    '#68A180',    '#968175',    '#E5D2DD'
  )

data_path <- "/mnt/lab-store/projects/scCLINIC/Reproduce"

create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
}

create_folder_if_not_exists(paste0(data_path,"/Fig3/"))

##### Figure 3

outdir <- paste0(data_path, "/scCLINIC_Figures/")
create_folder_if_not_exists(outdir)

#Load and save Robj
load(paste0(data_path, "/seu_kidney_codeocean.Robj")) #dataset obtained from https://codeocean.com/capsule/5650599/tree/v1
saveRDS(seu_kidney,paste0(data_path,"/Fig3/seu_kidney_codeocean.rds")) 
#Convert Old seurat (s4 object) to new seurat
obj <- CreateSeuratObject(counts = seu_kidney@raw.data)
obj@meta.data$CT.Park <- seu_kidney@meta.data$CT.Park
saveRDS(obj,paste0(data_path,"/Fig3/seu_kidney_codeocean_updated.rds"))

Name <- "kidneymouse_science"
Output <- paste0(data_path,"/Fig3/")
resol="Manual"
OverlapRatio="CT.Park"
ISThreshold = 0
gene_n=150
CELLANNOTATION = TRUE

obj <- readRDS(paste0(data_path,"/Fig3/seu_kidney_codeocean_updated.rds"))

obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)

obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold,CELLANNOTATION = TRUE)

saveRDS(obj,paste0(data_path,"/Fig3/step1d.rds"))

STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)

obj <- readRDS(paste0(data_path,"/Fig3/step1d.rds"))

obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)

obj <- readRDS(paste0(data_path,"/Fig3/kidneymouse_science_Step2/Output_annotation_index/scCLINICResult.rds"))

original_obj <- readRDS(paste0(data_path,"/Fig3/seu_kidney_codeocean_updated.rds"))

obj$CT.Park <- original_obj$CT.Park

PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)

#Follow DoubletFinder tutorial https://github.com/chris-mcginnis-ucsf/DoubletFinder 
###################################################
## Step 2: Perform and summarize parameter sweep ##
###################################################
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
library(DoubletFinder)               
seu_kidney <- readRDS(paste0(data_path,"/Fig3/seu_kidney_codeocean_updated.rds"))

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

saveRDS(seu_kidney,paste0(data_path,"/Fig3/seu_kidney_codeocean_updated_withDFScore.rds"))

#Plot Figure 3
scCLINICseuobj <-
  readRDS(paste0(
    data_path,"/Fig3/kidneymouse_science_Step2/Output_annotation_index/scCLINICResult.rds") ###SCCLINIC
  )

seu_kidney <- readRDS(paste0(data_path,"/Fig3/seu_kidney_codeocean_updated.rds"))

scCLINICseuobj@meta.data$CT.Park <- seu_kidney@meta.data$CT.Park

scCLINICseuobj$nCount_RNA <- colSums(x = scCLINICseuobj, slot = "counts")  # nCount_RNA
scCLINICseuobj$nFeature_RNA <- colSums(x = GetAssayData(object = scCLINICseuobj, slot = "counts") > 0)

DFseuobj <-
  readRDS(paste0(
    data_path,"/Fig3/seu_kidney_codeocean_updated_withDFScore.rds")
  )

# Copy DF metadata into scCLINICseuobj
scCLINICseuobj@meta.data$DFafqc <-
  DFseuobj@meta.data[, grep("DF.",  colnames(DFseuobj@meta.data), value = TRUE)]
rm(DFseuobj)
scCLINICseuobj$scCLINIC_Contaminated <-
  ifelse(test = as.numeric(levels(scCLINICseuobj$scCLINIC_Level))[scCLINICseuobj$scCLINIC_Level] < 4,
         yes =  "Artifact",
         no = "Singlet")
scCLINICseuobj$scCLINIC_contamclus <-
  ifelse(test = scCLINICseuobj$scCLINIC_Contaminated == "Artifact",
    yes = scCLINICseuobj$scCLINIC_ClusterID,
    no = "Singlet")

scCLINICseuobj$annotation_index_celltype <-
  dplyr::recode(
    scCLINICseuobj$annotation_index,
    "M1" = "CDIC_M1",
    "M2" = "CDPC_M2",
    "M3" = "CDTC_M3",
    "M4" = "DCT_M4",
    "M5" = "EC_M5",
    "M6" = "LH_M6",
    "M7" = "PT_M7",
    "M8" = "SC_M8"
  )

plot_A <-
  wrap_elements(
    DimPlot(
      scCLINICseuobj,
      group.by = "annotation_index_celltype",
      reduction = "umap",
      raster = FALSE,
      cols = manuscript_colors
    ) & ggtitle(NULL) & 
      theme(
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.75, 0.2),
        legend.title = element_text(size = text_size, family = font_type),
        legend.text = element_text(size = text_size, family = font_type),
        legend.key.size = unit(0.5, 'lines'),
        legend.spacing.x = unit(0.2, "lines"),
        legend.spacing.y = unit(0, "lines")
      ) & labs(color = "Subcluster")
  )

plot_SA <-
  DimPlot(
    scCLINICseuobj,
    group.by = "CT.Park",
    reduction = "umap",
    raster = FALSE,
    cols = manuscript_colors
  )

contamgeneinfo <-
  read.csv(
    paste0(
      data_path,
      "/Fig3/kidneymouse_science_Step2/Output_annotation_index/ArtifactsInfo.csv"
    ),
    na.strings = c("", "NA")
  )
#scCLINIC Subcluster ID
contamgeneinfo$MajorSub <-
  paste0(contamgeneinfo$Major_Cluster,
         "_",
         contamgeneinfo$Subcluster)
#Summarize ES score and their source of major cluster for each subclusters
result <- contamgeneinfo %>%
  group_by(MajorSub) %>%
  summarize(cluser_reflst = list(Source_of_Artifacts_SoA),
            dp1lst = list(Enrichment_Score)) %>%
  ungroup()

# Function to calculate the average ES score (dp1lst) for each source of artifacts (cluster_reflst)
tabulate_Cluster_Contribution_Score <-
  function(cluster_reflst, dp1lst) {
    components <- unlist(cluster_reflst)
    dp1_values <- unlist(dp1lst)
    CCS_Matrix <- tapply(dp1_values, components, mean, na.rm = TRUE)
    return(CCS_Matrix)
  }

# For each subclusters (each row in result), calculate the average ES score for each source of artifacts
tabulate_CCS <-
  mapply(
    tabulate_Cluster_Contribution_Score,
    result$cluser_reflst,
    result$dp1lst,
    SIMPLIFY = FALSE
  )

# List of all source of artifacts which contaminated major cluster X
lst_of_source_of_artifacts <-
  unique(unlist(lapply(tabulate_CCS, names)))

# Create a matrix, each row represent one subclusters and each column present each source of artifacts, to store the CCS for each major clusters
CCS_Matrix <-
  matrix(
    NA,
    nrow = length(tabulate_CCS),
    ncol = length(lst_of_source_of_artifacts),
    dimnames = list(NULL, lst_of_source_of_artifacts)
  )

# store the CCS for each major clusters in the matrix
for (i in seq_along(tabulate_CCS)) {
  CCS_Matrix[i, names(tabulate_CCS[[i]])] <- tabulate_CCS[[i]]
}

# Replace NA with 0
CCS_Matrix[is.na(CCS_Matrix)] <- 0

# Convert to dataframe for plotting
CCS_Matrix <- as.data.frame(CCS_Matrix)
CCS_Matrix$MajorSub <-
  result$MajorSub #Named each rows with their Subclusters ID
data <- CCS_Matrix %>%
  pivot_longer(cols = -MajorSub,
               names_to = "Component",
               values_to = "value") #dependency tidyr

#"M8" M4 M5
clusterx <- "M8"
file_name <- paste0("annotation_index_cluster_", clusterx, ".rds")
recluster <-
  readRDS(
    paste0(
      data_path,"/Fig3/kidneymouse_science_Step2/annotation_index_recluster/",
      file_name
    )
  )

recluster$scCLINIC_contamclus <- scCLINICseuobj$scCLINIC_contamclus
recluster$DFafqc <- scCLINICseuobj$DFafqc

B1 <-
  DimPlot(
    recluster,
    group.by = "scCLINIC_contamclus",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = c("#E59CC4", "#53A85F", "#9FA3A8")
  ) & labs(color = "scCLINIC") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

B2 <-
  DimPlot(
    recluster,
    reduction = "umap",
    group.by = "DFafqc",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = c("#E63863", "#9FA3A8")
  ) & labs(color = "DF") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.55, 0.95),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

filtered_data <- data %>%
  dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

heatmap_data <-
  dcast(filtered_data, MajorSub ~ Component, value.var = "value")
heatmap_data <-
  heatmap_data[,!colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

rownames(heatmap_data) <-
  heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL

# Summarize scCLINIC Score (CS) of each subclusters
CS_table <- scCLINICseuobj@meta.data %>%
  group_by(scCLINIC_ClusterID) %>% #scCLINIC_ClusterID
  summarize(scCLINICScore = mean(scCLINICScore)) %>% #scCLINICScore of each cells within each subclusters are same value)
  ungroup()

#Convert to plotting dataframe
heatmap_data <-
  as.data.frame(heatmap_data) %>% # Convert heatmap_data to a dataframe
  mutate(CS = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score #scCLINIC_ClusterID

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <- data.frame(heatmap_matrix)

#rownames(multicolplot) <- gsub(".*_", "", rownames(multicolplot))

rownames_sorted <-
  rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
multicolplot <- multicolplot[rownames_sorted, ]

colnames_sorted <-
  colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
multicolplot <- multicolplot[, colnames_sorted]

multicolplot$Category <-
  as.character(rownames(multicolplot))
multicolplot$Category <-
  factor(multicolplot$Category, levels = rownames_sorted)

x_limits <-
  c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <-
  seq(x_limits[1], x_limits[2], length.out = 3)

# Create individual plots without y-axis text
plot_list <- list()
for (i in seq_along(colnames_sorted)) {
  clus <- colnames_sorted[i]
  p <-
    ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(
      stat = "identity",
      fill = c(
        '#53A85F',
        '#F1BB72',
        '#57C3F3',
        '#D6E7A3',
        '#3A6963',
        '#E63863',
        '#E59CC4',
        '#AB3282'
      )[i]
    ) + 
    labs(title = paste0(
      dplyr::recode(
        clus,
        "M1" = "CDIC\nM1",
        "M2" = "CDPC\nM2",
        "M3" = "CDTC\nM3",
        "M4" = "DCT\nM4",
        "M5" = "EC\nM5",
        "M6" = "LH\nM6",
        "M7" = "PT\nM7",
        "M8" = "SC\nM8"
      )
    ),
    y = "",
    x = "Subcluster") +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "grey", size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = text_size, family = font_type),
      axis.title = element_text(size = text_size, family = font_type),
      plot.title = element_text(size = text_size, family = font_type),
      axis.text.x = element_blank(),
      plot.margin = unit(c(0, -1, 0, 0), "cm")
    ) +
    coord_flip() +
    scale_y_continuous(breaks = custom_breaks, limits = x_limits)
  
  plot_list[[i]] <- p
}

# Remove y-axis labels, title, and ticks from all but the left-most plot
plot_list[2:length(plot_list)] = plot_list[2:length(plot_list)] %>%
  map(
    ~ .x + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  )


# Combine the y-axis labels plot with the other plots
B3 <-
  wrap_elements(
    wrap_plots(plot_list, ncol = length(plot_list)) +
      plot_layout(guides = "collect") &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  )
#"M4"
clusterx <- "M4"
file_name <- paste0("annotation_index_cluster_", clusterx, ".rds")
recluster <-
  readRDS(
    paste0(
      data_path,"/Fig3/kidneymouse_science_Step2/annotation_index_recluster/",
      file_name
    )
  )

recluster$scCLINIC_contamclus <-
  scCLINICseuobj$scCLINIC_contamclus
recluster$DFafqc <- scCLINICseuobj$DFafqc

B4 <-
  DimPlot(
    recluster,
    group.by = "scCLINIC_contamclus",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = c("#E59CC4", "#9FA3A8")
  ) & labs(color = "scCLINIC") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.65, 0.9),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

B5 <-
  DimPlot(
    recluster,
    reduction = "umap",
    group.by = "DFafqc",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = c("#E63863", "#9FA3A8")
  ) & labs(color = "DF") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.6, 0.95),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

filtered_data <- data %>%
  dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

heatmap_data <-
  dcast(filtered_data, MajorSub ~ Component, value.var = "value")
heatmap_data <-
  heatmap_data[, !colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

rownames(heatmap_data) <-
  heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL

# Summarize scCLINIC Score (CS) of each subclusters
CS_table <- scCLINICseuobj@meta.data %>%
  group_by(scCLINIC_ClusterID) %>% #scCLINIC_ClusterID
  summarize(scCLINICScore = mean(scCLINICScore)) %>%  #scCLINICScore of each cells within each subclusters are same value... #scCLINICScore
  ungroup()

#Convert to plotting dataframe
heatmap_data <-
  as.data.frame(heatmap_data) %>% # Convert heatmap_data to a dataframe
  mutate(CS = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score #scCLINIC_ClusterID

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <- data.frame((heatmap_matrix))

#rownames(multicolplot) <- gsub(".*_", "", rownames(multicolplot))

rownames_sorted <-
  rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
multicolplot <-
  multicolplot[rownames_sorted,]

colnames_sorted <-
  colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
multicolplot <-
  multicolplot[, colnames_sorted]


multicolplot$Category <-
  as.character(rownames(multicolplot))
multicolplot$Category <-
  factor(multicolplot$Category, levels = rownames_sorted)

x_limits <-
  c(0, ceiling(max(multicolplot[,-ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <-
  seq(x_limits[1], x_limits[2], length.out = 3)

# Create individual plots without y-axis text
plot_list <- list()
for (i in seq_along(colnames_sorted)) {
  clus <- colnames_sorted[i]
  p <-
    ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(
      stat = "identity",
      fill = c(
        '#53A85F',
        '#F1BB72',
        '#57C3F3',
        '#D6E7A3',
        '#3A6963',
        '#E63863',
        '#E59CC4',
        '#AB3282'
      )[i]
    ) + # viridis(length(colnames_sorted))[i]) +
    #geom_bar(stat = "identity", fill = viridis(length(colnames_sorted))[i]) +
    labs(title = paste0(
      dplyr::recode(
        clus,
        "M1" = "CDIC\nM1",
        "M2" = "CDPC\nM2",
        "M3" = "CDTC\nM3",
        "M4" = "DCT\nM4",
        "M5" = "EC\nM5",
        "M6" = "LH\nM6",
        "M7" = "PT\nM7",
        "M8" = "SC\nM8"
      )
    ),
    y = "",
    x = "Subcluster") +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "grey", size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = text_size, family = font_type),
      axis.title = element_text(size = text_size, family = font_type),
      plot.title = element_text(size = text_size, family = font_type),
      axis.text.x = element_blank(),
      plot.margin = unit(c(0,-1, 0, 0), "cm")
    ) +
    coord_flip() +
    scale_y_continuous(breaks = custom_breaks, limits = x_limits)
  
  plot_list[[i]] <- p
}

# Remove y-axis labels, title, and ticks from all but the left-most plot
plot_list[2:length(plot_list)] = plot_list[2:length(plot_list)] %>%
  map(
    ~ .x + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  )


# Combine the y-axis labels plot with the other plots
B6 <-
  wrap_elements(
    wrap_plots(plot_list, ncol = length(plot_list)) +
      plot_layout(guides = "collect") &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  )
#"M5"
clusterx <- "M5"
file_name <- paste0("annotation_index_cluster_", clusterx, ".rds")
recluster <-
  readRDS(
    paste0(
      data_path,"/Fig3/kidneymouse_science_Step2/annotation_index_recluster/",
      file_name
    )
  )

recluster$scCLINIC_contamclus <-
  scCLINICseuobj$scCLINIC_contamclus
recluster$DFafqc <- scCLINICseuobj$DFafqc

B7 <-
  DimPlot(
    recluster,
    group.by = "scCLINIC_contamclus",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = c("#E59CC4", "#53A85F", '#D6E7A3', "#9FA3A8")
  ) & labs(color = "scCLINIC") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.75, 0.85),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

B8 <-
  DimPlot(
    recluster,
    reduction = "umap",
    group.by = "DFafqc",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = c("#E63863", "#9FA3A8")
  ) & labs(color = "DF") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.7, 0.95),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

filtered_data <- data %>%
  dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

heatmap_data <-
  dcast(filtered_data, MajorSub ~ Component, value.var = "value")
heatmap_data <-
  heatmap_data[,!colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

rownames(heatmap_data) <-
  heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL

# Summarize scCLINIC Score (CS) of each subclusters
CS_table <- scCLINICseuobj@meta.data %>%
  group_by(scCLINIC_ClusterID) %>% #scCLINIC_ClusterID
  summarize(scCLINICScore = mean(scCLINICScore)) %>%  #scCLINICScore of each cells within each subclusters are same value... #scCLINICScore
  ungroup()

#Convert to plotting dataframe
heatmap_data <-
  as.data.frame(heatmap_data) %>% # Convert heatmap_data to a dataframe
  mutate(CS = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score #scCLINIC_ClusterID

# Convert data to matrix
heatmap_matrix <-
  as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <-
  data.frame((heatmap_matrix))

#rownames(multicolplot) <- gsub(".*_", "", rownames(multicolplot))

rownames_sorted <-
  rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
multicolplot <-
  multicolplot[rownames_sorted, ]

colnames_sorted <-
  colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
multicolplot <-
  multicolplot[, colnames_sorted]

multicolplot$Category <-
  as.character(rownames(multicolplot))
multicolplot$Category <-
  factor(multicolplot$Category, levels = rownames_sorted)

x_limits <-
  c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <-
  seq(x_limits[1], x_limits[2], length.out = 3)

# Create individual plots without y-axis text
plot_list <- list()
for (i in seq_along(colnames_sorted)) {
  clus <- colnames_sorted[i]
  p <-
    ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(
      stat = "identity",
      fill = c(
        '#53A85F',
        '#F1BB72',
        '#57C3F3',
        '#D6E7A3',
        '#3A6963',
        '#E63863',
        '#E59CC4',
        '#AB3282'
      )[i]
    ) + # viridis(length(colnames_sorted))[i]) +
    #geom_bar(stat = "identity", fill = viridis(length(colnames_sorted))[i]) +
    labs(title = paste0(
      dplyr::recode(
        clus,
        "M1" = "CDIC\nM1",
        "M2" = "CDPC\nM2",
        "M3" = "CDTC\nM3",
        "M4" = "DCT\nM4",
        "M5" = "EC\nM5",
        "M6" = "LH\nM6",
        "M7" = "PT\nM7",
        "M8" = "SC\nM8"
      )
    ),
    y = "",
    x = "Subcluster") +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "grey", size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = text_size, family = font_type),
      axis.title = element_text(size = text_size, family = font_type),
      plot.title = element_text(size = text_size, family = font_type),
      axis.text.x = element_blank(),
      plot.margin = unit(c(0, -1, 0, 0), "cm")
    ) +
    coord_flip() +
    scale_y_continuous(breaks = custom_breaks, limits = x_limits)
  
  plot_list[[i]] <- p
}

# Remove y-axis labels, title, and ticks from all but the left-most plot
plot_list[2:length(plot_list)] = plot_list[2:length(plot_list)] %>%
  map(
    ~ .x + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  )


# Combine the y-axis labels plot with the other plots
B9 <-
  wrap_elements(
    wrap_plots(plot_list, ncol = length(plot_list)) +
      plot_layout(guides = "collect") &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  )

plot_M8 <-
  wrap_elements(B3 / (B1 + B2) +
                  plot_annotation(#title = "SC_M8",
                    theme = theme(
                      plot.title = element_text(
                        face = "bold",
                        size = text_size,
                        family = font_type,
                        hjust = 0.5
                      ),
                      plot.margin = unit(c(0, -0.1, 0, 0), "cm")
                    )))

#+draw_label("K562_M3 Artifact Enrichment Score", fontface='bold',y=1.1,fontfamily = font_type,size = text_size)

plot_M4 <-
  wrap_elements(B6 / (B4 + B5) +
                  plot_annotation(#title = "DCT_M4",
                    theme = theme(
                      plot.title = element_text(
                        face = "bold",
                        size = text_size,
                        family = font_type,
                        hjust = 0.5
                      ),
                      plot.margin = unit(c(0, -0.1, 0, 0), "cm")
                    )))
plot_M5 <-
  wrap_elements(B9 / (B7 + B8) +
                  plot_annotation(#title = "EC_M5",
                    theme = theme(
                      plot.title = element_text(
                        face = "bold",
                        size = text_size,
                        family = font_type,
                        hjust = 0.5
                      ),
                      plot.margin = unit(c(0, -0.1, 0, 0), "cm")
                    )))

plot_B <-
  wrap_elements(plot_M4 + plot_M5 + plot_M8)




plot_list <- list()

for (curcelltype in c("M7", "M4", "M5", "M8")) {
  recluster <-
    readRDS(
      paste0(
        data_path,"/Fig3/kidneymouse_science_Step2/annotation_index_recluster/annotation_index_cluster_",
        curcelltype,
        ".rds"
      )
    )
  plot_gene <-
    FeaturePlot(
      recluster,
      reduction = "umap",
      features = "Miox",
      pt.size = 0.1
    ) &
    scale_color_gradientn(colors = c("grey", "#E63863"),
                          limits = c(0, 5)) &
    theme(
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "None",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  plot_gene <-
    FeaturePlot(
      recluster,
      reduction = "umap",
      features = "Slc34a1",
      pt.size = 0.1
    ) &
    scale_color_gradientn(colors = c("grey", "#E63863"),
                          limits = c(0, 5)) &
    theme(
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "None",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  plot_list[[paste0("Miox_", curcelltype)]] <-
    plot_gene
  plot_list[[paste0("Slc34a1_", curcelltype)]] <-
    plot_gene
  
}

recluster <-
  readRDS(
    paste0(
      data_path,"/Fig3/kidneymouse_science_Step2/annotation_index_recluster/annotation_index_cluster_",
      "M8",
      ".rds"
    )
  )
plot_gene <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    features = "Slc34a1",
    pt.size = 0.1
  ) &
  scale_color_gradientn(
    colors = c("grey", "#E63863"),
    limits = c(0, 5),
    breaks = c(0, 2.5, 5)
  ) &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.2, 0.9),
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
plot_list[[paste0("Slc34a1_", "M8")]] <-
  plot_gene

# Create a separate plot to extract the legend
legend_plot <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    features = "Miox",
    pt.size = 0.1
  ) &
  scale_color_gradientn(colors = c("grey", "#E63863"), limits = c(0, 5)) &
  theme_void() &
  theme(legend.position = "right", plot.margin = unit(c(0, 0, 0, 0), "cm")) # Only keep the legend

# Extract the legend
legend <-
  cowplot::get_legend(legend_plot)

plot_C <-
  wrap_elements(wrap_plots(plot_list, nrow = 1, ncol = 8) &
                  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))

scCLINICcontam1 <-
  c("M1_S4", "M2_S6", "M4_S0", "M5_S0", "M6_S5", "M8_S2")#Stroma (M8), Endo (M5), Dist.Prox.Tubule (M4), Coll.Duct.IC (M1), Loop of Henle (M6)

Contam1Cell <-
  subset(scCLINICseuobj, subset = scCLINIC_ClusterID %in% scCLINICcontam1)
Contam1CellDF <-
  subset(Contam1Cell, subset = DFafqc == "Doublet") # scCLINIC + DF
Contam1CellNDF <-
  subset(Contam1Cell, subset = DFafqc == "Singlet") # scCLINIC

singlets <-
  scCLINICseuobj[, scCLINICseuobj$scCLINIC_Contaminated == "Singlet" &
                   scCLINICseuobj$DFafqc == "Singlet"]

scCLINICseuobj$grouping <- NA
scCLINICseuobj$grouping <-
  ifelse(
    colnames(scCLINICseuobj) %in% colnames(Contam1CellDF),
    "DF+scCLINIC",
    scCLINICseuobj$grouping
  )
scCLINICseuobj$grouping <-
  ifelse(
    colnames(scCLINICseuobj) %in% colnames(Contam1CellNDF),
    "scCLINIC",
    scCLINICseuobj$grouping
  )
scCLINICseuobj$grouping <-
  ifelse(
    colnames(scCLINICseuobj) %in% colnames(singlets),
    "Singlet",
    scCLINICseuobj$grouping
  )
ncountobjplot <-
  scCLINICseuobj[,!is.na(scCLINICseuobj$grouping)]

ncountobjplot$grouping <-
  factor(ncountobjplot$grouping,
         levels = c("Singlet", "scCLINIC", "DF+scCLINIC"))

my_comparisons <-
  list(c("Singlet", "scCLINIC"),
       c("Singlet", "DF+scCLINIC"),
       c("scCLINIC", "DF+scCLINIC"))

#not Doublet
plot_D <-
  wrap_elements(
    VlnPlot(ncountobjplot, features = "nFeature_RNA", group.by = "grouping") +
      geom_boxplot(
        width = 0.05,
        position = position_dodge(width = 0.1),
        alpha = 1,
        fill = "white",
        outlier.shape = NA
      ) + # Add white boxplot within each violin
      #stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", fill = "yellow", position = position_dodge(width = 0.75)) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 18,
        size = 3,
        color = "black",
        position = position_dodge(width = 0.75)
      ) +
      theme_minimal(base_size = 14) +
      scale_fill_manual(values = manuscript_colors) + # Specify fill colors
      scale_color_manual(values = c("black", "black", "black")) + # Specify scatter plot point colors
      labs(title = "",
           x = "",
           y = "nFeature_RNA") +
      theme(legend.position = "none") +
      stat_compare_means(
        test_sign = grouping,
        label = "p.signif",
        comparisons = my_comparisons,
        method = "t.test"
      ) +
      scale_y_continuous(limits = c(0, 4000))
  )

# Assuming 'ncountobjplot' is your Seurat object and 'grouping' is the grouping variable
# Extract the metadata
metadata <-
  ncountobjplot@meta.data

# Calculate the average nFeature_RNA for each group
average_nFeature_RNA <-
  metadata %>%
  group_by(grouping) %>%
  summarize(mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE))

# Calculate the average nCount_RNA for each group
average_nCount_RNA <-
  metadata %>%
  group_by(grouping) %>%
  summarize(mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE))

recluster <-
  readRDS(
    paste0(
      data_path,"/Fig3/kidneymouse_science_Step2/annotation_index_recluster/annotation_index_cluster_",
      "M5",
      ".rds"
    )
  )
recluster$Cluster <-
  ifelse(recluster$scCLINIC_subcluster == "S5", "M5_S5", "Others")

CDPCMarker <-
  read.csv(paste0(data_path,"/CollDuctPC.csv"),
           na.strings = c("", "NA"))
#M2
scCLINICseuobj <-
  AddModuleScore(
    object = scCLINICseuobj,
    features = list(CDPCMarker$genes),
    #Global_Top20[Global_Top20$cluster == "M2","gene"],
    ctrl = 5,
    name = "CDPC_M2_Marker"
  )
recluster$CDPC_M2_Marker <-
  scCLINICseuobj$CDPC_M2_Marker1
E1 <- FeaturePlot(
  scCLINICseuobj,
  reduction = "umap",
  features = "CDPC_M2_Marker1",
  pt.size = 0.001
) + scale_color_viridis(
  direction = -1,
  limits = c(-0.5, 1.5),
  breaks = c(-0.5, 1.5)
) &
  labs(title = "CDPC_M2_Module") &
  theme(
    plot.title = element_text(
      face = "bold",
      size = text_size,
      family = font_type
    ),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(1.5, 1),
    legend.direction = "horizontal"
  )
threshold = 0.5
filtered_df <- E1[[1]]$data %>%
  group_by(ident) %>%
  filter(
    umap_1 <= quantile(umap_1, 0.75) + threshold * IQR(umap_1) &
      umap_1 >= quantile(umap_1, 0.25) - threshold * IQR(umap_1) &
      umap_2 <= quantile(umap_2, 0.75) + threshold * IQR(umap_2) &
      umap_2 >= quantile(umap_2, 0.25) - threshold * IQR(umap_2)
  )
E1 <- E1 + geom_mark_hull(
  data = filtered_df,
  aes(
    x = umap_1,
    y = umap_2,
    label = ident,
    filter = ident == "M2"
  ),
  radius = unit(2.5, "mm"),
  expand = unit(2.5, "mm"),
  concavity = 1.5,
)
E2 <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    shape.by = "Cluster",
    features = "CDPC_M2_Marker",
    pt.size = 1
  ) +
  scale_color_viridis(
    direction = -1,
    limits = c(-0.5, 1.5),
    breaks = c(-0.5, 1.5)
  ) +
  scale_shape_manual(values = c(18, 16)) &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  )
filtered_df <- E2[[1]]$data %>%
  group_by(Cluster) %>%
  filter(
    umap_1 <= quantile(umap_1, 0.75) + threshold * IQR(umap_1) &
      umap_1 >= quantile(umap_1, 0.25) - threshold * IQR(umap_1) &
      umap_2 <= quantile(umap_2, 0.75) + threshold * IQR(umap_2) &
      umap_2 >= quantile(umap_2, 0.25) - threshold * IQR(umap_2)
  )
E2 <- E2 + geom_mark_hull(
  data = filtered_df,
  aes(
    x = umap_1,
    y = umap_2,
    label = Cluster,
    filter = Cluster == "M5_S5"
  ),
  radius = unit(4, "mm"),
  expand = unit(4, "mm"),
  concavity = 1.5,
)
DCTMarker <-
  read.csv(
    paste0(data_path,"/DistProxTubule.csv"),
    na.strings = c("", "NA")
  )
#M4
scCLINICseuobj <-
  AddModuleScore(
    object = scCLINICseuobj,
    features = list(DCTMarker$genes),
    #Global_Top20[Global_Top20$cluster == "M4","gene"],
    ctrl = 5,
    name = "DCT_M4_Marker"
  )
recluster$DCT_M4_Marker <-
  scCLINICseuobj$DCT_M4_Marker1
E3 <-
  FeaturePlot(
    scCLINICseuobj,
    reduction = "umap",
    features = "DCT_M4_Marker1",
    pt.size = 0.001
  ) + scale_color_viridis(
    direction = -1,
    limits = c(-0.5, 1.5),
    breaks = c(-0.5, 1.5)
  ) &
  labs(title = "DCT_M4_Marker") &
  theme(
    plot.title = element_text(
      face = "bold",
      size = text_size,
      family = font_type
    ),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(1.5, 1),
    legend.direction = "horizontal",
    #    legend.text = element_text(size = text_size, family = font_type)
  )
filtered_df <- E3[[1]]$data %>%
  group_by(ident) %>%
  filter(
    umap_1 <= quantile(umap_1, 0.75) + threshold * IQR(umap_1) &
      umap_1 >= quantile(umap_1, 0.25) - threshold * IQR(umap_1) &
      umap_2 <= quantile(umap_2, 0.75) + threshold * IQR(umap_2) &
      umap_2 >= quantile(umap_2, 0.25) - threshold * IQR(umap_2)
  )
E3 <- E3 + geom_mark_hull(
  data = filtered_df,
  aes(
    x = umap_1,
    y = umap_2,
    label = ident,
    filter = ident == "M4"
  ),
  radius = unit(2.5, "mm"),
  expand = unit(2.5, "mm"),
  concavity = 1.5,
)
E4 <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    shape.by = "Cluster",
    features = "DCT_M4_Marker",
    pt.size = 1
  ) +
  scale_color_viridis(
    direction = -1,
    limits = c(-0.5, 1.5),
    breaks = c(-0.5, 1.5)
  ) +
  scale_shape_manual(values = c(18, 16)) &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  )
filtered_df <- E4[[1]]$data %>%
  group_by(Cluster) %>%
  filter(
    umap_1 <= quantile(umap_1, 0.75) + threshold * IQR(umap_1) &
      umap_1 >= quantile(umap_1, 0.25) - threshold * IQR(umap_1) &
      umap_2 <= quantile(umap_2, 0.75) + threshold * IQR(umap_2) &
      umap_2 >= quantile(umap_2, 0.25) - threshold * IQR(umap_2)
  )
E4 <- E4 + geom_mark_hull(
  data = filtered_df,
  aes(
    x = umap_1,
    y = umap_2,
    label = Cluster,
    filter = Cluster == "M5_S5"
  ),
  radius = unit(2.5, "mm"),
  expand = unit(2.5, "mm"),
  concavity = 1.5,
)
LHMarker <-
  read.csv(
    paste0(data_path,"/LoopofHenle.csv"),
    na.strings = c("", "NA")
  )
#M6
scCLINICseuobj <-
  AddModuleScore(
    object = scCLINICseuobj,
    features = list(LHMarker$genes),
    #Global_Top20[Global_Top20$cluster == "M6","gene"],
    ctrl = 5,
    name = "LH_M6_Marker"
  )
recluster$LH_M6_Marker <-
  scCLINICseuobj$LH_M6_Marker1
E5 <-
  FeaturePlot(
    scCLINICseuobj,
    reduction = "umap",
    features = "LH_M6_Marker1",
    pt.size = 0.001
  ) + scale_color_viridis(
    direction = -1,
    limits = c(-0.5, 1.5),
    breaks = c(-0.5, 1.5)
  ) &
  labs(title = "LH_M6_Marker") &
  theme(
    plot.title = element_text(
      face = "bold",
      size = text_size,
      family = font_type
    ),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(1.5, 1),
    legend.direction = "horizontal"
  )
filtered_df <- E5[[1]]$data %>%
  group_by(ident) %>%
  filter(
    umap_1 <= quantile(umap_1, 0.75) + threshold * IQR(umap_1) &
      umap_1 >= quantile(umap_1, 0.25) - threshold * IQR(umap_1) &
      umap_2 <= quantile(umap_2, 0.75) + threshold * IQR(umap_2) &
      umap_2 >= quantile(umap_2, 0.25) - threshold * IQR(umap_2)
  )
E5 <- E5 + geom_mark_hull(
  data = filtered_df,
  aes(
    x = umap_1,
    y = umap_2,
    label = ident,
    filter = ident == "M6"
  ),
  radius = unit(2.5, "mm"),
  expand = unit(2.5, "mm"),
  concavity = 1.5,
)
E6 <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    shape.by = "Cluster",
    features = "LH_M6_Marker",
    pt.size = 1
  ) +
  scale_color_viridis(
    direction = -1,
    limits = c(-0.5, 1.5),
    breaks = c(-0.5, 1.5)
  ) +
  scale_shape_manual(values = c(18, 16)) &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None"
  )
filtered_df <- E6[[1]]$data %>%
  group_by(Cluster) %>%
  filter(
    umap_1 <= quantile(umap_1, 0.75) + threshold * IQR(umap_1) &
      umap_1 >= quantile(umap_1, 0.25) - threshold * IQR(umap_1) &
      umap_2 <= quantile(umap_2, 0.75) + threshold * IQR(umap_2) &
      umap_2 >= quantile(umap_2, 0.25) - threshold * IQR(umap_2)
  )
E6 <- E6 + geom_mark_hull(
  data = filtered_df,
  aes(
    x = umap_1,
    y = umap_2,
    label = Cluster,
    filter = Cluster == "M5_S5"
  ),
  radius = unit(2.5, "mm"),
  expand = unit(2.5, "mm"),
  concavity = 1.5,
)
plot_E <-
  wrap_elements(E1 | E2 | E3 | E4 | E5 | E6)

# Define a theme with all required elements
thm <- theme(
  text = element_text(
    size = text_size,
    family = font_type,
    face = "bold"
  ),
  plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
  # Remove all margins
)

# Create the final figure layout
figure_3 <-
  (((plot_A | plot_B) +
      plot_layout(widths = c(0.2, 0.8), ncol = 2) &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  ) /
    (plot_C +
       plot_layout() &
       theme(plot.margin = unit(c(
         0, 0, 0, 0
       ), "cm"))) /
    ((plot_D | plot_E) +
       plot_layout(widths = c(0.3, 0.7), ncol = 2) &
       theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    )) +
  plot_layout(heights = c(0.45, 0.25, 0.25), nrow = 3) &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

figure_3 <-
  figure_3 + plot_annotation(tag_levels = "A", theme = thm)
ggsave(
  filename = paste0(outdir, "Figure_3.png"),
  plot = figure_3,
  height = 12,
  width = 18,
  dpi = 300
)