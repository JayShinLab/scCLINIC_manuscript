#!/usr/bin/env Rscript
data_path = commandArgs(trailingOnly = TRUE)

if (length(data_path) == 0) {
  stop("Please specify path to manuscript data.", call. = FALSE)
}

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
library(patchwork)
library(ggpubr)
library(purrr)


text_size <- 12
font_type <- "Arial"
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

create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
}

create_folder_if_not_exists(paste0(data_path,"/Fig5/"))

#Please download blish_covid.seu.rds from the "Blish Lab Peripheral Blood Mononuclear Cells (PBMCs)" entry archive available at https://www.covid19cellatlas.org/index.patient.html.
#Reference https://www.nature.com/articles/s41591-020-0944-y

covid19natmed <- readRDS(paste0(data_path,"/blish_covid.seu.rds"))

covid19natmed.updated = UpdateSeuratObject(object = covid19natmed)

saveRDS(covid19natmed.updated,paste0(data_path,"/Fig5/covid19naturemed.rds"))

covid19obj <- readRDS(paste0(data_path,"/Fig5/covid19naturemed.rds"))

covidnat <- covid19obj@assays$RNA$counts

#apply DoubletFinder
library(DoubletFinder)

subobj <- CreateSeuratObject(covidnat)

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
covid_obj <- subobj
covid_obj <- NormalizeData(covid_obj)
covid_obj <- FindVariableFeatures(covid_obj, selection.method = "vst", nfeatures = 2000)
covid_obj <- ScaleData(covid_obj)
covid_obj <- RunPCA(covid_obj)
covid_obj <- RunUMAP(covid_obj, dims = 1:10)

sweep.vector <- DoubletFinder::paramSweep_v3(covid_obj, PCs = 1:10, sct = FALSE)
sweep.table <- DoubletFinder::summarizeSweep(sweep.vector, GT = FALSE)
bcmvn <- DoubletFinder::find.pK(sweep.table)
pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
pK <- as.numeric(levels(pK))[pK]
covid_obj <- DoubletFinder::doubletFinder_v3(covid_obj, PCs = 1:10, pN = 0.25, pK = pK,
                                              nExp = 0.1, reuse.pANN = FALSE, sct = FALSE)

metadata_names <- colnames(covid_obj@meta.data)
matching_column <- grep("pANN", metadata_names, value = TRUE)
score <- covid_obj@meta.data[, matching_column]

saveRDS(score, paste0(data_path,"/Fig5/covid19naturemedQC_DFscore.rds"))

covid19obj$DFScore <- score

saveRDS(covid19obj, paste0(data_path,"/Fig5/covid19naturemedQC.rds"))

#apply scCLINIC

Name <- "PBMCCOVIDNatMed"
Input <- paste0(data_path,"/Fig5/covid19naturemedQC.rds")
Output <- paste0(data_path,"/Fig5/")
resol=0.8
overlapRatioList=c(0.1,0.25,0.5,0.75,0.9)
OverlapRatio=0.5
ISThreshold=0
gene_n=150

obj <- STEP1A_GlobalMarkers(Input,Output,Name,resol)
obj <- STEP1B_MergingCluster(obj,Output,Name,resol,overlapRatioList,gene_n)
obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold)
saveRDS(obj, paste0(data_path,"/Fig5/step1d.rds"))
STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n)
#PlotContaminationPattern(obj,Output,Name,OverlapRatio)


##### Figure 5
outdir <- paste0(data_path, "/scCLINIC_Figures/")

create_folder_if_not_exists(outdir)

# Define a theme with all required elements
thm <- theme(
  text = element_text(
    size = text_size,
    family = font_type,
    face = "bold"
  ),
  plot.margin = unit(c(0, 0, 0, 0), "cm"),
  # Remove all margins
  plot.title = element_text(
    size = text_size,
    family = font_type,
    face = "bold"
  ),
  plot.subtitle = element_text(size = text_size, family = font_type),
  # Add subtitle element
  plot.tag = element_text(
    size = text_size,
    family = font_type,
    face = "bold"
  ) # Add tag element
)

singletlevel = 4
scCLINICseuobj <-
  readRDS(
    paste0(
      data_path,
      "/Fig5/PBMCCOVIDNatMed_Step2/Output_Overlap_Ratio_0.5/scCLINICResult.rds"
    )
  )

scCLINICseuobj$scCLINIC_artifact <-
  ifelse(as.numeric(levels(scCLINICseuobj$scCLINIC_Level))[scCLINICseuobj$scCLINIC_Level] < singletlevel,
         "Artifact",
         "Singlet")
scCLINICseuobj$scCLINIC_contamID <-
  ifelse(
    scCLINICseuobj$scCLINIC_artifact == "Artifact",
    scCLINICseuobj$scCLINIC_ClusterID,
    "Singlet")

doublet_afqc <-
  readRDS(paste0(data_path, "/Fig5/covid19naturemedQC.rds"))

metatable_unique <-
  readRDS(paste0(data_path, "/Fig5/covid19naturemed.rds"))

scCLINICseuobj@meta.data <-
  cbind(metatable_unique@meta.data, scCLINICseuobj@meta.data)

#7.5% of 44721 is 3354
top_indices <-
  order(doublet_afqc@meta.data$DFScore, decreasing = TRUE)[3354]
cutoff <- doublet_afqc@meta.data$DFScore[top_indices]

scCLINICseuobj$DFafqc_pANN <- doublet_afqc$DFScore
scCLINICseuobj$DFafqc <- NA
scCLINICseuobj$DFafqc <-
  ifelse(doublet_afqc$DFScore < cutoff, "Singlet", "Doublet")

# Replace elements in step1d dataset
scCLINICseuobj$Overlap_Ratio_0.5 <-
  dplyr::recode(
    scCLINICseuobj$Overlap_Ratio_0.5,
    "M1" = "M1_T_NK",
    "M2" = "M2_MoCD14",
    "M3" = "M3_B",
    "M4" = "M4_PB",
    "M5" = "M5_MoCD16",
    "M6" = "M6_CD8T",
    "M7" = "M7_RBC",
    "M8" = "M8_PLT",
    "M9" = "M9_DC",
    "M10" = "M10_Neu",
    "M11" = "M11_Eos",
    "M12" = "M12_pDC",
    "M13" = "M13_dNeu"
  )

scCLINICseuobj$Overlap_Ratio_0.5 <- factor(scCLINICseuobj$Overlap_Ratio_0.5,
                                           levels = unique(scCLINICseuobj$Overlap_Ratio_0.5)[order(as.numeric(gsub(
                                             "M(\\d+)_.*",
                                             "\\1",
                                             unique(scCLINICseuobj$Overlap_Ratio_0.5)
                                           )))])

plot_A0 <-
  wrap_elements(
    DimPlot(scCLINICseuobj, group.by = "Overlap_Ratio_0.5", cols = manuscript_colors) &
      ggtitle(NULL) & #CT.Park
      theme(
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size, family = font_type),
        legend.text = element_text(size = text_size, family = font_type),
        legend.key.size = unit(0.5, 'lines'),
        legend.spacing.x = unit(0.2, "lines"),
        legend.spacing.y = unit(0, "lines")
      ) & labs(color = "Cluster") &
      guides(
        color = guide_legend(
          title.position = "top",
          # Move legend title to the top
          title.hjust = 0,
          # Center the legend title
          ncol = 3,
          # Create three rows for the legend
          byrow = FALSE,
          override.aes = list(size = 3)
        )
      )
  )

#Azimuth annotation
obj <- readRDS(paste0(data_path,"/Fig5/step1d.rds"))
obj_azimuth <- Azimuth::RunAzimuth(obj,reference = "pbmcref")
saveRDS(obj_azimuth, paste0(data_path,"/Fig5/Fig5_Azimuth.rds"))

obj_azimuth <-
  readRDS(paste0(data_path,"/Fig5/Fig5_Azimuth.rds"))

#M8 Platelet

clusterx <- "M8"
file_name <- paste0("Overlap_Ratio_0.5_cluster_", clusterx, ".rds")
recluster <-
  readRDS(
    paste0(
      data_path,"/Fig5/PBMCCOVIDNatMed_Step2/Overlap_Ratio_0.5_recluster/",
      file_name
    )
  )

recluster$scCLINIC_contamID <- scCLINICseuobj$scCLINIC_contamID
recluster$scCLINICScore <- scCLINICseuobj$scCLINICScore

recluster$azimuthl1 <- obj_azimuth$predicted.celltype.l1
recluster$azimuthl2 <- obj_azimuth$predicted.celltype.l2

recluster$DFafqc <- scCLINICseuobj$DFafqc
recluster$DFafqc_pANN <- scCLINICseuobj$DFafqc_pANN

A1 <-
  DimPlot(
    recluster,
    group.by = "scCLINIC_contamID",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = c("#E59CC4", "#53A85F", '#57C3F3', "#9FA3A8")
  ) & labs(color = "scCLINIC") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.6, 0.2),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

A2 <-
  FeaturePlot(
    recluster,
    features = "scCLINICScore",
    raster = FALSE,
    cols = c("#9FA3A8", '#E63863')
  ) & labs(color = "scCLINIC") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.6, 0.25),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.75, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

A3 <-
  DimPlot(
    recluster,
    group.by = "azimuthl1",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = manuscript_colors
  ) & labs(color = "Azimuth") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.4, 0.2),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  ) &
  guides(
    color = guide_legend(
      title.position = "top",
      # Move legend title to the top
      title.hjust = 0,
      # Center the legend title
      ncol = 2,
      # Create three rows for the legend
      byrow = FALSE,
      override.aes = list(size = 3)
    )
  )

A4 <-
  DimPlot(
    recluster,
    group.by = "DFafqc",
    raster = FALSE,
    sizes.highlight = 0.1,
    cols = manuscript_colors
  ) & labs(color = "DF") &
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.6, 0.2),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

contamgeneinfo <-
  read.csv(
    paste0(
      data_path,"/Fig5/PBMCCOVIDNatMed_Step2/Output_Overlap_Ratio_0.5/ArtifactsInfo.csv"
    ),
    na.strings = c("", "NA")
  )
#scCLINIC Subcluster ID
contamgeneinfo$MajorSub <-
  paste0(contamgeneinfo$Major_Cluster, "_", contamgeneinfo$Subcluster)
#Summarize ES score and their source of major cluster for each subclusters
result <- contamgeneinfo %>%
  group_by(MajorSub) %>%
  summarize(
    cluser_reflst = list(Source_of_Artifacts_SoA),
    dp1lst = list(Enrichment_Score)
  ) %>%
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
  summarize(scCLINICScore = mean(scCLINICScore)) %>% #scCLINICScore of each cells within each subclusters are same value... #scCLINICScore
  ungroup()

#Convert to plotting dataframe
heatmap_data <-
  as.data.frame(heatmap_data) %>%  # Convert heatmap_data to a dataframe
  mutate(CS = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score #scCLINIC_ClusterID

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <- data.frame((heatmap_matrix))

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
    geom_bar(stat = "identity", fill = manuscript_colors[i]) + # viridis(length(colnames_sorted))[i]) +
    #geom_bar(stat = "identity", fill = viridis(length(colnames_sorted))[i]) +
    labs(title = paste0(
      dplyr::recode(
        clus,
        "M1" = "T_NK\nM1",
        "M2" = "MoCD14\nM2",
        "M3" = "B\nM3",
        "M4" = "PB\nM4",
        "M5" = "MoCD16\nM5",
        "M6" = "CD8T\nM6",
        "M7" = "RBC\nM7",
        "M8" = "PLT\nM8",
        "M9" = "DC\nM9",
        "M10" = "Neu\nM10",
        "M11" = "Eos\nM11",
        "M12" = "pDC\nM12",
        "M13" = "dNeu\nM13"
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
      plot.margin = unit(c(0, 0, 0, 0), "cm")
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
A5 <-
  wrap_elements(
    wrap_plots(plot_list, ncol = length(plot_list)) +
      plot_layout(guides = "collect") &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  )

plot_A <- wrap_elements(
  wrap_elements(wrap_plots(
    list(A1, A2, A3, A4), nrow = 2, ncol = 2
  ) & theme(plot.margin = unit(c(
    0, 0, 0, 0
  ), "cm"))) +
    plot_annotation(
      title = "M8_PLT_Cluster",
      theme = thm + theme(
        plot.title = element_text(
          face = "bold",
          size = text_size,
          family = font_type,
          hjust = 0.5
        )
      )
    )
)

###Pathway analysis for M2 and M8 (Comparing within COVID19 sample: particular cell type vs others)

library(ggrepel)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(pheatmap)
pacman::p_load(fgsea,tidyverse,pheatmap,openxlsx)
library(Seurat)

Hgmt<-gmtPathways(paste0(data_path,"/h.all.v2023.2.Hs.symbols.gmt"))
#Please download the h.all.v2023.2.Hs.symbols.gmt from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.v2023.2.Hs.symbols.gmt

GOgmt<-gmtPathways(paste0(data_path,"/c5.go.v2023.2.Hs.symbols.gmt"))
#Please download the c5.go.v2023.2.Hs.symbols.gmt from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt

GOgmt <- readRDS(paste0(data_path,"/GOgmt.rds"))

Hgmt <- readRDS(paste0(data_path,"/Hgmt.rds"))

singletlevel = 4

scclinicseuobj <- readRDS(paste0(data_path,"/Fig5/PBMCCOVIDNatMed_Step2/Output_Overlap_Ratio_0.5/scCLINICResult.rds"))

doublet_afqc <- readRDS(paste0(data_path,"/Fig5/covid19naturemedQC.rds"))

metatable_unique<- readRDS(paste0(data_path,"/Fig5/covid19naturemed.rds"))

scclinicseuobj@meta.data <- cbind(metatable_unique@meta.data,scclinicseuobj@meta.data)

#7.5% of 44721 is 3354
top_indices <- order(doublet_afqc@meta.data$DFScore, decreasing = TRUE)[3354]
cutoff <- doublet_afqc@meta.data$DFScore[top_indices]

scclinicseuobj$DFafqc_pANN <- doublet_afqc$DFScore
scclinicseuobj$DFafqc <- NA
scclinicseuobj$DFafqc <- ifelse(doublet_afqc$DFScore < cutoff, "Singlet", "Doublet")

qcsclst <- sort(unique(scclinicseuobj@meta.data[,"Overlap_Ratio_0.5"]))

OriginalALL <- scclinicseuobj

scclinicseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(scclinicseuobj$scCLINIC_Level))[scclinicseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")
scclinicSinglet <- subset(scclinicseuobj,subset = scCLINIC_artifact == "Singlet")

DFSinglet <- subset(scclinicseuobj,subset = DFafqc == "Singlet")



#subset Covid sample
CovidOri <- subset(OriginalALL, subset = Status == "COVID")
CovidDE <- subset(scclinicSinglet, subset = Status == "COVID")
CovidDF <- subset(DFSinglet, subset = Status == "COVID")

cell_type_list = c("M2","M8") #"M2","M8" Monocyte, Platelet

#Calculate marker gene of the particular cell types vs. others (within covid19 samples)
for (cell_type in cell_type_list){
#All COVID
obj <- CovidOri
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

saveRDS(markers,paste0(data_path,"/Fig5/GSEA_Marker_Ori_COVID_",cell_type,".rds"))

#DE COVID
obj <- CovidDE
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
saveRDS(markers,paste0(data_path,"/Fig5/GSEA_Marker_DE_COVID_",cell_type,".rds"))

#DF COVID
obj <- CovidDF
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
saveRDS(markers,paste0(data_path,"/Fig5/GSEA_Marker_DF_COVID_",cell_type,".rds"))

dbused <- "GOgmt"

#CELL TYPE SPECIFIC GSEA
markersscclinic <- readRDS(paste0(data_path,"/Fig5/GSEA_Marker_DE_COVID_",cell_type,".rds"))
markersDF <- readRDS(paste0(data_path,"/Fig5/GSEA_Marker_DF_COVID_",cell_type,".rds"))
markersOri <- readRDS(paste0(data_path,"/Fig5/GSEA_Marker_Ori_COVID_",cell_type,".rds"))

#markers$gene <- rownames(markers)
#markers <- markers[ markers$avg_log2FC >0,]
markersscclinic$stat=markersscclinic$avg_log2FC
markersDF$stat=markersDF$avg_log2FC
markersOri$stat=markersOri$avg_log2FC
#replace Inf with max finite value
markersscclinic$stat<-replace(markersscclinic$stat,markersscclinic$stat>max(markersscclinic$stat[is.finite(markersscclinic$stat)]),max(markersscclinic$stat[is.finite(markersscclinic$stat)])*10)
markersscclinic$stat<-replace(markersscclinic$stat,markersscclinic$stat<min(markersscclinic$stat[is.finite(markersscclinic$stat)]),min(markersscclinic$stat[is.finite(markersscclinic$stat)])*10)

markersDF$stat<-replace(markersDF$stat,markersDF$stat>max(markersDF$stat[is.finite(markersDF$stat)]),max(markersDF$stat[is.finite(markersDF$stat)])*10)
markersDF$stat<-replace(markersDF$stat,markersDF$stat<min(markersDF$stat[is.finite(markersDF$stat)]),min(markersDF$stat[is.finite(markersDF$stat)])*10)

markersOri$stat<-replace(markersOri$stat,markersOri$stat>max(markersOri$stat[is.finite(markersOri$stat)]),max(markersOri$stat[is.finite(markersOri$stat)])*10)
markersOri$stat<-replace(markersOri$stat,markersOri$stat<min(markersOri$stat[is.finite(markersOri$stat)]),min(markersOri$stat[is.finite(markersOri$stat)])*10)


Hlst <- list()
markerssub <- markersOri
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)

Hlst <- c(OriResult=H,Hlst)

markerssub <- markersDF
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)

Hlst <- c(DFResult=H,Hlst)

markerssub <- markersscclinic
rank<-setNames(markerssub$stat,rownames(markerssub))
H<-fgsea(GOgmt,rank,eps= 0.0,nproc=1,nPermSimple = 100000)

Hlst <- c(scclinicResult=H,Hlst)

saveRDS(Hlst,paste0(data_path,"/Fig5/COVID_Result_",dbused,"_",cell_type,".rds"))

}

###Platelet-monocyte aggregation
cell_type <- "M2"

dbused <- "GOgmt"

Hlst <-
  readRDS(
    paste0(
      data_path,"/Fig5/COVID_Result_",
      dbused,
      "_",
      cell_type,
      ".rds"
    )
  )

data <- data.frame(
  Group = c(rep("Ori", length(
    Hlst$OriResult.padj
  )), rep(
    "scclinic", length(Hlst$scclinicResult.padj)
  ), rep("DF", length(Hlst$DFResult.padj))),
  padj = c(
    Hlst$OriResult.padj,
    Hlst$scclinicResult.padj,
    Hlst$DFResult.padj
  ),
  #Value = c(OriH$NES,scclinicH$NES,DFH$NES),
  Gene = c(
    Hlst$OriResult.pathway,
    Hlst$scclinicResult.pathway,
    Hlst$DFResult.pathway
  ),
  NES = c(Hlst$OriResult.NES, Hlst$scclinicResult.NES, Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$scclinicResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$scclinicResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <-
  data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <-
  data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n = 50), "Gene"]

scCLINIC_p <-
  data[data$Group == "scclinic" &
         data$NES > 0 & data$padj < 0.05, ]
scCLINIC_p <-
  scCLINIC_p[head(order(scCLINIC_p$padj), n = 50), "Gene"]

Ori_p <-
  data[data$Group == "Ori" &
         data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n = 50), "Gene"]

pathwaydisplay2 <- unique(c(Ori_p, DF_p, scCLINIC_p))

cell_type <- "M8"
dbused <- "GOgmt"

Hlst <-
  readRDS(
    paste0(
      data_path,"/Fig5/COVID_Result_",
      dbused,
      "_",
      cell_type,
      ".rds"
    )
  )

data <- data.frame(
  Group = c(rep("Ori", length(
    Hlst$OriResult.padj
  )), rep(
    "scclinic", length(Hlst$scclinicResult.padj)
  ), rep("DF", length(Hlst$DFResult.padj))),
  padj = c(
    Hlst$OriResult.padj,
    Hlst$scclinicResult.padj,
    Hlst$DFResult.padj
  ),
  #Value = c(OriH$NES,scclinicH$NES,DFH$NES),
  Gene = c(
    Hlst$OriResult.pathway,
    Hlst$scclinicResult.pathway,
    Hlst$DFResult.pathway
  ),
  NES = c(Hlst$OriResult.NES, Hlst$scclinicResult.NES, Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$scclinicResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$scclinicResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)
data <-
  data[grepl("GOBP", data$Gene, ignore.case = TRUE), ]

DF_p <-
  data[data$Group == "DF" & data$NES > 0 & data$padj < 0.05, ]
DF_p <- DF_p[head(order(DF_p$padj), n = 50), "Gene"]

scCLINIC_p <-
  data[data$Group == "scclinic" &
         data$NES > 0 & data$padj < 0.05, ]
scCLINIC_p <-
  scCLINIC_p[head(order(scCLINIC_p$padj), n = 50), "Gene"]

Ori_p <-
  data[data$Group == "Ori" &
         data$NES > 0 & data$padj < 0.05, ]
Ori_p <- Ori_p[head(order(Ori_p$padj), n = 50), "Gene"]

pathwaydisplay3 <-
  unique(c(Ori_p, DF_p, scCLINIC_p))

overlap <- intersect(pathwaydisplay2, pathwaydisplay3)
print(overlap)

Hlst <-
  readRDS(
    paste0(
      data_path,"/Fig5/COVID_Result_",
      dbused,
      "_",
      cell_type,
      ".rds"
    )
  )

data <- data.frame(
  Group = c(rep("Ori", length(
    Hlst$OriResult.padj
  )), rep(
    "scclinic", length(Hlst$scclinicResult.padj)
  ), rep("DF", length(Hlst$DFResult.padj))),
  padj = c(
    Hlst$OriResult.padj,
    Hlst$scclinicResult.padj,
    Hlst$DFResult.padj
  ),
  #Value = c(OriH$NES,scclinicH$NES,DFH$NES),
  Gene = c(
    Hlst$OriResult.pathway,
    Hlst$scclinicResult.pathway,
    Hlst$DFResult.pathway
  ),
  NES = c(Hlst$OriResult.NES, Hlst$scclinicResult.NES, Hlst$DFResult.NES)
  
  #Lead = c(list(Hlst$OriResult.leadingEdge),list(Hlst$scclinicResult.leadingEdge),list(Hlst$DFResult.leadingEdge))
  #NGene = c(length(Hlst$OriResult.leadingEdge),length(Hlst$scclinicResult.leadingEdge),length(Hlst$DFResult.leadingEdge))
)

data <- data[data$Gene %in% unique(overlap),]
data$pstat <- 1  # Initialize the new column
data$pstat[data$padj <= 0.05] <- 2
data$pstat[data$padj <= 0.01] <- 3
data$pstat[data$padj <= 0.001] <- 4
data$pstat[is.na(data$padj)] <- NA

NES_break <- c(-4,-2, 0, 2, 4)
# Convert pstat to a factor
data$pstat <- factor(data$pstat)

# Modify the size breaks to match the levels of pstat
data$Gene_short <- data$Gene
data <- data[data$Group != "scclinic+",]
data$Gene_short <-
  factor(data$Gene_short, levels = overlap)
data <- data[order(data$Group, data$Gene_short), ]
substitute_underscores <- function(string) {
  string <- tolower(string)
  gsub("(([^_]*_){3})", "\\1\n", string, perl = TRUE)
}

# Apply the function to Gene_short
data$Gene_short <-
  substitute_underscores(data$Gene_short)
data$Gene_short <-
  factor(data$Gene_short, levels = unique(data$Gene_short))
C1 <-
  ggplot(data = data, aes(
    x = -log10(padj),
    y = Gene_short,
    fill = Group
  )) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = manuscript_colors,
                    breaks = c("scclinic", "DF", "Ori"),##added
                    labels = c("scCLINIC", "DoubletFinder", "All cells")) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = "dotted",
    color = "black",
    size = 1
  ) +  # Add dotted vertical line
  theme_bw() +
  ylab("") +
  xlab("-log10(padj)") +
  ggtitle("") +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = text_size, family = font_type),
    axis.text.x = element_text(size = text_size, family = font_type),
    axis.title = element_text(size = text_size, family = font_type),
    plot.title = element_text(size = text_size, family = font_type),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = c(0.7, 0.1),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

C2 <-
  ggplot(data = data, aes(x = NES, y = Gene_short, fill = Group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = manuscript_colors,
                    breaks = c("scclinic", "DF", "Ori"),##added
                    labels = c("scCLINIC", "DF", "All cells")) +
  theme_bw() +
  ylab("") +
  xlab("Normalized Enrichment Score") +
  scale_x_continuous(breaks = seq(0, max(data$NES), by = 1)) + # Set x-axis breaks
  ggtitle("") + # No title, can be adjusted if needed
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = text_size, family = font_type),
    axis.title = element_text(size = text_size, family = font_type),
    plot.title = element_text(size = text_size, family = font_type),
    axis.text.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "None"
  )

plot_C <- wrap_elements(
  wrap_elements(
    (C1 + C2) +
      plot_layout(widths = c(0.65, 0.35), ncol = 2) &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  ) +
    plot_annotation(
      title = "Pathway confounded by artifact in M8_PLT COVID-19 Sample",
      theme = thm + theme(
        plot.title = element_text(
          face = "bold",
          size = text_size,
          family = font_type,
          hjust = 0.5
        )
      )
    )
)

plot_list_2 <- list()
for (curgene in c("FCN1", "CD14", "AIF1", "MS4A1", "TRAC")) {
  plot_gene <-
    FeaturePlot(
      recluster,
      reduction = "umap",
      features = curgene,
      pt.size = 0.1
    ) &
    scale_color_gradientn(
      colors = c("grey", "#E63863"),
      limits = c(0, 5),
      oob = scales::squish
    ) &
    theme(
      plot.title = element_text(
        face = "italic",
        size = text_size,
        family = font_type
      ),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "None",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  
  plot_list_2[[paste0("M8_", curgene)]] <- plot_gene
}

for (curgene in c("FCN1", "CD14", "AIF1", "MS4A1", "TRAC")) {
  plot_gene <-
    FeaturePlot(
      scCLINICseuobj,
      reduction = "umap",
      features = curgene,
      pt.size = 0.1
    ) &
    scale_color_gradientn(
      colors = c("grey", "#E63863"),
      limits = c(0, 5),
      oob = scales::squish
    ) &
    theme(
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "None",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  
  plot_list_2[[paste0("Global_", curgene)]] <- plot_gene
}

# Create a separate plot to extract the legend
legend_plot <-
  FeaturePlot(
    recluster,
    reduction = "umap",
    features = curgene,
    pt.size = 0.1
  ) &
  scale_color_gradientn(
    colors = c("grey", "#E63863"),
    limits = c(0, 5),
    oob = scales::squish
  ) &
  theme_void() &
  theme(legend.position = "right", plot.margin = unit(c(0, 0, 0, 0), "cm")) # Only keep the legend

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

plot_B <- wrap_elements(
  wrap_elements(wrap_plots(
    list(
      plot_list_2$M8_FCN1,
      plot_list_2$M8_CD14,
      plot_list_2$M8_AIF1,
      plot_list_2$Global_FCN1,
      plot_list_2$Global_CD14,
      plot_list_2$Global_AIF1
    ),
    nrow = 2,
    ncol = 3
  ) &
    theme(plot.margin = unit(c(
      0, 0, 0.5, 0
    ), "cm"))) + wrap_elements(legend) + plot_layout(widths = c(1, 0.1)) +
    plot_annotation(
      title = "Monocyte Features in M8_PLT Cluster",
      theme = thm + theme(
        plot.title = element_text(
          face = "bold",
          size = text_size,
          family = font_type,
          hjust = 0.5
        )
      )
    )
)

plot_D <- wrap_elements(
  wrap_elements(wrap_plots(
    list(
      plot_list_2$M8_MS4A1,
      plot_list_2$M8_TRAC,
      plot_list_2$Global_MS4A1,
      plot_list_2$Global_TRAC
    ),
    nrow = 2,
    ncol = 2
  ) &
    theme(plot.margin = unit(c(
      0, 0, 0.5, 0
    ), "cm"))) + wrap_elements(legend) + plot_layout(widths = c(1, 0.1)) +
    plot_annotation(
      title = "Lymphocytes Features in M8_PLT Cluster",
      theme = thm + theme(
        plot.title = element_text(
          face = "bold",
          size = text_size,
          family = font_type,
          hjust = 0.5
        )
      )
    )
)

figure_5 <- ((((plot_A0 | plot_A | A5) +
                 plot_layout(widths = c(0.2, 0.25, 0.55)))
/ ((plot_C | plot_B | plot_D) +
     plot_layout(widths = c(0.5, 0.3, 0.2))
) +
  plot_layout(heights = c(0.5, 0.5))) + # Nest plot_B and plot_C/plot_D
  plot_layout(widths = c(0.2, 0.8))) + # Set widths for the columns
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

figure_5 <-
  figure_5 + plot_annotation(tag_levels = "A", theme = thm)
ggsave(
  filename = paste0(outdir, "Figure_5.png"),
  plot = figure_5,
  height = 11,
  width = 18,
  dpi = 300
)
