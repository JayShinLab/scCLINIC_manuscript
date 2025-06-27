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

create_folder_if_not_exists(paste0(data_path,"/Fig4/"))

#Please download the genes, barcodes, and matrix files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129363
#Reference https://www.nature.com/articles/s42255-019-0152-6 Nature Metabolism

library(ggrepel)
library(dplyr)
library(scCustomize)
library(enrichR)
library(fgsea)
library(tidyverse)
library(openxlsx)

Hgmt<-gmtPathways(paste0(data_path,"/h.all.v2023.2.Hs.symbols.gmt"))
#Please download the h.all.v2023.2.Hs.symbols.gmt from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.v2023.2.Hs.symbols.gmt

GOgmt<-gmtPathways(paste0(data_path,"/c5.go.v2023.2.Hs.symbols.gmt"))
#Please download the c5.go.v2023.2.Hs.symbols.gmt from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt

saveRDS(GOgmt,paste0(data_path,"/Fig4/GOgmt.rds"))
saveRDS(Hgmt,paste0(data_path,"/Fig4/Hgmt.rds"))

adinat <- Read10X(paste0(data_path)) 

adinatobj <- CreateSeuratObject(adinat)

adinatobj$SampleIDX <- sapply(strsplit(rownames(adinatobj@meta.data), "-"), "[[", 2)
adinatobj$nCount_RNA <- colSums(x = adinatobj, slot = "counts")  # nCount_RNA
adinatobj$nFeatureRNA <- colSums(x = GetAssayData(object = adinatobj, slot = "counts") > 0)
adinatobj[["percent.mt"]] <- PercentageFeatureSet(adinatobj, pattern = "^MT-")

adinatobjQC <- subset(adinatobj, 
                      nCount_RNA >= 200 & nFeatureRNA >= 200 & percent.mt <= 20)#following the QC step in original publication


library(DoubletFinder)

subobj <- adinatobjQC

#DF
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu_adi <- subobj
seu_adi <- NormalizeData(seu_adi)
seu_adi <- FindVariableFeatures(seu_adi, selection.method = "vst", nfeatures = 2000)
seu_adi <- ScaleData(seu_adi)
seu_adi <- RunPCA(seu_adi)
seu_adi <- RunUMAP(seu_adi, dims = 1:10)

# Remove columns with any NA values
seu_adi@meta.data <- seu_adi@meta.data[, colSums(is.na(seu_adi@meta.data)) == 0]
# Convert all factor columns to characters
seu_adi@meta.data[] <- lapply(seu_adi@meta.data, function(x) {
  if (is.factor(x)) as.character(x) else x
})

sweep.vector <- DoubletFinder::paramSweep_v3(seu_adi, PCs = 1:10, sct = FALSE)  #remove _v3 ##Seurat 4 need use v3
sweep.table <- DoubletFinder::summarizeSweep(sweep.vector, GT = FALSE)
bcmvn <- DoubletFinder::find.pK(sweep.table)

#sapply(seu_adi@meta.data, class) to see if any col is not factor, numeric, integer, character

pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
pK <- as.numeric(levels(pK))[pK]
#nExp <- round(ncol(seu_adi) * 0.1)
seu_adi <- DoubletFinder::doubletFinder_v3(seu_adi, PCs = 1:10, pN = 0.25, pK = pK,
                                              nExp = 0.1, reuse.pANN = FALSE, sct = FALSE)   #remove _v3 ##Seurat 4 need change reuse.pANN to FALSE and use v3

metadata_names <- colnames(seu_adi@meta.data)

matching_column <- grep("pANN", metadata_names, value = TRUE)

score <- seu_adi@meta.data[, matching_column]

saveRDS(score, paste0(data_path,"/Fig4/adinatobjQC_DFscore.rds"))

adinatobjQC$DFScore <- score

saveRDS(adinatobjQC, paste0(data_path,"/Fig4/adinatobjQC.rds"))


#scCLINIC

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

Name <- "AdiposeNatMet"
Input <- paste0(data_path,"/Fig4/adinatobjQC.rds")
Output <- paste0(data_path,"/Fig4/")
resol=0.8
overlapRatioList=c(0.1,0.25,0.5,0.75,0.9)
OverlapRatio=0.5
ISThreshold=0
gene_n=150

obj <- STEP1A_GlobalMarkers(Input,Output,Name,resol)
obj <- STEP1B_MergingCluster(obj,Output,Name,resol,overlapRatioList,gene_n)
obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold)
saveRDS(obj,paste0(data_path,"/Fig4/step1d.rds"))
STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n)
obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n)
PlotContaminationPattern(obj,Output,Name,OverlapRatio)


##### Figure 4
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

singletlevel = 6
scCLINICseuobj <-
  readRDS(
    paste0(
      data_path,"/Fig4/AdiposeNatMet_Step2/Output_Overlap_Ratio_0.5/scCLINICResult.rds"
    )
  )

scCLINICseuobj$scCLINIC_artifact <-
  ifelse(test = as.numeric(levels(scCLINICseuobj$scCLINIC_Level))[scCLINICseuobj$scCLINIC_Level] < singletlevel,
         yes =  "Artifact",
         no =  "Singlet")
scCLINICseuobj$scCLINIC_contamID <-
  ifelse(test = scCLINICseuobj$scCLINIC_artifact == "Artifact",
    yes = scCLINICseuobj$scCLINIC_ClusterID,
    no = "Singlet")

# Replace elements in step1d dataset
scCLINICseuobj$Overlap_Ratio_0.5 <-
  dplyr::recode(
    scCLINICseuobj$Overlap_Ratio_0.5,
    "M1" = "M1_ASPC",
    "M2" = "M2_ASPC",
    "M3" = "M3_Mp_DC",
    "M4" = "M4_NK_B_T",
    "M5" = "M5_EC",
    "M6" = "M6_SMC",
    "M7" = "M7_LEC",
    "M8" = "M8_NK_B_T"
  )

plot_A <-
  DimPlot(scCLINICseuobj, group.by = "Overlap_Ratio_0.5", cols = manuscript_colors) &
  ggtitle(NULL) & #CT.Park
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.05, 0.1),
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
      ncol = 2,
      # Create three rows for the legend
      byrow = FALSE,
      override.aes = list(size = 3)
    )
  )



# Create a separate plot to extract the legend
legend_plot <-
  DimPlot(scCLINICseuobj, group.by = "Overlap_Ratio_0.5", cols = manuscript_colors) &
  theme_void() &
  theme(
    legend.position = "right",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  ) & labs(color = "Cluster")# Only keep the legend

# Extract the legend
legend_A <- cowplot::get_legend(legend_plot)

#M4
#Azimuth
obj <- readRDS(paste0(data_path,"/Fig4/step1d.rds"))

obj_azimuth <- Azimuth::RunAzimuth(obj,reference = "adiposeref")
saveRDS(obj_azimuth, paste0(data_path,"/Fig4/Fig4_Azimuth.rds"))
obj_azimuth <-
  readRDS(paste0(data_path, "/Fig4/Fig4_Azimuth.rds"))

clusterx <- "M4"
file_name <- paste0("Overlap_Ratio_0.5_cluster_", clusterx, ".rds")
recluster <-
  readRDS(
    paste0(
      data_path,"/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/",
      file_name
    )
  )

recluster$scCLINIC_contamID <- scCLINICseuobj$scCLINIC_contamID

recluster$azimuthl1 <- obj_azimuth$predicted.celltype.l1
recluster$azimuthl2 <- obj_azimuth$predicted.celltype.l2


B1 <-
  DimPlot(
    recluster,
    group.by = "scCLINIC_contamID",
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
    legend.position = "right",
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )

contamgeneinfo <-
  read.csv(
    paste0(
      data_path,"/Fig4/AdiposeNatMet_Step2/Output_Overlap_Ratio_0.5/ArtifactsInfo.csv"
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
        "M1" = "ASPC\nM1",
        "M2" = "ASPC\nM2",
        "M3" = "Mp_DC\nM3",
        "M4" = "NK_B_T\nM4",
        "M5" = "EC\nM5",
        "M6" = "SMC\nM6",
        "M7" = "LEC\nM7",
        "M8" = "NK_B_T\nM8"
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
B3 <-
  wrap_elements(
    wrap_plots(plot_list, ncol = length(plot_list)) +
      plot_layout(guides = "collect") &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  )


plot_B <- wrap_elements(
  wrap_elements(
    (B1 + B2) +
      plot_layout(widths = c(0.5, 0.5), ncol = 2) &
      theme(plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
  ) +
    plot_annotation(
      title = "M4_NK_B_T_Cluster",
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

####Pathway analysis for M4
cell_type = "M4" 

library(ggrepel)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(pheatmap)
pacman::p_load(fgsea,tidyverse,pheatmap,openxlsx)
library(Seurat)

singletlevel = 6

scclinicseuobj <- readRDS(paste0(data_path,"/Fig4/AdiposeNatMet_Step2/Output_Overlap_Ratio_0.5/scCLINICResult.rds"))

doublet_afqc <- readRDS(paste0(data_path,"/Fig4/adinatobjQC.rds"))

metatable_unique <- read.table(file = paste0(data_path, "/ClinicalMetaData.csv"), sep = ",", header = TRUE) ####ClinicalMetaData.csv need to be in reproduce folder

scclinicseuobj$SampleIDX <- sapply(strsplit(rownames(scclinicseuobj@meta.data), "-"), "[[", 2)

indices <- match(scclinicseuobj@meta.data$SampleIDX, metatable_unique$x.SampleIDX)

# Assign the values from metatable_unique to the corresponding rows in scclinicseuobj@meta.data
scclinicseuobj@meta.data$ncells <- metatable_unique$ncells[indices]
scclinicseuobj@meta.data$Condition <- metatable_unique$y.Condition[indices]
scclinicseuobj@meta.data$Tissue <- metatable_unique$y.Tissue[indices]
scclinicseuobj@meta.data$SampleName <- metatable_unique$y.SampleName[indices]


#7.5% of 26540 is 1990
top_indices <- order(doublet_afqc@meta.data$DFScore, decreasing = TRUE)[1190]
cutoff <- doublet_afqc@meta.data$DFScore[top_indices]

scclinicseuobj$DFafqc_pANN <- doublet_afqc$DFScore
scclinicseuobj$DFafqc <- NA
scclinicseuobj$DFafqc <- ifelse(doublet_afqc$DFScore < cutoff, "Singlet", "Doublet")

qcsclst <- sort(unique(scclinicseuobj@meta.data[,"Overlap_Ratio_0.5"]))

OriginalALL <- scclinicseuobj

scclinicseuobj$scCLINIC_artifact <- ifelse(as.numeric(levels(scclinicseuobj$scCLINIC_Level))[scclinicseuobj$scCLINIC_Level] < singletlevel, "Artifact", "Singlet")

scclinicSinglet <- subset(scclinicseuobj,subset = scCLINIC_artifact == "Singlet")


DFSinglet <- subset(scclinicseuobj,subset = DFafqc == "Singlet")

OriginalALL$ConditionTissue <- paste0(OriginalALL$Condition, OriginalALL$Tissue)

#Diabetic + SAT
DiaSATOri <- subset(OriginalALL, subset = Condition == "Diabetic" & Tissue == "SAT")
DiaSATDE <- subset(scclinicSinglet, subset = Condition == "Diabetic" & Tissue == "SAT")
DiaSATDF <- subset(DFSinglet, subset = Condition == "Diabetic" & Tissue == "SAT")
#NonDiabetic + SAT
NonSATOri <- subset(OriginalALL, subset = Condition == "NonDiabetic" & Tissue == "SAT")
NonSATDE <- subset(scclinicSinglet, subset = Condition == "NonDiabetic" & Tissue == "SAT")
NonSATDF <- subset(DFSinglet, subset = Condition == "NonDiabetic" & Tissue == "SAT")
#Diabetic + VAT
DiaVATOri <- subset(OriginalALL, subset = Condition == "Diabetic" & Tissue == "VAT")
DiaVATDE <- subset(scclinicSinglet, subset = Condition == "Diabetic" & Tissue == "VAT")
DiaVATDF <- subset(DFSinglet, subset = Condition == "Diabetic" & Tissue == "VAT")
#NonDiabetic + VAT
NonVATOri <- subset(OriginalALL, subset = Condition == "NonDiabetic" & Tissue == "VAT")
NonVATDE <- subset(scclinicSinglet, subset = Condition == "NonDiabetic" & Tissue == "VAT")
NonVATDF <- subset(DFSinglet, subset = Condition == "NonDiabetic" & Tissue == "VAT")


#Compare DiaSAT vs NonSAT
#Batch Script start Fig4_published_batch1.R
subDiaSATOri <- subset(DiaSATOri,subset = Overlap_Ratio_0.5 == cell_type)
subDiaSATDE <- subset(DiaSATDE,subset = Overlap_Ratio_0.5 == cell_type)
subDiaSATDF <- subset(DiaSATDF,subset = Overlap_Ratio_0.5 == cell_type)
subNonSATOri <- subset(NonSATOri,subset = Overlap_Ratio_0.5 == cell_type)
subNonSATDE <- subset(NonSATDE,subset = Overlap_Ratio_0.5 == cell_type)
subNonSATDF <- subset(NonSATDF,subset = Overlap_Ratio_0.5 == cell_type)

#calculate marker gene for all cells diabetic SAT vs non-diabetic SAT
mergeseurat <- Merge_Seurat_List(c(subDiaSATOri,subNonSATOri),c("DiaSATOri","NonSATOri"))
mergeseurat <- DietSeurat(mergeseurat, counts = TRUE, data = TRUE, scale.data = FALSE)
mergeseurat$Group <- sub("_.*", "", rownames(mergeseurat@meta.data))
genes.percent.expression <- rowMeans(mergeseurat[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(object = mergeseurat) <- "Group"
#mergeseurat <- JoinLayers(mergeseurat) ###Seurat 4 remove this
markers<-FindMarkers(mergeseurat,ident.1 ="DiaSATOri", ident.2 = "NonSATOri", 
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)  ##Output and rank all filtered genes
saveRDS(markers,paste0(data_path,"/Fig4/GSEA_Marker_DiaSATOrivsNonSATOri_",cell_type,".rds"))

#calculate marker gene for scCLINIC Diabetic SAT vs Non Diabetic SAT
mergeseurat <- Merge_Seurat_List(c(subDiaSATDE,subNonSATDE),c("DiaSATDE","NonSATDE"))
mergeseurat <- DietSeurat(mergeseurat, counts = TRUE, data = TRUE, scale.data = FALSE)
mergeseurat$Group <- sub("_.*", "", rownames(mergeseurat@meta.data))
genes.percent.expression <- rowMeans(mergeseurat[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(object = mergeseurat) <- "Group"
#mergeseurat <- JoinLayers(mergeseurat) ###Seurat 4 remove this
markers<-FindMarkers(mergeseurat,ident.1 ="DiaSATDE", ident.2 = "NonSATDE", 
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)   ##Output and rank all filtered genes
saveRDS(markers,paste0(data_path,"/Fig4/GSEA_Marker_DiaSATDEvsNonSATDE_",cell_type,".rds"))

#calculate marker gene for DoubletFinder Diabetic SAT vs Non Diabetic SAT
mergeseurat <- Merge_Seurat_List(c(subDiaSATDF,subNonSATDF),c("DiaSATDF","NonSATDF"))
mergeseurat <- DietSeurat(mergeseurat, counts = TRUE, data = TRUE, scale.data = FALSE)
mergeseurat$Group <- sub("_.*", "", rownames(mergeseurat@meta.data))
genes.percent.expression <- rowMeans(mergeseurat[["RNA"]]$counts>0 )*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])
Idents(object = mergeseurat) <- "Group"
#mergeseurat <- JoinLayers(mergeseurat) ###Seurat 4 remove this
markers<-FindMarkers(mergeseurat,ident.1 ="DiaSATDF", ident.2 = "NonSATDF", 
                     min.cells.group = 0,
                     min.cells.feature = 0,
                     min.pct = 0,
                     logfc.threshold = 0,
                     return.thresh = Inf,
                     only.pos = FALSE,
                     features = genes.filter)   ##Output and rank all filtered genes
saveRDS(markers,paste0(data_path,"/Fig4/GSEA_Marker_DiaSATDFvsNonSATDF_",cell_type,".rds"))

dbused <- "GOgmt"

#CELL TYPE SPECIFIC GSEA
markersscclinic <- readRDS(paste0(data_path,"/Fig4/GSEA_Marker_DiaSATDEvsNonSATDE_",cell_type,".rds"))
markersDF <- readRDS(paste0(data_path,"/Fig4/GSEA_Marker_DiaSATDFvsNonSATDF_",cell_type,".rds"))
markersOri <- readRDS(paste0(data_path,"/Fig4/GSEA_Marker_DiaSATOrivsNonSATOri_",cell_type,".rds"))

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

saveRDS(Hlst,paste0(data_path,"/Fig4/DiaSATvsNonSAT_Result_",dbused,"_",cell_type,".rds"))


dbused <- "GOgmt"
cell_type <- "M4" #NK T B Cell
overlap <-
  c(
    "GOBP_RESPONSE_TO_ZINC_ION",
    "GOBP_RESPONSE_TO_COPPER_ION",
    "GOBP_INTRACELLULAR_ZINC_ION_HOMEOSTASIS"
  )

Hlst <-
  readRDS(
    paste0(
      data_path,"/Fig4/DiaSATvsNonSAT_Result_",
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
    "scclinic", length(Hlst$scclinicResult.padj)    #####scCLINIC
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
  scale_x_continuous(limits = c(0, 10)) +
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
    legend.position = c(0.8, 0.3),
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


plot_C1 <- wrap_elements((C1 + C2) +
                           plot_layout(widths = c(0.65, 0.35), ncol = 2) &
                           theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))

dbused <- "GOgmt"
cell_type <- "M4" #NK T B Cell
overlap <-
  c(
    "GOBP_LIPID_LOCALIZATION",
    "GOBP_FATTY_ACID_TRANSPORT",
    "GOBP_ERK1_AND_ERK2_CASCADE",
    "GOBP_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
    "GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM",
    "GOBP_VASCULATURE_DEVELOPMENT",
    "GOBP_TUBE_MORPHOGENESIS"
  )

Hlst <-
  readRDS(
    paste0(
      data_path,"/Fig4/DiaSATvsNonSAT_Result_",
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
C3 <-
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
    legend.position = c(0.8, 0.15),
    legend.title = element_text(size = text_size, family = font_type),
    legend.text = element_text(size = text_size, family = font_type),
    legend.key.size = unit(0.5, 'lines'),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0, "lines")
  )


C4 <-
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

plot_C2 <- wrap_elements((C3 + C4) +
                           plot_layout(widths = c(0.65, 0.35), ncol = 2) &
                           theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))

plot_C <-
  wrap_elements(plot_C2 / plot_C1 + plot_layout(heights = c(0.6, 0.4)) +
                  theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")))


clusterx <- cell_type
#Positive_regulation_of_immune_system_process_plot <- c("MT2A" ,  "MT1E"   ,"MT1M","APOD","PLA2G2A","DCN")
recluster <-
  readRDS(
    paste0(
      data_path,"/Fig4/AdiposeNatMet_Step2/Overlap_Ratio_0.5_recluster/Overlap_Ratio_0.5_cluster_",
      clusterx,
      ".rds"
    )
  )

plot_list_2 <- list()
for (curgene in c("APOD", "PLA2G2A", "DCN", "MT2A", "MT1E", "MT1M")) {
  plot_gene <-
    FeaturePlot(
      recluster,
      reduction = "umap",
      features = curgene,
      pt.size = 0.1
    ) &
    scale_color_gradientn(
      colors = c("grey", "#E63863"),
      limits = c(0, 6),
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
  
  plot_list_2[[paste0(curgene)]] <- plot_gene
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
    limits = c(0, 6),
    oob = scales::squish
  ) &
  theme_void() &
  theme(legend.position = "right", plot.margin = unit(c(0, 0, 0, 0), "cm")) # Only keep the legend

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

plot_D <- wrap_elements(
  wrap_elements(
    wrap_plots(plot_list_2, nrow = 2, ncol = 3) &
      theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
  ) + wrap_elements(legend) + plot_layout(widths = c(1, 0.1)) +
    plot_annotation(
      title = "M4_NK_B_T_Cluster",
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

metatable_unique <-
  read.table(paste0(data_path,"/ClinicalMetaData.csv"), ###in reproduce
             sep = ",")

scCLINICseuobj$SampleIDX <-
  sapply(strsplit(rownames(scCLINICseuobj@meta.data), "-"), "[[", 2)
indices <-
  match(scCLINICseuobj@meta.data$SampleIDX,
        metatable_unique$x.SampleIDX)

# Assign the values from metatable_unique to the corresponding rows in scCLINICseuobj@meta.data
scCLINICseuobj@meta.data$Condition <-
  metatable_unique$y.Condition[indices]
scCLINICseuobj@meta.data$Tissue <-
  metatable_unique$y.Tissue[indices]
scCLINICseuobj$ConditionTissue <-
  paste0(scCLINICseuobj$Condition, " (", scCLINICseuobj$Tissue, ")")
recluster$ConditionTissue <-
  scCLINICseuobj$ConditionTissue
plot_E <-
  wrap_elements(
    DimPlot(
      recluster,
      group.by = "ConditionTissue",
      raster = FALSE,
      sizes.highlight = 0.1,
      cols = manuscript_colors
    ) & labs(color = "Status & Distribution") &
      theme(
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.45, 0.9),
        legend.title = element_text(size = text_size, family = font_type),
        legend.text = element_text(size = text_size, family = font_type),
        legend.key.size = unit(0.5, 'lines'),
        legend.spacing.x = unit(0.2, "lines"),
        legend.spacing.y = unit(0, "lines"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
      )
  )

figure_4 <-
  (((plot_A / plot_E + plot_layout(heights = c(
    0.6, 0.4
  ))) |
    ((plot_B |
        B3) / (plot_C |
                 plot_D) + plot_layout(heights = c(0.3, 0.7))
    )) + # Nest plot_B and plot_C/plot_D
    plot_layout(widths = c(0.2, 0.8))) + # Set widths for the columns
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

figure_4 <-
  figure_4 + plot_annotation(tag_levels = "A", theme = thm)
ggsave(
  filename = paste0(outdir, "Figure_4.png"),
  plot = figure_4,
  height = 10,
  width = 18,
  dpi = 300
)
  