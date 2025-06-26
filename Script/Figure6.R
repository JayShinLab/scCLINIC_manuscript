library(hdf5r)

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
library(ggforce)
library(ggrepel)
library(grid) #for gpar() FigS6I
library(arrow, lib.loc = "/mnt/software/R/4.3.0/lib64/R/library")


localdir <- "~/testVisiumHDHumanCRC/human_coloncancer_ffpe/output_and_supplement_files/binned_outputs/square_016um/"
CRC <- Load10X_Spatial(data.dir = localdir, filename = "filtered_feature_bc_matrix.h5")

output_dir <- "/mnt/lab-store/projects/scCLINIC/Reproduce/TXY/Fig6/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)  # Create the directory if it doesn't exist
}
output_dir <- "/mnt/lab-store/projects/scCLINIC/Reproduce/TXY/Fig6/"
saveRDS(CRC,paste0(output_dir, "CRC016_raw.rds"))

#==================Preprocessing======================================================

CRC <- readRDS(paste0(output_dir, "CRC016_raw.rds"))
#Spatial Normalisation
CRC_SpatialQC <- NormalizeData(CRC)

#Unsupervised clustering 
# note that data is already normalized
DefaultAssay(CRC_SpatialQC) <- "Spatial"
CRC_SpatialQC <- FindVariableFeatures(CRC_SpatialQC)
CRC_SpatialQC <- ScaleData(CRC_SpatialQC)

# we select 50,0000 cells and create a new 'sketch' assay
CRC_SpatialQC <- SketchData(
  object = CRC_SpatialQC,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(CRC_SpatialQC) <- "sketch"

# perform clustering workflow
CRC_SpatialQC <- FindVariableFeatures(CRC_SpatialQC)
CRC_SpatialQC <- ScaleData(CRC_SpatialQC)
CRC_SpatialQC <- RunPCA(CRC_SpatialQC, assay = "sketch", reduction.name = "pca.sketch")
CRC_SpatialQC <- FindNeighbors(CRC_SpatialQC, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
CRC_SpatialQC <- FindClusters(CRC_SpatialQC, cluster.name = "seurat_cluster.sketched", resolution = 3)
CRC_SpatialQC <- RunUMAP(CRC_SpatialQC, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

options(future.globals.maxSize= 922746880) #; 880x1024x1024
#Projecting back to entire dataset
CRC_SpatialQCprojected <- ProjectData(
  object = CRC_SpatialQC,
  assay = "Spatial",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
saveRDS(CRC_SpatialQCprojected, paste0(output_dir, "CRC_projected.rds"))


#===============Cell type annotation=================================================
ref <- readRDS("~/testscC/OutputHumanCRC/online_scRNAseq_dataset/ValdeolivasA_ProfilingheterogeneityofCRCpaper/scRNAseq_dataset.rds")  
CRC <- readRDS(paste0(output_dir, "CRC_projected.rds"))


library(spacexr)
Idents(ref) <- "Cell_subtype"    

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$Cell_subtype)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- CRC[["Spatial"]]$counts
CRC_cells_hd <- colnames(CRC[["Spatial"]])
coords <- GetTissueCoordinates(CRC)[CRC_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
# Filter cell types with at least 25 cells
cell_counts <- table(reference@cell_types)
valid_cell_types <- names(cell_counts[cell_counts >= 25])
valid_cells <- which(reference@cell_types %in% valid_cell_types)
filtered_counts <- reference@counts[, valid_cells]
filtered_nUMI <- reference@nUMI[valid_cells]
filtered_cell_types <- reference@cell_types[valid_cells]
filtered_cell_types <- droplevels(filtered_cell_types)

# Re-create the Reference object with filtered data
reference_filtered <- Reference(counts = filtered_counts,
                                cell_types = filtered_cell_types,  # Dropped levels should ensure only valid types remain
                                nUMI = filtered_nUMI)

RCTD <- create.RCTD(query, reference_filtered, max_cores = 24)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD, file = paste0(output_dir, "myRCTD_FULL_CRC_withscCLINIC.rds"))

# add results back to Seurat object
CRC <- AddMetaData(CRC, metadata = RCTD@results$results_df)
#saveRDS(CRC,file = "~/testscC/OutputHumanCRC/Optimizingsubclusters/scCLINIC_subclusteroptimization_2Aresol0.8_Step2/postscClinic_downstream/FULL_CRC_withscCLINIC_withRCTD.rds")

# project RCTD labels from sketched cortical cells to all cortical cells
CRC$first_type <- as.character(CRC$first_type)
CRC$first_type[is.na(CRC$first_type)] <- "Unknown"

options(future.globals.maxSize= 891289600)
CRC <- ProjectData(
  object = CRC,
  assay = "Spatial",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)
saveRDS(CRC,file = paste0(output_dir, "CRC_projected_celltypelabelled.rds"))

#========================scCLINIC_userannotation=================================================
Name <- "scCLINIC_userannotation"
Output <- output_dir
resol="Manual"
OverlapRatio="merged_cell_type"
ISThreshold = 0
gene_n=150
CELLANNOTATION = TRUE

obj <- readRDS(paste0(output_dir, "CRC_projected_celltypelabelled.rds"))

merge_mapping <- c(
  "Unknown" = "Low-quality cells",
  "CMS2" = "CRC",
  "Stromal 2" = "Mesenchymal cells",
  "Stem-like/TA" = "Epithelial cells",
  "Stem-like-TA" = "Epithelial cells", # In case it still exists
  "Enteric glial cells" = "Enteric glial cells",
  "Mature Enterocytes type 2" = "Epithelial cells",
  "Mature Enterocytes type 1" = "Epithelial cells",
  "Pro-inflammatory" = "Immune cells",     #Pro-inflammatory Macrophages
  "Myofibroblasts" = "Mesenchymal cells",
  "Smooth muscle cells" = "Mesenchymal cells",
  "Goblet cells" = "Epithelial cells", 
  "Stalk-like ECs" = "Mesenchymal cells",
  "Proliferating" = "Immune cells",        #Proliferating Macrophages
  "CMS3"  = "CRC" ,
  "cDC" ="Immune cells",
  "Intermediate" = "Epithelial cells",
  "Tip-like ECs" = "Mesenchymal cells",
  "Proliferative ECs" = "Mesenchymal cells",
  "SPP1+" = "Immune cells", 
  "Lymphatic ECs" = "Mesenchymal cells",
  "Stromal 1" = "Mesenchymal cells",
  "Regulatory T cells" = "Immune cells", 
  "CMS1" = "CRC", 
  "Mast cells" = "Immune cells",
  "Pericytes" = "Mesenchymal cells",
  "IgG+ Plasma" = "Immune cells",
  "IgA+ Plasma" = "Immune cells",
  "Stromal 3" = "Mesenchymal cells",
  "CD8+ T cells" = "Immune cells",
  "CD4+ T cells" = "Immune cells",
  "NK cells" = "Immune cells",
  "T helper 17 cells" = "Immune cells",
  "T follicular helper cells" = "Immune cells",
  "gamma delta T cells" = "Immune cells",
  "CD19+CD20+ B" = "Immune cells",
  "Other" = "Low-quality cells"
)

# Map the existing cell types to the new merged names
obj@meta.data$merged_cell_type <- unname(
  merge_mapping[obj@meta.data$full_first_type]
)

# Verify the new merged cell types
unique(obj@meta.data$merged_cell_type)

# Set the merged cell types as the active identities
Idents(obj) <- obj@meta.data$merged_cell_type


obj[["RNA"]] <- NULL
obj[["RNA"]] <- obj[["Spatial"]]


obj@meta.data$nCount_RNA <- obj@meta.data$nCount_Spatial
obj@meta.data$nFeature_RNA <- obj@meta.data$nFeature_Spatial

obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)
saveRDS(obj, paste0(Output, Name, "_Step1/obj_Step1C.rds"))

obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold,CELLANNOTATION = TRUE)

STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)

#Calculating contamination scores
obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n,CELLANNOTATION = TRUE)

#PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)

#================Manuscript images=================================================
text_size <- 16
font_type <- "Arial"
line_width <- 1.2
bar_width <- 2.5

manuscript_colors <- c(
  '#E59CC4', '#3A6963', '#D6E7A3', '#FF6347', '#57C3F3', '#9FA3A8',
  '#E63863', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#8B4513',
  '#53A85F', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
  '#AA9A59', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#476D87',
  '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180',
  '#968175', '#E5D2DD', '#F1BB72', '#E95C59'
)

# Define a theme with all required elements
thm <- theme(
  text = element_text(size = text_size, family = font_type, face = "bold"),
  plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Remove all margins
  plot.title = element_text(size = text_size, family = font_type, face = "bold"),
  plot.subtitle = element_text(size = text_size, family = font_type), # Add subtitle element
  plot.tag = element_text(size = text_size, family = font_type, face = "bold") # Add tag element
)

# Create the output directory if it does not exist
CRC <- readRDS(paste0(output_dir, "CRC_projected_celltypelabelled.rds"))
scCLINICResult <- readRDS(paste0(output_dir, Name, "_Step2/Output_annotation_index/scCLINICResult.rds"))
M1 <- readRDS(paste0(output_dir, Name, "_Step2/annotation_index_recluster/annotation_index_cluster_M1.rds"))
M4 <- readRDS(paste0(output_dir, Name, "_Step2/annotation_index_recluster/annotation_index_cluster_M4.rds"))
M6 <- readRDS(paste0(output_dir, Name, "_Step2/annotation_index_recluster/annotation_index_cluster_M6.rds"))
M3 <- readRDS(paste0(output_dir, Name, "_Step2/annotation_index_recluster/annotation_index_cluster_M3.rds"))

#obj_Step1C <- readRDS("~/testscC/OutputHumanCRC/016/newUSERANNOTATION_LEVEL1mergedcelltypes_1B0252A15npeak7_ISfilter_Step1/obj_Step1C.rds")
obj_Step1C <- readRDS(paste0(output_dir, Name, "_Step1/obj_Step1C.rds"))

merge_mapping <- c(
  "Unknown" = "Low-quality cells",
  "CMS2" = "CRC",
  "Stromal 2" = "Mesenchymal cells",
  "Stem-like/TA" = "Epithelial cells",
  "Stem-like-TA" = "Epithelial cells", # In case it still exists
  "Enteric glial cells" = "Enteric glial cells",
  "Mature Enterocytes type 2" = "Epithelial cells",
  "Mature Enterocytes type 1" = "Epithelial cells",
  "Pro-inflammatory" = "Immune cells",     #Pro-inflammatory Macrophages
  "Myofibroblasts" = "Mesenchymal cells",
  "Smooth muscle cells" = "Mesenchymal cells",
  "Goblet cells" = "Epithelial cells", 
  "Stalk-like ECs" = "Mesenchymal cells",
  "Proliferating" = "Immune cells",        #Proliferating Macrophages
  "CMS3"  = "CRC" ,
  "cDC" ="Immune cells",
  "Intermediate" = "Epithelial cells",
  "Tip-like ECs" = "Mesenchymal cells",
  "Proliferative ECs" = "Mesenchymal cells",
  "SPP1+" = "Immune cells", 
  "Lymphatic ECs" = "Mesenchymal cells",
  "Stromal 1" = "Mesenchymal cells",
  "Regulatory T cells" = "Immune cells", 
  "CMS1" = "CRC", 
  "Mast cells" = "Immune cells",
  "Pericytes" = "Mesenchymal cells",
  "IgG+ Plasma" = "Immune cells",
  "IgA+ Plasma" = "Immune cells",
  "Stromal 3" = "Mesenchymal cells",
  "CD8+ T cells" = "Immune cells",
  "CD4+ T cells" = "Immune cells",
  "NK cells" = "Immune cells",
  "T helper 17 cells" = "Immune cells",
  "T follicular helper cells" = "Immune cells",
  "gamma delta T cells" = "Immune cells",
  "CD19+CD20+ B" = "Immune cells",
  "Other" = "Low-quality cells"
)

# Map the existing cell types to the new merged names
CRC@meta.data$merged_cell_type <- unname(
  merge_mapping[CRC@meta.data$full_first_type]
)

# Verify the new merged cell types
unique(CRC@meta.data$merged_cell_type)


# Add other relevant metadata to Breast
CRC$scCLINICScore <- scCLINICResult$scCLINICScore
CRC$scCLINIC_Level <- scCLINICResult$scCLINIC_Level
CRC$scCLINIC_ClusterID <- scCLINICResult$scCLINIC_ClusterID
CRC$annotation_index <- scCLINICResult$annotation_index
scCLINICResult$merged_cell_type <- CRC$merged_cell_type

DefaultAssay(scCLINICResult) <- "Spatial"

# Define new labels for your merged_cell_type
new_labels <- c("CRC" = "CRC (M1)", 
                "Enteric glial cells" = "Enteric glial cells (M2)",
                "Epithelial cells" = "Epithelial cells (M3)",
                "Immune cells" = "Immune cells (M4)",
                "Low-quality cells" = "Low-quality cells (M5)",
                "Mesenchymal cells" = "Mesenchymal cells (M6)")


# Add all old and new labels

# Update the factor levels in your Seurat object
scCLINICResult$merged_cell_type_newlabels <- factor(scCLINICResult$merged_cell_type, 
                                                    levels = names(new_labels), 
                                                    labels = new_labels)
# Update the factor levels in your Seurat object
obj_Step1C$merged_cell_type_newlabels <- factor(obj_Step1C$merged_cell_type, 
                                                levels = names(new_labels), 
                                                labels = new_labels)

#=================================Fig6A==========================================================
Idents(obj_Step1C) <- c("merged_cell_type_newlabels")
color_palette <- c("Low-quality cells (M5)" = "#9FA3A8", setNames(manuscript_colors, levels(obj_Step1C$merged_cell_type_newlabels)[levels(obj_Step1C$merged_cell_type_newlabels) != "Low-quality cells (M5)"]))
#names(color_palette) <- Idents(obj_Step1C) %>% levels()
Fig6A <- DimPlot(obj_Step1C, group.by = "merged_cell_type_newlabels", raster=FALSE,  sizes.highlight = 0.1,cols = color_palette, reduction = "umap", repel = TRUE, label = TRUE, label.size = 6, pt.size = 0.5)&labs(color = "scCLINIC")&
  theme( plot.title = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(), axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16),legend.position = c(0.02,0.13),legend.title = element_text(size = 18, family = font_type),legend.text = element_text(size = 16, family = font_type),legend.key.size = unit(0.5, 'lines'),legend.spacing.x = unit(0.2, "lines"),legend.spacing.y =unit(0, "lines"))

#ggsave(filename = paste0(output_dir, "Fig6A.png"), plot = Fig6A, width = 8, height = 8)

#=================================Fig6B==========================================================
CRC$merged_cell_type_newlabels <- scCLINICResult$merged_cell_type_newlabels
CRC$merged_cell_type_newlabels <- as.character(CRC$merged_cell_type_newlabels)  # Convert to character
CRC$merged_cell_type_newlabels[is.na(CRC$merged_cell_type_newlabels)] <- "Low-quality bins"  # Replace NA
CRC$merged_cell_type_newlabels <- as.factor(CRC$merged_cell_type_newlabels)  # Convert back to factor (if needed)

Idents(CRC) <- "merged_cell_type_newlabels"
names(manuscript_colors) <- Idents(CRC) %>% levels()

color_palette <- c("Low-quality bins" = "#9FA3A8", setNames(manuscript_colors, levels(CRC$merged_cell_type_newlabels)[levels(CRC$merged_cell_type_newlabels) != "Low-quality bins"]))

Fig6B <- SpatialDimPlot(CRC, group.by = "merged_cell_type_newlabels", 
                        cols = color_palette, image.alpha = 0,pt.size.factor = 5, 
                        label = FALSE) +
  theme(axis.text = element_text(size = text_size),
        legend.title         = element_blank(),
        legend.position      = c(0.45, 0.98),
        legend.justification = c("right", "top"), 
        legend.text = element_text(size = 18, family = font_type),
        legend.spacing = unit(0, "line"), 
        legend.key.height = unit(0, "line")
  ) +
  guides(
    fill = guide_legend(
      ncol  = 1,
      byrow = TRUE,
      override.aes = list(size = 4)        #to change legend point size
    )
  )


#ggsave(filename = paste0(output_dir, "Fig6B.png"), plot = Fig6B, width = 15, height = 10)

#=================================Fig6C==========================================================
Fig6C <-FeaturePlot(scCLINICResult, 
                    features = c("scCLINICScore"),
                    pt.size = 0.5,
                    reduction = "umap", 
                    raster = FALSE)+labs(color = "scCLINIC Score") +
  #ggtitle("ClusterID for Cells with Score > 0.1") +
  theme( plot.title = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),
         axis.text.x = element_blank(),axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
         legend.position = c(0.03,0.15),
         legend.title = element_text(size = 24, family = font_type),
         legend.text = element_text(size = 21, family = font_type),
         legend.key.size = unit(1.5, 'lines'),legend.spacing.x = unit(0.2, "lines"),
         legend.spacing.y =unit(0, "lines"))
#ggsave(filename = paste0(output_dir, "Fig6C.png"), plot = Fig6C, width = 8, height = 8)


#=================================Fig6D===Fig6E=======================================================
library(RColorBrewer)
area1 <- c("M1_S4", "M1_S7", "M4_S2",  "M6_S1")
area2 <- c("M1_S14", "M1_S19", "M3_S0" ,"M6_S17")
#area2 <- c("M1_S14", "M1_S19", "M6_S17", "M4_S4", "M3_S0")

scCLINICResult$level1_subclusters <- ifelse(
  scCLINICResult$scCLINIC_ClusterID %in% area1, scCLINICResult$scCLINIC_ClusterID,
  ifelse(
    scCLINICResult$scCLINIC_ClusterID %in% area2, scCLINICResult$scCLINIC_ClusterID,
    "Others"
  )
)

scCLINICResult$level1_subclusters <- factor(
  scCLINICResult$level1_subclusters,
  levels = c(area1, area2, "Others")
)

#red_cols  <- colorRampPalette(c("#fb6a4a", "#cb181d"))(length(area1))
#blue_cols <- colorRampPalette(c("#6baed6", "#08519c"))(length(area2))
red_cols_full  <- brewer.pal(5, "Reds")   # from very light to deep red
blue_cols_full <- brewer.pal(5, "Blues")  # from very light to deep blue
red_cols  <- red_cols_full[2:5]
blue_cols <- blue_cols_full[2:5]

custom_colors <- c(
  setNames(red_cols,  area1),
  setNames(blue_cols, area2),
  Others = "grey80"
)

Idents(scCLINICResult) <- "level1_subclusters"
levels_to_show <- c(area1, area2)

Fig6D <- DimPlot(
  scCLINICResult,
  reduction = "umap",
  group.by  = "level1_subclusters",
  cols      = custom_colors,
  pt.size = 0.5,
  label     = TRUE,
  label.size = 6,
  repel     = TRUE
) + 
  scale_color_manual(
    values       = custom_colors,       # full palette (including Others)
    breaks       = levels_to_show,      # only these levels in the legend
    labels       = levels_to_show      # same labels, in order
    
  )+
  theme(legend.position = c(0.02, 0.2), legend.text = element_text(size = 22) ,plot.title = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))


Fig6E <- SpatialDimPlot(scCLINICResult, 
                        group.by = "level1_subclusters", 
                        cols = custom_colors, 
                        image.alpha = 0, 
                        pt.size = 6) + 
  scale_fill_manual(
    values = custom_colors, 
    breaks = levels_to_show,
    labels = levels_to_show,
    drop   = TRUE
  ) +
  theme(
    legend.title         = element_blank(),
    legend.position      = c(0.4, 0.98),
    legend.justification = c("right", "top"), 
    legend.text = element_text(size = 18, family = font_type),
    legend.spacing = unit(0, "line"), 
    legend.key.height = unit(0, "line"),
    plot.margin = unit(c(0.5, 0, 0, 0), "lines")
  ) +
  guides(
    fill = guide_legend(
      ncol  = 2,
      byrow = FALSE,
      override.aes = list(size = 5)        #to change legend point size
      
    )
  )

#================================Fig6F========================================================

# Define new labels for your merged_cell_type
new_labels <- c(
  "Unknown" = "Unknown",
  "CMS2" = "CMS2",
  "Stromal 2" = "Stromal 2",
  "Stem-like/TA" = "Stem-like/TA",
  "Stem-like-TA" = "Stem-like-TA", # In case it still exists
  "Enteric glial cells" = "Enteric glial cells",
  "Mature Enterocytes type 2" = "Mature Enterocytes type 2",
  "Mature Enterocytes type 1" = "Mature Enterocytes type 1",
  "Pro-inflammatory" = "Pro-inflammatory Macrophages",     #Pro-inflammatory Macrophages
  "Myofibroblasts" = "Myofibroblasts",
  "Smooth muscle cells" = "Smooth muscle cells",
  "Goblet cells" = "Goblet cells", 
  "Stalk-like ECs" = "Stalk-like ECs",
  "Proliferating" = "Proliferating Macrophages",        #Proliferating Macrophages
  "CMS3"  = "CMS3" ,
  "cDC" ="cDC",
  "Intermediate" = "Intermediate",
  "Tip-like ECs" = "Tip-like ECs",
  "Proliferative ECs" = "Proliferative ECs",
  "SPP1+" = "SPP1+", 
  "Lymphatic ECs" = "Lymphatic ECs",
  "Stromal 1" = "Stromal 1",
  "Regulatory T cells" = "Regulatory T cells", 
  "CMS1" = "CMS1", 
  "Mast cells" = "Mast cells",
  "Pericytes" = "Pericytes",
  "IgG+ Plasma" = "IgG+ Plasma",
  "IgA+ Plasma" = "IgA+ Plasma",
  "Stromal 3" = "Stromal 3",
  "CD8+ T cells" = "CD8+ T cells",
  "CD4+ T cells" = "CD4+ T cells",
  "NK cells" = "NK cells",
  "T helper 17 cells" = "T helper 17 cells",
  "T follicular helper cells" = "T follicular helper cells",
  "gamma delta T cells" = "gamma delta T cells",
  "CD19+CD20+ B" = "CD19+CD20+ B",
  "Other" = "Other"
)

# Add all old and new labels

# Update the factor levels in your Seurat object
CRC$full_first_type_new <- factor(CRC$full_first_type, 
                                  levels = names(new_labels), 
                                  labels = new_labels)

my45colors <- c('#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', 
                '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
                '#F1BB72', '#4B0082', '#2E8B57', '#FFFF00', '#91D0BE', 
                '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', 
                '#625D9E', '#68A180', '#F3B1A0', '#968175', '#D6E7A3', 
                '#57C3F3','#476D87', '#E95C59', '#E59CC4', '#AB3282', 
                '#23452F', '#BD956A', '#8C549C', '#585658', '#E5D2DD', 
                '#000000', '#00FF00', '#00FFFF', '#FF0000', 
                '#FFD700', '#5F9EA0','#C1E6F3', '#6778AE', '#3A6963')

# Filter for major cluster 2 and subcluster M2S2
major_cluster4 <- CRC@meta.data %>% filter(annotation_index == "M4") %>% mutate(Group = "Major Cluster 4")
subcluster_M4S2 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M4_S2") %>% mutate(Group = "Subcluster M4S2")

major_cluster4_counts <- major_cluster4 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M4")

subcluster_M4S2_counts <- subcluster_M4S2 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M4_S2")

# Combine data
plot_data1 <- bind_rows(major_cluster4_counts, subcluster_M4S2_counts)

names(my45colors) <- NULL
# Stacked bar plot
A1 <- ggplot(plot_data1, aes(x = Group, y = count, fill = full_first_type_new)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +  # Use "fill" for proportion
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0, 1, by = 1)) +  # Convert to percentage
  scale_fill_manual(values = my45colors) +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  #  theme_minimal() +
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), axis.text.y = element_text(size = 25),axis.text = (element_text(size = 18)),axis.title.x = element_text(size = 18), 
        axis.title.y = element_blank(),legend.key.size = unit(15, "pt"), legend.text = element_text(size = 20), legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", size = 0.3),  # Minor gridlines
        panel.background = element_rect(fill = "white"))+  # White background
  guides(fill = guide_legend(ncol=2)) +
  scale_x_discrete(expand = c(0.3, 0.3))


major_cluster6 <- CRC@meta.data %>% filter(annotation_index == "M6") %>% mutate(Group = "Major Cluster 6")
subcluster_M6S1 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M6_S1") %>% mutate(Group = "Subcluster M6S1")
subcluster_M6S17 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M6_S17") %>% mutate(Group = "Subcluster M6S17")

major_cluster6_counts <- major_cluster6 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M6")

subcluster_M6S1_counts <- subcluster_M6S1 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M6_S1")
subcluster_M6S17_counts <- subcluster_M6S17 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M6_S17")

# Combine data
plot_data2 <- bind_rows(major_cluster6_counts, subcluster_M6S1_counts, subcluster_M6S17_counts)

my45colors <- c('#D6E7A3', '#57C3F3','#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', 
                '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', 
                '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', 
                '#AA9A59', '#E63863', '#E39A35', '#F1BB72', '#C1E6F3', '#6778AE', '#3A6963', 
                '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', 
                '#CCC9E6', '#625D9E', '#68A180', '#F3B1A0', '#968175', 
                '#E5D2DD', '#000000', '#00FF00', '#00FFFF', '#FF0000', 
                '#FFD700', '#5F9EA0', '#4B0082', '#2E8B57', '#FFFF00')

names(my45colors) <- NULL
# Stacked bar plot
A2 <- ggplot(plot_data2, aes(x = Group, y = count, fill = full_first_type_new)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +  # Use "fill" for proportion
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0, 1, by = 1)) +  # Convert to percentage
  scale_fill_manual(values = my45colors) +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  #  theme_minimal() +
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), axis.text.y = element_text(size = 25),axis.text = (element_text(size = 18)),axis.title.x = element_text(size = 18), 
        axis.title.y = element_blank(),legend.key.size = unit(15, "pt"), legend.text = element_text(size = 20), legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", size = 0.3),  # Minor gridlines
        panel.background = element_rect(fill = "white"))+  # White background
  guides(fill = guide_legend(ncol=2)) +
  scale_x_discrete(expand = c(0.3, 0.3))

major_cluster1 <- CRC@meta.data %>% filter(annotation_index == "M1") %>% mutate(Group = "Major Cluster 1")
subcluster_M1S7 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M1_S7") %>% mutate(Group = "Subcluster M1S7")
subcluster_M1S4 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M1_S4") %>% mutate(Group = "Subcluster M1S4")
subcluster_M1S14 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M1_S14") %>% mutate(Group = "Subcluster M1S14")
subcluster_M1S19 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M1_S19") %>% mutate(Group = "Subcluster M1S19")

major_cluster1_counts <- major_cluster1 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M1")

subcluster_M1S7_counts <- subcluster_M1S7 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M1_S7")
subcluster_M1S4_counts <- subcluster_M1S4 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M1_S4")
subcluster_M1S14_counts <- subcluster_M1S14 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M1_S14")
subcluster_M1S19_counts <- subcluster_M1S19 %>%
  group_by(full_first_type_new) %>%
  summarise(count = n()) %>%
  mutate(Group = "M1_S19")

# Combine data
plot_data3 <- bind_rows(major_cluster1_counts, subcluster_M1S7_counts, subcluster_M1S4_counts, subcluster_M1S14_counts, subcluster_M1S19_counts)
plot_data3$Group <- factor(plot_data3$Group, levels = c("M1", "M1_S4", "M1_S7", "M1_S14", "M1_S19"))
my45colors <- c('#C1E6F3', '#6778AE', '#3A6963', '#D6E7A3', '#57C3F3', 
                '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', 
                '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', 
                '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', 
                '#AA9A59', '#E63863', '#E39A35', '#F1BB72', '', 
                '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', 
                '#CCC9E6', '#625D9E', '#68A180', '#F3B1A0', '#968175', 
                '#E5D2DD', '#000000', '#00FF00', '#00FFFF', '#FF0000', 
                '#FFD700', '#5F9EA0', '#4B0082', '#2E8B57', '#FFFF00')
names(my45colors) <- NULL
# Stacked bar plot
A3 <- ggplot(plot_data3, aes(x = Group, y = count, fill = full_first_type_new)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +  # Use "fill" for proportion
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0, 1, by = 1)) +  # Convert to percentage
  scale_fill_manual(values = my45colors) +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  #  theme_minimal() +
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), axis.text.y = element_text(size = 25), axis.text = (element_text(size = 18)),axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 25),legend.key.size = unit(15, "pt"), legend.text = element_text(size = 20), legend.title = element_text(size = 20),
        legend.position = "bottom", legend.justification = c(0, 1),
        panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", size = 0.3),  # Minor gridlines
        panel.background = element_rect(fill = "white"))+  # White background
  guides(fill = guide_legend(ncol=1)) +
  scale_x_discrete(expand = c(0.3, 0.3))

Fig6F <- A3|A1|A2

#ggsave(filename = paste0(output_dir, "Fig6F.png"), plot = Fig6F, width = 28, height = 15)

#============================Fullplot==================================================
plot_A <- wrap_elements(
  wrap_elements(
    (Fig6A + Fig6B + Fig6C + Fig6D + plot_layout(ncol = 4, guides = "keep")) /
      (Fig6E + Fig6F + plot_layout(ncol = 2)) + 
      plot_layout(ncol = 1, heights = c(0.9, 0.8, 0.4, 1), widths = c(1)) &  
      theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
    
  ) +  
    plot_annotation(
      title = "Colorectal Cancer Spatial Dataset (16 um bins)",
      theme = thm + theme(
        plot.title = element_text(
          face = "bold",
          size = 20,
          family = font_type,
          hjust = 0.5
        )
      )
    )
)  


#saveplot for plot_A
ggsave(filename = paste0(output_dir, "Figure6.png"), plot = plot_A, width = 30, height = 30)



#===================================FigS6A===================================================
M1$scCLINIC_ClusterID <- scCLINICResult$scCLINIC_ClusterID

Idents(M1) <- c("scCLINIC_ClusterID")
names(manuscript_colors) <- Idents(M1) %>% levels() 
FigS6A <- DimPlot(M1, group.by = "scCLINIC_ClusterID",raster=FALSE, cols = manuscript_colors, pt.size = 0.5)&labs(color = "M1 subcluster")&
  theme( plot.title = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),legend.position = c(0.5,0.2),legend.title = element_text(size = 20, family = font_type),
         axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
         legend.text = element_text(size = 18, family = font_type),legend.key.size = unit(0.3, 'lines'),legend.spacing.x = unit(0.1, "lines"),legend.spacing.y =unit(0, "lines")
  )&
  guides(color = guide_legend(ncol=3, override.aes = list(size=3)) )

#ggsave(filename = "~/scCLINIC_Manuscript_figure6/Finalized23June2025/FigS6A.png", plot = FigS6A, width = 8, height = 8)

#===========================FigS6B=================================================================
#PLEASE SAVE A COPY OF ARTIFACTINFO.CSV FILE AND NAME IT "ArtifactsInfo_Bef.csv" BEFORE PROCEEDING
#Load artifacts information and scCLINIC result rds
contamgeneinfo <- read.csv(paste0(output_dir, Name, "_Step2/Output_annotation_index/ArtifactsInfo_Bef.csv"),	###SCCLINIC
                           na.strings = c("", "NA")
)

#scCLINIC Subcluster ID
contamgeneinfo$MajorSub <- paste0(contamgeneinfo$Major_Cluster,"_",contamgeneinfo$Subcluster)
#Summarize ES score and their source of major cluster for each subclusters
result <- contamgeneinfo %>%
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

# Convert to dataframe for plotting
CCS_Matrix <- as.data.frame(CCS_Matrix)
CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
data <- CCS_Matrix %>%
  pivot_longer(cols = -MajorSub, names_to = "Component", values_to = "value") #dependency tidyr

# Filter out other clusters, only keep cluster X which wish to display
clusterx <- "M1"
filtered_data <- data %>%
  dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

# Heatmap
# Convert dataframe
heatmap_data <- dcast(filtered_data, MajorSub ~ Component, value.var = "value")
heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

rownames(heatmap_data) <- heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL #remove MajorSub columns

# Summarize scCLINIC Score (CS) of each subclusters
CS_table <- scCLINICResult@meta.data %>%
  group_by(scCLINIC_ClusterID) %>%
  summarize(
    scCLINICScore = mean(scCLINICScore), #scCLINICScore of each cells within each subclusters are same value...
  ) %>%
  ungroup()

#Convert to plotting dataframe
heatmap_data <- as.data.frame(heatmap_data)  # Convert heatmap_data to a dataframe
heatmap_data <- heatmap_data %>%
  mutate(CS = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <- data.frame((heatmap_matrix))
colnames(multicolplot) <- dplyr::recode(
  colnames(multicolplot),
  "M1" = "CRC_M1",
  "M2" = "Entericglialcells_M2",
  "M3" = "Epithelial_M3",
  "M4" = "Immunecells_M4", 
  "M6" = "Mesenchymal_M6"
)

rownames_sorted <- rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
multicolplot <- multicolplot[rownames_sorted, ]

colnames_sorted <- colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
multicolplot <- multicolplot[, colnames_sorted]

multicolplot$Category <- as.character(rownames(multicolplot))
multicolplot$Category <- factor(multicolplot$Category, levels = rownames_sorted)

x_limits <- c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <- seq(x_limits[1], x_limits[2], length.out = 2)
text_size <- 16

# Create individual plots without y-axis text
plot_list <- list()
for (i in seq_along(colnames_sorted)) {
  clus <- colnames_sorted[i]
  p <- ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(stat = "identity", fill = c('#3A6963', '#D6E7A3', '#FF6347','#57C3F3', "#FFD300")[i]) +
    labs(title = paste(clus), y = "", x = "") +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "grey", size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 24),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(0,0.5,0,0), "cm")
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
y_axis_plot <- ggplot(multicolplot, aes_string(x = "Category", y = 1)) +
  geom_blank() +
  #geom_bar(stat = "identity") +
  labs(title = "", y = "", x = "Subcluster") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    #axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    plot.title = element_text(size = 24),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 22),
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  coord_flip() +
  scale_y_continuous(breaks = custom_breaks, limits = x_limits)


FigS6B <- y_axis_plot + 
  patchwork::wrap_plots(plot_list,nrow=1) + 
  patchwork::plot_annotation(title = "Artifact Source",theme = theme(plot.title = element_text(size=24,hjust = 0.5), axis.text = element_text(size = 20), axis.title = element_text(size = 22))) + 
  patchwork::plot_layout(nrow=1,widths = c(1,10000))

#ggsave(filename = "~/scCLINIC_Manuscript_figure6/Finalized23June2025/FigS6B.png", plot = FigS6B, width = 21, height = 10)


#===========================FigS6C==============================================================
M4$scCLINIC_ClusterID <- scCLINICResult$scCLINIC_ClusterID

Idents(M4) <- c("scCLINIC_ClusterID")
names(manuscript_colors) <- Idents(M4) %>% levels() 
FigS6C <- DimPlot(M4, group.by = "scCLINIC_ClusterID",raster=FALSE, cols = manuscript_colors, pt.size = 1)&labs(color = "M4 subcluster")&
  theme( plot.title = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),
         axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
         legend.position = c(0.02,0.15),legend.title = element_text(size = 20, family = font_type),legend.text = element_text(size = 20, family = font_type),legend.key.size = unit(0.3, 'lines'),legend.spacing.x = unit(0.1, "lines"),legend.spacing.y =unit(0, "lines")
  )&
  guides(color = guide_legend(ncol=3, override.aes = list(size=3)) )

#===========================FigS6D=================================================================
#Load artifacts information and scCLINIC result rds
contamgeneinfo <- read.csv(paste0(output_dir, Name, "_Step2/Output_annotation_index/ArtifactsInfo_Bef.csv"),	###SCCLINIC
                           na.strings = c("", "NA")
)

#scCLINIC Subcluster ID
contamgeneinfo$MajorSub <- paste0(contamgeneinfo$Major_Cluster,"_",contamgeneinfo$Subcluster)
#Summarize ES score and their source of major cluster for each subclusters
result <- contamgeneinfo %>%
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

# Convert to dataframe for plotting
CCS_Matrix <- as.data.frame(CCS_Matrix)
CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
data <- CCS_Matrix %>%
  pivot_longer(cols = -MajorSub, names_to = "Component", values_to = "value") #dependency tidyr

#for (clusterx in unique(obj@meta.data[,OverlapRatio])){
# Filter out other clusters, only keep cluster X which wish to display
clusterx <- "M4"
filtered_data <- data %>%
  dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

# Heatmap
# Convert dataframe
heatmap_data <- dcast(filtered_data, MajorSub ~ Component, value.var = "value")
heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

rownames(heatmap_data) <- heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL #remove MajorSub columns

# Summarize scCLINIC Score (CS) of each subclusters
CS_table <- scCLINICResult@meta.data %>%
  group_by(scCLINIC_ClusterID) %>%
  summarize(
    scCLINICScore = mean(scCLINICScore), #scCLINICScore of each cells within each subclusters are same value...
  ) %>%
  ungroup()

#Convert to plotting dataframe
heatmap_data <- as.data.frame(heatmap_data)  # Convert heatmap_data to a dataframe
heatmap_data <- heatmap_data %>%
  mutate(CS = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <- data.frame((heatmap_matrix))
colnames(multicolplot) <- dplyr::recode(
  colnames(multicolplot),
  "M1" = "CRC_M1",
  "M2" = "Entericglialcells_M2",
  "M3" = "Epithelial_M3",
  "M4" = "Immunecells_M4", 
  "M6" = "Mesenchymal_M6"
)

rownames_sorted <- rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
multicolplot <- multicolplot[rownames_sorted, ]

colnames_sorted <- colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
multicolplot <- multicolplot[, colnames_sorted]

multicolplot$Category <- as.character(rownames(multicolplot))
multicolplot$Category <- factor(multicolplot$Category, levels = rownames_sorted)

x_limits <- c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <- seq(x_limits[1], x_limits[2], length.out = 2)
text_size <- 20

# Create individual plots without y-axis text
plot_list <- list()
for (i in seq_along(colnames_sorted)) {
  clus <- colnames_sorted[i]
  p <- ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(stat = "identity", fill = c('#E59CC4', '#3A6963', '#D6E7A3', '#57C3F3', "#FFD300")[i]) +
    labs(title = paste(clus), y = "", x = "") +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "grey", size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 24),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(0,0.5,0,0), "cm")
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
y_axis_plot <- ggplot(multicolplot, aes_string(x = "Category", y = 1)) +
  geom_blank() +
  #geom_bar(stat = "identity") +
  labs(title = "", y = "", x = "Subcluster") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    plot.title = element_text(size = 25),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  coord_flip() +
  scale_y_continuous(breaks = custom_breaks, limits = x_limits)


FigS6D <- y_axis_plot + 
  patchwork::wrap_plots(plot_list,nrow=1) + 
  patchwork::plot_annotation(title = "Artifact Source",theme = theme(plot.title = element_text(size=24,hjust = 0.5),
                                                                     axis.title = element_text(size = 25))) + 
  patchwork::plot_layout(nrow=1,widths = c(1,10000))

#===========================FigS6E==============================================================
M6$scCLINIC_ClusterID <- scCLINICResult$scCLINIC_ClusterID

Idents(M6) <- c("scCLINIC_ClusterID")
names(manuscript_colors) <- Idents(M6) %>% levels() 
FigS6E <- DimPlot(M6, group.by = "scCLINIC_ClusterID",raster=FALSE, cols = manuscript_colors, pt.size = 0.5)&labs(color = "M6 subcluster")&
  theme(plot.title = element_blank(),axis.title = element_text(size = 18) ,axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.position = c(0.01,0.15),legend.title = element_text(size = 16, family = font_type),legend.text = element_text(size = 16, family = font_type),legend.key.size = unit(0.3, 'lines'),legend.spacing.x = unit(0.1, "lines"),legend.spacing.y =unit(0, "lines")
  )&
  guides(color = guide_legend(ncol=3, override.aes = list(size=3)) )

#===========================FigS6F=================================================================
#Load artifacts information and scCLINIC result rds
contamgeneinfo <- read.csv(paste0(output_dir, Name, "_Step2/Output_annotation_index/ArtifactsInfo_Bef.csv"),	###SCCLINIC
                           na.strings = c("", "NA")
)

#scCLINIC Subcluster ID
contamgeneinfo$MajorSub <- paste0(contamgeneinfo$Major_Cluster,"_",contamgeneinfo$Subcluster)
#Summarize ES score and their source of major cluster for each subclusters
result <- contamgeneinfo %>%
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

# Convert to dataframe for plotting
CCS_Matrix <- as.data.frame(CCS_Matrix)
CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
data <- CCS_Matrix %>%
  pivot_longer(cols = -MajorSub, names_to = "Component", values_to = "value") #dependency tidyr

# Filter out other clusters, only keep cluster X which wish to display
clusterx <- "M6"
filtered_data <- data %>%
  dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

# Heatmap
# Convert dataframe
heatmap_data <- dcast(filtered_data, MajorSub ~ Component, value.var = "value")
heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

rownames(heatmap_data) <- heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
heatmap_data$MajorSub <- NULL #remove MajorSub columns

# Summarize scCLINIC Score (CS) of each subclusters
CS_table <- scCLINICResult@meta.data %>%
  group_by(scCLINIC_ClusterID) %>%
  summarize(
    scCLINICScore = mean(scCLINICScore), #scCLINICScore of each cells within each subclusters are same value...
  ) %>%
  ungroup()

#Convert to plotting dataframe
heatmap_data <- as.data.frame(heatmap_data)  # Convert heatmap_data to a dataframe
heatmap_data <- heatmap_data %>%
  mutate(CS = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score

# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Reorder and prepare the data
multicolplot <- data.frame((heatmap_matrix))
colnames(multicolplot) <- dplyr::recode(
  colnames(multicolplot),
  "M1" = "CRC_M1",
  "M2" = "Entericglialcells_M2",
  "M3" = "Epithelial_M3",
  "M4" = "Immunecells_M4", 
  "M6" = "Mesenchymal_M6"
)


rownames_sorted <- rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
multicolplot <- multicolplot[rownames_sorted, ]

colnames_sorted <- colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
multicolplot <- multicolplot[, colnames_sorted]

multicolplot$Category <- as.character(rownames(multicolplot))
multicolplot$Category <- factor(multicolplot$Category, levels = rownames_sorted)

x_limits <- c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
custom_breaks <- seq(x_limits[1], x_limits[2], length.out = 2)
text_size <- 16

# Create individual plots without y-axis text
plot_list <- list()
for (i in seq_along(colnames_sorted)) {
  clus <- colnames_sorted[i]
  p <- ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
    geom_bar(stat = "identity", fill = c('#E59CC4', '#3A6963', '#D6E7A3', '#FF6347', "#FFD300")[i]) +
    labs(title = paste(clus), y = "", x = "") +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "grey", size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 24),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(0,0.5,0,0), "cm")
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
y_axis_plot <- ggplot(multicolplot, aes_string(x = "Category", y = 1)) +
  geom_blank() +
  #geom_bar(stat = "identity") +
  labs(title = "", y = "", x = "Subcluster") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    plot.title = element_text(size = 24),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  coord_flip() +
  scale_y_continuous(breaks = custom_breaks, limits = x_limits)



FigS6F <- y_axis_plot + 
  patchwork::wrap_plots(plot_list,nrow=1) + 
  patchwork::plot_annotation(title = "Artifact Source",theme = theme(plot.title = element_text(size=24,hjust = 0.5))) + 
  patchwork::plot_layout(nrow=1,widths = c(1,10000))


#==============================FigS6G=====================================================
M4$scCLINIC_ClusterID <- scCLINICResult$scCLINIC_ClusterID
M1$scCLINIC_ClusterID <- scCLINICResult$scCLINIC_ClusterID
M6$scCLINIC_ClusterID <- scCLINICResult$scCLINIC_ClusterID

M4$Contam_clus <- ifelse(scCLINICResult$scCLINIC_ClusterID == "M4_S2", "M4_S2", "Others")

Idents(M4) <- "Contam_clus"
custom_colors <- c("M4_S2" = "#E95C59")

levels(M4$Contam_clus)[levels(M4$Contam_clus) == "Others"] <- NA 
T2 <- DimPlot(M4, group.by = "Contam_clus",raster=FALSE, cols = custom_colors)&labs(color = "scCLINIC")&
  theme( plot.title = element_blank(), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),
         legend.position = c(0.2,0.2),legend.title = element_text(size = 20, family = font_type),legend.text = element_text(size = 18, family = font_type),legend.key.size = unit(0.5, 'lines'),legend.spacing.x = unit(0.2, "lines"),legend.spacing.y =unit(0, "lines"))


M1$Contam_clus <- ifelse(scCLINICResult$scCLINIC_ClusterID == "M1_S4", "M1_S4",
                         ifelse(scCLINICResult$scCLINIC_ClusterID == "M1_S7", "M1_S7", 
                                ifelse(scCLINICResult$scCLINIC_ClusterID == "M1_S14", "M1_S14",
                                       ifelse(scCLINICResult$scCLINIC_ClusterID == "M1_S19", "M1_S19", "Others"))))
Idents(M1) <- "Contam_clus"
custom_colors <- c("M1_S4" = "#E95C59", "M1_S7" = "#8C549C", "M1_S14" = "#1e88e5", "M1_S19" = "#71b84e")

levels(M1$Contam_clus)[levels(M1$Contam_clus) == "Others"] <- NA 
T1 <- DimPlot(M1, group.by = "Contam_clus",raster=FALSE, cols = custom_colors)&labs(color = "scCLINIC")&
  theme( plot.title = element_blank(),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),legend.position = c(0.8,0.2),legend.title = element_text(size = 20, family = font_type),legend.text = element_text(size = 18, family = font_type),legend.key.size = unit(0.5, 'lines'),legend.spacing.x = unit(0.2, "lines"),legend.spacing.y =unit(0, "lines"))


M6$Contam_clus <- ifelse(scCLINICResult$scCLINIC_ClusterID == "M6_S1", "M6_S1",
                         ifelse(scCLINICResult$scCLINIC_ClusterID == "M6_S17", "M6_S17", "Others"))

Idents(M6) <- "Contam_clus"
custom_colors <- c("M6_S1" = "#E95C59", "M6_S17" = "#8C549C")

levels(M6$Contam_clus)[levels(M6$Contam_clus) == "Others"] <- NA 
T3 <- DimPlot(M6, group.by = "Contam_clus",raster=FALSE, cols = custom_colors)&labs(color = "scCLINIC")&
  theme( plot.title = element_blank(),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),legend.position = c(0.8,0.2),legend.title = element_text(size = 20, family = font_type),legend.text = element_text(size = 18, family = font_type),legend.key.size = unit(0.5, 'lines'),legend.spacing.x = unit(0.2, "lines"),legend.spacing.y =unit(0, "lines"))


FigS6G <- T1/T2/T3

#==============================FigS6H=====================================================
R1 <- FeaturePlot(M4,
                  reduction = "umap",
                  features = c("CEACAM5" ,"CD68","RGS5", "SPARC", "AGR2"), ncol = 5, 
                  pt.size = 1.5) &     #CEACAM5 WAS USED PREVIOUSLY BUT GOBLET CLUSTER ALSO EXPRESS
  labs(title = NULL) &
  scale_color_gradientn(                                                 #C1QC IMMUNE GENE IS REPLACED WITH LYZ
    colors = c("grey", '#E63863'),
    limits = c(0, 4),
    oob = scales::squish
  )& 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
R1 <- (wrap_elements(grid::textGrob("M4", rot = 90, gp = gpar(fontsize = 19.3)))|R1) +  plot_layout(widths = c(0.01, 0.99))


R2 <- FeaturePlot(M1,
                  reduction = "umap",
                  features = c("CEACAM5" ,"CD68","RGS5", "SPARC", "AGR2"), ncol = 5, 
                  pt.size = 1.5) &     #CEACAM5 WAS USED PREVIOUSLY BUT GOBLET CLUSTER ALSO EXPRESS
  #labs(title = NULL) &
  scale_color_gradientn(                                                 #C1QC IMMUNE GENE IS REPLACED WITH LYZ
    colors = c("grey", '#E63863'),
    limits = c(0, 4),
    oob = scales::squish
  )& 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))

R2 <- (wrap_elements(grid::textGrob("M1", rot = 90, gp = gpar(fontsize = 19.3)))|R2) +  plot_layout(widths = c(0.01, 0.99))

R3 <- FeaturePlot(M6,
                  reduction = "umap",
                  features = c("CEACAM5" ,"CD68","RGS5", "SPARC", "AGR2"), ncol = 5, 
                  pt.size = 1.5) &     #CEACAM5 WAS USED PREVIOUSLY BUT GOBLET CLUSTER ALSO EXPRESS
  labs(title = NULL) &
  scale_color_gradientn(                                                 #C1QC IMMUNE GENE IS REPLACED WITH LYZ
    colors = c("grey", '#E63863'),
    limits = c(0, 4),
    oob = scales::squish
  )& 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
R3 <- (wrap_elements(grid::textGrob("M6", rot = 90, gp = gpar(fontsize = 19.3)))|R3) +  plot_layout(widths = c(0.01, 0.99))


FigS6H <- R2/R1/R3

FigS6H <- wrap_plots(FigS6H, guides = "collect")& 
  theme(legend.position = "right")

#===============================FigS6I======================================================
major_cluster3 <- CRC@meta.data %>% filter(annotation_index == "M3") %>% mutate(Group = "Major Cluster 3")
subcluster_M3S0 <- CRC@meta.data %>% filter(scCLINIC_ClusterID == "M3_S0") %>% mutate(Group = "Subcluster M3_S0")

major_cluster3_counts <- major_cluster3 %>%
  group_by(full_first_type) %>%
  summarise(count = n()) %>%
  mutate(Group = "M3")

subcluster_M3S0_counts <- subcluster_M3S0 %>%
  group_by(full_first_type) %>%
  summarise(count = n()) %>%
  mutate(Group = "M3_S0")


# Combine data
plot_data <- bind_rows(major_cluster3_counts, subcluster_M3S0_counts)
#plot_data <- bind_rows(major_cluster2_counts, subcluster_M2S2_counts)
#plot_data$full_first_type <- factor(plot_data$full_first_type, levels = unique(c(major_cluster2$full_first_type, subcluster_M2S2$full_first_type)))

names(manuscript_colors) <- NULL
# Stacked bar plot
FigS6I <- ggplot(plot_data, aes(x = Group, y = count, fill = full_first_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +  # Use "fill" for proportion
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0, 1, by = 1)) +  # Convert to percentage
  scale_fill_manual(values = manuscript_colors) +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  #  theme_minimal() +
  theme(axis.text.x = element_text(size = 24, angle = 45, hjust = 1), axis.text.y = element_text(size = 24),axis.text = (element_text(size = 24)),axis.title.x = element_text(size = 24), 
        axis.title.y = element_text(size = 25),legend.key.size = unit(15, "pt"), legend.text = element_text(size=22),
        legend.title = element_text(size = 24) ,legend.position = "bottom",legend.title.position = "top",
        panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", size = 0.3),  # Minor gridlines
        panel.background = element_rect(fill = "white"))+  # White background
  guides(fill = guide_legend(ncol=1)) +
  scale_x_discrete(expand = c(0.3, 0.3))
#ggsave(filename = "~/scCLINIC_Manuscript_figure6/Finalized23June2025/FigS6I.png", plot = FigS6I, width = 8, height = 15)


#============================Fullplot==================================================
plot_B <- wrap_elements(
    (FigS6A + wrap_elements(FigS6B)  + plot_layout(ncol = 2, widths = c(0.2,0.7),guides = "keep")) /
      (FigS6C + wrap_elements(FigS6D) + plot_layout(ncol = 2, widths = c(0.2,0.7),guides = "keep"))/
      (FigS6E + wrap_elements(FigS6F) + plot_layout(ncol = 2, widths = c(0.2,0.7),guides = "keep"))/ 
      (
        wrap_elements(FigS6G) +
          wrap_elements(FigS6H) +
          wrap_elements(FigS6I) +
          plot_layout(ncol = 3, widths = c(0.2, 0.6, 0.2), guides = "collect")
      )+ 
      plot_layout(ncol = 1, heights = c(0.2, 0.2, 0.2, 0.4), widths = c(1)) +  
      theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
    
  ) 


#saveplot for plot_A
ggsave(filename = paste0(output_dir, "ExtendedDataFig5.png"), plot = plot_B, width = 40, height = 40)


#Run this only after saving the ArtifactInfo.csv
#PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)
