# --------------------------------------------------- #
# Analysis of scRNA sequencing data as described in Tahri et al. 2023
# Script by Sabrin Tahri, optimized for speed by Mathijs Sanders and Gregory van Beek 
# Dept. of Hematology, Erasmus MC Cancer Institute, Rotterdam, the Netherlands
# --------------------------------------------------- #

# --------------------------------------------------- #
# READ-ME (IMPORTANT)
# This script is meant to generate the NK cell object by in silico selection of NK cells from the CD38+ immune dataset. 
# The CD38+ immune data from individual patients is first integrated per timepoint and per source (controls, newly diagnosed MM, C4D28, D100, W25, W52, W105 and relapse).
# The timepoints and sources are then integrated into one major CD38+ immune dataset. 
# NK cells are then subsetted from the major CD38+ immune object by selection on KLRD1, KLRF1, NKG7 and GNLY transcription (Yang et al.)
# Finally NK cells are reprocessed in to the final major NK cell object. 
# The script for dataset integration per timepoint is separately available on github, named 'integration_[timepoint].R'.
# --------------------------------------------------- #

# --------------------------------------------------- #
# Loading libraries
# --------------------------------------------------- #

print("Loading packages ...")
packages <- c("Seurat",
              "cowplot",
              "harmony",
              "readxl",
              "ggplot2",
              "readxl",
              "dittoSeq",
              "reticulate",
              "parallel",
              "future")

invisible(suppressMessages(suppressWarnings(lapply(packages, library, character.only = T))))

# --------------------------------------------------- #
# Function for defining number of cores 
# --------------------------------------------------- #

cores <- 10 
plan("multisession", workers = cores)
options(future.fork.enable = TRUE)
options(future.globals.maxSize = 80000 * 1024^2)
options(bitmapType = "cairo")

# --------------------------------------------------- #
# Function for speeding up the loading and saving of Seurat objects
# --------------------------------------------------- #

mcsaveRDS <- function(object,file,threads=parallel::detectCores()) {
  ### https://stackoverflow.com/questions/28927750/what-is-the-method-to-save-objects-in-a-compressed-form-using-multiple-threads-c
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

mcreadRDS <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}


# --------------------------------------------------- #
# Load object per timepoint
# --------------------------------------------------- #

print("Loading RDS files ...")
baseline_control <- mcreadRDS("/data/output/baseline_control.rds", threads = 4)
baseline_myeloma <- mcreadRDS("/data/output/baseline_myeloma.rds", threads = 4)
C4D28 <- mcreadRDS("/data/output/C4D28.rds", threads = 4)
D100 <- mcreadRDS("/data/output/D100.rds", threads = 4)
W25 <- mcreadRDS("/data/output/W25.rds", threads = 4)
W52 <- mcreadRDS("/data/output/W52.rds", threads = 4)
W105 <- mcreadRDS("/data/output/W105.rds", threads = 4)
relapse <- mcreadRDS("/data/output/relapse.rds", threads = 4)

# --------------------------------------------------- #
# Making and saving MM CD38+ immune object
# --------------------------------------------------- #

# Identification of integration anchors
print("Integrating first batch ...")
reference.list <- c(baseline_myeloma, C4D28, D100, W25, W52, W105, relapse)
total.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
total.integrated <- IntegrateData(anchorset = total.anchors, dims = 1:30)
myeloma.integrated <- total.integrated

# Post-processing merged data
print("Post-processing first batch ...")
DefaultAssay(object = myeloma.integrated) <- "integrated"
myeloma.integrated <- ScaleData(object = myeloma.integrated, verbose = F)
myeloma.integrated <- RunPCA(object = myeloma.integrated, verbose = F)
#ElbowPlot(object = myeloma.integrated, ndims = 40)
myeloma.integrated <- RunUMAP(object = myeloma.integrated, reduction = "pca", dims = 1:20)
myeloma.integrated <- FindNeighbors(object = myeloma.integrated, dims = 1:20)
myeloma.integrated <- FindClusters(object = myeloma.integrated, resolution = 0.3)

# Saving and loading RDS files
print("Saving integrated RDS file first batch ...")
mcsaveRDS(myeloma.integrated, file="/data/output/myeloma_alltimepoints.rds", threads = 4)
myeloma_alltimepoints <- myeloma.integrated 

# --------------------------------------------------- #
# Making total CD38+ immune object (MM + control)
# --------------------------------------------------- #

# Identification of integration anchors
print("Integrating second batch ...")
reference.list <- c(baseline_control, myeloma_alltimepoints)
total.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
total.integrated <- IntegrateData(anchorset = total.anchors, dims = 1:30)

# Post-processing merged data
print("Post-processing second batch ...")
DefaultAssay(object = total.integrated) <- "integrated"
total.integrated <- ScaleData(object = total.integrated, verbose = F)
total.integrated <- RunPCA(object = total.integrated, verbose = F)
ElbowPlot(object = total.integrated, ndims = 40)
total.integrated <- RunUMAP(object = total.integrated, reduction = "pca", dims = 1:20)
total.integrated <- FindNeighbors(object = total.integrated, dims = 1:20)
total.integrated <- FindClusters(object = total.integrated, resolution = 0.3)

# Saving and loading RDS files
print("Saving integrated RDS file second batch ...")
mcsaveRDS(total.integrated, file="/data/output/Total_CD38.rds", threads = 4)
Total_CD38 <- mcreadRDS("/data/output/Total_CD38.rds", threads = 4)

###
plan('sequential')
###

# --------------------------------------------------- #
# Post-processing
# --------------------------------------------------- #

# Removed unwanted admixture cells either by subsetting, or by selecting cells with Cellselector (vector called cells), and:
# cells_total <- myeloma.integrated@assays$RNA@counts@Dimnames
# cells_total <- cells_total[2]
# cells_total <- unlist(cells_total)
# cells_total <- cells_total[!cells_total %in% cells]
# myeloma.integrated <- subset(myeloma.integrated, cells = "cells_total")

DefaultAssay(object = Total_CD38) <- "integrated"
Total_CD38 <- ScaleData(object = Total_CD38, verbose = F)
Total_CD38 <- RunPCA(object = Total_CD38, verbose = F)
ElbowPlot(Total_CD38, ndims = 30, reduction = "pca")
Total_CD38 <- RunUMAP(object = Total_CD38, reduction = "pca", dims = 1:15, return.model = T)
Total_CD38 <- FindNeighbors(object = Total_CD38, dims = 1:15)
Total_CD38 <- FindClusters(object = Total_CD38, resolution = 0.3)


# --------------------------------------------------- #
# Inspect NK cell and other immune cell markers 
# --------------------------------------------------- #
DimPlot(Total_CD38)
#PC markers
FeaturePlot(Total_CD38, features = c("SDC1", "LAMP5", "SLAMF7"))
#T cell markers
FeaturePlot(Total_CD38, features = c("CD8A", "CD8B", "CD4"))
FeaturePlot(Total_CD38, features = c("CD3D", "CCR7"))
# B cell markers
FeaturePlot(Total_CD38, features = c("MS4A1", "CD19", "VPREB1"))
# Monocyte
FeaturePlot(Total_CD38, features = c("LYZ", "CD14", "FCGR3B"))
#Erythrocytes
FeaturePlot(Total_CD38, features = c("HBA1"))
#DCs
FeaturePlot(Total_CD38, features = c("FCER1A", "CD1C"))
#Neutrophils
FeaturePlot(Total_CD38, features = c("ELANE"))
#NK
FeaturePlot(Total_CD38, features = c("KLRD1", "KLRF1", "NKG7", "GNLY"))
FeaturePlot(Total_CD38, features = c("GZMK", "GZMB"))
FeaturePlot(Total_CD38, features = c("NCAM1", "FCGR3A"))
FeaturePlot(Total_CD38, features = c("NCR1", "NCR3"))

# --------------------------------------------------- #
# Inspect differential expressed genes between clusters to confirm NK cell identiy
# --------------------------------------------------- #

top10 <- head(VariableFeatures(Total_CD38), 10)
top10

Total_CD38_DEG <- FindAllMarkers(Total_CD38, only.pos = TRUE, min.pct = 0.25, logfc.treshold = 0.25)
Total_CD38_DEG

# --------------------------------------------------- #
# Subset and reprocess NK cell dataset
# --------------------------------------------------- #

Total_NK <- Total_CD38[,Idents(Total_CD38) %in% c(#selected_cluster,#selected_cluster,#selected_cluster)]

DefaultAssay(Total_NK) <- "integrated"

Total_NK <- FindVariableFeatures(Total_NK, selection.method = "vst", nfeatures = 2000)
Total_NK <- ScaleData(Total_NK, verbose = FALSE)
Total_NK <- RunPCA(Total_NK, npcs = 30, verbose = FALSE)
Total_NK <- RunUMAP(Total_NK, reduction = "pca", dims = 1:25)
Total_NK <- FindNeighbors(Total_NK, reduction = "pca", dims = 1:25)
Total_NK <- FindClusters(Total_NK, resolution = 0.4)
DimPlot(Total_NK, reduction = "umap", label = TRUE)
DimPlot(Total_NK, reduction = "umap", group.by = "source")
DimPlot(Total_NK, reduction = "umap", group.by = "patient")
DimPlot(Total_NK, reduction = "umap", group.by = "timepoint")

#repeat process if there is still contamination of other immune subsets

# --------------------------------------------------- #
# For convenience: saving and loading after preprocessing
# --------------------------------------------------- #
saveRDS(Total_NK, file = "~/total_NK.rds")
Total.NK <- readRDS(file = "~/total_NK.rds")
