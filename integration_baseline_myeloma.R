# --------------------------------------------------- #
# Analysis of scRNA sequencing data as described in Tahri et al. 2023
# Script by Sabrin Tahri, optimized for speed by Mathijs Sanders and Gregory van Beek 
# Dept. of Hematology, Erasmus MC Cancer Institute, Rotterdam, the Netherlands
# --------------------------------------------------- #

# --------------------------------------------------- #
# READ-ME (IMPORTANT)
# This script is meant to generate an integrated CD38+ immune object from all newly diagnosed MM patients (n=19). 
# A seperate script is available on github explaining the creation of the total NK cell object by in silico selection of NK cells from 
# the total CD38+ immune dataset, This script is named 'script total BM NK.R'. 

# --------------------------------------------------- #
# Loading libraries
# --------------------------------------------------- #

pacman::p_load('Seurat', 'future', 'parallel', 'cowplot', 'ggplot2', 'dplyr')

options(future.rng.onMisuse="ignore")
options(future.globals.maxSize = 1e4 * 1024^2)
options(future.fork.enable = TRUE)

# --------------------------------------------------- #
# Function for loading in objects
# --------------------------------------------------- #

loadSeuratObject <- function(idx, projectName, matLoc, bcLoc, selectionCriteria, details) {
  info <- structure(selectionCriteria[idx, ], names = colnames(selectionCriteria))
  currentDetails <- structure(details[idx, ], names = colnames(details))
  so <- Read10X(data.dir = matLoc[idx])
  so <- CreateSeuratObject(counts = so, min.cells = 3, min.features = 200, project = projectName)
  so <- AddMetaData(so, PercentageFeatureSet(so, pattern = "^MT-"), col.name = "percent.mito")
  so <- subset(x = so, subset = nFeature_RNA > info$nFeature_RNA_lower & nFeature_RNA < info$nFeature_RNA_higher & nCount_RNA > info$nCount_RNA_lower & nCount_RNA < info$nCount_RNA_higher & percent.mito < info$percent.mito)
  so <- NormalizeData(object = so, normalization.method = "LogNormalize", scale.factor = 1e4)
  so <- FindVariableFeatures(object = so, selection.method = "vst", nfeatures = 2000)
  so$source <- currentDetails$Source
  so$patient <- currentDetails$Patient
  so
}

# --------------------------------------------------- #
# Function for defining number of cores
# --------------------------------------------------- #

mcPredefined <- function(threads) {
  function(...) {
    mclapply(..., mc.cores = threads)
  }
}

# --------------------------------------------------- #
# Object information
# --------------------------------------------------- #

projectName <- 'NDMM'
threads <- 10
matLoc <- readLines('~/samples.txt')
selectionCriteria <- read.table('~/selection.txt', header = TRUE)
details <- read.table('~/details.txt', header = TRUE)

# --------------------------------------------------- #
# Making objects
# --------------------------------------------------- #

so.list <- ifelse(threads > 1, get('mcPredefined')(threads), get('lapply'))(seq(matLoc), loadSeuratObject, projectName, matLoc, bcLoc, selectionCriteria, details)

# --------------------------------------------------- #
# Integration 
# --------------------------------------------------- #

plan(strategy = 'multicore', workers = 15)

reference.list <- so.list[c(1:19)]
NDMM.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
NDMM.integrated <- IntegrateData(anchorset = NDMM.anchors, dims = 1:30)

# --------------------------------------------------- #
# Post-processing
# --------------------------------------------------- #

DefaultAssay(object = NDMM.integrated) <- "integrated"
NDMM.integrated <- ScaleData(object = NDMM.integrated, verbose = F)
NDMM.integrated <- RunPCA(object = NDMM.integrated, verbose = F)
ElbowPlot(NDMM.integrated, ndims = 30, reduction = "pca")
NDMM.integrated <- RunUMAP(object = NDMM.integrated, reduction = "pca", dims = 1:15, return.model = T)
NDMM.integrated <- FindNeighbors(object = NDMM.integrated, dims = 1:15)
NDMM.integrated <- FindClusters(object = NDMM.integrated, resolution = 0.3)

# --------------------------------------------------- #
# For convenience: saving and loading after preprocessing
# --------------------------------------------------- #
saveRDS(NDMM.integrated, file = "~/baseline_myeloma.rds")
baseline_myeloma <- readRDS(file = "~/baseline_myeloma.rds")