# Gene look-up in the sequenced data of Tahri et al. 

# --------------------------------------------------- #
# Installing Seurat
# --------------------------------------------------- #

# We provided you with a Seurat object, so you'll need the package 'Seurat' to look up genes.
# If already installed, you can skip this part.
# To install, the code is as follows:

install.packages('Seurat')

# --------------------------------------------------- #
# Loading Seurat
# --------------------------------------------------- #

# Everytime you restart R studio, you have to reload the Seurat package

library(Seurat)

# --------------------------------------------------- #
# Loading dataset
# --------------------------------------------------- #

all.integrated <- readRDS(file = "location_of_the_file_on_your_disk/Total_NK_CD38pos.rds")

# --------------------------------------------------- #
# Visualization 
# --------------------------------------------------- #

DimPlot(object = all.integrated, reduction = "umap", pt.size = 1.5, label=T, cols = c("#ef8a62", "#b2182b", "#2166ac", "#67a9cf", "#053061")) # This will visualize the plain UMAP
DimPlot(object = all.integrated, reduction = "umap", pt.size = 1.5, split.by = "patient", ncol=3) # This will visualise the UMAP per patient
DimPlot(object = all.integrated, reduction = "umap", pt.size = 1.5, split.by = "source", ncol=3) # This will visualise the UMAP according source (control vs myeloma)

# To look up genes, you have a few different option.
# Fill out your gene of interest after 'features'. Make sure to use "". 
FeaturePlot(object = all.integrated, features = "source", split.by = "source", reduction = "umap", min.cutoff = 0,  pt.size=1, sort.cell = T, cols = c("lightgrey", "darkred")) 
VlnPlot(object = all.integrated, features = "source", pt.size = 1, assay = "RNA")
DotPlot(object = all.integrated, features = c("GZMK", "GZMB", "FCGR3A", "NCAM1"), dot.scale = 8, col.min = 0, cols = c("lightgrey", "red"))

