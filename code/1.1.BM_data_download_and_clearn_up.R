# R version 4.3

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HCAData")

library(Seurat)
library(dplyr)
library(HCAData)

# 1. load sce_BM data
suppressPackageStartupMessages({
  library("ExperimentHub")
  library("SingleCellExperiment")
})
eh <- ExperimentHub()
bonemarrow_h5densematrix <- eh[["EH2047"]]
bonemarrow_coldata <- eh[["EH2048"]]
bonemarrow_rowdata <- eh[["EH2049"]]
sce_bonemarrow <- HCAData("ica_bone_marrow")
save(sce_bonemarrow, file = 'sce_bonemarrow.rdata')

# 2. Load data and seurat object generation
load('sce_bonemarrow.rdata')
# Modify cell barcodes
col <- gsub('-1', '', colData(sce_bonemarrow)$Barcode)
colData(sce_bonemarrow)$Barcode <- col
rownames(colData(sce_bonemarrow)) <- col

# Read cell type annotation data
celltype <- read.delim2('CensusImmune-BoneMarrow-10x_cell_type_2020-03-12.csv', fill = TRUE, stringsAsFactors = FALSE, sep = ',')
celltype$barcode1 <- paste0('Manton', gsub('._','',gsub('_cells','',celltype$cell_suspension.biomaterial_core.biomaterial_id)), '_HiSeq_', gsub('_.*','',celltype$cell_suspension.biomaterial_core.biomaterial_id),'-', celltype$barcode)

# Subset data
sce_bonemarrow_sub <- sce_bonemarrow[, celltype$barcode1[celltype$barcode1 %in% sce_bonemarrow$Barcode]]

# Convert to Seurat object
counts <- assay(sce_bonemarrow_sub, "counts")
counts_sparse <- as(counts, "sparseMatrix")
assay(sce_bonemarrow_sub, "counts") <- counts_sparse
seurat_bonemarrow_sub <- as.Seurat(sce_bonemarrow_sub, counts = "counts", data = NULL)

# 3. Metadata processing and cell type annotation
# Merge metadata
rownames(celltype) <- celltype$barcode1
seurat_bonemarrow_sub@meta.data <-left_join(seurat_bonemarrow_sub@meta.data, celltype, by = c("Barcode" = "barcode1"))

# Modify cell type format
pd <- seurat_bonemarrow_sub@meta.data
pd$celltype <- NULL

# Define cell type mapping
cell_type_mapping <- c( 'B cell T cell doublet' = 'doublet',  'CD14+ monocyte type 1' = 'monocyte',  'CD14+ monocyte type 2' = 'monocyte',  
                       'CD16+ monocyte' = 'monocyte',  'CD4+ naive T cell' = 'T',  'conventional dendritic cell' = 'DC',  'cytotoxic T cell type 1' = 'T',  
                       'cytotoxic T cell type 2' = 'T',  'erythroid cell type 1' = 'erythroid',  'erythroid cell type 2' = 'erythroid',  
                       'hematopoietic stem cell' = 'HSC',  'megakaryocyte' = 'megakaryocyte',  'memory B cell' = 'B',  'mesenchymal stem cell' = 'MSC',  
                       'naive B cell' = 'B',  'naive CD8+ T cell' = 'T',  'natural killer cell' = 'NK',  'plasma cell' = 'B',  'plasmacytoid dendritic cell' = 'DC',  
                       'precursor B cell' = 'B',  'pro-B cell' = 'B',  'T-helper cell' = 'T')
pd$celltype <- cell_type_mapping[pd$annotated_cell_identity.text]
seurat_bonemarrow_sub@meta.data <- pd

# 4. Gene name conversion
# Gene name conversion
counts <- seurat_bonemarrow_sub@assays$RNA@counts
a <- distinct(as.data.frame(rowData(sce_bonemarrow)), Symbol, .keep_all = TRUE)
a <- a[!a$Symbol %in% grep('_', a$Symbol, value = TRUE), ]
rownames(a) <- a$Symbol
counts <- counts[a$ID, ]

# Create new assay
Hum_gene_assay <- CreateAssayObject(counts)
seurat_bonemarrow_sub[['RNA']] <- Hum_gene_assay
seurat_bonemarrow_sub@assays$RNA@counts@Dimnames[[1]] <- a$Symbol
seurat_bonemarrow_sub@assays$RNA@data@Dimnames[[1]] <- a$Symbol
seurat_bonemarrow_sub@assays$RNA@meta.features <- a

# Save processed object
save(seurat_bonemarrow_sub, file = 'seurat_bonemarrow_sub.rdata')
