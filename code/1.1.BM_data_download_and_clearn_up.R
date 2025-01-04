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
query(eh, "HCAData")
AnnotationHub::query( eh, 'HCAData')
sce_bonemarrow <- HCAData::HCAData("ica_bone_marrow")
save(sce_bonemarrow, file = 'sce_bonemarrow.rdata')

# 2. Load data and seurat object generation
load('sce_bonemarrow.rdata')
# Modify cell barcodes
SummarizedExperiment::colData( sce_bonemarrow)[[ 'Barcode']] <- gsub( '-1$', '', SummarizedExperiment::colData( sce_bonemarrow)$Barcode)

# Read cell type annotation data
celltype <- read.delim2('CensusImmune-BoneMarrow-10x_cell_type_2020-03-12.csv', fill = TRUE, stringsAsFactors = FALSE, sep = ',')
celltype$barcode1 <- paste0('Manton', gsub('._','',gsub('_cells','',celltype$cell_suspension.biomaterial_core.biomaterial_id)), '_HiSeq_', gsub('_.*','',celltype$cell_suspension.biomaterial_core.biomaterial_id),'-', celltype$barcode)

# Subset data
counts <- assay( sce_bonemarrow, "counts")
barcode2use <- intersect( celltype$barcode1, sce_bonemarrow$Barcode)
column2use <- match( barcode2use, sce_bonemarrow$Barcode)
counts <- counts[, column2use]
dimnames( counts) <- list( rownames( rowData( sce_bonemarrow)), barcode2use)
counts <- counts[, barcode2use[ order( column2use)]]
seurat_bonemarrow_sub <- Seurat::CreateSeuratObject( counts = as( counts, 'dgCMatrix'))

# 3. Metadata processing and cell type annotation
# Merge metadata
i <- match( colnames( seurat_bonemarrow_sub), celltype[[ 'barcode1']])
cell_type_supplied <- celltype[ i, 'annotated_cell_identity.text']

# Define cell type mapping
cell_type_mapping <- c( 'B cell T cell doublet' = 'doublet',  'CD14+ monocyte type 1' = 'monocyte',  'CD14+ monocyte type 2' = 'monocyte',  
                       'CD16+ monocyte' = 'monocyte',  'CD4+ naive T cell' = 'T',  'conventional dendritic cell' = 'DC',  'cytotoxic T cell type 1' = 'T',  
                       'cytotoxic T cell type 2' = 'T',  'erythroid cell type 1' = 'erythroid',  'erythroid cell type 2' = 'erythroid',  
                       'hematopoietic stem cell' = 'HSC',  'megakaryocyte' = 'megakaryocyte',  'memory B cell' = 'B',  'mesenchymal stem cell' = 'MSC',  
                       'naive B cell' = 'B',  'naive CD8+ T cell' = 'T',  'natural killer cell' = 'NK',  'plasma cell' = 'B',  'plasmacytoid dendritic cell' = 'DC',  
                       'precursor B cell' = 'B',  'pro-B cell' = 'B',  'T-helper cell' = 'T')
# Modify cell type format
cell_type_used <- cell_type_mapping[ cell_type_supplied]
seurat_bonemarrow_sub@meta.data[, 'celltype'] <- cell_type_used

# 4. Gene name conversion
# Gene name conversion
counts <- seurat_bonemarrow_sub@assays$RNA@counts
row2rm <- which( duplicated( rowData( sce_bonemarrow)[, 'Symbol']))
seurat_bonemarrow_sub <- seurat_bonemarrow_sub[ -row2rm,]

#Replace gene IDs with gene symbols
row_name <- rownames( seurat_bonemarrow_sub)
row_name <- rowData( sce_bonemarrow)[ row_name, 'Symbol']
seurat_bonemarrow_sub@assays$RNA@counts@Dimnames[[ 1]] <- row_name
seurat_bonemarrow_sub@assays$RNA@data@Dimnames[[ 1]] <- row_name
dimnames( seurat_bonemarrow_sub@assays$RNA@meta.features)[[ 1]] <- row_name

#Filter out genes with undesirable symbols
row2remove <- grep( '_', rownames( seurat_bonemarrow_sub)) 
if( length( row2remove))
    seurat_bonemarrow_sub <- seurat_bonemarrow_sub[ -row2remove]
save( seurat_bonemarrow_sub, file = 'seurat_bonemarrow_sub.rdata') 

