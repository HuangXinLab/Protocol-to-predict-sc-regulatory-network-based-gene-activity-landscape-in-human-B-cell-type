# R version 4.3

library(Seurat)
library(ggplot2)

# load data
load('seurat_bonemarrow_sub.rdata')

# Quality control
seurat_bonemarrow_sub[[ 'percent.mt']] <- Seurat::PercentageFeatureSet( seurat_bonemarrow_sub, pattern = '^MT-')
seurat_bonemarrow_sub <- subset( seurat_bonemarrow_sub, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

# Gene expression normalization
seurat_bonemarrow_sub <- Seurat::NormalizeData( seurat_bonemarrow_sub)
seurat_bonemarrow_sub <- Seurat::FindVariableFeatures( seurat_bonemarrow_sub, selection.method = "vst", nfeatures = 2000)

# Cell-Cycle Scoring
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes
seurat_bonemarrow_sub <- Seurat::CellCycleScoring( seurat_bonemarrow_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_bonemarrow_sub$CC.Difference <- seurat_bonemarrow_sub$S.Score - seurat_bonemarrow_sub$G2M.Score

# Scaling the data
seurat_bonemarrow_sub <- Seurat::ScaleData( seurat_bonemarrow_sub, vars.to.regress = c( "nFeature_RNA",  "percent.mt", "CC.Difference"))

# Dimensionality Reduction
seurat_bonemarrow_sub <- Seurat::RunPCA( seurat_bonemarrow_sub, features = Seurat::VariableFeatures( object = seurat_bonemarrow_sub))
seurat_bonemarrow_sub <- Seurat::RunUMAP( seurat_bonemarrow_sub, dims = 1:30)

# Cell Clustering and annotation
seurat_bonemarrow_sub <- Seurat::FindNeighbors( seurat_bonemarrow_sub, dims = 1:30)
seurat_bonemarrow_sub <- Seurat::FindClusters( seurat_bonemarrow_sub, resolution = 0.4)

#Identify unique and highly expressed genes in each cluster
seurat_bonemarrow_sub.markers <- Seurat::FindAllMarkers(seurat_bonemarrow_sub, only.pos = TRUE)

# extract B cell and subclustering
B_sub <- seurat_bonemarrow_sub[, seurat_bonemarrow_sub[[ 'celltype']][[ 1]] %in% c( 'B', 'HSC')]

#NormalizeData
B_sub <- Seurat::NormalizeData( B_sub)
#FindVariableFeatures
B_sub <- Seurat::FindVariableFeatures( B_sub, selection.method = 'vst', nfeatures = 2000)
#scaleData
B_sub <- Seurat::ScaleData( B_sub, vars.to.regress = c( 'nFeature_RNA',  'percent.mt', 'CC.Difference')) # slow
#Dimensionality Reduction
B_sub <- Seurat::RunPCA( B_sub)
B_sub <- Seurat::RunUMAP( B_sub, dims = 1:30)
#Cell Clustering and annotation
B_sub <- Seurat::FindNeighbors( B_sub, dims = 1:30)
B_sub <- Seurat::FindClusters( B_sub, resolution = 0.4)
save( B_sub, file = 'B_sub.RData')

#Cell-Cycle scoring visualization
Seurat::DimPlot( B_sub, reduction = 'umap', group.by = 'Phase')

