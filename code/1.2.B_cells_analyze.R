# R version 4.3

library(Seurat)
library(ggplot2)

# load data
load('seurat_bonemarrow_sub.rdata')

# Quality control
seurat_bonemarrow_sub[["percent.mt"]] <- PercentageFeatureSet(seurat_bonemarrow_sub, pattern = "^MT-")
seurat_bonemarrow_sub <- subset(seurat_bonemarrow_sub, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

# Gene expression normalization
seurat_bonemarrow_sub <- NormalizeData( seurat_bonemarrow_sub)
seurat_bonemarrow_sub <- FindVariableFeatures( seurat_bonemarrow_sub, selection.method = "vst", nfeatures = 2000)

# Cell-Cycle Scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes 
seurat_bonemarrow_sub <- CellCycleScoring(seurat_bonemarrow_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_bonemarrow_sub $CC.Difference <- seurat_bonemarrow_sub $S.Score - seurat_bonemarrow_sub$G2M.Score

# Scaling the data
seurat_bonemarrow_sub <- ScaleData( seurat_bonemarrow_sub, vars.to.regress = c("nFeature_RNA", "percent.mit","CC.Difference"))

# Dimensionality Reduction
seurat_bonemarrow_sub <- RunPCA(seurat_bonemarrow_sub, features = VariableFeatures(object = seurat_bonemarrow_sub))
seurat_bonemarrow_sub <- RunUMAP(seurat_bonemarrow_sub, dims = 1:30)

# Cell Clustering and annotation
seurat_bonemarrow_sub <- FindNeighbors(seurat_bonemarrow_sub, dims = 1:30)
seurat_bonemarrow_sub <- FindClusters(seurat_bonemarrow_sub, resolution = 0.4)
seurat_bonemarrow_sub.markers <- FindAllMarkers(seurat_bonemarrow_sub, only.pos = TRUE)


# extract B cell and subclustering
B_temp <- seurat_bonemarrow_sub[, Idents(seurat_bonemarrow_sub) %in% c("B","HSC")]
B_sub <- CreateSeuratObject(counts = B_temp @assays$RNA@counts,
                             meta.data = B_temp@meta.data)
#NormalizeData
B_sub <- NormalizeData(B_sub)
#FindVariableFeatures
B_sub <- FindVariableFeatures(B_sub, selection.method = "vst", nfeatures = 2000)
#scaleData
all.genes <- rownames(B_sub)
B_sub <- ScaleData(B_sub, features = all.genes)
#Dimensionality Reduction
B_sub <- RunPCA(B_sub, features = VariableFeatures(B_sub))
B_sub <- RunUMAP(B_sub, dims = 1:30)
#Cell Clustering and annotation
B_sub <- FindNeighbors(B_sub, dims = 1:30)
B_sub <- FindClusters(B_sub, resolution = 0.4)
save(B_sub, file="B_sub.RData")


# plot
load('B_sub.RData')
pdf('cell_type_UMAP.pdf')
DimPlot(B_sub, reduction = "umap", label = TRUE)
dev.off()

pdf('cell_phase.pdf')
DimPlot(B_sub, reduction = "umap", group.by = 'Phase') + ggtitle('')
dev.off()
