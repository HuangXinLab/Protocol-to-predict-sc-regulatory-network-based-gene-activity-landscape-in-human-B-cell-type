##R version 4.3.2

library(Seurat)
library(ggplot2)

load('B_sub.rdata')
pdf('fig/cell_type_UMAP.pdf')
DimPlot(B_sub, reduction = "umap", label = TRUE)
dev.off()

pdf('fig/cell_phase.pdf')
DimPlot(B_sub, reduction = "umap", group.by = 'Phase') + ggtitle('')
dev.off()
