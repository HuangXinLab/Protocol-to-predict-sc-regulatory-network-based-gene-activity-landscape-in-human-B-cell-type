##R version 4.3.2

library(scMINER)
library(Seurat)


# 1_seurat2eset #######
load('./B_sub.rdata')

#expression_data
expression_data<- B_sub@assays$RNA@data 
#cell_metadata:
data <- B_sub@meta.data
i <- match( rownames( data), celltype[, 'barcode1'])
data[, 'celltype'] <- gsub( '[^a-zA-Z0-9_\\.]', '_', celltype[ i, 'annotated_cell_identity.text'])
data[, 'celltype_show'] <- celltype[ i, 'annotated_cell_identity.text']
cell_metadata <- new( 'AnnotatedDataFrame', data = data)

#Generate gene annotation and create sparse expression set (eSet)
#gene_annotation
gene_annotation <- with( list( x = row.names( expression_data)), data.frame( geneSymbol = x, row.names = x))
gene_annotation$nCells <- Matrix::rowSums( expression_data != 0)
#new eset
eset <- scMINER::createSparseEset( input_matrix = expression_data, cellData = cell_metadata@data, featureData = gene_annotation, addMetaData = FALSE)

# filter out genes expressed in less than 50cells
eset.filter.cell.log2 <- eset[ fData( eset)$nCells > 50,]
save( eset.filter.cell.log2, file = 'eset.filter.cell.log2.rdata')
# generateSJARACNeInput
load( 'eset.filter.cell.log2.rdata')
scMINER::generateSJARACNeInput( input_eset = eset.filter.cell.log2, 
												species_type = 'hg', 
												sjaracne_dir = 'SJARACNe', 
												group_name = 'celltype', 
												downSample_N = 1e2, 
												driver_type = 'TF_SIG')

