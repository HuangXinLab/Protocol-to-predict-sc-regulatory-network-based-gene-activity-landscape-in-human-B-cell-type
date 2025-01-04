##R version 4.3.2

library(scMINER)
library(NetBID2)
library(ggplot2)

#Heatmap visualization
genes_of_interest <-c("FLT3", "EBF1", "MS4A1","BCL2", "BCL2L1", "MCL1")
feature_heatmap(input_eset = acs_sc, 
				target = genes_of_interest, 
				group_name = "celltype",
				save_plot = FALSE, 
				width = 6, 
				height = 6, 
				name = "log2Exp")

#Feature plot visualization
feature_highlighting(input_eset = acs_sc, 
					target = genes_of_interest, 
					feature = "geneSymbol", 
					ylabel = "log2Exp", 
					x = "X", 
					y = "Y", 
					pct.size = 0.5)

#Violin plot visualization
feature_vlnplot(input_eset = acs_sc, 
				target = genes_of_interest, 
				feature = "geneSymbol", 
				group_by = "celltype", 
				ylabel = "log2Exp", 
				ncol = 4)
				
#Network visualization
DAG_result <- get.DA(input_eset = acs_sc[,acs_sc$scminer==celltype], group_name = "scminer")
# Calculate edge score for each target and celltype
edge_score <-SIG.table$target_list[[use_driver]]$MI*sign(SIG.table$target_list[[use_driver]]$spearman)
names(edge_score) <- SIG.table$target_list[[use_driver]]$target
# draw the targetNet
NetBID2::draw.targetNet(source_label=use_driver,
						source_z=DAG_result[DAG_result$geneSymbol==use_driver,
						sprintf('degree_%s',celltype)],
						edge_score = edge_score,
						pdf_file=sprintf('networks/%s_argetNet_out_%s_label_cex01.pdf',
						celltype,use_driver),
						label_cex = 0.5,
						n_layer=1,
						source_cex = 0.5, 
						alphabetical_order=FALSE)
