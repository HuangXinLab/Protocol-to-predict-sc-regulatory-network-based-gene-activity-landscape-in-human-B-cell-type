## R 4.3.2
# /mnt/sda/Public/Environment/miniconda3/envs/sc/bin/Rscript 4.2.network_visualization.r HSC

suppressMessages(library(NetBID2))
suppressMessages(library(scMINER))
args <- commandArgs(trailingOnly = TRUE)

celltype <- args[1]
load('/mnt/sda/Public/Project/B-ALL/Bdev/activity/acs_sc')
load(sprintf('/mnt/sda/Public/Project/B-ALL/Bdev/networks/%s.TF.network',celltype))
load(sprintf('/mnt/sda/Public/Project/B-ALL/Bdev/networks/%s.SIG.network',celltype))

fData(acs_sc)['geneSymbol'] <- fData(acs_sc)['gene_short_name']
capt <- capture.output(DAG_result <- get.DA(input_eset = acs_sc[,acs_sc$scminer==celltype], group_name = "scminer"))

### Calculate edge score for each target
if(celltype %in% c('PreB','ImmatureB','MatureB')) {
   
} else {
    use_driver='FLT3'
    edge_score <-SIG.table$target_list[[use_driver]]$MI*sign(SIG.table$target_list[[use_driver]]$spearman)
    names(edge_score) <- SIG.table$target_list[[use_driver]]$target
    ### DAG_result is the TF DA master table ##########
    NetBID2::draw.targetNet(source_label=use_driver,source_z=DAG_result[DAG_result$geneSymbol==use_driver,sprintf('degree_%s',celltype)], 
                        edge_score = edge_score,pdf_file=sprintf('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/networks/%s_argetNet_out_%s.pdf',celltype,use_driver),label_cex = 0.4,n_layer=1,source_cex = 0.5, alphabetical_order=FALSE)
}

use_driver='EBF1'
edge_score <-TF.table$target_list[[use_driver]]$MI*sign(TF.table$target_list[[use_driver]]$spearman)
names(edge_score) <- TF.table$target_list[[use_driver]]$target
### DAG_result is the TF DA master table ##########
NetBID2::draw.targetNet(source_label=use_driver,source_z=DAG_result[DAG_result$geneSymbol==use_driver,sprintf('degree_%s',celltype)], 
                    edge_score = edge_score,pdf_file=sprintf('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/networks/%s_argetNet_out_%s.pdf',celltype,use_driver),label_cex = 0.4,n_layer=1,source_cex = 0.5, alphabetical_order=FALSE)
use_driver='MS4A1'
edge_score <-SIG.table$target_list[[use_driver]]$MI*sign(SIG.table$target_list[[use_driver]]$spearman)
names(edge_score) <- SIG.table$target_list[[use_driver]]$target
### DAG_result is the TF DA master table ##########
NetBID2::draw.targetNet(source_label=use_driver,source_z=DAG_result[DAG_result$geneSymbol==use_driver,sprintf('degree_%s',celltype)], 
                        edge_score = edge_score,pdf_file=sprintf('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/networks/%s_argetNet_out_%s.pdf',celltype,use_driver),label_cex = 0.4,n_layer=1,source_cex = 0.5, alphabetical_order=FALSE)

use_driver='BCL2'
edge_score <-SIG.table$target_list[[use_driver]]$MI*sign(SIG.table$target_list[[use_driver]]$spearman)
names(edge_score) <- SIG.table$target_list[[use_driver]]$target
### DAG_result is the TF DA master table ##########
NetBID2::draw.targetNet(source_label=use_driver,source_z=DAG_result[DAG_result$geneSymbol==use_driver,sprintf('degree_%s',celltype)], 
                        edge_score = edge_score,pdf_file=sprintf('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/networks/%s_argetNet_out_%s.pdf',celltype,use_driver),label_cex = 0.4,n_layer=1,source_cex = 0.5, alphabetical_order=FALSE)

use_driver='BCL2L1'
edge_score <-SIG.table$target_list[[use_driver]]$MI*sign(SIG.table$target_list[[use_driver]]$spearman)
names(edge_score) <- SIG.table$target_list[[use_driver]]$target
### DAG_result is the TF DA master table ##########
NetBID2::draw.targetNet(source_label=use_driver,source_z=DAG_result[DAG_result$geneSymbol==use_driver,sprintf('degree_%s',celltype)], 
                        edge_score = edge_score,pdf_file=sprintf('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/networks/%s_argetNet_out_%s.pdf',celltype,use_driver),label_cex = 0.4,n_layer=1,source_cex = 0.5, alphabetical_order=FALSE)

use_driver='MCL1'
edge_score <-SIG.table$target_list[[use_driver]]$MI*sign(SIG.table$target_list[[use_driver]]$spearman)
names(edge_score) <- SIG.table$target_list[[use_driver]]$target
### DAG_result is the TF DA master table ##########
NetBID2::draw.targetNet(source_label=use_driver,source_z=DAG_result[DAG_result$geneSymbol==use_driver,sprintf('degree_%s',celltype)], 
                        edge_score = edge_score,pdf_file=sprintf('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/networks/%s_argetNet_out_%s.pdf',celltype,use_driver),label_cex = 0.4,n_layer=1,source_cex = 0.5, alphabetical_order=FALSE)
