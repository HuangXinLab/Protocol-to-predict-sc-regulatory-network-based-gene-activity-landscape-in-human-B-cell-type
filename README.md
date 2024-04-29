# Protocol-to-predict-sc-regulatory-network-based-gene-activity-landscape-during-human-B-cell-type
## 1.Single-cell classic plot: cell_subset_and_cell_cycle_scoring
we utilized the Seurat R packages to analyzis biggest healthy human BM single-cell RNA sequencing dataset in human cell atlas Atlas Census of Immune Cells project. See scripts: 1.cell_subset_and_cell_cycle_scoring.R
## 2.network building
we employed the SJARACNe algorithm to reconstruct cell type-specific interactomes for each B cell stage(HSC, CLP, Pre-pro-B, Pro-B, Pre-B, Immature B, Mature B, and Plasma cells), generating separate TF and SIG networks respectively. See scripts: 2.1.generateSJARACNeInput.R; 2.2.generateSJARACNenetwork.sh
## 3.infer activity
we used the NetBID2 to estimate the activity of each driver gene based on gene expression profiles and the corresponding network for each cell type. See scripts: 3.GetActivityFromSJARACNe.R
## 4.visualization
For the transcriptional activity of genes and the drivers activity inferred by NetBID2, we use different visualization methods for display. See scripts: 4.1.visualization.R; 4.2.network_visualization.R
