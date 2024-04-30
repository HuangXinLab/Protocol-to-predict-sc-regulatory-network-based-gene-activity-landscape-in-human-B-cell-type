# Protocol-to-predict-sc-regulatory-network-based-gene-activity-landscape-in-human-B-cell-type
## 1.Single-cell classic plot: cell_subset_and_cell_cycle_scoring
The human bone marrow (BM) data were sourced from the immune cell section of the Human Cell Atlas. We utilized Seurat for single-cell data processing, removal of low-quality cells, and standard single-cell data analysis, including B cell annotation. See scripts: 1.cell_subset_and_cell_cycle_scoring.R
## 2.network building
We used the SJARACNe algorithm to reconstruct cell-type-specific interactomes for each B cell type (HSC, CLP, Pre-pro-B, Pro-B, Pre-B, Immature B, Mature B, and Plasma cells) from the scRNA-seq profiles in the Human Cell Atlas Immune Cells study. This generated separate TF and SIG networks. See scripts: 2.1.generateSJARACNeInput.R; 2.2.generateSJARACNenetwork.sh
## 3.infer activity
We employed NetBID2 and scMINER to analyze the networks, identifying putative drivers of B cell developmental states. See scripts: 3.GetActivityFromSJARACNe.R
## 4.visualization
For the network-based activity of TF and SIG drivers inferred by NetBID2 and scMINER, we used different visualization methods for display. See scripts: 4.1.visualization.R; 4.2.network_visualization.R
