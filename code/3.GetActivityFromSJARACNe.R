##R version 4.3.2

library(NetBID2)
library(scMINER)

load('SJARACNe/Input.eset')

#Calculate activity
acs_sc <- scMINER::getActivity_inBatch(input_eset = eset.filter.cell.log2,
										sjaracne_dir = 'SJARACNe',
										group_name = 'celltype',
										driver_type = 'TF_SIG',
										activity_method = 'mean',
										do.z_normalization = FALSE)

#Driver estimation by differential activity analysis
DAG_result <- scMINER::getDA( input_eset = acs_sc, group_by = 'celltype')