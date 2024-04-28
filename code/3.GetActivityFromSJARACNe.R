##R version 4.3.2

library(NetBID2)
library(scMINER)

load('SJARACNe/Input.eset')
acs_sc <- GetActivityFromSJARACNe( 
SJARACNe_output_path ="./SJARACNe/",
SJARACNe_input_eset = input_eset,
activity.method="unweighted", # we highly recommend using 'unweighted' as activity calculation method
activity.norm=TRUE, 
group_name = "mergeCycling", # which group was used to partition expression profiles
save_network_file=T, # whether or not save network for each group
save_path="./SJARACNe/")
