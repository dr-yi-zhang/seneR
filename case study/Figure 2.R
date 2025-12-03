library(seneR)
env_Load() 

t2d_counts <- read.csv('t2d/RNAseq_counts.csv',row.names = 1)
meta <- read.csv('t2d/donor.csv')

t2d_sene <- SenCID(t2d_counts,binarize = T)
table(t2d_sene$recSID$RecSID)
sene_meta <- t2d_sene$pred_dict$SID4

sene_meta$diagnosis <- meta$diagnosis[match(rownames(sene_meta),meta$record_id)]
sene_meta$diagnosis_computed <- meta$diagnosis_computed[match(rownames(sene_meta),meta$record_id)]

#F2a
plot_group(sene_meta[,1:3], sene_meta[,4:5], 'diagnosis','SID_Score')

#F2d, F2e
gsva_res <- seneGSVA(t2d_counts)
plot_violin(gsva_res,meta,'age_group')
plot_violin(gsva_res,meta,'SID_group')
plot_violin(gsva_res,meta,'diagnosis')