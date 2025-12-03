library(seneR)
env_Load() 

example_expr <- read.csv(system.file("demo", "GSE246425_counts.csv", package = "seneR"),row.names = 1)
example_expr[1:4,1:4]
example_meta <- read.csv(system.file("demo", "GSE246425_meta.csv", package = "seneR"),row.names = 1)
head(example_meta)

SID_res <- SenCID(example_expr)
table(SID_res$recSID$RecSID)
score_res <- SID_res$pred_dict$SID3

#F1b
plot_group(score_res, example_meta, 'group','SID_Score',comparisons = list(c('Old','Young')))

#F1c
gsva_res <- seneGSVA(example_expr)
plot_violin(gsva_res,example_meta,'group', adjust_fdr = TRUE, p_threshold = 0.05)
