library(randomForest)
#library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

ptv_tab = read.table(args[1], sep='\t', header=T)
ptv_tab$rescue_conserved = ptv_tab$rescue_loeuf_min <= ptv_tab$min_LOEUF

load(args[2])
load(args[3])

predictions_A = predict(my_forest_c, ptv_tab, type='vote')
predictions_B = predict(my_forest_d, ptv_tab, type='vote')

ptv_tab$Model_GN_A_prob = predictions_A[, 2]
ptv_tab$GN_A_LoF_confidence = ifelse(ptv_tab$Model_GN_A_prob > 0.575, 'LOW',  
				     ifelse(ptv_tab$Model_GN_A_prob > 0.25, 'MEDIUM', 'HIGH'))

ptv_tab$Model_CLV_prob = predictions_B[, 2]
ptv_tab$CLV_LoF_confidence = ifelse(ptv_tab$Model_CLV_prob > 0.59, 'LOW',
                           ifelse(ptv_tab$Model_CLV_prob > 0.25, 'MEDIUM', 'HIGH'))

write.table(ptv_tab, file='results.tsv', sep='\t', row.names=F, quote=F)
