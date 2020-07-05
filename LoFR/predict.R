library(randomForest)
#library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

load(args[2])

ptv_tab = read.table(args[1], sep='\t', header=T)
ptv_tab$rescue_conserved = ptv_tab$rescue_loeuf_min <= ptv_tab$min_LOEUF

predictions = predict(my_forest_b, ptv_tab, type='vote')
ptv_tab$predict_prob = predictions[, 2]
ptv_tab$LoF_confidence = ifelse(ptv_tab$predict_prob > 0.5, 'LOW',
			   ifelse(ptv_tab$predict_prob > 0.25, 'MEDIUM', 'HIGH'))

write.table(ptv_tab, file='results.tsv', sep='\t', row.names=F, quote=F)
