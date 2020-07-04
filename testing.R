#library(STAR)
library(ggplot2)
#library(MCMCpack)
library(reshape2)
library(cowplot)
library(colorRamps)
library(Hmisc)
library(pROC)
library(randomForest)
mypal = matlab.like2(10000)

setwd("/media/barbitoff/DATA/Working issues/BI_Science/Natural Selection/Remastered")

########Analysis workflow################################
#
# 1. PTV sites are selected using run_ptv.py
# 2. Indel sites are removed by comparing ref and alt lengths
# 2*. CCR values are added using add_ccr.py ran on intersection of CCR and variant data
# 3. Data is aggregated per-gene
# 4. Gene-level data is inputted into the first chunk in this script
# 5. Filtering and annotation with U/ExAC_shet is done, missing genes are removed
# 6. Descriptive AC/AN statistics are collected in the second chunk
# 7. IG_params.py is ran on the full_out.tsv file
# 8. estimate_S.py is ran on the results of IG_params.py
# 9. Output of IG_params and estimate_S are used in the third and fourth chunk in this script
# 10. Output of estimate_S.py is forwarded to CPL.py to infer likelihoods
# 11. Output of CPL.py is used in the fifth chunk to generate plots.
# 12.1. ALoFT scores are calculated and added to data with add_aloft.py
# 12.2. Variant and gene-level data is merged.
# 13. Functional evidence for pLoF variants is investigated in sixth and seventh chunk.
# 14. Most common variants in each gene are selected in each gene
# 15. In shell, most common variants and exons with these variants are purged from data.
# 16. Effects of such removal are evaluated in the eighth chunk,
# 17. Constraint scores are added by add_constraint.py from ptv_with_vep.vcf
# 18. Relationships between constraint and violation are investigated in the tenth chunk.
# 19. Gene-level pLI/LOEUF are appended to gene-level data frame in the eleventh chunk.
# 20. Gene-level relationships are investigated.
# 21. Gene-level isoform counts and transcript counts are aggregated in shell from GENCODE v.19
# 22. Isoform and exon counts are added to gene data and compared.
#

#### Gene-level data aggregation and coverage filtering ########

raw_data = read.table('full_data_constr_pext.tsv', sep='\t', header=T)
cov_data = read.table('./covered_genes.tsv', header=T, sep='\t')

raw_data = raw_data[as.character(raw_data$gene) %in% as.character(cov_data$name), ]
raw_data$gene = droplevels(raw_data$gene)

ac_cols = grepl('AC_', colnames(raw_data))
an_cols = grepl('AN_', colnames(raw_data))

gene_data = data.frame(gene = sort(unique(raw_data$gene)))
for (i in colnames(raw_data)[ac_cols]){
  gene_data[, i] = as.numeric(by(raw_data[, i], raw_data$gene, sum))
}
gene_data$AC_main = gene_data$AC_nfe + gene_data$AC_afr + gene_data$AC_amr + 
  gene_data$AC_eas + gene_data$AC_sas
for (i in colnames(raw_data)[an_cols]){
  gene_data[, i] = as.numeric(by(raw_data[, i], raw_data$gene, mean))
}
gene_data$AN_main = gene_data$AN_nfe + gene_data$AN_afr + gene_data$AN_amr + 
  gene_data$AN_eas + gene_data$AN_sas

# Add covered genes with no pLoF
dummy_row = c(rep(0, 8), rep(colMeans(gene_data[, 10:17])))
zero_genes = unique(as.character(cov_data$name)[!(as.character(cov_data$name) 
                             %in% as.character(gene_data$gene))])
remainder = as.data.frame(t(sapply(zero_genes, function(x) c(x, dummy_row))))
rownames(remainder) = c()
colnames(remainder) = colnames(gene_data)

final_gene_data = as.data.frame(rbind(gene_data, remainder))
num_cols = colnames(final_gene_data)[grepl('A[CN]_', colnames(final_gene_data))]
for (col in num_cols){
  final_gene_data[, col] = as.numeric(as.character(final_gene_data[, col]))
}
final_gene_data$gene = as.factor(as.character(final_gene_data$gene))
final_gene_data = 
  final_gene_data[final_gene_data$AC_main/final_gene_data$AN_main <= 0.001, ]


# Cassa comparison
cassa = read.table('cassa_table.csv', header=T, sep='\t')
rownames(cassa) = cassa$gene_symbol
final_gene_data$U = sapply(final_gene_data$gene, 
                           function(x) cassa[as.character(x), 'U'])
final_gene_data = final_gene_data[!(is.na(final_gene_data$U)), ]
final_gene_data[is.na(final_gene_data)] = 1

perc = ecdf(final_gene_data$U)
final_gene_data$tercile = floor(perc(final_gene_data$U) * 3)
final_gene_data$tercile = ifelse(final_gene_data$tercile < 3,
                                 final_gene_data$tercile + 1,
                                 final_gene_data$tercile)
table(final_gene_data$tercile)
final_gene_data$shet_cassa = sapply(final_gene_data$gene, 
                                    function(x) cassa[as.character(x), 's_het'])

table(is.na(final_gene_data$shet_cassa))

write.table(final_gene_data, file='full_aggregated.tsv', sep='\t',
            row.names=F, quote=F)


##############Gene summary statistics#############
final_gene_data = read.table('full_aggregated.tsv', header=T, sep='\t')

tpl = melt(final_gene_data[, 2:8])
ggplot(tpl, aes(x=variable, y=value, fill=variable)) + geom_violin(adjust=5) + 
  geom_boxplot(width = 0.075, outlier.shape = NA, fill='white') +
  scale_y_continuous(limits=c(0, 25)) + theme_bw()

ans = data.frame(AN = colMeans(final_gene_data[, 10:16]), 
                 pops = colnames(final_gene_data)[10:16])
ggplot(ans, aes(x=pops, y=AN)) + geom_bar(stat='identity')

rownames(final_gene_data) = final_gene_data$gene
inheritance = read.table('genes_inheritance.txt', header=F, sep='\t', 
                         stringsAsFactors = F, row.names=1)

final_gene_data$disease = sapply(rownames(final_gene_data), 
                          function(x) ifelse(x%in% rownames(inheritance), 
                                             inheritance[x, 'V2'], 'none'))


final_gene_data$dbin = sapply(final_gene_data$disease, 
                    function(x) ifelse(x == "none", 0, 1))
final_gene_data$whole = apply(final_gene_data,
                              1,
                              function(elt)
                                chisq.test(matrix(c(as.numeric(elt[2:8]),
                                                    as.numeric(elt[10:16])),
                                                  byrow=T, nrow=2))$p.value)

final_gene_data$no_fin_asj = apply(final_gene_data,
                              1,
                              function(elt)
                                chisq.test(matrix(c(as.numeric(elt[2:6]),
                                                    as.numeric(elt[10:14])),
                                                  byrow=T, nrow=2))$p.value)

for (i in 2:6) { 
  final_gene_data[, ncol(final_gene_data) + 1] = apply(final_gene_data, 
                               1, 
                               function(elt) 
                                 chisq.test(matrix(c(as.numeric(elt[2:6][c(2:6) != i]),
                                                     as.numeric(elt[10:14][c(10:14) != 8+i])),
                                                   byrow=T, nrow=2))$p.value)
}

colnames(final_gene_data) = c(colnames(final_gene_data)[1:24], 'no_afr',
                                     'no_sas', 'no_amr', 'no_eas', 'no_nfe')
tpl = melt(final_gene_data[, c(1, 21:29)], id.vars=c('gene', 'disease', 'dbin'))
tpl = na.omit(tpl[tpl$value < 1e-5, ])

sum(tpl$variable == 'whole')
sum(tpl$variable == 'whole' & tpl$dbin == 1)

sum(tpl$variable == 'no_fin_asj')
sum(tpl$variable == 'no_fin_asj' & tpl$dbin == 1)

ggplot(tpl, aes(x=variable, fill=as.factor(dbin))) + geom_bar() +
  theme_bw() + coord_flip()

write.table(final_gene_data, file='full_out.tsv', sep='\t',
            row.names=F, quote=F)



### IG distribution fits #########################

dIG = function(x, alpha, beta) {
  fc = (beta/(2*pi*(x^3)))^(1/2)
  sc = exp(-((beta * (x - alpha)^2)/(2*alpha^2*x)))
  return(fc * sc)
}

p1 = sapply(seq(0.0001, 1, by=0.0001), 
            function(x) dIG(x, beta=0.0091, alpha=0.288))
p2 = sapply(seq(0.0001, 1, by=0.0001), 
            function(x) dIG(x, beta=0.0296, alpha=0.945))
p3 = sapply(seq(0.0001, 1, by=0.0001), 
            function(x) dIG(x, beta=0.0176, alpha=0.249))

#[(array([0.91184562, 0.00710335]), 10478.644733658062),
#(array([0.91368573, 0.00788945]), 16651.85653419709), 
#(array([0.93319117, 0.01253176]), 15636.08337229499)]

p_g1 = sapply(seq(0.0001, 1, by=0.0001), 
            function(x) dIG(x, beta=0.0103, alpha=0.089))
p_g2 = sapply(seq(0.0001, 1, by=0.0001), 
            function(x) dIG(x, beta=0.0104, alpha=0.908))
p_g3 = sapply(seq(0.0001, 1, by=0.0001), 
            function(x) dIG(x, beta=0.0193, alpha=0.079))


dist = data.frame(s = seq(0.0001, 1, 0.0001), ExAC = (p1 + p2 + p3)/3,
                  gnomAD = (p_g1 + p_g2 + p_g3)/3)

tpl = melt(dist, id.vars='s')

ggplot(tpl, aes(x=s, y=value, col=variable, fill=variable)) + geom_line(lwd=1) + scale_x_log10() + 
  theme_bw() + ylab('Probability density') + xlab('Selection coefficient')# + facet_wrap(~variable, nrow=5)



########gnomAD vs ExAC####################################

shet_gn_data = read.table('./full_corr_wshet.tsv', sep=',', header=T)
# Filtering out lowest values
# shet_gn_data = shet_gn_data[shet_gn_data$s_main > 0.0001, ]
head(shet_gn_data)

ggplot(shet_gn_data, aes(x=shet_cassa, y=s_main)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  scale_x_log10(limits=c(0.0001, 1)) + 
  scale_y_log10(limits=c(0.0001, 1)) +
  scale_fill_gradientn(colours = mypal) +
  xlab('Selection coefficient (Cassa et al., 2017)') + 
  ylab('Selection coefficient (gnomAD)')

shet_gn_data$log_s_cassa = -log10(shet_gn_data$shet_cassa)
shet_gn_data$log_s = -log10(shet_gn_data$s_main)

fit_1 <- lm(log_s ~ log_s_cassa, shet_gn_data)
summary(fit_1)

# Correlation = 0.740, CI = [0.733, 0.747]
cor(shet_gn_data$log_s, shet_gn_data$log_s_cassa)
cor.test(shet_gn_data$log_s, shet_gn_data$log_s_cassa)

s_comp = melt(shet_gn_data, id.vars='gene', 
           measure.vars = c('s_main', 'shet_cassa'))
s_comp$variable = ifelse(s_comp$variable == 's_main', 
                         'gnomAD', 'ExAC')
a <- ggplot(s_comp, aes(x=value, fill=variable)) + 
  geom_histogram(col='black') +
  theme_bw() + facet_wrap(~variable, nrow=2) +
  scale_x_log10() +
  xlab('Selection coefficient') + ylab('Count') +
  guides(fill=F)

b <- ggplot(s_comp, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(col='black', outlier.shape = NA) +
  theme_bw() + scale_y_log10() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab('Selection coefficient') + xlab('Dataset') +
  guides(fill=F)

plot_grid(a, b, ncol=2, rel_widths = c(1, 0.45),
          labels = c('a', 'b'))

aggregate(value~variable, s_comp, median)
aggregate(value~variable, s_comp, IQR)


#########Cross-population likelihood######################

#dat = read.table('out_wshet.tsv', sep=',', header=T)
#ggplot(dat, aes(x='main', y=s_main)) + geom_violin() + scale_y_log10()

cpl_head = read.table('full_out_boot_stats.csv', row.names=1, header=T, sep=',')
cpl_head = na.omit(cpl_head)
cpl_head = cpl_head[cpl_head$s_main > 0.0001 & cpl_head$L != Inf, ]

L_comp = melt(cpl_head, id.vars='gene', measure.vars = c('L', 'L_replicate'))

ggplot(L_comp, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(col='black') +
  theme_bw()# + facet_wrap(~variable, nrow=2)

# Observed vs expected
ggplot(cpl_head, aes(x=L, y=L_replicate)) + geom_hex(aes(fill=log10(..count..))) + 
  scale_y_continuous(limits=c(-5, 280)) + theme_bw() + 
  xlab("Observed likelihood") + ylab("Expected likelihood") + 
  scale_fill_gradientn(colours = mypal)

ggplot(cpl_head, aes(x=as.numeric(L))) + geom_histogram(col='black', fill='red') +
  theme_bw()

# Q-Q
quant = data.frame(observed = sort(cpl_head$L), 
                   expected = sort(cpl_head$L_replicate))
ggplot(quant, aes(x=expected, y=observed)) + geom_point() +
  theme_bw() + scale_y_continuous(limits=c(0, 160)) + scale_x_continuous(limits=c(0, 20)) +
  xlab('Expected likelihhod') + ylab('Observed likelihood') + 
  geom_abline(slope=1, col='red', lwd=1)

# Q-Q in large shet
quant = data.frame(observed = sort(cpl_head[cpl_head$s_main > 0.02, ]$L), 
                   expected = sort(cpl_head[cpl_head$s_main > 0.02, ]$L_replicate))
ggplot(quant, aes(x=expected, y=observed)) + geom_point() +
  theme_bw() +
  xlab('Expected likelihhod') + ylab('Observed likelihood') + 
  geom_abline(slope=1, col='red', lwd=1)

# CPL vs chisq
rownames(cpl_head) = cpl_head$gene
cpl_head$log_chi_p = -log10(cpl_head$without_fin_asj)
cpl_head = na.omit(cpl_head)
cpl_head = cpl_head[cpl_head$log_chi_p != Inf, ]

ggplot(cpl_head, aes(log_chi_p, L)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() +
  xlab('chi-squared p-value') + ylab('Observed likelihood') + 
  scale_fill_gradientn(colours = mypal)


#OLS (z-score transform)

cpl_head$OLS = abs(cpl_head$L_mean - cpl_head$L)/cpl_head$L_sd
ggplot(cpl_head, aes(x=OLS)) + geom_histogram(fill='red', col='black') + theme_bw()

ggplot(cpl_head, aes(log_chi_p, OLS)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() +
  xlab('chi-squared p-value') + ylab('Likelihood Z-score (OLS)') + 
  scale_fill_gradientn(colours = mypal)

cor(as.numeric(cpl_head$OLS), as.numeric(cpl_head$log_chi_p))

# CPL vs selection

ggplot(cpl_head, aes(x=s_main, y=OLS)) + geom_point() + scale_x_log10() +
  geom_hline(yintercept = 14, col='red', lwd=1) + theme_bw()

inheritance = read.table('genes_inheritance.txt', header=F, sep='\t', 
                         stringsAsFactors = F, row.names=1)

cpl_head$disease = sapply(rownames(cpl_head), 
                          function(x) ifelse(x%in% rownames(inheritance), 
                                             inheritance[x, 'V2'], 'none'))

table(cpl_head$disease)

ggplot(cpl_head, aes(x=s_main, y=OLS, col=disease)) + 
  geom_point() + scale_x_log10() +
  geom_hline(yintercept = 6, col='red', lwd=1) + theme_bw() + 
  facet_wrap(~disease, nrow=3) + guides(col=F) + xlab('Selection coefficient') + 
  ylab('Observed likelihod')

# OLS vs Shet aspect ration 5:4.5
ggplot(cpl_head, aes(x=s_main, y=OLS)) + 
  geom_hex(aes(fill=log10(..count..))) + scale_x_log10() +
  scale_fill_gradientn(colours=mypal) + 
  theme_bw() + 
  geom_vline(xintercept = 0.006, col='red', lwd=1, lty=2) +
  facet_wrap(~disease, nrow=3) + guides(col=F) + xlab('Selection coefficient') + 
  ylab('Likelihood Z-score (OLS)')

cpl_head$dbin = ifelse(cpl_head$disease %in% c('AD', 'AR'), 1, 0)

summary(glm(dbin ~ OLS, cpl_head, family="binomial"))
summary(glm(dbin ~ s_main + OLS, cpl_head, family="binomial"))

# Caret training
library(caret)
cpl_head$d_class = as.factor(ifelse(cpl_head$dbin == 1, 
                                 'implicated',
                                 'non_implicated'))

train_control <- trainControl(method="cv", number=5, classProbs=T,
                              summaryFunction=twoClassSummary)
model_s = train(d_class~s_main, data=cpl_head, trControl=train_control, 
                method="rpart")
model_z = train(d_class~OLS, data=cpl_head, trControl=train_control, 
                method="rpart")
model_s_z = train(d_class~s_main + OLS, data=cpl_head, trControl=train_control, 
                  method="rpart")

print(model_s)
print(model_s_z)
print(model_z)


write.table(cpl_head, file='final_full_gene_stats.tsv', 
            sep='\t', row.names=F, quote=F)




######Functional evidence#############################

var_data = read.table('full_data_constr_pext.tsv', sep='\t', header=T)

var_data$OLS = sapply(var_data$gene, 
                      function(x) as.numeric(cpl_head[as.character(x), 'OLS']))

var_data = var_data[!is.na(var_data$OLS), ]
var_data$OLS = as.numeric(var_data$OLS)
var_data$s_main = sapply(var_data$gene, 
                function(x) as.numeric(cpl_head[as.character(x), 's_main']))

# Exploring ALoFT score vs shet
ggplot(var_data, aes(x=s_main, y=aloft_dom)) + 
  geom_hex(aes(fill=log10(..count..))) + scale_x_log10() +
  scale_fill_gradientn(colours=mypal) + 
  theme_bw()

ggplot(var_data, aes(x=s_main, y=aloft_tol)) + 
  geom_hex(aes(fill=log10(..count..))) + scale_x_log10() +
  scale_fill_gradientn(colours=mypal) + 
  theme_bw()

cor.test(var_data$s_main, var_data$aloft_tol, use = "complete.obs")
# -0.15 very significant (expectedly)



var_data = var_data[var_data$s_main > 0.02, ]

var_data$AC_main = var_data$AC_afr + var_data$AC_amr + 
                   var_data$AC_eas + var_data$AC_sas + var_data$AC_nfe
var_data$violator = ifelse(var_data$OLS > 6, 'yes', 'no')

var_data$splice = grepl('splice', var_data$conseq)


#Confidence intervals
violators = var_data[var_data$violator == 'yes', ]
non_violators = var_data[var_data$violator == 'no', ]

# var_data$is_mean = as.factor(var_data$is_mean)
# var_data$is_max = as.factor(var_data$is_max)
# var_data$is_5pct_mean = as.factor(var_data$is_5pct_mean)


#No separation of top variants
# Raw CCR pct
ggplot(var_data, aes(x=ccr.exon_pct, y=OLS)) + 
  geom_hex(aes(fill=log10(..count..))) +
  scale_fill_gradientn(colours=mypal) + 
  theme_bw()

ggplot(var_data, aes(x=ccr.exon_pct, fill=violator)) + 
  geom_density(alpha=0.5) + theme_bw()

ggplot(var_data, aes(x=violator, y=ccr.exon_pct, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
          outlier.shape=NA)

wilcox.test(violators$ccr.exon_pct, non_violators$ccr.exon_pct)

aggregate(ccr.exon_pct ~ violator, var_data, median)
aggregate(ccr.exon_pct ~ violator, var_data, IQR)


# CCR pct gene-based Z
ggplot(var_data, aes(x=ccr_exon_z, y=OLS)) + 
  geom_hex(aes(fill=log10(..count..))) +
  scale_fill_gradientn(colours=mypal) + 
  theme_bw()

ggplot(var_data, aes(x=ccr_exon_z, fill=violator)) + 
  geom_density(alpha=0.5) + theme_bw()



# CCR pct gene-based fraction (pct max)
ggplot(var_data, aes(x=ccr_exon_frac, y=OLS)) + 
  geom_hex(aes(fill=log10(..count..))) +
  scale_fill_gradientn(colours=mypal) + 
  theme_bw() + 
  scale_y_continuous(limits=c(0, 15))

ggplot(var_data, aes(x=ccr_exon_frac, fill=violator)) + 
  geom_density(alpha=0.5) + theme_bw()

a <- ggplot(var_data, aes(x=violator, y=ccr_exon_frac, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
               outlier.shape=NA) +
  guides(fill=F) + ylab('Exon constraint (%)')

ggplot(var_data, aes(x=violator, y=ccr_exon_frac, fill=violator)) + 
  theme_bw() + geom_boxplot(col='black', lwd=0.5, outlier.shape = NA)


pdf('CCR_pct_max.pdf', width=1.8, height=2.4)
print(a)
dev.off()

wilcox.test(violators$ccr_exon_frac, non_violators$ccr_exon_frac)



top_cons = as.data.frame(rbind(binconf(sum(na.omit(violators)$ccr_exon_frac == 1), 
                            nrow(na.omit(violators))),
                    binconf(sum(na.omit(non_violators)$ccr_exon_frac == 1), 
                            nrow(na.omit(non_violators)))))
top_cons$violator = factor(c('yes', 'no'))
b_1 <- ggplot(top_cons, aes(x=violator, y=PointEst, 
                  fill=violator)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(width=0.4, lwd=0.5, aes(ymin=Lower, ymax=Upper)) + 
  theme_bw() + guides(fill=F) +
  ylab('pLoF in most\nconserved exon, %')

#Expression

ggplot(var_data, aes(x=violator, fill=is_mean)) +
  geom_bar(position = "fill")

ggplot(var_data, aes(x=violator, fill=is_max)) +
  geom_bar(position = "fill")

# ggplot(var_data, aes(x=violator, fill=is_5pct_mean)) +
#   geom_bar(position = "fill")


means = as.data.frame(rbind(binconf(sum(violators$is_mean), nrow(violators)),
              binconf(sum(non_violators$is_mean), nrow(non_violators))))
means$violator = c('yes', 'no')
b_2 <- ggplot(means, aes(x=violator, y=PointEst, 
                  fill=violator)) +
  geom_bar(col='black', stat='identity')  + 
  geom_errorbar(width=0.4, lwd=0.5, aes(ymin=Lower, ymax=Upper)) + 
  theme_bw() + guides(fill=F) +
  ylab('pLoF in isoform\nwith high TPM, %')


maxes = as.data.frame(rbind(binconf(sum(violators$is_max), nrow(violators)),
                                binconf(sum(non_violators$is_max), nrow(non_violators))))
maxes$violator = c('yes', 'no')
b_3 <- ggplot(maxes, aes(x=violator, y=PointEst, 
                  fill=violator)) +
  geom_bar(col='black', stat='identity')  + 
  geom_errorbar(width=0.4, lwd=0.5, aes(ymin=Lower, ymax=Upper)) + 
  theme_bw() + guides(fill=F) +
  ylab('pLoF in isoform\nwith max. TPM, %')

pdf('func_evidence.pdf', width=6.3, height=1.85)
plot_grid(a, b_1, b_2, b_3, nrow=1)
dev.off()

# binconf(nrow(violators) - sum(violators$is_5pct_mean), nrow(violators))
# binconf(nrow(non_violators) - sum(non_violators$is_max), nrow(non_violators))

binconf(sum(violators$splice), nrow(violators))
binconf(sum(non_violators$splice), nrow(non_violators))

#ALoFT scores

ad <- ggplot(var_data, aes(x=violator, y=aloft_dom, fill=violator)) +
  geom_boxplot(outlier.shape=NA) +
  theme_bw() + guides(fill=F) +
  scale_y_continuous(limits = c(0, 1))

# Only recessive shows difference (strangely), dominant - in the opposite direction
ar <- ggplot(var_data, aes(x=violator, y=aloft_rec, fill=violator)) +
  geom_boxplot(outlier.shape=NA) +
  theme_bw() + guides(fill=F) +
  scale_y_continuous(limits = c(0, 1))

at <- ggplot(var_data, aes(x=violator, y=aloft_tol, fill=violator)) +
  geom_boxplot(outlier.shape=NA) +
  theme_bw() + guides(fill=F) +
  scale_y_continuous(limits = c(0, 1))

plot_grid(ad, ar, at, nrow=1)

table(violators$aloft_tol > 0.25)
table(non_violators$aloft_tol > 0.25)


# Newer constraint stats

rescuer = as.data.frame(rbind(binconf(sum(na.omit(violators)$rescue == "True"), 
                                       nrow(na.omit(violators))),
                               binconf(sum(na.omit(non_violators)$rescue == "True"), 
                                       nrow(na.omit(non_violators)))))
rescuer$violator = factor(c('yes', 'no'))
r_1 <- ggplot(rescuer, aes(x=violator, y=PointEst, 
                            fill=violator)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(width=0.4, lwd=0.5, aes(ymin=Lower, ymax=Upper)) + 
  theme_bw() + guides(fill=F) +
  ylab('% variants with a rescue isoform')
print(r_1)


r_2 <- ggplot(var_data, aes(x=violator, y=rescue_loeuf_min, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
               outlier.shape=NA) +
  guides(fill=F) + ylab('Min. rescue LOEUF')
print(r_2)

resc_cons = as.data.frame(rbind(binconf(sum(na.omit(violators)$rescue_loeuf_min < 1), 
                                       nrow(na.omit(violators))),
                               binconf(sum(na.omit(non_violators)$rescue_loeuf_min < 1), 
                                       nrow(na.omit(non_violators)))))
resc_cons$violator = factor(c('yes', 'no'))
r_3 <- ggplot(resc_cons, aes(x=violator, y=PointEst, 
                            fill=violator)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(width=0.4, lwd=0.5, aes(ymin=Lower, ymax=Upper)) + 
  theme_bw() + guides(fill=F) +
  ylab('Rescue isoform conserved, %')
print(r_3)

r_4 <- ggplot(var_data, aes(x=violator, y=rescue_oe_mis, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
               outlier.shape=NA) +
  guides(fill=F) + ylab('Min. rescue LOEUF')
print(r_4)


# ML
var_data$violator = as.factor(var_data$violator)
var_data$loeuf_diff = as.numeric(var_data$rescue_loeuf_min) - 
  as.numeric(var_data$min_LOEUF)
var_data$loeuf_diff[var_data$loeuf_diff < 0] = 0
ggplot(var_data, aes(x=violator, y=loeuf_diff, fill=violator)) + 
  geom_violin(lwd=0.5, alpha=0.5) + 
  geom_boxplot(lwd=0.5, width=0.4, fill='white') + 
  theme_bw() + ylab('unaff./aff. LOEUF')
var_data$diff_small = var_data$loeuf_diff < 0.25
aggregate(diff_small ~ violator, var_data, mean)

var_data$rescue_conserved = var_data$rescue_loeuf_min < 0.35 | 
                        (var_data$rescue_loeuf_min <= var_data$min_LOEUF)
all_vars$rescue_conserved = all_vars$rescue_loeuf_min < 0.35 | 
  (all_vars$rescue_loeuf_min <= all_vars$min_LOEUF)

aggregate(rescue_conserved~violator, var_data, mean)
var_data$rescue_logical = var_data$rescue == 'True'
all_vars$rescue_logical = all_vars$rescue == 'True'
aggregate(rescue_logical~violator, var_data, mean, na.action = "na.omit")


train_ind = sample(1:nrow(var_data), 0.8*nrow(var_data))
trainset = var_data[train_ind, ]
include = sapply(as.character(trainset$violator),
                 function(x) ifelse(x == "yes", TRUE,
                                    runif(1) < 0.015))
trainset_balanced = trainset[include, ]
testset = var_data[-train_ind, ]

my_fit <- glm(violator~ccr_exon_frac + is_mean + 
               is_max + rescue_loeuf_min + pext_avg, 
             trainset, family="binomial", na.action = "na.exclude")
summary(my_fit)

testset$prediction = predict(my_fit, testset)
table(testset$prediction)

pROC_obj <- roc(testset$violator, testset$prediction,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")
#plot(sens.ci, type="bars")
#dev.off()

violatoreg_res = data.frame(prob=pROC_obj$thresholds,
                           sens=pROC_obj$sensitivities,
                           spec=pROC_obj$specificities)

testset$glm_class = testset$prediction > 4.506798



# Trying random forest classifier
my_forest = randomForest(violator~ccr_exon_frac + is_mean + 
                           is_max + pext_avg + loeuf_norm, 
                         trainset, na.action = NULL)
prediction_forest = predict(my_forest, testset, type="vote")
testset$prediction_forest = prediction_forest[, 2]

pROC_obj <- roc(testset$violator, testset$prediction_forest,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj, specificities = pROC_obj$specificities,
                 sensitivities = pROC_obj$sensitivities)
plot(sens.ci, type="shape", col="lightblue")
#plot(sens.ci, type="bars")
#dev.off()

violatoRF_res = data.frame(prob=pROC_obj$thresholds,
                           sens=pROC_obj$sensitivities,
                           spec=pROC_obj$specificities)

testset$forest_class = testset$prediction_forest > 0.003

# Forest with balanced classes

my_forest_b = randomForest(violator~ccr_exon_frac + is_mean + 
                           is_max + loeuf_norm + pext_avg, 
                         trainset_balanced, na.action = NULL)
prediction_forest_b = predict(my_forest_b, testset, type="vote")
testset$prediction_forest_b = prediction_forest_b[, 2]

pROC_obj <- roc(testset$violator, testset$prediction_forest_b,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj, specificities = pROC_obj$specificities,
                 sensitivities = pROC_obj$sensitivities)
plot(sens.ci, type="shape", col="lightblue")
#plot(sens.ci, type="bars")
#dev.off()

violatoRF_b_res = data.frame(prob=pROC_obj$thresholds,
                           sens=pROC_obj$sensitivities,
                           spec=pROC_obj$specificities)

testset$forest_b_class = testset$prediction_forest > 0.003






#####Leave splice out####################################
# This section is to proof that FE is not only for splice variants


var_data = var_data[var_data$splice == F, ]
violators = var_data[var_data$violator == 'yes', ]
non_violators = var_data[var_data$violator == 'no', ]


a <- ggplot(var_data, aes(x=violator, y=ccr_exon_frac, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
               outlier.shape=NA) +
  guides(fill=F) + ylab('Exon constraint (%)')

print(a)


top_cons = as.data.frame(rbind(binconf(sum(na.omit(violators)$ccr_exon_frac == 1), 
                                       nrow(na.omit(violators))),
                               binconf(sum(na.omit(non_violators)$ccr_exon_frac == 1), 
                                       nrow(na.omit(non_violators)))))
top_cons$violator = c('yes', 'no')
b_1 <- ggplot(top_cons, aes(x=violator, y=PointEst, 
                            ymin=Lower, ymax=Upper, fill=violator)) +
  geom_bar(col='black', stat='identity') + geom_errorbar(width=0.4, lwd=0.5) + 
  theme_bw() + guides(fill=F) +
  ylab('pLoF in most\nconserved exon, %')

means = as.data.frame(rbind(binconf(sum(violators$is_mean), nrow(violators)),
                            binconf(sum(non_violators$is_mean), nrow(non_violators))))
means$violator = c('yes', 'no')
b_2 <- ggplot(means, aes(x=violator, y=PointEst, 
                         ymin=Lower, ymax=Upper, fill=violator)) +
  geom_bar(col='black', stat='identity') + geom_errorbar(width=0.4, lwd=0.5) + 
  theme_bw() + guides(fill=F) +
  ylab('pLoF in isoform\nwith high TPM, %')


maxes = as.data.frame(rbind(binconf(sum(violators$is_max), nrow(violators)),
                            binconf(sum(non_violators$is_max), nrow(non_violators))))
maxes$violator = c('yes', 'no')
b_3 <- ggplot(maxes, aes(x=violator, y=PointEst, 
                         ymin=Lower, ymax=Upper, fill=violator)) +
  geom_bar(col='black', stat='identity') + geom_errorbar(width=0.4, lwd=0.5) + 
  theme_bw() + guides(fill=F) +
  ylab('pLoF in isoform\nwith max. TPM, %')

plot_grid(a, b_1, b_2, b_3, nrow=1)




###############Top variant removed########################

generate_plots <- function(threshold, filename) {
  var_data = read.table(filename, sep='\t', header=T)
  
  var_data$OLS = sapply(var_data$gene, 
                        function(x) as.numeric(cpl_head[as.character(x), 'OLS']))
  
  var_data = var_data[!is.na(var_data$OLS), ]
  var_data$OLS = as.numeric(var_data$OLS)
  var_data$s_main = sapply(var_data$gene, 
                           function(x) as.numeric(cpl_head[as.character(x), 's_main']))
  
  
  var_data = var_data[var_data$s_main > threshold, ]
  
  var_data$AC_main = var_data$AC_afr + var_data$AC_amr + 
    var_data$AC_eas + var_data$AC_sas + var_data$AC_nfe
  var_data$violator = ifelse(var_data$OLS > 6, 'yes', 'no')
  
  var_data$splice = grepl('splice', var_data$conseq)
  
  violators = var_data[var_data$violator == 'yes', ]
  non_violators = var_data[var_data$violator == 'no', ]
  
  
  a <- ggplot(var_data, aes(x=violator, y=ccr_exon_frac, fill=violator)) + 
    geom_violin(alpha=0.5) + theme_bw() + 
    geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
                 outlier.shape=NA) +
    guides(fill=F) + ylab('Exon constraint (%)')
  
  print(a)
  
  
  top_cons = as.data.frame(rbind(binconf(sum(na.omit(violators)$ccr_exon_frac == 1), 
                                         nrow(na.omit(violators))),
                                 binconf(sum(na.omit(non_violators)$ccr_exon_frac == 1), 
                                         nrow(na.omit(non_violators)))))
  top_cons$violator = c('yes', 'no')
  b_1 <- ggplot(top_cons, aes(x=violator, y=PointEst, 
                              ymin=Lower, ymax=Upper, fill=violator)) +
    geom_bar(col='black', stat='identity') + geom_errorbar(width=0.4, lwd=0.5) + 
    theme_bw() + guides(fill=F) +
    ylab('pLoF in most\nconserved exon, %')
  
  means = as.data.frame(rbind(binconf(sum(violators$is_mean), nrow(violators)),
                              binconf(sum(non_violators$is_mean), nrow(non_violators))))
  means$violator = c('yes', 'no')
  b_2 <- ggplot(means, aes(x=violator, y=PointEst, 
                           ymin=Lower, ymax=Upper, fill=violator)) +
    geom_bar(col='black', stat='identity') + geom_errorbar(width=0.4, lwd=0.5) + 
    theme_bw() + guides(fill=F) +
    ylab('pLoF in isoform\nwith high TPM, %')
  
  
  maxes = as.data.frame(rbind(binconf(sum(violators$is_max), nrow(violators)),
                              binconf(sum(non_violators$is_max), nrow(non_violators))))
  maxes$violator = c('yes', 'no')
  b_3 <- ggplot(maxes, aes(x=violator, y=PointEst, 
                           ymin=Lower, ymax=Upper, fill=violator)) +
    geom_bar(col='black', stat='identity') + geom_errorbar(width=0.4, lwd=0.5) + 
    theme_bw() + guides(fill=F) +
    ylab('pLoF in isoform\nwith max. TPM, %')
  
  pls <- plot_grid(a, b_1, b_2, b_3, nrow=1)
  return(pls)
}

threshold = 0.02

print(generate_plots(0.02, 'selected_top1_vars.tsv'))

nr = generate_plots(threshold, 'full_data_aloft.tsv')
print(nr)

vr = generate_plots(threshold, 'full_data_aloft_varremoved.tsv')
print(vr)

er = generate_plots(threshold, 'full_data_aloft_exonremoved.tsv')
print(er)

plot_grid(nr, vr, er, nrow=3)


############Functional evidence in genes with different shet######

x = lapply(c(0.001, 0.006, 0.02), function(x) generate_plots(x, 'full_data_aloft.tsv'))

plot_grid(x[[1]], x[[2]], x[[3]], nrow=3)


############pLI scores############################################

var_data = read.table('full_data_pLI.tsv', sep='\t', header=T)

var_data$OLS = sapply(var_data$gene, 
                      function(x) as.numeric(cpl_head[as.character(x), 'OLS']))

var_data = var_data[!is.na(var_data$OLS), ]
var_data$OLS = as.numeric(var_data$OLS)
var_data$s_main = sapply(var_data$gene, 
                         function(x) as.numeric(cpl_head[as.character(x), 's_main']))

var_data$AC_main = var_data$AC_afr + var_data$AC_amr + 
  var_data$AC_eas + var_data$AC_sas + var_data$AC_nfe
var_data$violator = ifelse(var_data$OLS > 6, 'yes', 'no')

var_data$splice = grepl('splice', var_data$conseq)

# s vs pLI

#var_data = read.table('comprehensive_variant_annotation.tsv', header=T, sep='\t')

loeuf_s <- ggplot(var_data, aes(x=s_main, y=min_LOEUF, col=violator)) + 
  geom_smooth(lwd=1, method="glm", fullrange=T, 
              method.args = list(family = "binomial")) +
  theme_bw() + scale_x_log10()

pLI_s <- ggplot(var_data, aes(x=s_main, y=max_pLI, col=violator)) + 
  geom_smooth(lwd=1, method="glm", fullrange=T,
              method.args = list(family = "binomial")) +
  theme_bw() + scale_x_log10()

plot_grid(loeuf_s, pLI_s, nrow=2)

loeuf_hex <- ggplot(var_data, aes(x=s_main, y=min_LOEUF)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() + scale_fill_gradientn(colours = mypal) + 
  facet_wrap(~violator, nrow=2) +
  scale_x_log10() + guides(fill=F)

pli_hex <- ggplot(var_data, aes(x=s_main, y=max_pLI)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() + scale_fill_gradientn(colours = mypal) + 
  facet_wrap(~violator, nrow=2) +
  scale_x_log10()

plot_grid(loeuf_hex, pli_hex, nrow=1, rel_widths = c(0.75, 1))

# Selection and split
var_data = var_data[var_data$s_main > 0.02, ]
violators = var_data[var_data$violator == 'yes', ]
non_violators = var_data[var_data$violator == 'no', ]

#Confidence intervals
#violators = var_data[var_data$violator == 'yes' & var_data$s_main > 0.006, ]
#non_violators = var_data[var_data$violator == 'no' & var_data$s_main > 0.02, ]


exc <- ggplot(var_data, aes(x=violator, y=ccr_exon_frac, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
               outlier.shape=NA) +
  guides(fill=F) + ylab('Exon constraint (%)')

print(exc)

pli <- ggplot(var_data, aes(x=violator, y=max_pLI, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
               outlier.shape=NA) +
  guides(fill=F) + ylab('Maximum trascript pLI')

print(pli)

loeuf <- ggplot(var_data, aes(x=violator, y=min_LOEUF, fill=violator)) + 
  geom_violin(alpha=0.5) + theme_bw() + 
  geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
               outlier.shape=NA) +
  guides(fill=F) + ylab('Minimum transcript LOEUF')

print(loeuf)
plot_grid(exc, pli, loeuf, nrow=1)


low_cons_tr = as.data.frame(rbind(binconf(sum(na.omit(violators)$min_LOEUF > 0.75), 
                                       nrow(na.omit(violators))),
                               binconf(sum(na.omit(non_violators)$min_LOEUF > 0.75), 
                                       nrow(na.omit(non_violators)))))
low_cons_tr$violator = c('yes', 'no')
oe_1 <- ggplot(low_cons_tr, aes(x=violator, y=PointEst, 
                            ymin=Lower, ymax=Upper, fill=violator)) +
  geom_bar(col='black', stat='identity') + geom_errorbar(width=0.4, lwd=0.5) + 
  theme_bw() + guides(fill=F) +
  ylab('pLoF in poorly\nconserved transcript, %')
print(oe_1)

write.table(var_data, 'comprehensive_variant_annotation.tsv', sep='\t', row.names=F, quote=F)


generate_plots_pli <- function(threshold, filename) {
  var_data = read.table(filename, sep='\t', header=T)
  
  var_data$OLS = sapply(var_data$gene, 
                        function(x) as.numeric(cpl_head[as.character(x), 'OLS']))
  
  var_data = var_data[!is.na(var_data$OLS), ]
  var_data$OLS = as.numeric(var_data$OLS)
  var_data$s_main = sapply(var_data$gene, 
                           function(x) as.numeric(cpl_head[as.character(x), 's_main']))
  
  
  var_data = var_data[var_data$s_main > threshold, ]
  
  var_data$AC_main = var_data$AC_afr + var_data$AC_amr + 
    var_data$AC_eas + var_data$AC_sas + var_data$AC_nfe
  var_data$violator = ifelse(var_data$OLS > 6, 'yes', 'no')
  
  var_data$splice = grepl('splice', var_data$conseq)
  
  exc <- ggplot(var_data, aes(x=violator, y=ccr_exon_frac, fill=violator)) + 
    geom_violin(alpha=0.5) + theme_bw() + 
    geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
                 outlier.shape=NA) +
    guides(fill=F) + ylab('Exon constraint (%)')
  
  pli <- ggplot(var_data, aes(x=violator, y=max_pLI, fill=violator)) + 
    geom_violin(alpha=0.5) + theme_bw() + 
    geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
                 outlier.shape=NA) +
    guides(fill=F) + ylab('Maximum trascript pLI')
  
  loeuf <- ggplot(var_data, aes(x=violator, y=min_LOEUF, fill=violator)) + 
    geom_violin(alpha=0.5) + theme_bw() + 
    geom_boxplot(width=0.4, col='black', lwd=0.5, fill='white',
                 outlier.shape=NA) +
    guides(fill=F) + ylab('Minimum LOEUF')
  
  pls <- plot_grid(exc, pli, loeuf, nrow=1)
  return(pls)
}

x = lapply(c(0.001, 0.006, 0.02), function(x) generate_plots_pli(x, 'full_data_pLI.tsv'))

plot_grid(x[[1]], x[[2]], x[[3]], nrow=3)

#########pLI vs shet on gene level#################################

gene_cons = read.table('gnomad.v2.1.1.lof_metrics.by_gene.txt', header=T, sep='\t')

cpl_head$pLI = as.numeric(sapply(cpl_head$gene, function(x)
  as.numeric(gene_cons[gene_cons$gene == as.character(x), 'pLI'][1])))

cpl_head$loeuf = as.numeric(sapply(cpl_head$gene, function(x)
  as.numeric(gene_cons[gene_cons$gene == as.character(x), 'oe_lof_upper'][1])))
cpl_head$loeuf = ifelse(cpl_head$loeuf > 1, 1, cpl_head$loeuf)

cpl_head$violator = ifelse(cpl_head$OLS > 5, 'yes', 'no')

#cpl_head = read.table('comprehensive_gene_annotation.tsv', sep='\t', header=T)

loeuf_s <- ggplot(cpl_head, aes(x=s_main, y=loeuf, col=violator)) + 
  geom_smooth(lwd=1, method="glm", fullrange=T, 
              method.args = list(family = "binomial")) +
  theme_bw() + scale_x_log10() +
  xlab('Selection coefficient') +
  ylab('LOEUF')

pLI_s <- ggplot(cpl_head, aes(x=s_main, y=pLI, col=violator)) + 
  geom_smooth(lwd=1, method="glm", fullrange=T,
              method.args = list(family = "binomial")) +
  theme_bw() + scale_x_log10() +
  xlab('Selection coefficient') +
  ylab('pLI')

plot_grid(loeuf_s, pLI_s, nrow=2)

write.table(cpl_head, 'comprehensive_gene_annotation.tsv', 
            sep='\t', row.names=F, quote=F)

ggplot(cpl_head, aes(x=loeuf, y=pLI)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() + scale_fill_gradientn(colours = mypal)

loeuf <- ggplot(cpl_head, aes(x=s_main, y=loeuf)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() + scale_x_log10() +
  scale_fill_gradientn(colours = mypal) +
  facet_wrap(~violator, nrow=2)

pli <- ggplot(cpl_head, aes(x=s_main, y=pLI)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() + scale_x_log10() +
  scale_fill_gradientn(colours = mypal) +
  facet_wrap(~violator, nrow=2)

plot_grid(pli, loeuf)


############### Investigation of isoform counts ############

gene_data = read.table('./comprehensive_gene_annotation.tsv', sep='\t',
                       header=T)

iso_data = read.table('merged_counts.tsv', sep='\t', header=F)
colnames(iso_data) = c('isoforms', 'gene', 'exons')
head(iso_data)
ggplot(iso_data, aes(x=exons, y=isoforms)) + geom_hex() + 
  scale_x_continuous(limits=c(0,50))
cor(iso_data$isoforms, iso_data$exons, use='complete.obs')
cor.test(iso_data$isoforms, iso_data$exons, use='complete.obs')
rownames(iso_data) = iso_data$gene

gene_data$isoforms = sapply(as.character(gene_data$gene), 
                            function(x) iso_data[x, 'isoforms'])

gene_data$exons = sapply(as.character(gene_data$gene), 
                            function(x) iso_data[x, 'exons'])

iso <- ggplot(gene_data[gene_data$s_main > 0.02, ], 
              aes(x=violator, y=isoforms, fill=violator)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(outlier.shape=NA, width=0.3, fill='white') + theme_bw() + 
  scale_y_continuous(limits=c(0, 30)) + guides(fill=F)

aggregate(isoforms~violator, gene_data[gene_data$s_main > 0.02, ], mean)

exo <- ggplot(gene_data[gene_data$s_main > 0.02, ], 
              aes(x=violator, y=exons, fill=violator)) + 
  geom_violin(alpha=0.5) +
  geom_boxplot(outlier.shape=NA, width=0.3, fill='white') + theme_bw() + 
  scale_y_continuous(limits=c(0, 75)) + guides(fill=F)

aggregate(exons~violator, gene_data[gene_data$s_main > 0.02, ], mean)

plot_grid(iso, exo, nrow=1)


########### Variant pexts ##############################################

var_pext = read.table('vars_with_pexts.tsv', sep='\t', header=T)
var_pext_sub = var_pext[var_pext$s_main > 0.02, ]
pext <- ggplot(var_pext_sub, aes(x=violator, y=pext_avg, fill=violator)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.5) + theme_bw() + guides(fill=F)
print(pext)

plot_grid(iso, exo, pext, nrow=1)

var_pext_sub$low_pext = var_pext_sub$pext_avg < 0.1
aggregate(low_pext ~ violator, var_pext_sub,  mean)
aggregate(low_pext ~ violator, var_pext_sub,  length)
aggregate(low_pext ~ violator, var_pext_sub,  sum)
wilcox.test(pext_avg~violator, var_pext_sub)

flt_counts = matrix(c(23318, 972, 252, 26), nrow=2)
chisq.test(flt_counts)

var_pext_sub$exons = as.numeric(sapply(as.character(var_pext_sub$gene),
                            function(x) iso_data[x, 'exons']))

ggplot(var_pext_sub, aes(x=exons, y=pext_avg, col=violator)) + geom_point() +
  facet_wrap(~violator, nrow=2)
cor.test(as.numeric(var_pext_sub$exons),
    as.numeric(var_pext_sub$pext_avg),
    use='complete.obs')

pext_vs_exon <- lm(pext_avg~exons+violator, var_pext_sub)
summary(pext_vs_exon)


# Filter pext > 0.1 and look into the changes

var_expressed = var_pext[var_pext$pext_avg >= 0.1, ]
gene_pext_flt <- aggregate(AC_main~gene, var_expressed, sum)
rownames(gene_pext_flt) = gene_pext_flt$gene
gene_data$pext_AC = sapply(as.character(gene_data$gene),
                           function(x) gene_pext_flt[x, 'AC_main'])

vars_non_true = all_vars[!all_vars$pred_class, ]
gene_level = aggregate(AC_main ~ gene, vars_non_true, sum)
gene_level$source_ac = sapply(as.character(gene_level$gene), function(x) cpl_head[x, 'AC_main'])

ggplot(gene_level, aes(source_ac, AC_main)) + 
  geom_hex(aes(fill=log10(..count..))) + 
  theme_bw() + scale_fill_gradientn(colours = mypal)



########### Here we will do some testing with predictions #######

prediction_forest_all = predict(my_forest, all_vars, type="vote")
all_vars$prediction_prob = prediction_forest_all[, 2]
all_vars$pred_class = all_vars$prediction_prob > 0.020
all_vars$violator_class = all_vars$violator == 'yes'
av_pred = aggregate(pred_class~violator, all_vars, 
                    function(x) binconf(sum(x), length(x)))

prediction_forest_cons = predict(my_forest, var_data, type="vote")
var_data$prediction_prob = prediction_forest_cons[, 2]
var_data$pred_class = var_data$prediction_prob > 0.020
cv_pred = aggregate(pred_class~violator, var_data, 
                    function(x) binconf(sum(x), length(x)))

pred_all = as.data.frame(as.matrix(rbind(av_pred, cv_pred)))
colnames(pred_all) = c('violator', 'PointEst', 'Lower', 'Upper')
pred_all$PointEst = as.numeric(as.character(pred_all$PointEst))
pred_all$Lower = as.numeric(as.character(pred_all$Lower))
pred_all$Upper = as.numeric(as.character(pred_all$Upper))
pred_all$gene_group = rep(c('All', 'Constrained'), each=2)

ggplot(pred_all, aes(x=violator, y=PointEst, fill=violator)) + 
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), lwd=0.5, width=0.4) +
  theme_bw() + guides(fill=F) + xlab('Model violation') +
  ylab('% predicted positive') + facet_wrap(~gene_group, nrow=1)
