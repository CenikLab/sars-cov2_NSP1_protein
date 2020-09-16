library(ribor)
library(ggplot2)
library(reshape2)
library(edgeR)
library(erccdashboard)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(EnhancedVolcano)
library(MatchIt)

color.palette0 = colorRampPalette(c("#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"), space="Lab")

strip_extension = function(genenames) {
  return(sapply ( strsplit( genenames, split = "-") , "[[", 1 ) )
}

nsp.ribo = Ribo('./all.ribo', rename = rename_default)
                        )
final_set = c("20200717-WT-HEK-1-setB","20200717-WT-HEK-2-setB","20200717-WT-HEK-3-setB",
              "20200717-NSP1-HEK-1-setB","20200717-NSP1-HEK-2-setB","20200717-NSP1-HEK-3-setB",
              "20200717-NSP2-HEK-1-setB", "20200717-NSP2-HEK-2-setB", "20200717-NSP2-HEK-3-setB")

final_wt = c("20200717-WT-HEK-1-setB","20200717-WT-HEK-2-setB","20200717-WT-HEK-3-setB")
             
final_nsp1 = c("20200717-NSP1-HEK-1-setB","20200717-NSP1-HEK-2-setB","20200717-NSP1-HEK-3-setB")
              
final_nsp2 = c("20200717-NSP2-HEK-1-setB", "20200717-NSP2-HEK-2-setB", "20200717-NSP2-HEK-3-setB")

rc <- get_region_counts(nsp.ribo,
                        range.lower = 28,
                        range.upper = 35,
                        experiment = final_set,
                        length      = TRUE,
                        transcript  = FALSE,
                        tidy = F,
                        alias       = TRUE,
                        region      = c("CDS"), 
                        compact = F)

rc_5 <- get_region_counts(nsp.ribo,
                          range.lower = 28,
                          range.upper = 35,
                          experiment = final_set,
                          length      = TRUE,
                          transcript  = FALSE,
                          tidy = F,
                          alias       = TRUE,
                          region      = c("UTR5"),
                          compact = F)

rc_3 <- get_region_counts(nsp.ribo,
                          range.lower = 28,
                          range.upper = 35,
                          experiment = final_set,
                          length      = TRUE,
                          transcript  = FALSE,
                          tidy = F,
                          alias       = TRUE,
                          region = c("UTR3"), 
                          compact = F)

rnaseq <- get_rnaseq(ribo.object = nsp.ribo,
                     # experiment = final_set,
                     tidy        = F,
                     compact = F, 
                     region = "CDS",
                     alias       = TRUE)



rcw = dcast(rc, transcript ~ experiment)
rc_3w = dcast(rc_3, transcript ~ experiment)
rc_5w = dcast(rc_5, transcript ~ experiment)
rnaseq_w = dcast(rnaseq, transcript ~ experiment)

region_lengths = get_region_lengths(nsp.ribo, alias = T)
rcw_wlength = merge(rcw, region_lengths, by = "transcript")

ercc <- read.csv('./ERCC_combined_counts.csv')
rnaseq_w_order = c(1, 11, 6, 12, 16, 5, 4, 10, 9, 13, 14, 8, 3, 15, 17, 7)
ercc = ercc[, rnaseq_w_order]
colnames(ercc) = colnames(rnaseq_w)
ercc = rbind(ercc, rnaseq_w)

plot_length_distribution(x           = nsp.ribo,
                         region      = "CDS",
                         range.lower = 15,
                         range.upper = 40,
                         fraction    = TRUE)

## Split length distribution by Experiment

q1 = plot_length_distribution(x           = nsp.ribo,
                         experiment = final_wt, 
                         title = "WT",
                         region      = "CDS",
                         range.lower = 15,
                         range.upper = 40,
                         fraction    = TRUE)

q2 = plot_length_distribution(x           = nsp.ribo,
                         experiment = final_nsp1, 
                         title = "NSP1",
                         region      = "CDS",
                         range.lower = 15,
                         range.upper = 40,
                         fraction    = TRUE)

q3= plot_length_distribution(x           = nsp.ribo,
                         experiment = final_nsp2,
                         title = "NSP2",
                         region      = "CDS",
                         range.lower = 15,
                         range.upper = 40,
                         fraction    = TRUE)
pdf("./Figures/Length_Distribution_v2.pdf", height = 6, width = 6)
plot_grid(q1 + ylim(0,0.125), 
          q2 + ylim(0,0.125),
          q3 + ylim(0,0.125),
          align = "vh", labels = "AUTO", ncol = 1)
dev.off()

for(seq_len in 15:40) {
  p1 = plot_region_counts(x           = nsp.ribo,
                          range.lower = seq_len,
                          range.upper = seq_len, 
                          experiment = initial_experiments,
                          title= seq_len)
  print (p1)
}

for(seq_len in 15:40) {
  p1 = plot_region_counts(x           = nsp.ribo,
                          range.lower = seq_len,
                          range.upper = seq_len, 
                          experiment = final_set,
                          title= seq_len)
  print (p1)
}

plot_region_counts(x           = nsp.ribo,
                   experiment = final_set,
                   range.lower = 28,
                   range.upper = 35)

plot_metagene(nsp.ribo,
              site        = "stop",
              normalize   = TRUE,
              experiment = final_set,
              title       = "Stop Site Coverage",
              range.lower = 28,
              range.upper = 35)

plot_metagene(nsp.ribo,
              site        = "stop",
              normalize   = TRUE,
              experiment = final_nsp1,
              title       = "Stop Site Coverage",
              range.lower = 28,
              range.upper = 35)

plot_metagene(nsp.ribo,
              site        = "stop",
              normalize   = TRUE,
              experiment = final_wt,
              title       = "Stop Site Coverage",
              range.lower = 28,
              range.upper = 35)

plot_metagene(nsp.ribo,
              site        = "stop",
              normalize   = TRUE,
              experiment = final_nsp2,
              title       = "Stop Site Coverage",
              range.lower = 28,
              range.upper = 35)


plot_metagene(nsp.ribo,
              site        = "start",
              normalize   = TRUE,
              experiment = final_nsp1,
              title       = "Start Site Coverage",
              range.lower = 28,
              range.upper = 35)

plot_metagene(nsp.ribo,
              site        = "start",
              normalize   = TRUE,
              experiment = final_wt,
              title       = "Start Site Coverage",
              range.lower = 28,
              range.upper = 35)

plot_metagene(nsp.ribo,
              site        = "start",
              normalize   = TRUE,
              experiment = final_nsp2,
              title       = "Start Site Coverage",
              range.lower = 28,
              range.upper = 35)



expressed_ribo = rowSums ( cpm(rcw[,-1]) > 1) > 2
c1 = cor(rcw[expressed_ribo,-1], method = "spearman")
pdf('./Figures/riboseq_replicate_correlations.pdf', height = 5, width = 5)
pheatmap(c1, main  = "Ribosome Profiling")
dev.off()

expressed_rna = rowSums ( cpm(rnaseq_w[,-c(1:7)]) > 1) > 2
c2 = cor(rnaseq_w[expressed_rna,-c(1:7)], method = "spearman")
pdf('./Figures/rnaseq_replicate_correlations.pdf', height = 5, width = 5)
pheatmap(c2, main  = "RNA-Seq")
dev.off()
                    

## Representative pairwise comparisons
## Fix the x - y axis range and remove p-value
plot_pairwise_relationships <- function (counts_w, id1, id2, xlab = "", ylab= "") { 
  sp <- ggscatter(counts_w, x = id1, 
                  y = id2,
                  #                add = "reg.line", conf.int = FALSE,     
                  #                add.params = list(color = "blue", size = 0.5),
                  font.family = "Helvetica", size = 0.5,
                  color = "gray", alpha = 0.3, ggtheme = theme_bw()) 
formatted =   sp +   
    scale_x_log10(labels = scales::label_number_si(), limits = c(0.3, 100000)) +   
    scale_y_log10(labels = scales::label_number_si(), limits = c(0.3, 100000)) + 
    labs (x=xlab, y = ylab) +
    stat_cor(method = "spearman", 
             aes(label = ..r.label..), 
             label.x = 0.1, label.y = 4.5, 
             cor.coef.name = "rho", digits = 3)  + 
    geom_hex(bins= 58, aes(alpha=log10(..count..) ), fill="#440154FF" )   
return (formatted)  
}

p1 = plot_pairwise_relationships(rcw, "20200717-NSP1-HEK-1-setB", 
                                 "20200717-NSP1-HEK-2-setB", 
                                 xlab = "NSP1_Rep1", ylab = "NSP1_Rep2")
p2 = plot_pairwise_relationships(rcw, "20200717-NSP1-HEK-1-setB", 
                                 "20200717-NSP1-HEK-3-setB", 
                                 xlab = "NSP1_Rep1", ylab = "NSP1_Rep3")
p3 = plot_pairwise_relationships(rcw, "20200717-NSP1-HEK-2-setB",
                                 "20200717-NSP1-HEK-3-setB", 
                                 xlab = "NSP1_Rep2", ylab = "NSP1_Rep3")
p4 = plot_pairwise_relationships(rcw, "20200717-NSP2-HEK-1-setB", 
                                 "20200717-NSP2-HEK-2-setB", 
                                 xlab = "NSP2_Rep1", ylab = "NSP2_Rep2")
p5 = plot_pairwise_relationships(rcw, "20200717-NSP2-HEK-1-setB", 
                                 "20200717-NSP2-HEK-3-setB", 
                                 xlab = "NSP2_Rep1", ylab = "NSP2_Rep3")
p6 = plot_pairwise_relationships(rcw, "20200717-NSP2-HEK-2-setB", 
                                 "20200717-NSP2-HEK-3-setB", 
                                 xlab = "NSP2_Rep2", ylab = "NSP2_Rep3")
p7 = plot_pairwise_relationships(rcw, "20200717-WT-HEK-1-setB", 
                                 "20200717-WT-HEK-2-setB", 
                                 xlab = "WT_Rep1", ylab = "WT_Rep2")
p8 = plot_pairwise_relationships(rcw, "20200717-WT-HEK-1-setB", 
                                 "20200717-WT-HEK-3-setB", 
                                 xlab = "WT_Rep1", ylab = "WT_Rep3")
p9 = plot_pairwise_relationships(rcw, "20200717-WT-HEK-2-setB", 
                                 "20200717-WT-HEK-3-setB", 
                                 xlab = "WT_Rep2", ylab = "WT_Rep3")

pdf('./Figures/riboseq_correlation.pdf', height = 3, width = 8)

row1 = plot_grid(p1 +  theme(legend.position="none"),
          p2 +  theme(legend.position="none"),
          p3 +  theme(legend.position="none"),
          label_size = 12, labels = c("A", "B", "C"),
          label_fontfamily = "Helvetica", ncol =3, align = "vh" )
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(row1, legend, rel_widths = c(3, .4))

row2 = plot_grid( p4 +  theme(legend.position="none"),
          p5 + theme(legend.position="none"),
          p6 +  theme(legend.position="none"),
          label_size = 12, labels = c("D", "E", "F"),
          label_fontfamily = "Helvetica", ncol =3, align = "vh" )
legend <- get_legend(
  # create some space to the left of the legend
  p4 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(row2, legend, rel_widths = c(3, .4))

row3 = plot_grid (p7 + theme(legend.position="none"),
          p8 + theme(legend.position="none"),
          p9 + theme(legend.position="none"), 
          label_size = 12, labels = c("G", "H", "I"),
          label_fontfamily = "Helvetica", ncol =3, align = "vh")
legend <- get_legend(
  # create some space to the left of the legend
  p7 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(row3, legend, rel_widths = c(3, .4))

dev.off()

## RNA-Seq correlations
p1 = plot_pairwise_relationships(rnaseq_w, "20200717-NSP1-HEK-1-setB", 
                                 "20200717-NSP1-HEK-2-setB", 
                                 xlab = "NSP1_Rep1", ylab = "NSP1_Rep2")
p2 = plot_pairwise_relationships(rnaseq_w, "20200717-NSP1-HEK-1-setB", 
                                 "20200717-NSP1-HEK-3-setB", 
                                 xlab = "NSP1_Rep1", ylab = "NSP1_Rep3")
p3 = plot_pairwise_relationships(rnaseq_w, "20200717-NSP1-HEK-2-setB",
                                 "20200717-NSP1-HEK-3-setB", 
                                 xlab = "NSP1_Rep2", ylab = "NSP1_Rep3")
p4 = plot_pairwise_relationships(rnaseq_w, "20200717-NSP2-HEK-1-setB", 
                                 "20200717-NSP2-HEK-2-setB", 
                                 xlab = "NSP2_Rep1", ylab = "NSP2_Rep2")
p5 = plot_pairwise_relationships(rnaseq_w, "20200717-NSP2-HEK-1-setB", 
                                 "20200717-NSP2-HEK-3-setB", 
                                 xlab = "NSP2_Rep1", ylab = "NSP2_Rep3")
p6 = plot_pairwise_relationships(rnaseq_w, "20200717-NSP2-HEK-2-setB", 
                                 "20200717-NSP2-HEK-3-setB", 
                                 xlab = "NSP2_Rep2", ylab = "NSP2_Rep3")
p7 = plot_pairwise_relationships(rnaseq_w, "20200717-WT-HEK-1-setB", 
                                 "20200717-WT-HEK-2-setB", 
                                 xlab = "WT_Rep1", ylab = "WT_Rep2")
p8 = plot_pairwise_relationships(rnaseq_w, "20200717-WT-HEK-1-setB", 
                                 "20200717-WT-HEK-3-setB", 
                                 xlab = "WT_Rep1", ylab = "WT_Rep3")
p9 = plot_pairwise_relationships(rnaseq_w, "20200717-WT-HEK-2-setB", 
                                 "20200717-WT-HEK-3-setB", 
                                 xlab = "WT_Rep2", ylab = "WT_Rep3")

# pdf('./Figures/rnaseq_correlation.pdf', height = 3, width = 8)
tiff ('./Figures/rnaseq_correlation_row1.tiff', 
      height = 3, width = 8, 
      units = 'in', res = 300,
      compression = "none")
row1 = plot_grid(p1 +  theme(legend.position="none"),
                 p2 +  theme(legend.position="none"),
                 p3 +  theme(legend.position="none"),
                 label_size = 12, labels = c("A", "B", "C"),
                 label_fontfamily = "Helvetica", ncol =3, align = "vh" )
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(row1, legend, rel_widths = c(3, .4))
dev.off()

tiff ('./Figures/rnaseq_correlation_row2.tiff', 
      height = 3, width = 8, 
      units = 'in', res = 300,
      compression = "none")

row2 = plot_grid( p4 +  theme(legend.position="none"),
                  p5 + theme(legend.position="none"),
                  p6 +  theme(legend.position="none"),
                  label_size = 12, labels = c("D", "E", "F"),
                  label_fontfamily = "Helvetica", ncol =3, align = "vh" )
legend <- get_legend(
  # create some space to the left of the legend
  p4 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(row2, legend, rel_widths = c(3, .4))
dev.off()

tiff ('./Figures/rnaseq_correlation_row3.tiff', 
      height = 3, width = 8, 
      units = 'in', res = 300,
      compression = "none")
row3 = plot_grid (p7 + theme(legend.position="none"),
                  p8 + theme(legend.position="none"),
                  p9 + theme(legend.position="none"), 
                  label_size = 12, labels = c("G", "H", "I"),
                  label_fontfamily = "Helvetica", ncol =3, align = "vh")
legend <- get_legend(
  # create some space to the left of the legend
  p7 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(row3, legend, rel_widths = c(3, .4))

dev.off()

## ERCC analyses
## This should be run to compare pairwise comparisons.
## Mix 1 -> NSP1; Mix2 NSP2; Mix1 - WT

second_set = ercc[, c(1, 8:13)]
colnames(second_set) = c("Feature", "NSP1_1", "NSP1_2", "NSP1_3", 
                        "NSP2_1", "NSP2_2", "NSP2_3")


exDat_secondset <- runDashboard(datType="count", isNorm = FALSE,
                               exTable=second_set,
                               filenameRoot="SecondSet", 
                               sample1Name="NSP1",
                               sample2Name="NSP2", 
                               erccmix="RatioPair",
                               erccdilution=1/10, spikeVol=1,
                               totalRNAmass=0.100, choseFDR=0.1)
summary(exDat_secondset)
grid.arrange(exDat_secondset$Figures$dynRangePlot)
grid.arrange(exDat_secondset$Figures$rocPlot)
grid.arrange(exDat_secondset$Figures$lodrERCCPlot)
grid.arrange(exDat_secondset$Figures$maPlot)
grid.arrange(exDat_secondset$Figures$r_mPlot)

# 
exDat_secondset$erccInfo$idColsSRM
exDat_secondset$Results$r_mDat

# total_rnaseq_counts = c(20592143,21453893, 21681715,
#                         21458498, 21278366,25666116, 
#                         34840423, 40997938,25135784, 
#                         35003078, 42398224, 36914462, 
#                         30927366, 25935593, 46271445)
counts_per_million <- function (count_matrix, total_counts) { 
  return (t(apply(count_matrix, 1, function(x){1000000*x/total_counts})))
}

ercc_only = ercc[1:92,]
ercc_cpm = counts_per_million(ercc_only[,-1], colSums(rnaseq_w[,-1]))
row.names(ercc_cpm) = ercc[1:92,1]

## Jacknife estimate of the ratio assuming replicate pairing
# Durbin J (1959) A note on the application of Quenouille's method of bias reduction to estimation of ratios. Biometrika 46: 477-480
# Bias is at most O(n^-2)
jacknife_ratio = function (vec1, vec2) {
  n = length(vec1)
  r = mean(vec1) / mean(vec2)
  sub_r = 0
  for (pair in 1:n) { 
    sub_r = sub_r + mean(vec1[pair]) / mean(vec2[pair])    
  }
  r_cor = n*r - (n-1)/n * sub_r
  return(r_cor)
}

nsp1_idx = 7:9
nsp2_idx = 10:12
wt_idx = 13:15

jacknife_ratios = data.frame(NSP1_NSP2 = c(), NSP1_WT = c(), WT_NSP2 = c() )
for (ercc_spike in 1:nrow(ercc_cpm)) { 
  jacknife_ratios[ercc_spike, 1] = jacknife_ratio(ercc_cpm[ercc_spike, nsp1_idx], ercc_cpm[ercc_spike, nsp2_idx])
  jacknife_ratios[ercc_spike, 2] = jacknife_ratio(ercc_cpm[ercc_spike, wt_idx], ercc_cpm[ercc_spike, nsp2_idx])
  jacknife_ratios[ercc_spike, 3] = jacknife_ratio(ercc_cpm[ercc_spike, nsp1_idx], ercc_cpm[ercc_spike, wt_idx])
  
}
colnames(jacknife_ratios) = c("NSP1 to NSP2 Ratio", "WT to NSP2 Ratio", "NSP1 to WT Ratio")

expected_observed = jacknife_ratios[expressed_ercc,]
expected_observed$Feature =  row.names(ercc_cpm)[expressed_ercc]
exDat_secondset$erccInfo$idColsSRM$Feature = as.character(exDat_secondset$erccInfo$idColsSRM$Feature)

expected_observed2 = merge(expected_observed, exDat_secondset$erccInfo$idColsSRM, 
                          by= "Feature")

half = expected_observed2$Ratio == "1:2"
equal = expected_observed2$Ratio == "1:1"
twothirds = expected_observed2$Ratio == "1:1.5"
four = expected_observed2$Ratio == "4:1"

boxplot(expected_observed2$`NSP1 to NSP2 Ratio`[half])
abline(h = 0.5)
boxplot(expected_observed2$`NSP1 to NSP2 Ratio`[equal])
abline(h = 1)
boxplot(expected_observed2$`NSP1 to NSP2 Ratio`[twothirds])
abline(h = 1/1.5)
boxplot(expected_observed2$`NSP1 to NSP2 Ratio`[four])
abline(h = 4)

plot(jacknife_ratios[,1], jacknife_ratios[,2])
abline (a = 0, b = 1)
expressed_ercc = rowSums(ercc_cpm[,7:15] > 10) > 3

jacknife_ratios_tidy = melt(jacknife_ratios[expressed_ercc,-2])

ggplot(jacknife_ratios_tidy, 
       aes(x= variable , y = value) ) +
  geom_point(alpha = 0.5, cex = 0.75, position = "jitter") +
#  geom_boxplot(color = "black", alpha=0, varwidth = T, outlier.shape = NA) + 
  geom_hline(yintercept = 1, color = "darkblue", lty = 3) + 
  geom_hline(yintercept = 1/1.5, color = "darkblue", lty = 3) + 
  geom_hline(yintercept = 1/2, color = "darkblue", lty = 3) + 
  geom_hline(yintercept = 4, color = "darkblue", lty = 3) +
  coord_flip() + 
  ylab("Jackknife Ratio Estimate") + xlab("")  +  theme_bw() 


## We can also use the below data to simplify the explanation
str(exDat_secondset$Results$r_mDat)


## Differential Expression 
## For RNA-Seq analyses we can include the ERCC counts in DE. 


## Differential RNA-Seq w/ERCC
nsp_rna <- factor(c("NSP1","NSP1","NSP1",
                    "NSP2","NSP2","NSP2",
                    "WT","WT","WT" )  ) 
r <- DGEList(counts=ercc[,8:16], group= nsp_rna, genes = ercc[,1])
keep <- filterByExpr(r)
r <- r[keep,,keep.lib.sizes=FALSE]
r <- calcNormFactors(r, method= "TMM")
design <- model.matrix(~0+nsp_rna)
r <- estimateDisp(r,design)

plotBCV(r, cex = 0.3)
plotMDS(r, labels = c("NSP1_SetB", "NSP1_SetB", "NSP1_SetB", 
                      "NSP2_SetB", "NSP2_SetB", "NSP2_SetB",
                      "WT_SetB", "WT_SetB", "WT_SetB"))

fit <- glmQLFit(r,design)
## NSP1 - NSP2
qlf_rna <- glmQLFTest(fit, contrast=c(1, -1, 0 ))
## NSP1 - WT ~ No ERCC controls are DE. 
# qlf <- glmQLFTest(fit, contrast=c(1, 0, -1 ))
## NSP2 - WT
# qlf <- glmQLFTest(fit, contrast=c(0, 1, -1 ))



# We might want to remove ERCC from Volcano  
ercc_idx  = grep ('ERCC' , qlf_rna$genes$genes )
qlf_rna$table = qlf_rna$table[-ercc_idx,]
qlf_rna$genes = qlf_rna$genes[-ercc_idx,]

topTags(qlf_rna)
summary(decideTests(qlf_rna, p.value = 0.05, adjust.method = "fdr"))
summary(decideTests(qlf_rna, p.value = 0.05, lfc = 1, adjust.method = "fdr"))

p1 = plotMD(qlf_rna, p.value = 0.05, hl.cex = 0.5, main = "", 
            hl.col = color.palette0(2), legend = F)

nsp1_low_rna = qlf_rna$genes$genes [ decideTests(qlf_rna, p.value = 0.05, adjust.method = "fdr") == -1 ]
nsp1_high_rna  = qlf_rna$genes$genes [ decideTests(qlf_rna, p.value = 0.05, adjust.method = "fdr") == 1 ]
fdr5_nominalp_RNA  = 0.004415391

# histone H3-K4 methylation, mRNA processing/splicing - RNA high
# RNA low ~ ribosome; high logFC is driven by lowCPM genes
selected_genes_rna = c("ASH1L", "KMT2B", "KMT2C", "KMT2D",
                       "OGT", "PAXIP1", "SETD1A", "SETD1B", "TET3",
                       "BUD13","CASC3","CDK13","CHERP","CSTF2T",
                       "DDX17","DDX5","DHX9",
                       "HNRNPA0","HNRNPH3","HSPA8","PNN",
                       "PRPF3","PRPF4B",  "UPF3B", 
                       "SF1","SFPQ","SFSWAP", "SRRT", "SRSF4", "SRSF6", 
                       "MALT1", "ATF3"
                      )
# Removing logFC to improve visibility
selected_genes_rna = c("ASH1L", "KMT2B", 
                       "OGT", "SETD1A", "SETD1B", "TET3",
                       "BUD13","CASC3","CDK13","CHERP","CSTF2T",
                       "DDX17","DDX5",
                       "HNRNPA0","HSPA8","PNN",
                       "PRPF3","PRPF4B",
                       "SF1", "SRSF4",  
                       "MALT1", "ATF3"
)
keyvals = ifelse(qlf_rna$table$logFC < 0 & qlf_rna$table$PValue < fdr5_nominalp_RNA, '#6e005f',
                 ifelse(qlf_rna$table$logFC > 0 & qlf_rna$table$PValue < fdr5_nominalp_RNA, '#045275', 'grey60'))
names(keyvals)[keyvals == '#6e005f'] <- 'lowRNA'
names(keyvals)[keyvals == '#045275'] <- 'highRNA'
names(keyvals)[keyvals == 'grey60'] <- 'NS'
 

pdf("./Figures/Volcano_RNA.pdf", width = 5, height = 6)
EnhancedVolcano(qlf_rna$table,
                xlim = c(-3.3, 2.8),
                ylim = c(0, 23),
                lab = strip_extension(as.character(qlf_rna$genes)),
                x = 'logFC',
                y = 'PValue',
                title = "RNA Expression", 
                subtitle = "NSP1 vs NSP2 expression",
                axisLabSize = 14,
                titleLabSize = 16,
                subtitleLabSize = 14,
                legendIconSize = 2.0,
                legendLabSize = 12,
                pCutoff = fdr5_nominalp_RNA,
                FCcutoff = 0, 
                cutoffLineType = 'blank',
                vline = c(0),
                #                vlineCol = c('grey50', 'grey0','grey50'),
                vlineType = 'solid',
                vlineWidth = 0.5,
                drawConnectors = TRUE,
                selectLab = selected_genes_rna,
                pointSize = c(ifelse(qlf_rna$table$PValue< fdr5_nominalp_RNA, 1, 0.2)),
                caption = "",
                legendLabels = c("NS", "NS", "NS", "5% FDR"),
                widthConnectors = 0.5,
                colConnectors = 'black',
                colCustom = keyvals, 
                #                col=c('grey60', 'grey60', 'grey60', color.palette0(2)[2]),
                colAlpha = 0.9,
                gridlines.minor = FALSE,
                gridlines.major = T)
dev.off()



write.table(strip_extension(as.character(nsp1_low_rna) ) , './Figures/NSP1_NSP2_RNA_Low.txt', 
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(strip_extension(as.character(nsp1_high_rna) ) , './Figures/NSP1_NSP2_RNA_High.txt', 
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(strip_extension(as.character(qlf_rna$genes$genes) ) , './Figures/AllExpressed_Genes_RNA.txt', 
            row.names = F, col.names = F, quote = F, sep = "\t")

write.csv(topTags(qlf_rna, sort.by = "logFC", p.value = 0.1, 447) , './Figures/NSP1_NSP2_RNA_LogFC.csv', 
          row.names = F,  quote = F)


## Differential Expression - RiboSeq Only
nsp <- factor(c("NSP1","NSP1","NSP1",
                  "NSP2","NSP2","NSP2",
                  "WT","WT","WT" )  )
y <- DGEList(counts=rcw[,-1],group=nsp, genes = rcw[,1])
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method= "TMM")
design <- model.matrix(~0+nsp)

## Might be best to use a common dispersion for UTR analyses
y <- estimateCommonDisp(y,design)
y <- estimateDisp(y,design)


plotBCV(y, cex = 0.3)
plotMDS(y, labels = c("NSP1", "NSP1", "NSP1",
                      "NSP2", "NSP2", "NSP2",
                      "WT", "WT", "WT"))

fit <- glmQLFit(y,design)
# Most important comparison is NSP1 vs NSP2
qlf <- glmQLFTest(fit, contrast=c(1, -1, 0 ))

# The additional comparisons reveal NSP expression leads to major changes compared to WT cells
qlf <- glmQLFTest(fit, contrast=c(1, 0, -1 ))
qlf <- glmQLFTest(fit, contrast=c(0, 1, -1 ))
qlf <- glmQLFTest(fit, contrast=c(1, -1/2, -1/2 ))

topTags(qlf, sort.by = "logFC", p.value = 0.05)
summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))
plotMD(qlf, p.value = 0.05, hl.cex = 0.5, main = "",
       hl.col = color.palette0(2))

# nsp1_low = qlf$genes$genes [ decideTests(qlf, p.value = 0.05, adjust.method = "fdr") == -1 ]
# nsp1_high  = qlf$genes$genes [ decideTests(qlf, p.value = 0.05, adjust.method = "fdr") == 1 ]


## Differential Expression - TE
all_counts = merge(rnaseq_w[,-c(2:7)], rcw, by= "transcript")
nsp_exptype <- factor(c("NSP1.RNA","NSP1.RNA","NSP1.RNA",
                "NSP2.RNA","NSP2.RNA","NSP2.RNA",
                "WT.RNA","WT.RNA","WT.RNA", 
                "NSP1.Ribo","NSP1.Ribo","NSP1.Ribo",
                "NSP2.Ribo","NSP2.Ribo","NSP2.Ribo",
                "WT.Ribo","WT.Ribo","WT.Ribo"  )  ) 

y <- DGEList(counts=all_counts[,-1],group=nsp_exptype, genes = all_counts[,1])
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0+nsp_exptype)
colnames(design) <- levels(nsp_exptype)
y <- estimateDisp(y,design)

plotBCV(y)
plotMDS(y)

fit <- glmQLFit(y,design)

my.contrasts <- makeContrasts(
  TE_NSP1vsNSP2 = (NSP1.Ribo - NSP1.RNA) - (NSP2.Ribo - NSP2.RNA), 
  RiboNSP1vsNSP2 = NSP1.Ribo - NSP2.Ribo, 
  RiboNSP1vsWT = NSP1.Ribo - WT.Ribo,
  RiboNSP2vsWT = NSP2.Ribo - WT.Ribo, 
  levels = design
)

# qlf <- glmQLFTest(fit, contrast=my.contrasts[,"RiboNSP1vsNSP2"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"TE_NSP1vsNSP2"])    

summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))

summary(decideTests(qlf, p.value = 0.1, adjust.method = "fdr"))
summary(decideTests(qlf, p.value = 0.2, adjust.method = "fdr"))
fdr20_nominalp = 0.005972371
fdr10_nominalp = 0.001927183
fdr5_nominalp = 0.0007377652


plotMD(qlf, p.value = 0.05, hl.cex = 0.75, main = "",
       hl.col = color.palette0(2), legend = F)

# Note change of xlim   
selected_genes = c("RPS12", "RPS19", "RPS27", "RPL12", 
                   "EEF1A1",  "EEF1A2", "EEF1B2", "EEF1D", 
                   "EIF2A", "EIF3E", "EIF3F", "EIF3G", 
                   "EIF3H", "EIF3L", "EIF4B",
                   "EEF2", "MARS", "NOB1", "NHP2", "BOP1",
                   "CSDE1", "PCBP2", "STAU1", "CAPRIN1", 
                   "ATF4", "FAXDC2", "ELOVL7")
keyvals = ifelse(qlf$table$logFC < 0 & qlf$table$PValue < fdr5_nominalp, '#6e005f',
                 ifelse(qlf$table$logFC > 0 & qlf$table$PValue < fdr5_nominalp, '#045275', 'grey60'))
names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60'] <- 'NS'

pdf("./Figures/Volcano_TE.pdf", width = 5, height = 6)
EnhancedVolcano(qlf$table,
                xlim = c(-4,3),
                ylim = c(0, 18),
                lab = strip_extension(as.character(qlf$genes$genes)),
                x = 'logFC',
                y = 'PValue',
                title = "Translation Efficiency", 
                subtitle = "NSP1 vs NSP2 expression",
                axisLabSize = 14,
                titleLabSize = 16,
                subtitleLabSize = 14,
                legendIconSize = 2.0,
                legendLabSize = 12,
                pCutoff = fdr5_nominalp,
                FCcutoff = 0, 
                cutoffLineType = 'blank',
                vline = c(0),
#                vlineCol = c('grey50', 'grey0','grey50'),
                vlineType = 'solid',
                vlineWidth = 0.5,
                drawConnectors = TRUE,
                selectLab = selected_genes,
                pointSize = c(ifelse(qlf$table$PValue< fdr5_nominalp, 1, 0.2)),
                caption = "",
                legendLabels = c("NS", "NS", "NS", "5% FDR"),
                widthConnectors = 0.5,
                colConnectors = 'black',
                colCustom = keyvals, 
#                col=c('grey60', 'grey60', 'grey60', color.palette0(2)[2]),
                colAlpha = 0.9,
                gridlines.minor = FALSE,
                gridlines.major = T)
dev.off()


nsp1_low = qlf$genes$genes [ decideTests(qlf, p.value = 0.05, adjust.method = "fdr") == -1 ]
nsp1_high  = qlf$genes$genes [ decideTests(qlf, p.value = 0.05, adjust.method = "fdr") == 1 ]

write.table(strip_extension(as.character(nsp1_low) ) , './Figures/NSP1_NSP2_TE_Low.txt', 
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(strip_extension(as.character(nsp1_high) ) , './Figures/NSP1_NSP2_TE_High.txt', 
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(strip_extension(as.character(qlf$genes$genes) ) , './Figures/AllExpressed_Genes_TE.txt', 
            row.names = F, col.names = F, quote = F, sep = "\t")
write.csv(topTags(qlf, sort.by = "logFC", p.value = 0.1, 231) , './Figures/NSP1_NSP2_TE_LogFC.csv', 
            row.names = F,  quote = F)

topTags(qlf, sort.by = "logFC", p.value = 0.05)
topTags(qlf, sort.by = "PValue", p.value = 0.05)
o <- order(qlf$table$PValue)
cpm(y)[o[1:5],]
qlf$genes[o[1:5],]

## Calculate CPM Ratios and generate dotplots for selected genes

selected_genes_rp = c("RPS5", "RPS12", "RPS19", "RPS10", "RPS18", "RPS25", "RPS27", 
                      "RPL3", "RPL10", "RPL12", "RPL13A", "RPL34", "RPL3", "RPL18A") 
selected_genes_trans = c("EEF1A1",  "EEF1A2", "EEF1B2", "EEF1D", 
                         "EIF2A", "EIF3E", "EIF3F", "EIF3G", 
                         "EIF3H", "EIF3L", "EIF4B","EEF2", "MARS")
selected_genes_other = c("IPO5", "IPO7", "THOC5")
selected_genes_biogenesis = c( "NOB1", "NHP2", "BOP1")
selected_genes_chaperone = c ("HSPA8", "HSP90AB1", "BAG6")
selected_genes_rbp = c ("CSDE1", "PCBP2", "STAU1", "CAPRIN1", "FXR1")

cpm_ratios = data.frame (genes = y$genes,
                         NSP1_Rep1 = cpm(y)[,10]/cpm(y)[,1],
                         NSP1_Rep2 = cpm(y)[,11]/cpm(y)[,2],
                         NSP1_Rep3 = cpm(y)[,12]/cpm(y)[,3],
                         NSP2_Rep1 = cpm(y)[,13]/cpm(y)[,4],
                         NSP2_Rep2 = cpm(y)[,14]/cpm(y)[,5],
                         NSP2_Rep3 = cpm(y)[,15]/cpm(y)[,6]
                         )
subset_selected_genes = function (cpms, gene_list) { 
  tidy = melt(cpms[ strip_extension(as.character(cpms$genes)) %in% gene_list,  ] )
  tidy$variable = sapply(strsplit(as.character(tidy$variable), split = "_"), "[[", 1)
  return(tidy)
}

rp_cpm_ratio_tidy  = subset_selected_genes(cpm_ratios, selected_genes_rp)
trn_cpm_ratio_tidy  = subset_selected_genes(cpm_ratios, selected_genes_trans)
other_cpm_ratio_tidy  = subset_selected_genes(cpm_ratios, selected_genes_other)
biogenesis_cpm_ratio_tidy =   subset_selected_genes(cpm_ratios, selected_genes_biogenesis)
chaperone_cpm_ratio_tidy = subset_selected_genes(cpm_ratios, selected_genes_chaperone)
rbp_cpm_ratio_tidy = subset_selected_genes(cpm_ratios, selected_genes_rbp)

s1 = ggplot(rp_cpm_ratio_tidy, 
  aes( x =  strip_extension(as.character(genes)), y = value) ) + 
  geom_jitter(aes(color = variable),
      position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2),
      size = 1 ) +
  # stat_summary( aes(color = variable),
  #  geom="line", lwd=2, fun = mean, position = position_dodge(0.4) ) + 
  scale_color_manual(values =  c('#089099', '#d12959')) + theme_bw()+
  ylab("Ratio") + xlab("") + ylim (0, 3.6)+ theme(legend.position="none")

s2 = ggplot(trn_cpm_ratio_tidy, 
            aes( x = strip_extension(as.character(genes)), y = value) ) + 
  geom_jitter(aes(color = variable),
              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2),
              size = 1 ) +
  # stat_summary( aes(color = variable),
  #  geom="line", lwd=2, fun = mean, position = position_dodge(0.4) ) + 
  scale_color_manual(values =  c('#089099', '#d12959')) + theme_bw() +
  ylab("Ratio") + xlab("") + ylim (0, 3.6) + theme(legend.position="none")

s3 = ggplot(other_cpm_ratio_tidy, 
            aes( x =  strip_extension(as.character(genes)), y = value) ) + 
  geom_jitter(aes(color = variable),
              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2),
              size = 1 ) +
  scale_color_manual(values =  c('#089099', '#d12959')) + theme_bw()+
  ylab("Ratio") + xlab("") + ylim (0, 3.6) +  theme(legend.position="none")


s4 = ggplot(biogenesis_cpm_ratio_tidy, 
            aes( x =  strip_extension(as.character(genes)), y = value) ) + 
  geom_jitter(aes(color = variable),
              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2),
              size = 1 ) +
  scale_color_manual(values =  c('#089099', '#d12959')) + theme_bw()+
  ylab("Ratio") + xlab("") + ylim (0, 3.6)  +  theme(legend.position="none")


s5 = ggplot(chaperone_cpm_ratio_tidy, 
            aes( x =  strip_extension(as.character(genes)), y = value) ) + 
  geom_jitter(aes(color = variable),
              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2),
              size = 1 ) +
  scale_color_manual(values =  c('#089099', '#d12959')) + theme_bw()+
  ylab("Ratio") + xlab("") + ylim (0, 3.6) +  theme(legend.position="none")


s6 = ggplot(rbp_cpm_ratio_tidy, 
            aes( x =  strip_extension(as.character(genes)), y = value) ) + 
  geom_jitter(aes(color = variable),
              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2),
              size = 1 ) +
  scale_color_manual(values =  c('#089099', '#d12959')) + theme_bw()+
  ylab("Ratio") + xlab("") + ylim (0, 3.6)  +  theme(legend.position="none")

pdf('./Figures/TE_Dotplot.pdf', height = 9, width = 9 )
bottom_row = plot_grid(s3, s4, s5, s6,
                       align = "vh", labels = c("E", "F", "G", "H"), nrow = 1)
plot_grid(s1, s2, bottom_row,
          align = "h", labels = c("C", "D", ""), ncol = 1)
dev.off()  

## OAZ1 - frameshift
## NOB1, NHP2, BOP1 -> rRNA processing
## viral process genes: go:0016032
## IPO5 - IPO7 - THOC5 influenza involvement
## HNRNPA1 
## HSPA8 - chaperone!
## Nucleophosmin - NPM1 
## PCBP2
## VIM - https://pubmed.ncbi.nlm.nih.gov/27423069/
## slc25a5- slc25a6 and others


## ILF3

go <- goana(qlf, species="Hs")
topGO(go, sort="down")
topGO(go, sort="up")

keg <- kegga(qlf, species="Hs")
topKEGG(keg, sort="up")
topKEGG(keg, sort="down")

## Compare with QTI data
qti = read.csv('./QTI_Seq.csv')
nsp1_high_qti = strip_extension(as.character(nsp1_high) )
nsp1_low_qti = strip_extension(as.character(nsp1_low) )

sum ( qti$Gene %in% nsp1_high_qti   ) / length(nsp1_high_qti)
length ( unique (qti$Gene[qti$Gene %in% nsp1_high_qti] ) ) 
length (  qti$Gene[qti$Gene %in% nsp1_high_qti]   ) 

qti_nsp_high = qti[qti$Gene %in% nsp1_high_qti , ]
qti_nsp_low = qti[qti$Gene %in% nsp1_low_qti , ]


## Given the higher overall FPKM, higher LTM is expected. Need to correct for this.
qti_annotated = qti[qti_nsp_high$Relative.to.aTIS == 0,]
qti_annotated$color = "#D3D3D3"
qti_annotated$color[qti_annotated$Gene %in% nsp1_high_qti] = "#6e0062"
ggplot(qti_annotated, 
       aes(x= Gene %in% nsp1_high_qti, y = FPKM_control) ) +
  geom_point(alpha = 0.5, cex = 0.75, position = "jitter", color = qti_annotated$color) +
  geom_boxplot(color = "black", alpha=0, varwidth = T, outlier.shape = NA) + 
  coord_flip() + 
  ylab("FPKM") + xlab("High TE")  +  theme_minimal() + theme(legend.position="none") 

## 

## Nearest neighbor matching by FPKM; the result still holds 
## We see that nsp1_highTE genes have high LTM values in these cells even when controlling for expression
## They likely to have lower normalized fold upon starvation. In other words, their initiation is preferentially slowed?
set.seed(3)
qti_TIS = qti[qti$Relative.to.aTIS ==0,]
qti_TIS$treat = rep(0)
qti_TIS$treat[ qti_TIS$Gene %in% nsp1_high_qti] = 1
m.out <- matchit(treat ~ FPKM_control , data = qti_TIS, 
                 method = "nearest" ) 
summary(m.out)
matched_qti_TIS = match.data(m.out)

matched_qti_TIS$color = "#D3D3D3"
matched_qti_TIS$color[matched_qti_TIS$Gene %in% nsp1_high_qti] = "#6e0062"

g1 = ggplot(matched_qti_TIS, 
       aes(x= Gene %in% nsp1_high_qti, y = FPKM_control) ) +
  geom_point(alpha = 0.5, cex = 0.75, position = "jitter", color = matched_qti_TIS$color) +
  geom_boxplot(color = "black", alpha=0, varwidth = T, outlier.shape = NA) + 
  coord_flip() + 
  ylab("FPKM") + xlab("High TE")  +  theme_minimal() + theme(legend.position="none") 

g2 = ggplot(matched_qti_TIS, 
       aes(x= Gene %in% nsp1_high_qti, y = Normalized_LTM_value_control.UQ.) ) +
  geom_point(alpha = 0.5, cex = 0.75, position = "jitter", color = matched_qti_TIS$color) +
  geom_boxplot(color = "black", alpha=0, varwidth = T, outlier.shape = NA) + 
  coord_flip() + 
  ylab("HEK293 Normalized Read Counts (LTM + Puro)") + xlab("High TE")  +  theme_minimal() + theme(legend.position="none") 

g3 = ggplot(matched_qti_TIS, 
       aes(x= Gene %in% nsp1_high_qti, y = Normalized_Fold_change) ) +
  geom_point(alpha = 0.5, cex = 0.75, position = "jitter", color = matched_qti_TIS$color) +
  geom_boxplot(color = "black", alpha=0, varwidth = T, outlier.shape = NA) + 
  geom_hline(yintercept = 1, color = "darkblue", lty = 3) + coord_flip() + 
  ylab("Starvation Fold Change") + xlab("High TE")  +  theme_minimal() + theme(legend.position="none") 


pdf("./Figures/QTI_Seq.pdf", height = 6, width = 6)
plot_grid(g1 , 
          g2 ,
          g3 ,
          align = "vh", labels = "AUTO", ncol = 1)
dev.off()


wilcox.test ( matched_qti_TIS$Normalized_LTM_value_control.UQ.[matched_qti_TIS$Gene %in% nsp1_high_qti], 
              matched_qti_TIS$Normalized_LTM_value_control.UQ.[!matched_qti_TIS$Gene %in% nsp1_high_qti] ) 
wilcox.test ( matched_qti_TIS$Normalized_Fold_change[matched_qti_TIS$Gene %in% nsp1_high_qti], 
              matched_qti_TIS$Normalized_Fold_change[!matched_qti_TIS$Gene %in% nsp1_high_qti] ) 

# ### 


## 5TOP
known_top = read.csv('./LARP_TOP_Philippe_ThoreenPNAS/KnownTop_S1.csv', header= F)
core_top = read.csv('./LARP_TOP_Philippe_ThoreenPNAS/CoreTOP_S7.csv')
hek_topscore = read.csv('./LARP_TOP_Philippe_ThoreenPNAS/HEK_TopScore.csv')

known_top$V1 %in% strip_extension(as.character(nsp1_high))
core_top$gene %in% strip_extension(as.character(nsp1_high))

new_top = core_top[!core_top$gene %in% known_top$V1, ]
str(new_top)

new_top$gene %in% strip_extension(as.character(nsp1_high))
sum (new_top$gene %in% strip_extension(as.character(nsp1_high)) ) 

sum ( hek_topscore$gene[hek_topscore$topscore > 3] %in% strip_extension(as.character(nsp1_high))  ) 

# What is the HEK TopScore distribution of known_top and core_top is all very high
# There is more variability for TopScore of high TE

hek_topscore$KnownCoreTop = hek_topscore$gene %in% known_top$V1 | hek_topscore$gene %in% new_top$gene
hek_topscore$highTE = hek_topscore$gene %in% strip_extension(as.character(nsp1_high)) 

boxplot(hek_topscore$topscore ~ hek_topscore$highTE)

# https://www.nature.com/articles/s41418-019-0350-5 - UBAP2L - stress granule assembly including RPS
# BAG6 - https://www.genecards.org/cgi-bin/carddisp.pl?gene=BAG6 => Another chaperone


# One Idea is to remove the "False category" 
t1 = ggplot(hek_topscore[hek_topscore$highTE,], 
       aes(x= highTE , y = topscore, color = KnownCoreTop) ) +
  geom_point(alpha = 0.5, cex = 0.75, position = "jitter") +
#  geom_boxplot(color = "black", alpha=0, varwidth = T, outlier.shape = NA) + 
#  coord_flip() + 
  ylab("HEK293 TopScore") + xlab("High TE")  +  
  theme_bw() + theme(legend.position = "none")

t2 = ggplot(hek_topscore[hek_topscore$highTE,]) +
  geom_density(aes(x = topscore) ) + 
  theme_minimal_hgrid() + coord_flip() +
  xlab("") +   
  geom_vline(xintercept = 0.70801050, color = "darkblue", lty = 3) + 
  geom_vline(xintercept = 1.83393814, color = "darkblue", lty = 3) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

quantile (hek_topscore$topscore, seq(0,1,.05))

pdf("./Figures/TopScoreDistribution.pdf", height = 3, width = 3)
plot_grid(t1, t2,
          align = "vh", nrow = 1)
dev.off()



# We can also plot TE as a function of TopScore
# Note that there are handful of duplicate genes 
qlf$table$gene = strip_extension(as.character(qlf$genes$genes))
te_top = merge(hek_topscore, qlf$table, by = "gene")

## Limit to high CPM; change color ramp

t2 = ggplot(te_top, 
       aes(y= logFC , x = logCPM, color = topscore ) ) +
  geom_point(alpha = 0.9, cex = 0.75) +
  scale_color_gradient2(midpoint= quantile(te_top$topscore)[4], low=color.palette0(2)[1], mid="white",
                        high= color.palette0(2)[2]) + 
  xlim(1,max(te_top$logCPM)) + ylim (-2.5, 2.5) + 
  geom_hline(yintercept = median(te_top$logFC), color = "darkblue", lty = 1) + 
  ylab("logFC TE") + xlab("logCPM")  +  theme_bw() 

pdf("./Figures/TopRNA.pdf", height = 4.5, width = 4.5)
plot_grid(t1, t2, rel_heights = c(1,2),
          align = "vh", labels = "AUTO", ncol = 1)
dev.off()




