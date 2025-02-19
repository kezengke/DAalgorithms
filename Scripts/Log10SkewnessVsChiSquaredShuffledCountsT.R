#shuffle entire counts table 100x and test all pval distribution
rm(list = ls())
library(MetagenomeTools)
library(coin)
library(DESeq2)
library(edgeR)
library(moments)
library(ggExtra)
library(patchwork)
library(ggplot2)


makePlot <- function(SkewChi, test, dataType) {
  p<-ggplot(SkewChi, aes(x = skew, y = chi)) +
    geom_point(alpha = 0.7, color = "lightskyblue1") +
    geom_hline(yintercept = (log10(0.05) * -1), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
    labs(title = paste0("(", dataType, "-", name, ") ", test),
         x = "Skewness",
         y = "log10(Chi-squared Test p-value)") +
    theme_minimal() +
    xlim(-2, 2) +
    ylim(0, 31)

  combinedP<-ggMarginal(p, type = "histogram", margins = "x", size = 5, fill = "lightskyblue1", alpha = 0.7)
  return(combinedP)
}

# RDP
combined_plots<-NULL
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RDP/Log10SkewVsChiShuffleCountsT(RDP).png"), width=4800, height=3600, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # KS distribution
  ttestChi<-log10(as.numeric(apply(ttestpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
  wilcoxChi<-log10(as.numeric(apply(wilcoxpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
  deseq2Chi<-log10(as.numeric(apply(deseq2pvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
  edgerChi<-log10(as.numeric(apply(edgerpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1

  # Skewness
  ttestskew<-apply(ttestpvals, 2, skewness)
  wilcoxskew<-apply(wilcoxpvals, 2, skewness)
  deseq2skew<-apply(deseq2pvals, 2, skewness)
  edgerskew<-apply(edgerpvals, 2, skewness)

  colN<-c("skew", "chi")

  ttest<-data.frame(ttestskew, ttestChi)
  wilcox<-data.frame(wilcoxskew, wilcoxChi)
  deseq2<-data.frame(deseq2skew, deseq2Chi)
  edger<-data.frame(edgerskew, edgerChi)

  names(ttest)<-colN
  names(wilcox)<-colN
  names(deseq2)<-colN
  names(edger)<-colN

  p1<-makePlot(ttest, "t-test", "RDP")
  p2<-makePlot(wilcox, "Wilcoxon", "RDP")
  p3<-makePlot(deseq2, "DESeq2", "RDP")
  p4<-makePlot(edger, "edgeR", "RDP")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
  } else {
    combined_plots <- combined_plots + wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
  }

}

print(combined_plots + plot_layout(ncol = 4))
dev.off()

# # dada2
# combined_plots<-NULL
# all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)
#
# png(paste0("Plots/dada2/Log10SkewVsChiShuffleCountsT(dada2).png"), width=4800, height=3600, res = 300)
#
# for (file in all_files) {
#   name <- gsub(basename(file), pattern=".txt$", replacement="")
#
#   ttestpvals<-read.table(paste0("ShuffleCountsTDump/dada2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
#   wilcoxpvals<-read.table(paste0("ShuffleCountsTDump/dada2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
#   deseq2pvals<-read.table(paste0("ShuffleCountsTDump/dada2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
#   edgerpvals<-read.table(paste0("ShuffleCountsTDump/dada2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")
#
#   # KS distribution
#   ttestChi<-log10(as.numeric(apply(ttestpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
#   wilcoxChi<-log10(as.numeric(apply(wilcoxpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
#   deseq2Chi<-log10(as.numeric(apply(deseq2pvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
#   edgerChi<-log10(as.numeric(apply(edgerpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
#
#   # Skewness
#   ttestskew<-apply(ttestpvals, 2, skewness)
#   wilcoxskew<-apply(wilcoxpvals, 2, skewness)
#   deseq2skew<-apply(deseq2pvals, 2, skewness)
#   edgerskew<-apply(edgerpvals, 2, skewness)
#
#   colN<-c("skew", "chi")
#
#   ttest<-data.frame(ttestskew, ttestChi)
#   wilcox<-data.frame(wilcoxskew, wilcoxChi)
#   deseq2<-data.frame(deseq2skew, deseq2Chi)
#   edger<-data.frame(edgerskew, edgerChi)
#
#   names(ttest)<-colN
#   names(wilcox)<-colN
#   names(deseq2)<-colN
#   names(edger)<-colN
#
#   p1<-makePlot(ttest, "t-test", "dada2")
#   p2<-makePlot(wilcox, "Wilcoxon", "dada2")
#   p3<-makePlot(deseq2, "DESeq2", "dada2")
#   p4<-makePlot(edger, "edgeR", "dada2")
#
#   if (is.null(combined_plots)) {
#     combined_plots <- wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
#   } else {
#     combined_plots <- combined_plots + wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
#   }
#
# }
#
# print(combined_plots + plot_layout(ncol = 4))
# dev.off()

# WGS
combined_plots<-NULL
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/WGS/Log10SkewVsChiShuffleCountsT(WGS).png"), width=4800, height=2400, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # KS distribution
  ttestChi<-log10(as.numeric(apply(ttestpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
  wilcoxChi<-log10(as.numeric(apply(wilcoxpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
  deseq2Chi<-log10(as.numeric(apply(deseq2pvals, 2, CalcChiSquaredPval)) + 1e-30) * -1
  edgerChi<-log10(as.numeric(apply(edgerpvals, 2, CalcChiSquaredPval)) + 1e-30) * -1

  # Skewness
  ttestskew<-apply(ttestpvals, 2, skewness)
  wilcoxskew<-apply(wilcoxpvals, 2, skewness)
  deseq2skew<-apply(deseq2pvals, 2, skewness)
  edgerskew<-apply(edgerpvals, 2, skewness)

  colN<-c("skew", "chi")

  ttest<-data.frame(ttestskew, ttestChi)
  wilcox<-data.frame(wilcoxskew, wilcoxChi)
  deseq2<-data.frame(deseq2skew, deseq2Chi)
  edger<-data.frame(edgerskew, edgerChi)

  names(ttest)<-colN
  names(wilcox)<-colN
  names(deseq2)<-colN
  names(edger)<-colN

  p1<-makePlot(ttest, "t-test", "WGS")
  p2<-makePlot(wilcox, "Wilcoxon", "WGS")
  p3<-makePlot(deseq2, "DESeq2", "WGS")
  p4<-makePlot(edger, "edgeR", "WGS")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
  } else {
    combined_plots <- combined_plots + wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
  }

}

print(combined_plots + plot_layout(ncol = 4))
dev.off()


