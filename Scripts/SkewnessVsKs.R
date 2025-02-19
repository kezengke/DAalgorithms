#shuffle sample 100x and test all pval distribution
rm(list = ls())
library(MetagenomeTools)
library(coin)
library(DESeq2)
library(edgeR)
library(moments)
library(ggExtra)
library(patchwork)
library(ggplot2)


makePlot <- function(SkewKs, test, dataType) {
  p<-ggplot(SkewKs, aes(x = skew, y = ks)) +
    geom_point(alpha = 0.7, color = "orange") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
    labs(title = paste0("(", dataType, "-", name, ") ", test),
         x = "Skewness",
         y = "KS Test p-value") +
    theme_minimal()
  combinedP<-ggMarginal(p, type = "histogram", margins = "x", size = 5, fill = "orange", alpha = 0.7)
  return(combinedP)
}

# RDP
combined_plots<-NULL
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RDP/SkewVsKS(RDP).png"), width=4800, height=3600, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/RDPRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  NormcountsT<-read.table(paste0("CountsTables/RDPNorm/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  NormcountsT<-NormcountsT[, intersect(colnames(NormcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    ttestresults<-calcTtest(NormcountsT, shuffleMeta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormcountsT, shuffleMeta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(RawcountsT, shuffleMeta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(RawcountsT, shuffleMeta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  # KS distribution
  ttestKS<-as.numeric(apply(ttestpvals, 2, CalcKSpval))
  wilcoxKS<-as.numeric(apply(wilcoxpvals, 2, CalcKSpval))
  deseq2KS<-as.numeric(apply(deseq2pvals, 2, CalcKSpval))
  edgerKS<-as.numeric(apply(edgerpvals, 2, CalcKSpval))

  # Skewness
  ttestskew<-apply(ttestpvals, 2, skewness)
  wilcoxskew<-apply(wilcoxpvals, 2, skewness)
  deseq2skew<-apply(deseq2pvals, 2, skewness)
  edgerskew<-apply(edgerpvals, 2, skewness)

  colN<-c("skew", "ks")

  ttest<-data.frame(ttestskew, ttestKS)
  wilcox<-data.frame(wilcoxskew, wilcoxKS)
  deseq2<-data.frame(deseq2skew, deseq2KS)
  edger<-data.frame(edgerskew, edgerKS)

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

# dada2
combined_plots<-NULL
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/dada2/SkewVsKS(dada2).png"), width=4800, height=3600, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/dada2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  NormcountsT<-read.table(paste0("CountsTables/dada2Norm/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  NormcountsT<-NormcountsT[, intersect(colnames(NormcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    ttestresults<-calcTtest(NormcountsT, shuffleMeta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormcountsT, shuffleMeta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(RawcountsT, shuffleMeta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(RawcountsT, shuffleMeta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  # KS distribution
  ttestKS<-as.numeric(apply(ttestpvals, 2, CalcKSpval))
  wilcoxKS<-as.numeric(apply(wilcoxpvals, 2, CalcKSpval))
  deseq2KS<-as.numeric(apply(deseq2pvals, 2, CalcKSpval))
  edgerKS<-as.numeric(apply(edgerpvals, 2, CalcKSpval))

  # Skewness
  ttestskew<-apply(ttestpvals, 2, skewness)
  wilcoxskew<-apply(wilcoxpvals, 2, skewness)
  deseq2skew<-apply(deseq2pvals, 2, skewness)
  edgerskew<-apply(edgerpvals, 2, skewness)

  colN<-c("skew", "ks")

  ttest<-data.frame(ttestskew, ttestKS)
  wilcox<-data.frame(wilcoxskew, wilcoxKS)
  deseq2<-data.frame(deseq2skew, deseq2KS)
  edger<-data.frame(edgerskew, edgerKS)

  names(ttest)<-colN
  names(wilcox)<-colN
  names(deseq2)<-colN
  names(edger)<-colN

  p1<-makePlot(ttest, "t-test", "dada2")
  p2<-makePlot(wilcox, "Wilcoxon", "dada2")
  p3<-makePlot(deseq2, "DESeq2", "dada2")
  p4<-makePlot(edger, "edgeR", "dada2")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
  } else {
    combined_plots <- combined_plots + wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + wrap_elements(p4)
  }

}

print(combined_plots + plot_layout(ncol = 4))
dev.off()

# WGS
combined_plots<-NULL
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/WGS/SkewVsKS(WGS).png"), width=4800, height=2400, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/WGSRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  NormcountsT<-read.table(paste0("CountsTables/WGSNorm/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  NormcountsT<-NormcountsT[, intersect(colnames(NormcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    ttestresults<-calcTtest(NormcountsT, shuffleMeta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormcountsT, shuffleMeta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(RawcountsT, shuffleMeta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(RawcountsT, shuffleMeta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  # KS distribution
  ttestKS<-as.numeric(apply(ttestpvals, 2, CalcKSpval))
  wilcoxKS<-as.numeric(apply(wilcoxpvals, 2, CalcKSpval))
  deseq2KS<-as.numeric(apply(deseq2pvals, 2, CalcKSpval))
  edgerKS<-as.numeric(apply(edgerpvals, 2, CalcKSpval))

  # Skewness
  ttestskew<-apply(ttestpvals, 2, skewness)
  wilcoxskew<-apply(wilcoxpvals, 2, skewness)
  deseq2skew<-apply(deseq2pvals, 2, skewness)
  edgerskew<-apply(edgerpvals, 2, skewness)

  colN<-c("skew", "ks")

  ttest<-data.frame(ttestskew, ttestKS)
  wilcox<-data.frame(wilcoxskew, wilcoxKS)
  deseq2<-data.frame(deseq2skew, deseq2KS)
  edger<-data.frame(edgerskew, edgerKS)

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


