#ploting pvalue histograms
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(patchwork)

# RDP
png(paste0("Plots/RDP/OriginalPvalueHistograms.png"), width= 3600, height=4800, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL

all_files<-list.files("PkgResults/RDP/ttest", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_t.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(RDP-", name, ") t-test Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/RDP/Wilcoxon", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_wilcox.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(RDP-", name, ") Wilcoxon Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/RDP/DESeq2", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_deseq2.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(RDP-", name, ") DESeq2 Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/RDP/edgeR", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_edgeR.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(RDP-", name, ") edgeR Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots + plot_layout(ncol = 3))
dev.off()

# dada2
png(paste0("Plots/dada2/OriginalPvalueHistograms.png"), width= 3600, height=4800, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL

all_files<-list.files("PkgResults/dada2/ttest", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_t.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(dada2-", name, ") t-test Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/dada2/Wilcoxon", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_wilcox.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(dada2-", name, ") Wilcoxon Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/dada2/DESeq2", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_deseq2.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(dada2-", name, ") DESeq2 Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/dada2/edgeR", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_edgeR.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(dada2-", name, ") edgeR Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots + plot_layout(ncol = 3))
dev.off()

# WGS
png(paste0("Plots/WGS/OriginalPvalueHistograms.png"), width= 2400, height=4800, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL

all_files<-list.files("PkgResults/WGS/ttest", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_t.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(WGS-", name, ") t-test Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/WGS/Wilcoxon", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_wilcox.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(WGS-", name, ") Wilcoxon Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/WGS/DESeq2", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_deseq2.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(WGS-", name, ") DESeq2 Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/WGS/edgeR", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_edgeR.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "thistle1", paste0("(WGS-", name, ") edgeR Original p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots + plot_layout(ncol = 2))
dev.off()
