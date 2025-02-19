#ploting pvalue histograms
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(patchwork)

# RDP
png(paste0("Plots/RDP/ShuffledPvalueHistograms.png"), width= 3600, height=4800, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL

all_files<-list.files("PkgResults/RDP/shuffledttest", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_t.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(RDP-", name, ") t-test Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/RDP/shuffledWilcoxon", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_wilcox.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(RDP-", name, ") Wilcoxon Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/RDP/shuffledDESeq2", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_deseq2.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(RDP-", name, ") DESeq2 Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/RDP/shufflededgeR", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_edgeR.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(RDP-", name, ") edgeR Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots + plot_layout(ncol = 3))
dev.off()

# dada2
png(paste0("Plots/dada2/ShuffledPvalueHistograms.png"), width= 3600, height=4800, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL

all_files<-list.files("PkgResults/dada2/shuffledttest", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_t.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(dada2-", name, ") t-test Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/dada2/shuffledWilcoxon", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_wilcox.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(dada2-", name, ") Wilcoxon Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/dada2/shuffledDESeq2", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_deseq2.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(dada2-", name, ") DESeq2 Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/dada2/shufflededgeR", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_edgeR.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(dada2-", name, ") edgeR Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots + plot_layout(ncol = 3))
dev.off()

# WGS
png(paste0("Plots/WGS/ShuffledPvalueHistograms.png"), width= 2400, height=4800, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL

all_files<-list.files("PkgResults/WGS/shuffledttest", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_t.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(WGS-", name, ") t-test Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/WGS/shuffledWilcoxon", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_wilcox.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(WGS-", name, ") Wilcoxon Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/WGS/shuffledDESeq2", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_deseq2.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(WGS-", name, ") DESeq2 Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

all_files<-list.files("PkgResults/WGS/shufflededgeR", pattern = "*.txt",full.names = TRUE)
for (file in all_files) {
  name <- gsub(basename(file), pattern="_edgeR.txt$", replacement="")
  results<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-PvalHistogram(results, "powderblue", paste0("(WGS-", name, ") edgeR Shuffeled p-value"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots + plot_layout(ncol = 2))
dev.off()
