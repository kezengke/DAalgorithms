# plotting taxa correlation histogram of datasets
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(patchwork)

CorHistogram <- function(table, plotCol) {
  spearCor <- cor(t(table), method = "spearman")
  corVal <- spearCor[lower.tri(spearCor, diag = FALSE)]
  pval<-as.numeric(CalcKSpval(corVal))
  hist(corVal,
       breaks = seq(-1, 1, by = 0.05),  # Adjust number of bins
       col = plotCol,
       border = "black",
       main = NULL,
       xlab = "Spearman Correlation",
       ylab = "Frequency")
  return(pval)
}

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RDP/TaxaCorrelationHistograms.png"), width=1200 * length(all_files), height=1200, res = 300)
par(mfrow = c(1, length(all_files)), mar = c(5, 6, 4, 1) + 0.1)

for (file in all_files) {
  name <- gsub(basename(file), pattern = ".txt$", replacement = "")
  countsT <- read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

  pVal<-CorHistogram(countsT, "olivedrab3")
  title(paste0("(RDP-", name, ")\nTaxa Correlation\nKS p-val(uniform): ", round(pVal, 4)))
}

dev.off()

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/dada2/TaxaCorrelationHistograms.png"), width=1200 * length(all_files), height=1200, res = 300)
par(mfrow = c(1, length(all_files)), mar = c(5, 6, 4, 1) + 0.1)

for (file in all_files) {
  name <- gsub(basename(file), pattern = ".txt$", replacement = "")
  countsT <- read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

  pVal<-CorHistogram(countsT, "olivedrab3")
  title(paste0("(RDP-", name, ")\nTaxa Correlation\nKS p-val(uniform): ", round(pVal, 4)))
}

dev.off()

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/WGS/TaxaCorrelationHistograms.png"), width=1200 * length(all_files), height=1200, res = 300)
par(mfrow = c(1, length(all_files)), mar = c(5, 6, 4, 1) + 0.1)

for (file in all_files) {
  name <- gsub(basename(file), pattern = ".txt$", replacement = "")
  countsT <- read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

  pVal<-CorHistogram(countsT, "olivedrab3")
  title(paste0("(RDP-", name, ")\nTaxa Correlation\nKS p-val(uniform): ", round(pVal, 4)))
}

dev.off()

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RNAseq/TaxaCorrelationHistograms.png"), width=1200 * length(all_files), height=1200, res = 300)
par(mfrow = c(1, length(all_files)), mar = c(5, 6, 4, 1) + 0.1)

for (file in all_files) {
  name <- gsub(basename(file), pattern = ".txt$", replacement = "")
  countsT <- read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

  pVal<-CorHistogram(countsT, "olivedrab3")
  title(paste0("(RDP-", name, ")\nTaxa Correlation\nKS p-val(uniform): ", round(pVal, 4)))
}

dev.off()
