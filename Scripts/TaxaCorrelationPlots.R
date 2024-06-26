# plotting taxa correlation histogram of datasets
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(patchwork)

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RDP/TaxaCorrelationHistograms.png"), width= 1200 * length(all_files), height=1200, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-TaxaCorHistogram(countsT, paste0("(RDP-", name, ") Taxa Correlation"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots)
dev.off()

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/dada2/TaxaCorrelationHistograms.png"), width=1200 * length(all_files), height=1200, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-TaxaCorHistogram(countsT, paste0("(dada2-", name, ") Taxa Correlation"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots)
dev.off()

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/WGS/TaxaCorrelationHistograms.png"), width=1200 * length(all_files), height=1200, res = 300)
par(mar=c(5,6,4,1)+.1)

combined_plots<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)

  p<-TaxaCorHistogram(countsT, paste0("(WGS-", name, ") Taxa Correlation"))

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }
}

print(combined_plots)
dev.off()
