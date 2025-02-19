#shuffle taxa counts and make taxa correlation histograms
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(ggplot2)
library(patchwork)


# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RDP/ShuffledTaxaCorrelationHistograms.png"), width= 1200 * length(all_files), height=1200, res = 300)

combined_plots<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/RDPRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  newCountsT<-apply(RawcountsT, 2, sample)

  rownames(newCountsT)<-rownames(RawcountsT)
  colnames(newCountsT)<-colnames(RawcountsT)
  newCountsT<-data.frame(newCountsT)

  p<-TaxaCorHistogram(newCountsT, paste0("(RDP-", name, ") Shuffled Taxa Correlation"), "lightseagreen")

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

png(paste0("Plots/dada2/ShuffledTaxaCorrelationHistograms.png"), width= 1200 * length(all_files), height=1200, res = 300)

combined_plots<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/dada2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  newCountsT<-apply(RawcountsT, 2, sample)

  rownames(newCountsT)<-rownames(RawcountsT)
  colnames(newCountsT)<-colnames(RawcountsT)
  newCountsT<-data.frame(newCountsT)

  p<-TaxaCorHistogram(newCountsT, paste0("(dada2-", name, ") Shuffled Taxa Correlation"), "lightseagreen")

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

png(paste0("Plots/WGS/ShuffledTaxaCorrelationHistograms.png"), width= 1200 * length(all_files), height=1200, res = 300)

combined_plots<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/WGSRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  newCountsT<-apply(RawcountsT, 2, sample)

  rownames(newCountsT)<-rownames(RawcountsT)
  colnames(newCountsT)<-colnames(RawcountsT)
  newCountsT<-data.frame(newCountsT)

  p<-TaxaCorHistogram(newCountsT, paste0("(WGS-", name, ") Shuffled Taxa Correlation"), "lightseagreen")

  if (is.null(combined_plots)) {
    combined_plots <- p
  } else {
    combined_plots <- combined_plots + p
  }

}

print(combined_plots)
dev.off()
