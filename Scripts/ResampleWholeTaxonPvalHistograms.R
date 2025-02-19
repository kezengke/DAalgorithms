#shuffle taxa counts and compare pvalues
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(coin)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(patchwork)

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/RDP/(", name, ")ResampleWholeTaxonMethodsPvalueHistograms.png"), width= 4800, height=1200, res = 300)
  RawcountsT<-read.table(paste0("CountsTables/RDPRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  newRawT<-resampleWholeTaxonRNORM(RawcountsT, meta, 1)
  newNormT<-normFun(newRawT)

  ttestRES<-calcTtest(newNormT, meta)
  wilcoxRES<-calcWilcox(newNormT, meta)
  deseq2RES<-calcDESeq2(newRawT, meta)
  edgerRES<-calcEdgeR(newRawT, meta)

  p1<-PvalHistogram(ttestRES, "red", paste0("(RDP-", name, ") t-test"))
  p2<-PvalHistogram(wilcoxRES, "tan2", paste0("(RDP-", name, ") Wilcoxon"))
  p3<-PvalHistogram(deseq2RES, "purple", paste0("(RDP-", name, ") DESeq2"))
  p4<-PvalHistogram(edgerRES, "cornflowerblue", paste0("(RDP-", name, ") edgeR"))

  combined_plots <- p1 + p2 + p3 + p4
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/dada2/(", name, ")ResampleWholeTaxonMethodsPvalueHistograms.png"), width= 4800, height=1200, res = 300)
  RawcountsT<-read.table(paste0("CountsTables/dada2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  newRawT<-resampleWholeTaxonRNORM(RawcountsT, meta, 1)
  newNormT<-normFun(newRawT)

  ttestRES<-calcTtest(newNormT, meta)
  wilcoxRES<-calcWilcox(newNormT, meta)
  deseq2RES<-calcDESeq2(newRawT, meta)
  edgerRES<-calcEdgeR(newRawT, meta)

  p1<-PvalHistogram(ttestRES, "red", paste0("(dada2-", name, ") t-test"))
  p2<-PvalHistogram(wilcoxRES, "tan2", paste0("(dada2-", name, ") Wilcoxon"))
  p3<-PvalHistogram(deseq2RES, "purple", paste0("(dada2-", name, ") DESeq2"))
  p4<-PvalHistogram(edgerRES, "cornflowerblue", paste0("(dada2-", name, ") edgeR"))

  combined_plots <- p1 + p2 + p3 + p4
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/WGS/(", name, ")ResampleWholeTaxonMethodsPvalueHistograms.png"), width= 4800, height=1200, res = 300)
  RawcountsT<-read.table(paste0("CountsTables/WGSRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  newRawT<-resampleWholeTaxonRNORM(RawcountsT, meta, 1)
  newNormT<-normFun(newRawT)

  ttestRES<-calcTtest(newNormT, meta)
  wilcoxRES<-calcWilcox(newNormT, meta)
  deseq2RES<-calcDESeq2(newRawT, meta)
  edgerRES<-calcEdgeR(newRawT, meta)

  p1<-PvalHistogram(ttestRES, "red", paste0("(WGS-", name, ") t-test"))
  p2<-PvalHistogram(wilcoxRES, "tan2", paste0("(WGS-", name, ") Wilcoxon"))
  p3<-PvalHistogram(deseq2RES, "purple", paste0("(WGS-", name, ") DESeq2"))
  p4<-PvalHistogram(edgerRES, "cornflowerblue", paste0("(WGS-", name, ") edgeR"))

  combined_plots <- p1 + p2 + p3 + p4
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}


