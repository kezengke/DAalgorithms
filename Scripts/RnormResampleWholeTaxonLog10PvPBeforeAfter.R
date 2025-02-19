rm(list = ls())
library(MetagenomeTools)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(patchwork)
set.seed(9527)

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-LoadCountsT(file)
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/RDP/RnormWholeTaxonBeforeAfterLog10PvPPlots(", name, ").png"), width=3400, height=2400, res = 300)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  newT<-resampleWholeTaxonRNORM(countsT, meta, 1)

  Ottest<-read.table(paste0("PkgResults/RDP/ttest/", name, "_t.txt"), header = T, row.names = 1)
  ODESeq2<-read.table(paste0("PkgResults/RDP/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)

  newTnorm<-normFun(newT)
  Nttest<-calcTtest(newTnorm, meta)
  NDESeq2<-calcDESeq2(newT, meta)

  Ottest<-processTtestRes(Ottest)
  ODESeq2<-processDESeq2Res(ODESeq2)
  Nttest<-processTtestRes(Nttest)
  NDESeq2<-processDESeq2Res(NDESeq2)

  # OtvNt
  p1<-Log10PvPPlot(Ottest, Nttest, "t-test-Before", "t-test-After", paste0("(RDP-", name, ") t-test Before vs. t-test After"))
  # OtvOd
  p2<-Log10PvPPlot(ODESeq2, NDESeq2, "DESeq2-Before", "DESeq2-After", paste0("(RDP-", name, ") DESeq2 Before vs. DESeq2 After"))
  # OtvNd
  p3<-Log10PvPPlot(Ottest, ODESeq2, "t-test-Before", "DESeq2-Before", paste0("(RDP-", name, ") t-test Before vs. DESeq2 Before"))
  # NtvNd
  p4<-Log10PvPPlot(Nttest, NDESeq2, "t-test-After", "DESeq2-After", paste0("(RDP-", name, ") t-test After vs. DESeq2 After"))

  combined_plots <- p1 + p2 + p3 + p4
  print(combined_plots + plot_layout(ncol = 2))

  dev.off()
}

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-LoadCountsT(file)
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/dada2/RnormWholeTaxonBeforeAfterLog10PvPPlots(", name, ").png"), width=3400, height=2400, res = 300)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  newT<-resampleWholeTaxonRNORM(countsT, meta, 1)

  Ottest<-read.table(paste0("PkgResults/dada2/ttest/", name, "_t.txt"), header = T, row.names = 1)
  ODESeq2<-read.table(paste0("PkgResults/dada2/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)

  newTnorm<-normFun(newT)
  Nttest<-calcTtest(newTnorm, meta)
  NDESeq2<-calcDESeq2(newT, meta)

  Ottest<-processTtestRes(Ottest)
  ODESeq2<-processDESeq2Res(ODESeq2)
  Nttest<-processTtestRes(Nttest)
  NDESeq2<-processDESeq2Res(NDESeq2)

  # OtvNt
  p1<-Log10PvPPlot(Ottest, Nttest, "t-test-Before", "t-test-After", paste0("(dada2-", name, ") t-test Before vs. t-test After"))
  # OtvOd
  p2<-Log10PvPPlot(ODESeq2, NDESeq2, "DESeq2-Before", "DESeq2-After", paste0("(dada2-", name, ") DESeq2 Before vs. DESeq2 After"))
  # OtvNd
  p3<-Log10PvPPlot(Ottest, ODESeq2, "t-test-Before", "DESeq2-Before", paste0("(dada2-", name, ") t-test Before vs. DESeq2 Before"))
  # NtvNd
  p4<-Log10PvPPlot(Nttest, NDESeq2, "t-test-After", "DESeq2-After", paste0("(dada2-", name, ") t-test After vs. DESeq2 After"))

  combined_plots <- p1 + p2 + p3 + p4
  print(combined_plots + plot_layout(ncol = 2))

  dev.off()
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-LoadCountsT(file)
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/WGS/RnormWholeTaxonBeforeAfterLog10PvPPlots(", name, ").png"), width=3400, height=2400, res = 300)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  newT<-resampleWholeTaxonRNORM(countsT, meta, 1)

  Ottest<-read.table(paste0("PkgResults/WGS/ttest/", name, "_t.txt"), header = T, row.names = 1)
  ODESeq2<-read.table(paste0("PkgResults/WGS/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)

  newTnorm<-normFun(newT)
  Nttest<-calcTtest(newT, meta)
  NDESeq2<-calcDESeq2(newT, meta)

  Ottest<-processTtestRes(Ottest)
  ODESeq2<-processDESeq2Res(ODESeq2)
  Nttest<-processTtestRes(Nttest)
  NDESeq2<-processDESeq2Res(NDESeq2)

  # OtvNt
  p1<-Log10PvPPlot(Ottest, Nttest, "t-test-Before", "t-test-After", paste0("(WGS-", name, ") t-test Before vs. t-test After"))
  # OtvOd
  p2<-Log10PvPPlot(ODESeq2, NDESeq2, "DESeq2-Before", "DESeq2-After", paste0("(WGS-", name, ") DESeq2 Before vs. DESeq2 After"))
  # OtvNd
  p3<-Log10PvPPlot(Ottest, ODESeq2, "t-test-Before", "DESeq2-Before", paste0("(WGS-", name, ") t-test Before vs. DESeq2 Before"))
  # NtvNd
  p4<-Log10PvPPlot(Nttest, NDESeq2, "t-test-After", "DESeq2-After", paste0("(WGS-", name, ") t-test After vs. DESeq2 After"))

  combined_plots <- p1 + p2 + p3 + p4
  print(combined_plots + plot_layout(ncol = 2))

  dev.off()
}


