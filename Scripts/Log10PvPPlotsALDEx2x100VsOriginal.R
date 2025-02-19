rm(list = ls())
set.seed(9527)
library(ggplot2)
library(patchwork)
library(ALDEx2)
library(dplyr)
library(MetagenomeTools)

# RDP
all_names<-gsub(basename(list.files("CountsTables/RDPRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for (name in all_names) {

  Tpvals<-read.table(paste0("ALDEx2RerunDump/RDP/ALDEx2ttest/", name, "_ALDEx2T.txt"))
  Wpvals<-read.table(paste0("ALDEx2RerunDump/RDP/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"))

  avgTPvals<-apply(Tpvals, 1, mean)
  avgWPvals<-apply(Wpvals, 1, mean)

  originalTRes<-read.table(paste0("PkgResults/RDP/ttest/", name, "_t.txt"), row.names = 1, header = T, sep = "\t")
  avgTRes<-data.frame(originalTRes$stats, avgTPvals)
  colnames(avgTRes)<-colnames(originalTRes)

  originalWRes<-read.table(paste0("PkgResults/RDP/Wilcoxon/", name, "_wilcox.txt"), row.names = 1, header = T, sep = "\t")
  avgWRes<-data.frame(originalWRes$stats, avgWPvals)
  colnames(avgWRes)<-colnames(originalWRes)

  oTRES<-processTtestRes(originalTRes)
  aTRES<-processALDEx2Res(avgTRes)
  aTRES$p_directed<-aTRES$p_directed * (-1)

  oWRES<-processWilcoxonRes(originalWRes)
  aWRES<-processALDEx2Res(avgWRes)
  aWRES$p_directed<-aWRES$p_directed * (-1)

  p1<-Log10PvPPlot(oTRES, aTRES, "t-test", "ALDEx2t-test", paste0("(RDP-", name, ") \nOriginal t-test vs. 100x avg. ALDEx2 t-test"))
  p2<-Log10PvPPlot(oWRES, aWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") \nOriginal Wilcoxon vs. 100x avg. ALDEx2 Wilcoxon"))

  png(paste0("Plots/RDP/100xALDEx2OriginalTtestWilcoxonLog10PvPPlots(", name, ").png"), width= 3400, height=1200, res = 300)

  combined_plots <- p1 + p2
  print(combined_plots + plot_layout(ncol = 2))

  dev.off()

}

# dada2
all_names<-gsub(basename(list.files("CountsTables/dada2Raw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for (name in all_names) {

  Tpvals<-read.table(paste0("ALDEx2RerunDump/dada2/ALDEx2ttest/", name, "_ALDEx2T.txt"))
  Wpvals<-read.table(paste0("ALDEx2RerunDump/dada2/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"))

  avgTPvals<-apply(Tpvals, 1, mean)
  avgWPvals<-apply(Wpvals, 1, mean)

  originalTRes<-read.table(paste0("PkgResults/dada2/ttest/", name, "_t.txt"), row.names = 1, header = T, sep = "\t")
  avgTRes<-data.frame(originalTRes$stats, avgTPvals)
  colnames(avgTRes)<-colnames(originalTRes)

  originalWRes<-read.table(paste0("PkgResults/dada2/Wilcoxon/", name, "_wilcox.txt"), row.names = 1, header = T, sep = "\t")
  avgWRes<-data.frame(originalWRes$stats, avgWPvals)
  colnames(avgWRes)<-colnames(originalWRes)

  oTRES<-processTtestRes(originalTRes)
  aTRES<-processALDEx2Res(avgTRes)
  aTRES$p_directed<-aTRES$p_directed * (-1)

  oWRES<-processWilcoxonRes(originalWRes)
  aWRES<-processALDEx2Res(avgWRes)
  aWRES$p_directed<-aWRES$p_directed * (-1)

  p1<-Log10PvPPlot(oTRES, aTRES, "t-test", "ALDEx2t-test", paste0("(dada2-", name, ") \nOriginal t-test vs. 100x avg. ALDEx2 t-test"))
  p2<-Log10PvPPlot(oWRES, aWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(dada2-", name, ") \nOriginal Wilcoxon vs. 100x avg. ALDEx2 Wilcoxon"))

  png(paste0("Plots/dada2/100xALDEx2OriginalTtestWilcoxonLog10PvPPlots(", name, ").png"), width= 3400, height=1200, res = 300)

  combined_plots <- p1 + p2
  print(combined_plots + plot_layout(ncol = 2))

  dev.off()

}

# WGS
all_names<-gsub(basename(list.files("CountsTables/WGSRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for (name in all_names) {

  Tpvals<-read.table(paste0("ALDEx2RerunDump/WGS/ALDEx2ttest/", name, "_ALDEx2T.txt"))
  Wpvals<-read.table(paste0("ALDEx2RerunDump/WGS/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"))

  avgTPvals<-apply(Tpvals, 1, mean)
  avgWPvals<-apply(Wpvals, 1, mean)

  originalTRes<-read.table(paste0("PkgResults/WGS/ttest/", name, "_t.txt"), row.names = 1, header = T, sep = "\t")
  avgTRes<-data.frame(originalTRes$stats, avgTPvals)
  colnames(avgTRes)<-colnames(originalTRes)

  originalWRes<-read.table(paste0("PkgResults/WGS/Wilcoxon/", name, "_wilcox.txt"), row.names = 1, header = T, sep = "\t")
  avgWRes<-data.frame(originalWRes$stats, avgWPvals)
  colnames(avgWRes)<-colnames(originalWRes)

  oTRES<-processTtestRes(originalTRes)
  aTRES<-processALDEx2Res(avgTRes)
  aTRES$p_directed<-aTRES$p_directed * (-1)

  oWRES<-processWilcoxonRes(originalWRes)
  aWRES<-processALDEx2Res(avgWRes)
  aWRES$p_directed<-aWRES$p_directed * (-1)

  p1<-Log10PvPPlot(oTRES, aTRES, "t-test", "ALDEx2t-test", paste0("(WGS-", name, ") \nOriginal t-test vs. 100x avg. ALDEx2 t-test"))
  p2<-Log10PvPPlot(oWRES, aWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") \nOriginal Wilcoxon vs. 100x avg. ALDEx2 Wilcoxon"))

  png(paste0("Plots/WGS/100xALDEx2OriginalTtestWilcoxonLog10PvPPlots(", name, ").png"), width= 3400, height=1200, res = 300)

  combined_plots <- p1 + p2
  print(combined_plots + plot_layout(ncol = 2))

  dev.off()

}

# RNAseq
all_names<-gsub(basename(list.files("CountsTables/RNAseqRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for (name in all_names) {

  Tpvals<-read.table(paste0("ALDEx2RerunDump/RNAseq/ALDEx2ttest/", name, "_ALDEx2T.txt"))
  Wpvals<-read.table(paste0("ALDEx2RerunDump/RNAseq/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"))

  avgTPvals<-apply(Tpvals, 1, mean)
  avgWPvals<-apply(Wpvals, 1, mean)

  originalTRes<-read.table(paste0("PkgResults/RNAseq/ttest/", name, "_t.txt"), row.names = 1, header = T, sep = "\t")
  avgTRes<-data.frame(originalTRes$stats, avgTPvals)
  colnames(avgTRes)<-colnames(originalTRes)

  originalWRes<-read.table(paste0("PkgResults/RNAseq/Wilcoxon/", name, "_wilcox.txt"), row.names = 1, header = T, sep = "\t")
  avgWRes<-data.frame(originalWRes$stats, avgWPvals)
  colnames(avgWRes)<-colnames(originalWRes)

  oTRES<-processTtestRes(originalTRes)
  aTRES<-processALDEx2Res(avgTRes)
  aTRES$p_directed<-aTRES$p_directed * (-1)

  oWRES<-processWilcoxonRes(originalWRes)
  aWRES<-processALDEx2Res(avgWRes)
  aWRES$p_directed<-aWRES$p_directed * (-1)

  p1<-Log10PvPPlot(oTRES, aTRES, "t-test", "ALDEx2t-test", paste0("(RNAseq-", name, ") \nOriginal t-test vs. 100x avg. ALDEx2 t-test"))
  p2<-Log10PvPPlot(oWRES, aWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") \nOriginal Wilcoxon vs. 100x avg. ALDEx2 Wilcoxon"))

  png(paste0("Plots/RNAseq/100xALDEx2OriginalTtestWilcoxonLog10PvPPlots(", name, ").png"), width= 3400, height=1200, res = 300)

  combined_plots <- p1 + p2
  print(combined_plots + plot_layout(ncol = 2))

  dev.off()

}
