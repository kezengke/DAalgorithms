# Plot log10 pvalue plots
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(dplyr)
library(patchwork)

# RDP
all_names<-gsub(basename(list.files("CountsTables/RDPRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RDP/ANCOMBC2Log10PvPPlots(", name, ").png"), width=5025, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RDP/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RDP/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RDP/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RDP/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/RDP/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/RDP/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/RDP/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)


  # tvx
  p1<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(RDP-", name, ") t-test vs. ANCOMBC2"))
  # wvx
  p2<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(RDP-", name, ") Wilcoxon vs. ANCOMBC2"))
  # dvx
  p3<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(RDP-", name, ") DESeq2 vs. ANCOMBC2"))
  # evx
  p4<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(RDP-", name, ") edgeR vs. ANCOMBC2"))


  # xva
  p5<-Log10PvPPlot(xTRES, aRES, "ALDEx2t-test", "ANCOMBC2", paste0("(RDP-", name, ") ALDEx2t-test vs. ANCOMBC2"))
  # xva
  p6<-Log10PvPPlot(xWRES, aRES, "ALDEx2Wilcoxon", "ANCOMBC2", paste0("(RDP-", name, ") ALDEX2Wilcoxon vs. ANCOMBC2"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}

# dada2
all_names<-gsub(basename(list.files("CountsTables/dada2Raw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/dada2/ANCOMBC2Log10PvPPlots(", name, ").png"), width=5025, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/dada2/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/dada2/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/dada2/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/dada2/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/dada2/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/dada2/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/dada2/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)


  # tvx
  p1<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(dada2-", name, ") t-test vs. ANCOMBC2"))
  # wvx
  p2<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(dada2-", name, ") Wilcoxon vs. ANCOMBC2"))
  # dvx
  p3<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(dada2-", name, ") DESeq2 vs. ANCOMBC2"))
  # evx
  p4<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(dada2-", name, ") edgeR vs. ANCOMBC2"))


  # xva
  p5<-Log10PvPPlot(xTRES, aRES, "ALDEx2t-test", "ANCOMBC2", paste0("(dada2-", name, ") ALDEx2t-test vs. ANCOMBC2"))
  # xva
  p6<-Log10PvPPlot(xWRES, aRES, "ALDEx2Wilcoxon", "ANCOMBC2", paste0("(dada2-", name, ") ALDEX2Wilcoxon vs. ANCOMBC2"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}

# WGS
all_names<-gsub(basename(list.files("CountsTables/WGSRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/WGS/ANCOMBC2Log10PvPPlots(", name, ").png"), width=5025, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/WGS/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/WGS/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/WGS/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/WGS/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/WGS/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/WGS/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/WGS/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)


  # tvx
  p1<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(WGS-", name, ") t-test vs. ANCOMBC2"))
  # wvx
  p2<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(WGS-", name, ") Wilcoxon vs. ANCOMBC2"))
  # dvx
  p3<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(WGS-", name, ") DESeq2 vs. ANCOMBC2"))
  # evx
  p4<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(WGS-", name, ") edgeR vs. ANCOMBC2"))


  # xva
  p5<-Log10PvPPlot(xTRES, aRES, "ALDEx2t-test", "ANCOMBC2", paste0("(WGS-", name, ") ALDEx2t-test vs. ANCOMBC2"))
  # xva
  p6<-Log10PvPPlot(xWRES, aRES, "ALDEx2Wilcoxon", "ANCOMBC2", paste0("(WGS-", name, ") ALDEX2Wilcoxon vs. ANCOMBC2"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}

# RNAseq
all_names<-gsub(basename(list.files("CountsTables/RNAseqRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RNAseq/ANCOMBC2Log10PvPPlots(", name, ").png"), width=5025, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RNAseq/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RNAseq/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RNAseq/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RNAseq/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/RNAseq/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/RNAseq/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/RNAseq/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)


  # tvx
  p1<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(RNAseq-", name, ") t-test vs. ANCOMBC2"))
  # wvx
  p2<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(RNAseq-", name, ") Wilcoxon vs. ANCOMBC2"))
  # dvx
  p3<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(RNAseq-", name, ") DESeq2 vs. ANCOMBC2"))
  # evx
  p4<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(RNAseq-", name, ") edgeR vs. ANCOMBC2"))


  # xva
  p5<-Log10PvPPlot(xTRES, aRES, "ALDEx2t-test", "ANCOMBC2", paste0("(RNAseq-", name, ") ALDEx2t-test vs. ANCOMBC2"))
  # xva
  p6<-Log10PvPPlot(xWRES, aRES, "ALDEx2Wilcoxon", "ANCOMBC2", paste0("(RNAseq-", name, ") ALDEX2Wilcoxon vs. ANCOMBC2"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}
