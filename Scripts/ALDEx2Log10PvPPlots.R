# Plot log10 pvalue plots
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(dplyr)
library(patchwork)

# RDP
all_names<-gsub(basename(list.files("CountsTables/RDPRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RDP/ALDEx2Log10PvPPlots(", name, ").png"), width=6700, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RDP/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RDP/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RDP/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RDP/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/RDP/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/RDP/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)


  # tvx
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(RDP-", name, ") t-test vs. ALDEx2t-test"))
  # wvx
  p2<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(RDP-", name, ") Wilcoxon vs. ALDEx2t-test"))
  # dvx
  p3<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(RDP-", name, ") DESeq2 vs. ALDEx2t-test"))
  # evx
  p4<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(RDP-", name, ") edgeR vs. ALDEx2t-test"))


  # tvx
  p5<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") t-test vs. ALDEx2Wilcoxon"))
  # wvx
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  # dvx
  p7<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  # evx
  p8<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") edgeR vs. ALDEx2Wilcoxon"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# DADA2
all_names<-gsub(basename(list.files("CountsTables/DADA2Raw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/DADA2/ALDEx2Log10PvPPlots(", name, ").png"), width=6700, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/DADA2/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/DADA2/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/DADA2/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/DADA2/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/DADA2/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/DADA2/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)


  # tvx
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(DADA2-", name, ") t-test vs. ALDEx2t-test"))
  # wvx
  p2<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(DADA2-", name, ") Wilcoxon vs. ALDEx2t-test"))
  # dvx
  p3<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(DADA2-", name, ") DESeq2 vs. ALDEx2t-test"))
  # evx
  p4<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(DADA2-", name, ") edgeR vs. ALDEx2t-test"))

  # tvx
  p5<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(DADA2-", name, ") t-test vs. ALDEx2Wilcoxon"))
  # wvx
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(DADA2-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  # dvx
  p7<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(DADA2-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  # evx
  p8<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(DADA2-", name, ") edgeR vs. ALDEx2Wilcoxon"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# WGS
all_names<-gsub(basename(list.files("CountsTables/WGSRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/WGS/ALDEx2Log10PvPPlots(", name, ").png"), width=6700, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/WGS/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/WGS/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/WGS/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/WGS/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/WGS/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/WGS/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)


  # tvx
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(WGS-", name, ") t-test vs. ALDEx2t-test"))
  # wvx
  p2<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(WGS-", name, ") Wilcoxon vs. ALDEx2t-test"))
  # dvx
  p3<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(WGS-", name, ") DESeq2 vs. ALDEx2t-test"))
  # evx
  p4<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(WGS-", name, ") edgeR vs. ALDEx2t-test"))

  # tvx
  p5<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") t-test vs. ALDEx2Wilcoxon"))
  # wvx
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  # dvx
  p7<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  # evx
  p8<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") edgeR vs. ALDEx2Wilcoxon"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# RNAseq
all_names<-gsub(basename(list.files("CountsTables/RNAseqRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RNAseq/ALDEx2Log10PvPPlots(", name, ").png"), width=6700, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RNAseq/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RNAseq/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RNAseq/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RNAseq/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/RNAseq/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/RNAseq/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)


  # tvx
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(RNAseq-", name, ") t-test vs. ALDEx2t-test"))
  # wvx
  p2<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(RNAseq-", name, ") Wilcoxon vs. ALDEx2t-test"))
  # dvx
  p3<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(RNAseq-", name, ") DESeq2 vs. ALDEx2t-test"))
  # evx
  p4<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(RNAseq-", name, ") edgeR vs. ALDEx2t-test"))

  # tvx
  p5<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") t-test vs. ALDEx2Wilcoxon"))
  # wvx
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  # dvx
  p7<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  # evx
  p8<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") edgeR vs. ALDEx2Wilcoxon"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

