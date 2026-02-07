# Plot log10 pvalue plots
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(dplyr)
library(patchwork)

# RDP
all_names<-gsub(basename(list.files("CountsTables/RDPRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RDP/RestLog10PvPPlots(", name, ").png"), width=6700, height=4800, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RDP/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RDP/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RDP/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RDP/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/RDP/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/RDP/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/RDP/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)
  mRES<-read.table(paste0("PkgResults/RDP/metagenomeseq/", name, "_metagenomeseq.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)
  mRES<-processMetagenomeSeqRes(mRES)


  # t vs all
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(RDP-", name, ") t-test vs. ALDEx2t-test"))
  p2<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") t-test vs. ALDEx2Wilcoxon"))
  p3<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(RDP-", name, ") t-test vs. ANCOMBC2"))
  p4<-Log10PvPPlot(tRES, mRES, "t-test", "metagenomeSeq", paste0("(RDP-", name, ") t-test vs. metagenomeSeq"))

  # w vs all
  p5<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(RDP-", name, ") Wilcoxon vs. ALDEx2t-test"))
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  p7<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(RDP-", name, ") Wilcoxon vs. ANCOMBC2"))
  p8<-Log10PvPPlot(wRES, mRES, "Wilcoxon", "metagenomeSeq", paste0("(RDP-", name, ") Wilcoxon vs. metagenomeSeq"))

  # d vs all
  p9<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(RDP-", name, ") DESeq2 vs. ALDEx2t-test"))
  p10<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  p11<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(RDP-", name, ") DESeq2 vs. ANCOMBC2"))
  p12<-Log10PvPPlot(dRES, mRES, "DESeq2", "metagenomeSeq", paste0("(RDP-", name, ") DESeq2 vs. metagenomeSeq"))

  # e vs all
  p13<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(RDP-", name, ") edgeR vs. ALDEx2t-test"))
  p14<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(RDP-", name, ") edgeR vs. ALDEx2Wilcoxon"))
  p15<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(RDP-", name, ") edgeR vs. ANCOMBC2"))
  p16<-Log10PvPPlot(eRES, mRES, "edgeR", "metagenomeSeq", paste0("(RDP-", name, ") edgeR vs. metagenomeSeq"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + p14 + p15 + p16
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# dada2
all_names<-gsub(basename(list.files("CountsTables/dada2Raw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/dada2/RestLog10PvPPlots(", name, ").png"), width=6700, height=4800, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/dada2/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/dada2/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/dada2/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/dada2/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/dada2/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/dada2/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/dada2/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)
  mRES<-read.table(paste0("PkgResults/dada2/metagenomeseq/", name, "_metagenomeseq.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)
  mRES<-processMetagenomeSeqRes(mRES)


  # t vs all
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(dada2-", name, ") t-test vs. ALDEx2t-test"))
  p2<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(dada2-", name, ") t-test vs. ALDEx2Wilcoxon"))
  p3<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(dada2-", name, ") t-test vs. ANCOMBC2"))
  p4<-Log10PvPPlot(tRES, mRES, "t-test", "metagenomeSeq", paste0("(dada2-", name, ") t-test vs. metagenomeSeq"))

  # w vs all
  p5<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(dada2-", name, ") Wilcoxon vs. ALDEx2t-test"))
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(dada2-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  p7<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(dada2-", name, ") Wilcoxon vs. ANCOMBC2"))
  p8<-Log10PvPPlot(wRES, mRES, "Wilcoxon", "metagenomeSeq", paste0("(dada2-", name, ") Wilcoxon vs. metagenomeSeq"))

  # d vs all
  p9<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(dada2-", name, ") DESeq2 vs. ALDEx2t-test"))
  p10<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(dada2-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  p11<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(dada2-", name, ") DESeq2 vs. ANCOMBC2"))
  p12<-Log10PvPPlot(dRES, mRES, "DESeq2", "metagenomeSeq", paste0("(dada2-", name, ") DESeq2 vs. metagenomeSeq"))

  # e vs all
  p13<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(dada2-", name, ") edgeR vs. ALDEx2t-test"))
  p14<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(dada2-", name, ") edgeR vs. ALDEx2Wilcoxon"))
  p15<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(dada2-", name, ") edgeR vs. ANCOMBC2"))
  p16<-Log10PvPPlot(eRES, mRES, "edgeR", "metagenomeSeq", paste0("(dada2-", name, ") edgeR vs. metagenomeSeq"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + p14 + p15 + p16
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# WGS
all_names<-gsub(basename(list.files("CountsTables/WGSRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/WGS/RestLog10PvPPlots(", name, ").png"), width=6700, height=4800, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/WGS/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/WGS/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/WGS/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/WGS/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/WGS/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/WGS/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/WGS/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)
  mRES<-read.table(paste0("PkgResults/WGS/metagenomeseq/", name, "_metagenomeseq.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)
  mRES<-processMetagenomeSeqRes(mRES)


  # t vs all
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(WGS-", name, ") t-test vs. ALDEx2t-test"))
  p2<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") t-test vs. ALDEx2Wilcoxon"))
  p3<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(WGS-", name, ") t-test vs. ANCOMBC2"))
  p4<-Log10PvPPlot(tRES, mRES, "t-test", "metagenomeSeq", paste0("(WGS-", name, ") t-test vs. metagenomeSeq"))

  # w vs all
  p5<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(WGS-", name, ") Wilcoxon vs. ALDEx2t-test"))
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  p7<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(WGS-", name, ") Wilcoxon vs. ANCOMBC2"))
  p8<-Log10PvPPlot(wRES, mRES, "Wilcoxon", "metagenomeSeq", paste0("(WGS-", name, ") Wilcoxon vs. metagenomeSeq"))

  # d vs all
  p9<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(WGS-", name, ") DESeq2 vs. ALDEx2t-test"))
  p10<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  p11<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(WGS-", name, ") DESeq2 vs. ANCOMBC2"))
  p12<-Log10PvPPlot(dRES, mRES, "DESeq2", "metagenomeSeq", paste0("(WGS-", name, ") DESeq2 vs. metagenomeSeq"))

  # e vs all
  p13<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(WGS-", name, ") edgeR vs. ALDEx2t-test"))
  p14<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(WGS-", name, ") edgeR vs. ALDEx2Wilcoxon"))
  p15<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(WGS-", name, ") edgeR vs. ANCOMBC2"))
  p16<-Log10PvPPlot(eRES, mRES, "edgeR", "metagenomeSeq", paste0("(WGS-", name, ") edgeR vs. metagenomeSeq"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + p14 + p15 + p16
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}

# RNAseq
all_names<-gsub(basename(list.files("CountsTables/RNAseqRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RNAseq/RestLog10PvPPlots(", name, ").png"), width=6700, height=4800, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RNAseq/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RNAseq/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RNAseq/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RNAseq/edgeR/", name, "_edger.txt"), header = T, row.names = 1)
  xTRES<-read.table(paste0("PkgResults/RNAseq/ALDEx2ttest/", name, "_ALDEx2T.txt"), header = T, row.names = 1)
  xWRES<-read.table(paste0("PkgResults/RNAseq/ALDEx2Wilcoxon/", name, "_ALDEx2W.txt"), header = T, row.names = 1)
  aRES<-read.table(paste0("PkgResults/RNAseq/ancombc/", name, "_ancombc2.txt"), header = T, row.names = 1)
  mRES<-read.table(paste0("PkgResults/RNAseq/metagenomeseq/", name, "_metagenomeseq.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)
  xTRES<-processALDEx2Res(xTRES)
  xWRES<-processALDEx2Res(xWRES)
  aRES<-processANCOMBC2Res(aRES)
  mRES<-processMetagenomeSeqRes(mRES)


  # t vs all
  p1<-Log10PvPPlot(tRES, xTRES, "t-test", "ALDEx2t-test", paste0("(RNAseq-", name, ") t-test vs. ALDEx2t-test"))
  p2<-Log10PvPPlot(tRES, xWRES, "t-test", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") t-test vs. ALDEx2Wilcoxon"))
  p3<-Log10PvPPlot(tRES, aRES, "t-test", "ANCOMBC2", paste0("(RNAseq-", name, ") t-test vs. ANCOMBC2"))
  p4<-Log10PvPPlot(tRES, mRES, "t-test", "metagenomeSeq", paste0("(RNAseq-", name, ") t-test vs. metagenomeSeq"))

  # w vs all
  p5<-Log10PvPPlot(wRES, xTRES, "Wilcoxon", "ALDEx2t-test", paste0("(RNAseq-", name, ") Wilcoxon vs. ALDEx2t-test"))
  p6<-Log10PvPPlot(wRES, xWRES, "Wilcoxon", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") Wilcoxon vs. ALDEx2Wilcoxon"))
  p7<-Log10PvPPlot(wRES, aRES, "Wilcoxon", "ANCOMBC2", paste0("(RNAseq-", name, ") Wilcoxon vs. ANCOMBC2"))
  p8<-Log10PvPPlot(wRES, mRES, "Wilcoxon", "metagenomeSeq", paste0("(RNAseq-", name, ") Wilcoxon vs. metagenomeSeq"))

  # d vs all
  p9<-Log10PvPPlot(dRES, xTRES, "DESeq2", "ALDEx2t-test", paste0("(RNAseq-", name, ") DESeq2 vs. ALDEx2t-test"))
  p10<-Log10PvPPlot(dRES, xWRES, "DESeq2", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") DESeq2 vs. ALDEx2Wilcoxon"))
  p11<-Log10PvPPlot(dRES, aRES, "DESeq2", "ANCOMBC2", paste0("(RNAseq-", name, ") DESeq2 vs. ANCOMBC2"))
  p12<-Log10PvPPlot(dRES, mRES, "DESeq2", "metagenomeSeq", paste0("(RNAseq-", name, ") DESeq2 vs. metagenomeSeq"))

  # e vs all
  p13<-Log10PvPPlot(eRES, xTRES, "edgeR", "ALDEx2t-test", paste0("(RNAseq-", name, ") edgeR vs. ALDEx2t-test"))
  p14<-Log10PvPPlot(eRES, xWRES, "edgeR", "ALDEx2Wilcoxon", paste0("(RNAseq-", name, ") edgeR vs. ALDEx2Wilcoxon"))
  p15<-Log10PvPPlot(eRES, aRES, "edgeR", "ANCOMBC2", paste0("(RNAseq-", name, ") edgeR vs. ANCOMBC2"))
  p16<-Log10PvPPlot(eRES, mRES, "edgeR", "metagenomeSeq", paste0("(RNAseq-", name, ") edgeR vs. metagenomeSeq"))


  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + p14 + p15 + p16
  print(combined_plots + plot_layout(ncol = 4))

  dev.off()
}
