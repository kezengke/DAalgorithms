# Plot log10 pvalue plots for shuffled results
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(dplyr)
library(patchwork)

# RDP
all_names<-gsub(basename(list.files("CountsTables/RDPRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RDP/ShuffledLog10PvPPlots(", name, ").png"), width=5000, height=2400, res = 300)
  # par(mfrow=c(2,3))
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RDP/shuffledttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RDP/shuffledWilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RDP/shuffledDESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RDP/shufflededgeR/", name, "_edger.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)

  # tvw
  p1<-Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", paste0("(RDP-", name, ") shuffled t-test vs. Wilcoxon"))
  # tvd
  p2<-Log10PvPPlot(tRES, dRES, "t-test", "DESeq2", paste0("(RDP-", name, ") shuffled t-test vs. DESeq2"))
  # tve
  p3<-Log10PvPPlot(tRES, eRES, "t-test", "edgeR", paste0("(RDP-", name, ") shuffled t-test vs. edgeR"))

  # wvd
  p4<-Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2", paste0("(RDP-", name, ") shuffled Wilcoxon vs. DESeq2"))
  # wve
  p5<-Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR", paste0("(RDP-", name, ") shuffled Wilcoxon vs. edgeR"))

  # dve
  p6<-Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR", paste0("(RDP-", name, ") shuffled DESeq2 vs. edgeR"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}


# DADA2
all_names<-gsub(basename(list.files("CountsTables/dada2Raw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/dada2/ShuffledLog10PvPPlots(", name, ").png"), width=5000, height=2400, res = 300)
  # par(mfrow=c(2,3))
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/dada2/shuffledttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/dada2/shuffledWilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/dada2/shuffledDESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/dada2/shufflededgeR/", name, "_edger.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)

  # tvw
  p1<-Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", paste0("(dada2-", name, ") shuffled t-test vs. Wilcoxon"))
  # tvd
  p2<-Log10PvPPlot(tRES, dRES, "t-test", "DESeq2", paste0("(dada2-", name, ") shuffled t-test vs. DESeq2"))
  # tve
  p3<-Log10PvPPlot(tRES, eRES, "t-test", "edgeR", paste0("(dada2-", name, ") shuffled t-test vs. edgeR"))

  # wvd
  p4<-Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2", paste0("(dada2-", name, ") shuffled Wilcoxon vs. DESeq2"))
  # wve
  p5<-Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR", paste0("(dada2-", name, ") shuffled Wilcoxon vs. edgeR"))

  # dve
  p6<-Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR", paste0("(dada2-", name, ") shuffled DESeq2 vs. edgeR"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}

# WGS
all_names<-gsub(basename(list.files("CountsTables/WGSRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/WGS/ShuffledLog10PvPPlots(", name, ").png"), width=5000, height=2400, res = 300)
  # par(mfrow=c(2,3))
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/WGS/shuffledttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/WGS/shuffledWilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/WGS/shuffledDESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/WGS/shufflededgeR/", name, "_edger.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)

  # tvw
  p1<-Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", paste0("(WGS-", name, ") shuffled t-test vs. Wilcoxon"))
  # tvd
  p2<-Log10PvPPlot(tRES, dRES, "t-test", "DESeq2", paste0("(WGS-", name, ") shuffled t-test vs. DESeq2"))
  # tve
  p3<-Log10PvPPlot(tRES, eRES, "t-test", "edgeR", paste0("(WGS-", name, ") shuffled t-test vs. edgeR"))

  # wvd
  p4<-Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2", paste0("(WGS-", name, ") shuffled Wilcoxon vs. DESeq2"))
  # wve
  p5<-Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR", paste0("(WGS-", name, ") shuffled Wilcoxon vs. edgeR"))

  # dve
  p6<-Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR", paste0("(WGS-", name, ") shuffled DESeq2 vs. edgeR"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}
