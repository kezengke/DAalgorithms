rm(list = ls())
set.seed(9527)
library(ggplot2)
library(patchwork)
library(ALDEx2)
library(dplyr)
library(MetagenomeTools)

rawT<-read.table("CountsTables/RDPRaw/rhizo.txt", header = T, row.names = 1, sep = "\t")
meta<-read.table("MetaData/metadata_rhizo.txt", header = T, row.names = 1, sep = "\t")

rawT<-rawT[, intersect(colnames(rawT), rownames(meta)), drop = F]
meta<-meta[intersect(colnames(rawT), rownames(meta)), , drop = F]

rownames(meta)<-colnames(rawT)
colnames(meta)<-"conditions"

allT_pvals<-vector()
allW_pvals<-vector()
for (i in 1:100) {
  ALDEx2RES<-calcALDEx2Ttest(rawT, meta)
  pval<-ALDEx2RES[, 2, drop = F]
  allT_pvals[i]<-pval

  ALDEx2RES<-calcALDEx2Wilcoxon(rawT, meta)
  pval<-ALDEx2RES[, 2, drop = F]
  allW_pvals[i]<-pval
}

Tpvals<-do.call(cbind, allT_pvals)
Tpvals<-data.frame(Tpvals)
rownames(Tpvals)<-rownames(rawT)
avgTPvals<-apply(Tpvals, 1, mean)

Wpvals<-do.call(cbind, allW_pvals)
Wpvals<-data.frame(Wpvals)
rownames(Wpvals)<-rownames(rawT)
avgWPvals<-apply(Wpvals, 1, mean)

write.table(Tpvals, "ALDEx2RerunDump/Rhizo_ALDEx2_t.txt", sep = "\t")
write.table(Wpvals, "ALDEx2RerunDump/Rhizo_ALDEx2_wilcoxon.txt", sep = "\t")

originalTRes<-read.table("PkgResults/RDP/ttest/rhizo_t.txt", row.names = 1, header = T, sep = "\t")
avgTRes<-data.frame(originalTRes$stats, avgTPvals)
colnames(avgTRes)<-colnames(originalTRes)

originalWRes<-read.table("PkgResults/RDP/Wilcoxon/rhizo_wilcox.txt", row.names = 1, header = T, sep = "\t")
avgWRes<-data.frame(originalWRes$stats, avgWPvals)
colnames(avgWRes)<-colnames(originalWRes)

oTRES<-processTtestRes(originalTRes)
aTRES<-processALDEx2Res(avgTRes)
aTRES$p_directed<-aTRES$p_directed * (-1)

oWRES<-processWilcoxonRes(originalWRes)
aWRES<-processALDEx2Res(avgWRes)
aWRES$p_directed<-aWRES$p_directed * (-1)

sortedAT<-aTRES[order((aTRES$p_directed)), , drop = F]


TopSeven<-c(rownames(sortedAT)[1:6] , rownames(sortedAT)[nrow(sortedAT)])

png("Plots/RDP/TopSeven100xALDEx2PvalHistogram.png", width= 3000, height=1200 * (length(TopSeven)), res = 300)
combined_plots<-NULL
for (taxon in TopSeven) {
  all_Tpvals<-as.numeric(Tpvals[taxon, ])
  all_Wpvals<-as.numeric(Wpvals[taxon, ])

  norm_Tpvals<-log10(all_Tpvals)
  norm_Tpvals<-norm_Tpvals * -1

  norm_Wpvals<-log10(all_Wpvals)
  norm_Wpvals<-norm_Wpvals * -1

  p1<-ggplot(data = data.frame(x = norm_Tpvals), aes(x = x)) +
    geom_histogram(binwidth = 0.05, fill = "pink", color = "black", alpha = 0.7) +
    labs(x = "-log10(p-value)", y = "Frequency", title = paste(taxon, "\n100x ALDEx2 t-test P-values")) +
    theme_minimal()

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p1)
  } else {
    combined_plots <- combined_plots + wrap_elements(p1)
  }

  p2<-ggplot(data = data.frame(x = norm_Wpvals), aes(x = x)) +
    geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(x = "-log10(p-value)", y = "Frequency", title = paste(taxon, "\n100x ALDEx2 Wilcoxon P-values")) +
    theme_minimal()

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p2)
  } else {
    combined_plots <- combined_plots + wrap_elements(p2)
  }

}
print(combined_plots + plot_layout(ncol = 2))
dev.off()

