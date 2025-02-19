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

ATres<-read.table("ALDEx2RerunDump/Rhizo_ALDEx2_t.txt", sep = "\t", row.names = 1, header = T)
AWres<-read.table("ALDEx2RerunDump/Rhizo_ALDEx2_wilcoxon.txt", sep = "\t", row.names = 1, header = T)

avgTPvals<-apply(ATres, 1, mean)
avgWPvals<-apply(AWres, 1, mean)

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



plot_list <- list()
counter <- 1
pdf("Plots/RDP/SevenTaxaCountsBoxPlots.pdf", width=12, height=13)
par(mfrow=c(4, 3))
par(mar=c(5, 6, 4, 1)+.1)
for (i in 1:length(TopSeven)) {
  mainText<-paste0(TopSeven[i],
                   "\nt-test P=", signif(oTRES[TopSeven[i], 2], 3),
                   "\nALDEx2 t-test P=", signif(aTRES[TopSeven[i], 2], 3),
                   "\nWilcoxon P=", signif(oWRES[TopSeven[i], 2], 3),
                   "\nALDEx2 Wilcoxon P=", signif(aWRES[TopSeven[i], 2], 3))

  taxon_data <- data.frame(
    Count = as.numeric(unlist(rawT[TopSeven[i], ])),
    Condition = meta$conditions
  )

  plot_list[[counter]]<-ggplot(taxon_data, aes(x = Condition, y = Count, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "coral3") +
    ggtitle(mainText) +
    xlab("Conditions") +
    ylab("Counts") +
    theme_minimal() +
    scale_fill_manual(values = c("olivedrab3", "goldenrod1"))

  if (counter == 12 || i == length(TopSeven)) {
    combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 5)
    print(combined_plot)

    plot_list <- list()
    counter <- 0
  }

  counter <- counter + 1
}
dev.off()

