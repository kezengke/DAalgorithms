#counts boxplot with pvalues
rm(list = ls())
library(ggplot2)
library(patchwork)
wilcoxP<-read.table("PkgResults/RDP/Wilcoxon/rhizo_wilcox.txt", header = T, row.names = 1, sep = "\t")
aldexWP<-read.table("PkgResults/RDP/ALDEx2Wilcoxon/rhizo_ALDEx2W.txt", header = T, row.names = 1, sep = "\t")
rawT<-read.table("CountsTables/RDPRaw/rhizo.txt", header = T, row.names = 1, sep = "\t")
normT<-read.table("CountsTables/RDPNorm/rhizo.txt", header = T, row.names = 1, sep = "\t")
meta<-read.table("MetaData/metadata_rhizo.txt", header = T, row.names = 1, sep = "\t")

# taxaOrder<-rownames(aldexWP)[order(aldexWP$pval)]
taxaOrder<-rownames(wilcoxP)[order(wilcoxP$pval)]

plot_list <- list()
counter <- 1
pdf("Plots/RDP/SortedByAldex2WilcoxonCountsBoxPlots.pdf", width=12, height=13)
par(mfrow=c(4, 3))
par(mar=c(5, 6, 4, 1)+.1)
for (i in 1:length(taxaOrder)) {
  taxon<-taxaOrder[i]
  mainText<-paste0(taxon, "\nWilcoxon P=", signif(wilcoxP[taxon, 2], 3),
                   "\nALDEx2 Wilcoxon P=", signif(aldexWP[taxon, 2], 3))

  taxon_data <- data.frame(
    Count = as.numeric(unlist(rawT[taxon, ])),
    Condition = meta$characters
  )

  plot_list[[counter]]<-ggplot(taxon_data, aes(x = Condition, y = Count, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "coral3") +
    ggtitle(mainText) +
    xlab("Conditions") +
    ylab("Counts") +
    theme_minimal() +
    scale_fill_manual(values = c("olivedrab3", "goldenrod1"))

  if (counter == 12 || i == length(taxaOrder)) {
    combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 5)
    print(combined_plot)

    plot_list <- list()
    counter <- 0
  }

  counter <- counter + 1
}
dev.off()

plot_list <- list()
counter <- 1
pdf("Plots/RDP/SortedByAldex2WilcoxonCountsBoxPlots(Norm).pdf", width=12, height=13)
par(mfrow=c(4, 3))
par(mar=c(5, 6, 4, 1)+.1)
for (i in 1:length(taxaOrder)) {
  taxon<-taxaOrder[i]
  mainText<-paste0(taxon, "\nWilcoxon P=", signif(wilcoxP[taxon, 2], 3),
                   "\nALDEx2 Wilcoxon P=", signif(aldexWP[taxon, 2], 3))

  taxon_data <- data.frame(
    Count = as.numeric(unlist(normT[taxon, ])),
    Condition = meta$characters
  )

  plot_list[[counter]]<-ggplot(taxon_data, aes(x = Condition, y = Count, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "coral3") +
    ggtitle(mainText) +
    xlab("Conditions") +
    ylab("Normalized Counts") +
    theme_minimal() +
    scale_fill_manual(values = c("skyblue", "pink"))

  if (counter == 12 || i == length(taxaOrder)) {
    combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 5)
    print(combined_plot)

    plot_list <- list()
    counter <- 0
  }

  counter <- counter + 1
}
dev.off()

