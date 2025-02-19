#rnorm resample up, log10 pval R^2 plots
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(reshape2)
library(patchwork)

times<-c(0.5, 2, seq(5, 45, 5))

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RDP/PvalAbsDiffChange(RDP).png"), width=3000, height=1200, res = 300)
combined_plts<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  tRES<-read.table(paste0("PkgResults/RDP/ttest/", name, "_t.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RDP/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  tRES<-processTtestRes(tRES)
  dRES<-processDESeq2Res(dRES)

  testingTaxa<-rownames(dRES[order(dRES$pval), , drop = F])[1:50]

  tRES<-tRES[testingTaxa, , drop = F]
  dRES<-dRES[testingTaxa, , drop = F]

  oRES<-cbind(tRES$p_directed, dRES$p_directed)
  rownames(oRES)<-testingTaxa
  colnames(oRES)<-c("X1", "X2")
  oRES<-data.frame(oRES)
  oDiff<-abs(oRES$X1 - oRES$X2)

  allDiff<-c()
  Rho2<-c()
  for (m in times) {
    p<-read.table(paste0("ResampleDump/RDP/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", header = T, row.names = 1)
    newRho2<-(cor(p$X1, p$X2, method = "spearman"))^2
    Rho2<-c(Rho2, newRho2)
    absDiff<-abs(p$X1 - p$X2)
    allDiff<-cbind(allDiff, absDiff)
  }
  data <- data.frame(allDiff)
  colnames(data)<-times
  rownames(data)<-rownames(p)

  nDiff<-data[testingTaxa, which.max(Rho2), drop = F]

  Diffs<-cbind(oDiff, nDiff)
  colnames(Diffs)<-c("Original", "Resampled")

  data_melted <- melt(as.matrix(Diffs))
  data_melted$Var2 <- factor(data_melted$Var2, levels = colnames(Diffs))

  plt<-ggplot(data_melted, aes(x = Var2, y = value)) +
    geom_boxplot(color = "cornflowerblue") +
    labs(title = paste0("(RDP-", name, ")"),
         x = "P-value difference type",
         y = "Absolute Difference") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  if (is.null(combined_plts)) {
    combined_plts <- plt
  } else {
    combined_plts <- combined_plts + plt
  }

}
print(combined_plts + plot_layout(ncol = 3))

dev.off()

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/dada2/PvalAbsDiffChange(dada2).png"), width=3000, height=1200, res = 300)
combined_plts<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  tRES<-read.table(paste0("PkgResults/dada2/ttest/", name, "_t.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/dada2/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  tRES<-processTtestRes(tRES)
  dRES<-processDESeq2Res(dRES)

  testingTaxa<-rownames(dRES[order(dRES$pval), , drop = F])[1:50]

  tRES<-tRES[testingTaxa, , drop = F]
  dRES<-dRES[testingTaxa, , drop = F]

  oRES<-cbind(tRES$p_directed, dRES$p_directed)
  rownames(oRES)<-testingTaxa
  colnames(oRES)<-c("X1", "X2")
  oRES<-data.frame(oRES)
  oDiff<-abs(oRES$X1 - oRES$X2)

  allDiff<-c()
  Rho2<-c()
  for (m in times) {
    p<-read.table(paste0("ResampleDump/dada2/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", header = T, row.names = 1)
    newRho2<-(cor(p$X1, p$X2, method = "spearman"))^2
    Rho2<-c(Rho2, newRho2)
    absDiff<-abs(p$X1 - p$X2)
    allDiff<-cbind(allDiff, absDiff)
  }
  data <- data.frame(allDiff)
  colnames(data)<-times
  rownames(data)<-rownames(p)

  nDiff<-data[testingTaxa, which.max(Rho2), drop = F]

  Diffs<-cbind(oDiff, nDiff)
  colnames(Diffs)<-c("Original", "Resampled")

  data_melted <- melt(as.matrix(Diffs))
  data_melted$Var2 <- factor(data_melted$Var2, levels = colnames(Diffs))

  plt<-ggplot(data_melted, aes(x = Var2, y = value)) +
    geom_boxplot(color = "cornflowerblue") +
    labs(title = paste0("(dada2-", name, ")"),
         x = "P-value difference type",
         y = "Absolute Difference") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  if (is.null(combined_plts)) {
    combined_plts <- plt
  } else {
    combined_plts <- combined_plts + plt
  }

}
print(combined_plts + plot_layout(ncol = 3))

dev.off()

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/WGS/PvalAbsDiffChange(WGS).png"), width=2000, height=1200, res = 300)
combined_plts<-NULL
for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  tRES<-read.table(paste0("PkgResults/WGS/ttest/", name, "_t.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/WGS/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  tRES<-processTtestRes(tRES)
  dRES<-processDESeq2Res(dRES)

  testingTaxa<-rownames(dRES[order(dRES$pval), , drop = F])[1:50]

  tRES<-tRES[testingTaxa, , drop = F]
  dRES<-dRES[testingTaxa, , drop = F]

  oRES<-cbind(tRES$p_directed, dRES$p_directed)
  rownames(oRES)<-testingTaxa
  colnames(oRES)<-c("X1", "X2")
  oRES<-data.frame(oRES)
  oDiff<-abs(oRES$X1 - oRES$X2)

  allDiff<-c()
  Rho2<-c()
  for (m in times) {
    p<-read.table(paste0("ResampleDump/WGS/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", header = T, row.names = 1)
    newRho2<-(cor(p$X1, p$X2, method = "spearman"))^2
    Rho2<-c(Rho2, newRho2)
    absDiff<-abs(p$X1 - p$X2)
    allDiff<-cbind(allDiff, absDiff)
  }
  data <- data.frame(allDiff)
  colnames(data)<-times
  rownames(data)<-rownames(p)

  nDiff<-data[testingTaxa, which.max(Rho2), drop = F]

  Diffs<-cbind(oDiff, nDiff)
  colnames(Diffs)<-c("Original", "Resampled")

  data_melted <- melt(as.matrix(Diffs))
  data_melted$Var2 <- factor(data_melted$Var2, levels = colnames(Diffs))

  plt<-ggplot(data_melted, aes(x = Var2, y = value)) +
    geom_boxplot(color = "cornflowerblue") +
    labs(title = paste0("(WGS-", name, ")"),
         x = "P-value difference type",
         y = "Absolute Difference") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  if (is.null(combined_plts)) {
    combined_plts <- plt
  } else {
    combined_plts <- combined_plts + plt
  }

}
print(combined_plts + plot_layout(ncol = 2))

dev.off()
