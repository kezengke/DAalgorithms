#rnorm resample down, log10 pval rho^2 plots
rm(list = ls())
library(MetagenomeTools)
library(DESeq2)
library(dplyr)
library(ggplot2)
set.seed(9527)

times<-seq(0.005, 1.5, 0.005)

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/RDP/RDPDecreaseRnormResampleRhosquared(", name, ").png"), width=1800, height=1200, res = 300)
  par(mar=c(5,6,4,1)+.1)
  countsT<-LoadCountsT(file)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  rho2<-c()
  for (m in times) {
    newT<-resampleRNORM(countsT, meta, m)
    dRES<-read.table(paste0("PkgResults/RDP/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)

    tRES<-calcTtest(newT, meta)
    tRES<-na.omit(tRES)
    dRES<-dRES[rownames(tRES), , drop = F]

    p<-directPFun(tRES, dRES)
    colnames(p)<-c("X1", "X2")
    newRho2<-(cor(p$X1, p$X2, method = "spearman"))^2
    rho2<-c(rho2, newRho2)
  }
  data <- data.frame(times, rho2)

  p<-ggplot(data, aes(x = times, y = rho2)) +
    geom_point(color = "coral3", size = 3) +
    geom_smooth(method = "loess", alpha = 0.3, color = "gold2", fill = "gold2", se = TRUE) +
    labs(x = "Multiple of variance", y = "rho-squared values", title = paste0("(RDP-", name, ") Resample of multiple var. vs. Rho-squared of Log10 p-value plots")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = rel(0.8))
    )

  print(p)

  dev.off()
}

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/dada2/dada2DecreaseRnormResampleRhosquared(", name, ").png"), width=1800, height=1200, res = 300)
  par(mar=c(5,6,4,1)+.1)
  countsT<-LoadCountsT(file)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  rho2<-c()
  for (m in times) {
    newT<-resampleRNORM(countsT, meta, m)
    dRES<-read.table(paste0("PkgResults/dada2/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)

    tRES<-calcTtest(newT, meta)
    tRES<-na.omit(tRES)
    dRES<-dRES[rownames(tRES), , drop = F]

    p<-directPFun(tRES, dRES)
    colnames(p)<-c("X1", "X2")
    newRho2<-(cor(p$X1, p$X2, method = "spearman"))^2
    rho2<-c(rho2, newRho2)
  }
  data <- data.frame(times, rho2)

  p<-ggplot(data, aes(x = times, y = rho2)) +
    geom_point(color = "coral3", size = 3) +
    geom_smooth(method = "loess", alpha = 0.3, color = "gold2", fill = "gold2", se = TRUE) +
    labs(x = "Multiple of variance", y = "rho-squared values", title = paste0("(dada2-", name, ") Resample of multiple var. vs. Rho-squared of Log10 p-value plots")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = rel(0.8))
    )

  print(p)

  dev.off()
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/WGS/WGSDecreaseRnormResampleRhosquared(", name, ").png"), width=1800, height=1200, res = 300)
  par(mar=c(5,6,4,1)+.1)
  countsT<-LoadCountsT(file)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  rho2<-c()
  for (m in times) {
    newT<-resampleRNORM(countsT, meta, m)
    dRES<-read.table(paste0("PkgResults/WGS/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)

    tRES<-calcTtest(newT, meta)
    tRES<-na.omit(tRES)
    dRES<-dRES[rownames(tRES), , drop = F]

    p<-directPFun(tRES, dRES)
    colnames(p)<-c("X1", "X2")
    newRho2<-(cor(p$X1, p$X2, method = "spearman"))^2
    rho2<-c(rho2, newRho2)
  }
  data <- data.frame(times, rho2)

  p<-ggplot(data, aes(x = times, y = rho2)) +
    geom_point(color = "coral3", size = 3) +
    geom_smooth(method = "loess", alpha = 0.3, color = "gold2", fill = "gold2", se = TRUE) +
    labs(x = "Multiple of variance", y = "rho-squared values", title = paste0("(WGS-", name, ") Resample of multiple var. vs. Rho-squared of Log10 p-value plots")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = rel(0.8))
    )

  print(p)

  dev.off()
}
