#rnorm resample up, log10 pval R^2 plots
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(patchwork)

times<-c(0.5, 2, seq(5, 45, 5))

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/RDP/RDPIncreaseRnormResampleRSquared(", name, ").png"), width=1800, height=1200, res = 300)
  par(mar=c(5,6,4,1)+.1)

  R2<-c()
  for (m in times) {
    p<-read.table(paste0("ResampleDump/RDP/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", header = T, row.names = 1)
    newR2<-(cor(p$X1, p$X2, method = "pearson"))^2
    R2<-c(R2, newR2)
  }

  data <- data.frame(times, R2)

  # data <- data.frame(times = factor(times, levels = as.character(times)), R2 = R2)
  # data$times <- paste0(data$times, "x")
  # data$times<-factor(data$times, levels = data$times)

  p<-ggplot(data, aes(x = times, y = R2)) +
    geom_point(color = "coral3", size = 3) +
    geom_smooth(method = "loess", alpha = 0.3, color = "gold2", fill = "gold2", se = TRUE) +
    labs(x = "Multiple of variance", y = "R-squared values", title = paste0("(RDP-", name, ") Resample of multiple var. vs. R-squared of \nLog10 p-value plots")) +
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
  png(paste0("Plots/dada2/dada2IncreaseRnormResampleRSquared(", name, ").png"), width=1800, height=1200, res = 300)
  par(mar=c(5,6,4,1)+.1)

  R2<-c()
  for (m in times) {
    p<-read.table(paste0("ResampleDump/dada2/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", header = T, row.names = 1)
    newR2<-(cor(p$X1, p$X2, method = "pearson"))^2
    R2<-c(R2, newR2)
  }

  data <- data.frame(times, R2)

  # data <- data.frame(times = factor(times, levels = as.character(times)), R2 = R2)
  # data$times <- paste0(data$times, "x")
  # data$times<-factor(data$times, levels = data$times)

  p<-ggplot(data, aes(x = times, y = R2)) +
    geom_point(color = "coral3", size = 3) +
    geom_smooth(method = "loess", alpha = 0.3, color = "gold2", fill = "gold2", se = TRUE) +
    labs(x = "Multiple of variance", y = "R-squared values", title = paste0("(dada2-", name, ") Resample of multiple var. vs. R-squared of \nLog10 p-value plots")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = rel(1))
    )

  print(p)

  dev.off()
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/WGS/WGSIncreaseRnormResampleRSquared(", name, ").png"), width=1800, height=1200, res = 300)
  par(mar=c(5,6,4,1)+.1)

  R2<-c()
  for (m in times) {
    p<-read.table(paste0("ResampleDump/WGS/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", header = T, row.names = 1)
    newR2<-(cor(p$X1, p$X2, method = "pearson"))^2
    R2<-c(R2, newR2)
  }

  data <- data.frame(times, R2)

  # data <- data.frame(times = factor(times, levels = as.character(times)), R2 = R2)
  # data$times <- paste0(data$times, "x")
  # data$times<-factor(data$times, levels = data$times)

  p<-ggplot(data, aes(x = times, y = R2)) +
    geom_point(color = "coral3", size = 3) +
    geom_smooth(method = "loess", alpha = 0.3, color = "gold2", fill = "gold2", se = TRUE) +
    labs(x = "Times", y = "R-squared values", title = paste0("(WGS-", name, ") Resample of multiple var. vs. R-squared of \nLog10 p-value plots")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = rel(1))
    )

  print(p)

  dev.off()
}
