# density plot of resampled two conditions
library(ggplot2)
library(ggExtra)
library(patchwork)
library(MetagenomeTools)
rm(list = ls())
file<-"CountsTables/RDPRaw/assal.txt"
countsT<-read.table(file, sep = "\t", header = T, row.names = 1)
meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                 header = T, row.names = 1)
meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

rownames(meta)<-colnames(countsT)
colnames(meta)<-"conditions"

newT<-resampleRNORM(countsT, meta, 1)

pdf("Plots/Density(Assal).pdf", width=14, height=7)
for (i in 1:nrow(countsT)) {
  combined_plots<-NULL
  taxon<-rownames(countsT)[i]

  # original
  group1<-unlist(countsT[taxon, rownames(meta)[meta$conditions == "Pre"]])
  group2<-unlist(countsT[taxon, rownames(meta)[meta$conditions == "Post"]])

  data <- data.frame(
    value = c(group1, group2),
    group =c(rep("Group 1", length(group1)), rep("Group 2", length(group2)))
  )

  p1<-ggplot(data, aes(x = value, fill = group)) +
    geom_density(alpha = 0.5) +  # Add density plot with transparency
    scale_fill_manual(values = c("blue", "red")) +
    labs(title = "Density Distribution of Group 1 and Group 2",
         x = "Value",
         y = "Density") +
    theme_classic()

  # resampled
  group1<-unlist(newT[taxon, rownames(meta)[meta$conditions == "Pre"]])
  group2<-unlist(newT[taxon, rownames(meta)[meta$conditions == "Post"]])

  data <- data.frame(
    value = c(group1, group2),
    group =c(rep("Group 1", length(group1)), rep("Group 2", length(group2)))
  )

  p2<-ggplot(data, aes(x = value, fill = group)) +
    geom_density(alpha = 0.5) +  # Add density plot with transparency
    scale_fill_manual(values = c("tan3", "green")) +
    labs(title = "Density Distribution of Group 1 and Group 2 Resampled",
         x = "Value",
         y = "Density") +
    theme_classic()

  combined_plots<- p1 + p2
  print(combined_plots + plot_layout(ncol = 2))
}
dev.off()
