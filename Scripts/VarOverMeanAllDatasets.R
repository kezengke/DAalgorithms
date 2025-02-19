#variance over mean for all datasets
rm(list = ls())
library(reshape2)
library(ggplot2)

all_VarOverMean<-list()
all_Names<-c()
# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  name<-paste0("RDP-", name)
  all_Names<-c(all_Names, name)
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]

  allMean<-apply(countsT, 1, mean)
  allVar<-apply(countsT, 1, var)

  allMean<-log10(allMean)
  allVar<-log10(allVar)

  VarOverMean<-allVar/allMean

  all_VarOverMean[[name]]<-VarOverMean
}

# DADA2
all_files<-list.files("CountsTables/DADA2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  name<-paste0("DADA2-", name)
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]

  allMean<-apply(countsT, 1, mean)
  allVar<-apply(countsT, 1, var)

  allMean<-log10(allMean)
  allVar<-log10(allVar)

  VarOverMean<-allVar/allMean

  all_VarOverMean[[name]]<-VarOverMean
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  name<-paste0("WGS-", name)
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]

  allMean<-apply(countsT, 1, mean)
  allVar<-apply(countsT, 1, var)

  allMean<-log10(allMean)
  allVar<-log10(allVar)

  VarOverMean<-allVar/allMean

  all_VarOverMean[[name]]<-VarOverMean
}

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  name<-paste0("RNAseq-", name)
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]

  allMean<-apply(countsT, 1, mean)
  allVar<-apply(countsT, 1, var)

  allMean<-log10(allMean)
  allVar<-log10(allVar)

  VarOverMean<-allVar/allMean

  all_VarOverMean[[name]]<-VarOverMean
}

long_data <- melt(all_VarOverMean, id.vars = "original_order")
long_data$L1 <- factor(long_data$L1, levels = unique(long_data$L1))
long_data$DataType <- ifelse(grepl("RDP", long_data$L1), "RDP",
                             ifelse(grepl("DADA2", long_data$L1), "DADA2",
                                    ifelse(grepl("WGS", long_data$L1), "WGS", "RNAseq")))
long_data$DataType <- factor(long_data$DataType, levels = c("RDP", "DADA2", "WGS", "RNAseq"))
category_colors <- c("RDP" = "#FCBB44",
                     "DADA2" = "#F1766D",
                     "WGS" = "olivedrab3",
                     "RNAseq" = "#839DD1")

png(paste0("Plots/VarOverMeanAllDatasets.png"), width=1200, height=1200, res = 300)
p<-ggplot(long_data, aes(x = L1, y = value, fill = DataType)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  labs(title = "All Datasets Variance Over Mean", x = "Datasets", y = "Variance/Mean") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

dev.off()
