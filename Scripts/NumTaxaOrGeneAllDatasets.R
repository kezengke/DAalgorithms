#variance over mean for all datasets
rm(list = ls())
library(reshape2)
library(ggplot2)

all_NumTaxa<-list()
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

  numTaxa<-nrow(countsT)

  all_NumTaxa[[name]]<-numTaxa
}

# DADA2
all_files<-list.files("CountsTables/DADA2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  name<-paste0("DADA2-", name)
  all_Names<-c(all_Names, name)
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]

  numTaxa<-nrow(countsT)

  all_NumTaxa[[name]]<-numTaxa
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  name<-paste0("WGS-", name)
  all_Names<-c(all_Names, name)
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]

  numTaxa<-nrow(countsT)

  all_NumTaxa[[name]]<-numTaxa
}

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  name<-paste0("RNAseq-", name)
  all_Names<-c(all_Names, name)
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]

  numTaxa<-nrow(countsT)

  all_NumTaxa[[name]]<-numTaxa
}


long_data <- melt(all_NumTaxa, id.vars = "original_order")
long_data$L1 <- factor(long_data$L1, levels = unique(long_data$L1))
long_data$DataType <- ifelse(grepl("RDP", long_data$L1), "RDP",
                             ifelse(grepl("DADA2", long_data$L1), "DADA2",
                                    ifelse(grepl("WGS", long_data$L1), "WGS", "RNAseq")))
long_data$DataType <- factor(long_data$DataType, levels = c("RDP", "DADA2", "WGS", "RNAseq"))
category_colors <- c("RDP" = "#076685",
                     "DADA2" = "#2caaa4",
                     "WGS" = "#c8a67d",
                     "RNAseq" = "#f2a016")

png(paste0("Plots/NumberOfTaxaOrGenesAllDatasets.png"), width=1200, height=1200, res = 300)
p<-ggplot(long_data, aes(x = L1, y = value, fill = DataType)) +
  geom_col() +
  geom_text(aes(label = value),
            vjust = -0.2,
            size = 2) +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  labs(title = "All Datasets Number of Taxa/Genes", x = "Datasets", y = "Number of Taxa/Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

dev.off()
