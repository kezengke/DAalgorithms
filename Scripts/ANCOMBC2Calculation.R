# Calculating ANCOMBC2 results
rm(list = ls())
library(MetagenomeTools)
library(phyloseq)
library(ANCOMBC)

# RDP
rm(list = ls())
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcANCOMBC2(countsT, meta)

  write.table(results, paste0("PkgResults/RDP/ancombc/", gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
}

# dada2
rm(list = ls())
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcANCOMBC2(countsT, meta)

  write.table(results, paste0("PkgResults/dada2/ancombc/", gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
}

# WGS
rm(list = ls())
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcANCOMBC2(countsT, meta)

  write.table(results, paste0("PkgResults/WGS/ancombc/", gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
}

# RNAseq
rm(list = ls())
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcANCOMBC2(countsT, meta)

  write.table(results, paste0("PkgResults/RNAseq/ancombc/", gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
}
