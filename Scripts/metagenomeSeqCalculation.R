# Calculating metagenomeSeq results
rm(list = ls())
library(MetagenomeTools)
library(metagenomeSeq)

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

  results<-calcMetagenomeSeq(countsT, meta)

  write.table(results, paste0("PkgResults/RDP/metagenomeSeq/", gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
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

  results<-calcMetagenomeSeq(countsT, meta)

  write.table(results, paste0("PkgResults/dada2/metagenomeSeq/", gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
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

  results<-calcMetagenomeSeq(countsT, meta)

  write.table(results, paste0("PkgResults/WGS/metagenomeSeq/", gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
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

  results<-calcMetagenomeSeq(countsT, meta)

  write.table(results, paste0("PkgResults/RNAseq/metagenomeSeq/", gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
}


