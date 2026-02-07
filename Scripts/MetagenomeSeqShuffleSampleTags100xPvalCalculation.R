#shuffle sample tags 100x and test all pval distribution
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(metagenomeSeq)

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/RDPRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  metagenomeseqpvals<-c()

  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    metagenomeseqresults<-calcMetagenomeSeq(RawcountsT, shuffleMeta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  }
}

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/dada2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  metagenomeseqpvals<-c()

  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    metagenomeseqresults<-calcMetagenomeSeq(RawcountsT, shuffleMeta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  }
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/WGSRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  metagenomeseqpvals<-c()

  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    metagenomeseqresults<-calcMetagenomeSeq(RawcountsT, shuffleMeta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  }
}

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/RNAseqRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  metagenomeseqpvals<-c()

  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    metagenomeseqresults<-calcMetagenomeSeq(RawcountsT, shuffleMeta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  }
}

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/RNAseqRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  metagenomeseqpvals<-c()

  for (i in 1:100) {
    shuffleMeta<-meta
    shuffleMeta$conditions<-sample(shuffleMeta$conditions)

    metagenomeseqresults<-calcMetagenomeSeq(RawcountsT, shuffleMeta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  }
}
