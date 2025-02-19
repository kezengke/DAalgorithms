rm(list = ls())
set.seed(9527)
library(ggplot2)
library(patchwork)
library(ALDEx2)
library(dplyr)
library(MetagenomeTools)

# # RDP
# all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)
#
# for (file in all_files) {
#   countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
#   meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
#                    header = T, row.names = 1)
#
#   countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
#   meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]
#
#   rownames(meta)<-colnames(countsT)
#   colnames(meta)<-"conditions"
#
#
#   allT_pvals<-vector()
#   allW_pvals<-vector()
#   for (i in 1:100) {
#     ALDEx2RES<-calcALDEx2Ttest(countsT, meta)
#     pval<-ALDEx2RES[, 2, drop = F]
#     allT_pvals[i]<-pval
#
#     ALDEx2RES<-calcALDEx2Wilcoxon(countsT, meta)
#     pval<-ALDEx2RES[, 2, drop = F]
#     allW_pvals[i]<-pval
#   }
#
#   Tpvals<-do.call(cbind, allT_pvals)
#   Tpvals<-data.frame(Tpvals)
#   rownames(Tpvals)<-rownames(countsT)
#   avgTPvals<-apply(Tpvals, 1, mean)
#
#   Wpvals<-do.call(cbind, allW_pvals)
#   Wpvals<-data.frame(Wpvals)
#   rownames(Wpvals)<-rownames(countsT)
#   avgWPvals<-apply(Wpvals, 1, mean)
#
#   save_dir<-"ALDEx2RerunDump/RDP/ALDEx2ttest/"
#   if (file.exists(file.path(save_dir))){
#     write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ALDEx2RerunDump/RDP/ALDEx2Wilcoxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
#   }
#
# }

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"


  allT_pvals<-vector()
  allW_pvals<-vector()
  for (i in 1:100) {
    ALDEx2RES<-calcALDEx2Ttest(countsT, meta)
    pval<-ALDEx2RES[, 2, drop = F]
    allT_pvals[i]<-pval

    ALDEx2RES<-calcALDEx2Wilcoxon(countsT, meta)
    pval<-ALDEx2RES[, 2, drop = F]
    allW_pvals[i]<-pval
  }

  Tpvals<-do.call(cbind, allT_pvals)
  Tpvals<-data.frame(Tpvals)
  rownames(Tpvals)<-rownames(countsT)
  avgTPvals<-apply(Tpvals, 1, mean)

  Wpvals<-do.call(cbind, allW_pvals)
  Wpvals<-data.frame(Wpvals)
  rownames(Wpvals)<-rownames(countsT)
  avgWPvals<-apply(Wpvals, 1, mean)

  save_dir<-"ALDEx2RerunDump/dada2/ALDEx2ttest/"
  if (file.exists(file.path(save_dir))){
    write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ALDEx2RerunDump/dada2/ALDEx2Wilcoxon/"
  if (file.exists(file.path(save_dir))){
    write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  }

}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"


  allT_pvals<-vector()
  allW_pvals<-vector()
  for (i in 1:100) {
    ALDEx2RES<-calcALDEx2Ttest(countsT, meta)
    pval<-ALDEx2RES[, 2, drop = F]
    allT_pvals[i]<-pval

    ALDEx2RES<-calcALDEx2Wilcoxon(countsT, meta)
    pval<-ALDEx2RES[, 2, drop = F]
    allW_pvals[i]<-pval
  }

  Tpvals<-do.call(cbind, allT_pvals)
  Tpvals<-data.frame(Tpvals)
  rownames(Tpvals)<-rownames(countsT)
  avgTPvals<-apply(Tpvals, 1, mean)

  Wpvals<-do.call(cbind, allW_pvals)
  Wpvals<-data.frame(Wpvals)
  rownames(Wpvals)<-rownames(countsT)
  avgWPvals<-apply(Wpvals, 1, mean)

  save_dir<-"ALDEx2RerunDump/WGS/ALDEx2ttest/"
  if (file.exists(file.path(save_dir))){
    write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ALDEx2RerunDump/WGS/ALDEx2Wilcoxon/"
  if (file.exists(file.path(save_dir))){
    write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  }

}

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"


  allT_pvals<-vector()
  allW_pvals<-vector()
  for (i in 1:100) {
    ALDEx2RES<-calcALDEx2Ttest(countsT, meta)
    pval<-ALDEx2RES[, 2, drop = F]
    allT_pvals[i]<-pval

    ALDEx2RES<-calcALDEx2Wilcoxon(countsT, meta)
    pval<-ALDEx2RES[, 2, drop = F]
    allW_pvals[i]<-pval
  }

  Tpvals<-do.call(cbind, allT_pvals)
  Tpvals<-data.frame(Tpvals)
  rownames(Tpvals)<-rownames(countsT)
  avgTPvals<-apply(Tpvals, 1, mean)

  Wpvals<-do.call(cbind, allW_pvals)
  Wpvals<-data.frame(Wpvals)
  rownames(Wpvals)<-rownames(countsT)
  avgWPvals<-apply(Wpvals, 1, mean)

  save_dir<-"ALDEx2RerunDump/RNAseq/ALDEx2ttest/"
  if (file.exists(file.path(save_dir))){
    write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(Tpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ALDEx2RerunDump/RNAseq/ALDEx2Wilcoxon/"
  if (file.exists(file.path(save_dir))){
    write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(Wpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  }

}
