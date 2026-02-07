#shuffle counts in every taxon 100x and test all pval distribution
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(ALDEx2)

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

  aldextpvals<-c()
  aldexwpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    colnames(shuffledCountsT)<-colnames(RawcountsT)

    aldextresults<-calcALDEx2Ttest(shuffledCountsT, meta)
    aldextpvals<-cbind(aldextpvals, aldextresults$pval)

    aldexwresults<-calcALDEx2Wilcoxon(shuffledCountsT, meta)
    aldexwpvals<-cbind(aldexwpvals, aldexwresults$pval)
  }

  rownames(aldextpvals)<-rownames(RawcountsT)
  rownames(aldexwpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(aldextpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(aldextpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(aldexwpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(aldexwpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  }
}

# DADA2
all_files<-list.files("CountsTables/DADA2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/DADA2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  aldextpvals<-c()
  aldexwpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    colnames(shuffledCountsT)<-colnames(RawcountsT)

    aldextresults<-calcALDEx2Ttest(shuffledCountsT, meta)
    aldextpvals<-cbind(aldextpvals, aldextresults$pval)

    aldexwresults<-calcALDEx2Wilcoxon(shuffledCountsT, meta)
    aldexwpvals<-cbind(aldexwpvals, aldexwresults$pval)
  }

  rownames(aldextpvals)<-rownames(RawcountsT)
  rownames(aldexwpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleInTaxonDump/DADA2/"
  if (file.exists(file.path(save_dir))){
    write.table(aldextpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(aldextpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ShuffleInTaxonDump/DADA2/"
  if (file.exists(file.path(save_dir))){
    write.table(aldexwpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(aldexwpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
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

  aldextpvals<-c()
  aldexwpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    colnames(shuffledCountsT)<-colnames(RawcountsT)

    aldextresults<-calcALDEx2Ttest(shuffledCountsT, meta)
    aldextpvals<-cbind(aldextpvals, aldextresults$pval)

    aldexwresults<-calcALDEx2Wilcoxon(shuffledCountsT, meta)
    aldexwpvals<-cbind(aldexwpvals, aldexwresults$pval)
  }

  rownames(aldextpvals)<-rownames(RawcountsT)
  rownames(aldexwpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(aldextpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(aldextpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2T.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(aldexwpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(aldexwpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ALDEx2W.txt"), sep = "\t", row.names = T)
  }
}
