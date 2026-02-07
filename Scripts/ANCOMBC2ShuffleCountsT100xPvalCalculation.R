#shuffle counts in entire counts table 100x and test all pval distribution
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(phyloseq)
library(ANCOMBC)

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

  ancombc2pvals<-c()
  for (i in 1:100) {
    all_values <- as.vector(as.matrix(RawcountsT))
    repeat {
      shuffled_values <- sample(all_values)  # Shuffle the countsT
      shuffledCountsT <- matrix(shuffled_values, nrow = nrow(RawcountsT), ncol = ncol(RawcountsT))  # Reshape to original size
      shuffledCountsT <- as.data.frame(shuffledCountsT)  # Convert to data frame

      # Check if any row is entirely zeros
      if (all(rowSums(shuffledCountsT) > 0)) {
        break  # Exit the loop if no rows are all zeros
      }
    }

    rownames(shuffledCountsT) <- rownames(RawcountsT)
    colnames(shuffledCountsT) <- colnames(RawcountsT)

    ancombc2results<-calcANCOMBC2(shuffledCountsT, meta)
    ancombc2pvals<-cbind(ancombc2pvals, ancombc2results$pval)

  }

  rownames(ancombc2pvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ancombc2pvals[is.na(ancombc2pvals)]<-1

  save_dir<-"ShuffleCountsTDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(ancombc2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ancombc2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
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

  ancombc2pvals<-c()
  for (i in 1:100) {
    all_values <- as.vector(as.matrix(RawcountsT))
    repeat {
      shuffled_values <- sample(all_values)  # Shuffle the countsT
      shuffledCountsT <- matrix(shuffled_values, nrow = nrow(RawcountsT), ncol = ncol(RawcountsT))  # Reshape to original size
      shuffledCountsT <- as.data.frame(shuffledCountsT)  # Convert to data frame

      # Check if any row is entirely zeros
      if (all(rowSums(shuffledCountsT) > 0)) {
        break  # Exit the loop if no rows are all zeros
      }
    }

    rownames(shuffledCountsT) <- rownames(RawcountsT)
    colnames(shuffledCountsT) <- colnames(RawcountsT)

    ancombc2results<-calcANCOMBC2(shuffledCountsT, meta)
    ancombc2pvals<-cbind(ancombc2pvals, ancombc2results$pval)

  }

  rownames(ancombc2pvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ancombc2pvals[is.na(ancombc2pvals)]<-1

  save_dir<-"ShuffleCountsTDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(ancombc2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ancombc2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
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

  ancombc2pvals<-c()
  for (i in 1:100) {
    all_values <- as.vector(as.matrix(RawcountsT))
    repeat {
      shuffled_values <- sample(all_values)  # Shuffle the countsT
      shuffledCountsT <- matrix(shuffled_values, nrow = nrow(RawcountsT), ncol = ncol(RawcountsT))  # Reshape to original size
      shuffledCountsT <- as.data.frame(shuffledCountsT)  # Convert to data frame

      # Check if any row is entirely zeros
      if (all(rowSums(shuffledCountsT) > 0)) {
        break  # Exit the loop if no rows are all zeros
      }
    }

    rownames(shuffledCountsT) <- rownames(RawcountsT)
    colnames(shuffledCountsT) <- colnames(RawcountsT)

    ancombc2results<-calcANCOMBC2(shuffledCountsT, meta)
    ancombc2pvals<-cbind(ancombc2pvals, ancombc2results$pval)

  }

  rownames(ancombc2pvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ancombc2pvals[is.na(ancombc2pvals)]<-1

  save_dir<-"ShuffleCountsTDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(ancombc2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ancombc2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_ancombc2.txt"), sep = "\t", row.names = T)
  }
}
