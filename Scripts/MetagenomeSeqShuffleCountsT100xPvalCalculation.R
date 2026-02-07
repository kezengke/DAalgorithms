#shuffle counts in entire counts table 100x and test all pval distribution
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

    nr <- nrow(RawcountsT)
    nc <- ncol(RawcountsT)
    all_values <- as.vector(as.matrix(RawcountsT))

    max_iter <- 5000
    iter <- 0

    repeat {
      iter <- iter + 1
      shuffled_values <- sample(all_values, length(all_values), replace = FALSE)
      shuffledCountsM <- matrix(shuffled_values, nrow = nr, ncol = nc)  # fills by column

      # require no all-zero rows and no all-zero columns
      ok <- all(rowSums(shuffledCountsM != 0) > 0) && all(colSums(shuffledCountsM != 0) >= 2)

      if (ok) break
      if (iter >= max_iter) {
        stop("Could not get a shuffle without all-zero rows/cols after ", max_iter,
             " tries. The table is likely too sparse for this constraint.")
      }
    }

    rownames(shuffledCountsM) <- rownames(RawcountsT)
    colnames(shuffledCountsM) <- colnames(RawcountsT)

    metagenomeseqresults<-calcMetagenomeSeq(shuffledCountsM, meta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)

  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  metagenomeseqpvals[is.na(metagenomeseqpvals)]<-1

  save_dir<-"ShuffleCountsTDump/RDP/"
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

    nr <- nrow(RawcountsT)
    nc <- ncol(RawcountsT)
    all_values <- as.vector(as.matrix(RawcountsT))

    max_iter <- 5000
    iter <- 0

    repeat {
      iter <- iter + 1
      shuffled_values <- sample(all_values, length(all_values), replace = FALSE)
      shuffledCountsM <- matrix(shuffled_values, nrow = nr, ncol = nc)  # fills by column

      # require no all-zero rows and no all-zero columns
      ok <- all(rowSums(shuffledCountsM != 0) > 0) && all(colSums(shuffledCountsM != 0) >= 2)

      if (ok) break
      if (iter >= max_iter) {
        stop("Could not get a shuffle without all-zero rows/cols after ", max_iter,
             " tries. The table is likely too sparse for this constraint.")
      }
    }

    rownames(shuffledCountsM) <- rownames(RawcountsT)
    colnames(shuffledCountsM) <- colnames(RawcountsT)

    metagenomeseqresults<-calcMetagenomeSeq(shuffledCountsM, meta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)

  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  metagenomeseqpvals[is.na(metagenomeseqpvals)]<-1

  save_dir<-"ShuffleCountsTDump/dada2/"
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

    nr <- nrow(RawcountsT)
    nc <- ncol(RawcountsT)
    all_values <- as.vector(as.matrix(RawcountsT))

    max_iter <- 5000
    iter <- 0

    repeat {
      iter <- iter + 1
      shuffled_values <- sample(all_values, length(all_values), replace = FALSE)
      shuffledCountsM <- matrix(shuffled_values, nrow = nr, ncol = nc)  # fills by column

      # require no all-zero rows and no all-zero columns
      ok <- all(rowSums(shuffledCountsM != 0) > 0) && all(colSums(shuffledCountsM != 0)  >= 2)

      if (ok) break
      if (iter >= max_iter) {
        stop("Could not get a shuffle without all-zero rows/cols after ", max_iter,
             " tries. The table is likely too sparse for this constraint.")
      }
    }

    rownames(shuffledCountsM) <- rownames(RawcountsT)
    colnames(shuffledCountsM) <- colnames(RawcountsT)

    metagenomeseqresults<-calcMetagenomeSeq(shuffledCountsM, meta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)

  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  metagenomeseqpvals[is.na(metagenomeseqpvals)]<-1

  save_dir<-"ShuffleCountsTDump/WGS/"
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

    nr <- nrow(RawcountsT)
    nc <- ncol(RawcountsT)
    all_values <- as.vector(as.matrix(RawcountsT))

    max_iter <- 5000
    iter <- 0

    repeat {
      iter <- iter + 1
      shuffled_values <- sample(all_values, length(all_values), replace = FALSE)
      shuffledCountsM <- matrix(shuffled_values, nrow = nr, ncol = nc)  # fills by column

      # require no all-zero rows and no all-zero columns
      ok <- all(rowSums(shuffledCountsM != 0) > 0) && all(colSums(shuffledCountsM != 0)  >= 2)

      if (ok) break
      if (iter >= max_iter) {
        stop("Could not get a shuffle without all-zero rows/cols after ", max_iter,
             " tries. The table is likely too sparse for this constraint.")
      }
    }

    rownames(shuffledCountsM) <- rownames(RawcountsT)
    colnames(shuffledCountsM) <- colnames(RawcountsT)

    metagenomeseqresults<-calcMetagenomeSeq(shuffledCountsM, meta)
    metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)

  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  metagenomeseqpvals[is.na(metagenomeseqpvals)]<-1

  save_dir<-"ShuffleCountsTDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  }
}
