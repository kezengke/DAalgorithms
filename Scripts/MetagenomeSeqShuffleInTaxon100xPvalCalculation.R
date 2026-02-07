#shuffle counts in every taxon 100x and test all pval distribution
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

    max_iter <- 5000
    iter <- 0
    success <- FALSE

    repeat {
      iter <- iter + 1
      shuffledCountsT <- t(apply(RawcountsT, 1, sample))  # transpose back

      # require no all-zero rows and no all-zero columns (>= 2 nonzeros for cols)
      ok <- all(rowSums(shuffledCountsT != 0) > 0) &&
        all(colSums(shuffledCountsT != 0) >= 2)

      if (ok) {
        success <- TRUE
        break
      }

      if (iter >= max_iter) {
        # ---- PATCH: minimally repair instead of filling with NA ----
        M <- shuffledCountsT

        # Fix all-zero rows
        zero_rows <- which(rowSums(M != 0) == 0)
        if (length(zero_rows)) {
          col_order <- order(colSums(M != 0), decreasing = FALSE)
          for (r in zero_rows) {
            take <- head(col_order, 2)
            M[r, take] <- 1
          }
        }

        # Fix columns with <2 nonzeros
        bad_cols <- which(colSums(M != 0) < 2)
        for (cj in bad_cols) {
          need <- 2 - sum(M[, cj] != 0)
          if (need > 0) {
            candidates <- which(M[, cj] == 0)
            if (length(candidates) > 0) {
              pick <- head(sample(candidates), need)
              M[pick, cj] <- 1
            }
          }
        }

        shuffledCountsT <- M
        dimnames(shuffledCountsT) <- dimnames(RawcountsT)

        success <- FALSE  # still mark as fail (unless you want to analyze patched)
        break
      }
    }

    shuffledCountsT <- data.frame(shuffledCountsT)
    colnames(shuffledCountsT) <- colnames(RawcountsT)

    if (success) {
      metagenomeseqresults <- calcMetagenomeSeq(shuffledCountsT, meta)
      metagenomeseqpvals <- cbind(metagenomeseqpvals, metagenomeseqresults$pval)
    } else {
      # If fail: append NA vector of same length as pvals
      na_pvals <- rep(NA, nrow(RawcountsT))
      metagenomeseqpvals <- cbind(metagenomeseqpvals, na_pvals)
    }
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleInTaxonDump/RDP/"
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

    max_iter <- 5000
    iter <- 0
    success <- FALSE

    repeat {
      iter <- iter + 1
      shuffledCountsT <- t(apply(RawcountsT, 1, sample))  # transpose back

      # require no all-zero rows and no all-zero columns (>= 2 nonzeros for cols)
      ok <- all(rowSums(shuffledCountsT != 0) > 0) &&
        all(colSums(shuffledCountsT != 0) >= 2)

      if (ok) {
        success <- TRUE
        break
      }

      if (iter >= max_iter) {
        # ---- PATCH: minimally repair instead of filling with NA ----
        M <- shuffledCountsT

        # Fix all-zero rows
        zero_rows <- which(rowSums(M != 0) == 0)
        if (length(zero_rows)) {
          col_order <- order(colSums(M != 0), decreasing = FALSE)
          for (r in zero_rows) {
            take <- head(col_order, 2)
            M[r, take] <- 1
          }
        }

        # Fix columns with <2 nonzeros
        bad_cols <- which(colSums(M != 0) < 2)
        for (cj in bad_cols) {
          need <- 2 - sum(M[, cj] != 0)
          if (need > 0) {
            candidates <- which(M[, cj] == 0)
            if (length(candidates) > 0) {
              pick <- head(sample(candidates), need)
              M[pick, cj] <- 1
            }
          }
        }

        shuffledCountsT <- M
        dimnames(shuffledCountsT) <- dimnames(RawcountsT)

        success <- FALSE  # still mark as fail (unless you want to analyze patched)
        break
      }
    }

    shuffledCountsT <- data.frame(shuffledCountsT)
    colnames(shuffledCountsT) <- colnames(RawcountsT)

    if (success) {
      metagenomeseqresults <- calcMetagenomeSeq(shuffledCountsT, meta)
      metagenomeseqpvals <- cbind(metagenomeseqpvals, metagenomeseqresults$pval)
    } else {
      # If fail: append NA vector of same length as pvals
      na_pvals <- rep(NA, nrow(RawcountsT))
      metagenomeseqpvals <- cbind(metagenomeseqpvals, na_pvals)
    }
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleInTaxonDump/dada2/"
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

    max_iter <- 5000
    iter <- 0
    success <- FALSE

    repeat {
      iter <- iter + 1
      shuffledCountsT <- t(apply(RawcountsT, 1, sample))  # transpose back

      # require no all-zero rows and no all-zero columns (>= 2 nonzeros for cols)
      ok <- all(rowSums(shuffledCountsT != 0) > 0) &&
        all(colSums(shuffledCountsT != 0) >= 2)

      if (ok) {
        success <- TRUE
        break
      }

      if (iter >= max_iter) {
        # ---- PATCH: minimally repair instead of filling with NA ----
        M <- shuffledCountsT

        # Fix all-zero rows
        zero_rows <- which(rowSums(M != 0) == 0)
        if (length(zero_rows)) {
          col_order <- order(colSums(M != 0), decreasing = FALSE)
          for (r in zero_rows) {
            take <- head(col_order, 2)
            M[r, take] <- 1
          }
        }

        # Fix columns with <2 nonzeros
        bad_cols <- which(colSums(M != 0) < 2)
        for (cj in bad_cols) {
          need <- 2 - sum(M[, cj] != 0)
          if (need > 0) {
            candidates <- which(M[, cj] == 0)
            if (length(candidates) > 0) {
              pick <- head(sample(candidates), need)
              M[pick, cj] <- 1
            }
          }
        }

        shuffledCountsT <- M
        dimnames(shuffledCountsT) <- dimnames(RawcountsT)

        success <- FALSE  # still mark as fail (unless you want to analyze patched)
        break
      }
    }

    shuffledCountsT <- data.frame(shuffledCountsT)
    colnames(shuffledCountsT) <- colnames(RawcountsT)

    if (success) {
      metagenomeseqresults <- calcMetagenomeSeq(shuffledCountsT, meta)
      metagenomeseqpvals <- cbind(metagenomeseqpvals, metagenomeseqresults$pval)
    } else {
      # If fail: append NA vector of same length as pvals
      na_pvals <- rep(NA, nrow(RawcountsT))
      metagenomeseqpvals <- cbind(metagenomeseqpvals, na_pvals)
    }
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleInTaxonDump/WGS/"
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

    max_iter <- 5000
    iter <- 0
    success <- FALSE

    repeat {
      iter <- iter + 1
      shuffledCountsT <- t(apply(RawcountsT, 1, sample))  # transpose back

      # require no all-zero rows and no all-zero columns (>= 2 nonzeros for cols)
      ok <- all(rowSums(shuffledCountsT != 0) > 0) &&
        all(colSums(shuffledCountsT != 0) >= 2)

      if (ok) {
        success <- TRUE
        break
      }

      if (iter >= max_iter) {
        # ---- PATCH: minimally repair instead of filling with NA ----
        M <- shuffledCountsT

        # Fix all-zero rows
        zero_rows <- which(rowSums(M != 0) == 0)
        if (length(zero_rows)) {
          col_order <- order(colSums(M != 0), decreasing = FALSE)
          for (r in zero_rows) {
            take <- head(col_order, 2)
            M[r, take] <- 1
          }
        }

        # Fix columns with <2 nonzeros
        bad_cols <- which(colSums(M != 0) < 2)
        for (cj in bad_cols) {
          need <- 2 - sum(M[, cj] != 0)
          if (need > 0) {
            candidates <- which(M[, cj] == 0)
            if (length(candidates) > 0) {
              pick <- head(sample(candidates), need)
              M[pick, cj] <- 1
            }
          }
        }

        shuffledCountsT <- M
        dimnames(shuffledCountsT) <- dimnames(RawcountsT)

        success <- FALSE  # still mark as fail (unless you want to analyze patched)
        break
      }
    }

    shuffledCountsT <- data.frame(shuffledCountsT)
    colnames(shuffledCountsT) <- colnames(RawcountsT)

    if (success) {
      metagenomeseqresults <- calcMetagenomeSeq(shuffledCountsT, meta)
      metagenomeseqpvals <- cbind(metagenomeseqpvals, metagenomeseqresults$pval)
    } else {
      # If fail: append NA vector of same length as pvals
      na_pvals <- rep(NA, nrow(RawcountsT))
      metagenomeseqpvals <- cbind(metagenomeseqpvals, na_pvals)
    }
  }

  rownames(metagenomeseqpvals)<-rownames(RawcountsT)

  save_dir<-"ShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(metagenomeseqpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_metagenomeseq.txt"), sep = "\t", row.names = T)
  }

}
