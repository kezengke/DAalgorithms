#testing normality before shuffling in sample
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(coin)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(patchwork)

NormalityKS <- function(x) {
  # Check if the standard deviation is zero
  if (sd(x) == 0) {
    return(NA)  # Can't standardize if all values are the same
  }
  # Standardize the data
  standardized_x <- (x - mean(x)) / sd(x)
  # Perform the KS test against a standard normal distribution
  ks_result <- ks.test(standardized_x, "pnorm", mean = 0, sd = 1)
  return(ks_result$p.value)
}

# Function to calculate mean and sd for each group
calculateMeanSd <- function(z) {
  g1 <- unlist(z[meta$conditions == group1])
  g2 <- unlist(z[meta$conditions == group2])
  c(mean1 = mean(g1), mean2 = mean(g2), sd1 = sd(g1), sd2 = sd(g2))
}

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  combined_plots<-NULL
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/RDP/ResampledShuffledInSampleKSNormalityTest(", name,").png"), width=1200*20, height=1200*10, res = 500)
  RawcountsT<-read.table(paste0("CountsTables/RDPRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # calculating mean var
  groups <- unique(meta$conditions)

  group1 <- groups[1]
  group2 <- groups[2]
  sample1<-rownames(meta)[meta$conditions == group1]
  sample2<-rownames(meta)[meta$conditions == group2]

  # Apply the function to each row of the table and combine results into a data frame
  MeanSd_table <- t(apply(RawcountsT, 1, calculateMeanSd))
  MeanSd_table <- as.data.frame(MeanSd_table)
  rownames(MeanSd_table) <- rownames(table)

  # resample
  RawcountsT<-resampleRNORM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 2, sample)

    shuffledCountsT1<-shuffledCountsT[, meta$conditions == group1, drop = F]
    shuffledCountsT2<-shuffledCountsT[, meta$conditions == group2, drop = F]

    shuffledCountsT1<-shuffledCountsT1[apply(shuffledCountsT1, 1, function(row) any(row != 0)), ]
    shuffledCountsT2<-shuffledCountsT2[apply(shuffledCountsT2, 1, function(row) any(row != 0)), ]

    pvals1<-data.frame(apply(shuffledCountsT1, 1, NormalityKS))
    pvals2<-data.frame(apply(shuffledCountsT2, 1, NormalityKS))

    colnames(pvals1)<-"pval"
    colnames(pvals2)<-"pval"

    p1<-PvalHistogram(pvals1, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 1)")
    p2<-PvalHistogram(pvals2, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 2)")

    if (is.null(combined_plots)) {
      combined_plots <- wrap_elements(p1)
    } else {
      combined_plots <- combined_plots + wrap_elements(p1)
      combined_plots <- combined_plots + wrap_elements(p2)
    }

  }
  print(combined_plots + plot_layout(ncol = 20))
  dev.off()
}

# DADA2
all_files<-list.files("CountsTables/DADA2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  combined_plots<-NULL
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/DADA2/ResampledShuffledInSampleKSNormalityTest(", name,").png"), width=1200*20, height=1200*10, res = 500)
  RawcountsT<-read.table(paste0("CountsTables/DADA2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # calculating mean var
  groups <- unique(meta$conditions)

  group1 <- groups[1]
  group2 <- groups[2]
  sample1<-rownames(meta)[meta$conditions == group1]
  sample2<-rownames(meta)[meta$conditions == group2]

  # Apply the function to each row of the table and combine results into a data frame
  MeanSd_table <- t(apply(RawcountsT, 1, calculateMeanSd))
  MeanSd_table <- as.data.frame(MeanSd_table)
  rownames(MeanSd_table) <- rownames(table)

  # resample
  RawcountsT<-resampleRNORM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 2, sample)

    shuffledCountsT1<-shuffledCountsT[, meta$conditions == group1, drop = F]
    shuffledCountsT2<-shuffledCountsT[, meta$conditions == group2, drop = F]

    shuffledCountsT1<-shuffledCountsT1[apply(shuffledCountsT1, 1, function(row) any(row != 0)), ]
    shuffledCountsT2<-shuffledCountsT2[apply(shuffledCountsT2, 1, function(row) any(row != 0)), ]

    pvals1<-data.frame(apply(shuffledCountsT1, 1, NormalityKS))
    pvals2<-data.frame(apply(shuffledCountsT2, 1, NormalityKS))

    colnames(pvals1)<-"pval"
    colnames(pvals2)<-"pval"

    p1<-PvalHistogram(pvals1, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 1)")
    p2<-PvalHistogram(pvals2, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 2)")

    if (is.null(combined_plots)) {
      combined_plots <- wrap_elements(p1)
    } else {
      combined_plots <- combined_plots + wrap_elements(p1)
      combined_plots <- combined_plots + wrap_elements(p2)
    }

  }
  print(combined_plots + plot_layout(ncol = 20))
  dev.off()
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  combined_plots<-NULL
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/WGS/ResampledShuffledInSampleKSNormalityTest(", name,").png"), width=1200*20, height=1200*10, res = 500)
  RawcountsT<-read.table(paste0("CountsTables/WGSRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # calculating mean var
  groups <- unique(meta$conditions)

  group1 <- groups[1]
  group2 <- groups[2]
  sample1<-rownames(meta)[meta$conditions == group1]
  sample2<-rownames(meta)[meta$conditions == group2]

  # Apply the function to each row of the table and combine results into a data frame
  MeanSd_table <- t(apply(RawcountsT, 1, calculateMeanSd))
  MeanSd_table <- as.data.frame(MeanSd_table)
  rownames(MeanSd_table) <- rownames(table)

  # resample
  RawcountsT<-resampleRNORM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 2, sample)

    shuffledCountsT1<-shuffledCountsT[, meta$conditions == group1, drop = F]
    shuffledCountsT2<-shuffledCountsT[, meta$conditions == group2, drop = F]

    shuffledCountsT1<-shuffledCountsT1[apply(shuffledCountsT1, 1, function(row) any(row != 0)), ]
    shuffledCountsT2<-shuffledCountsT2[apply(shuffledCountsT2, 1, function(row) any(row != 0)), ]

    pvals1<-data.frame(apply(shuffledCountsT1, 1, NormalityKS))
    pvals2<-data.frame(apply(shuffledCountsT2, 1, NormalityKS))

    colnames(pvals1)<-"pval"
    colnames(pvals2)<-"pval"

    p1<-PvalHistogram(pvals1, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 1)")
    p2<-PvalHistogram(pvals2, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 2)")

    if (is.null(combined_plots)) {
      combined_plots <- wrap_elements(p1)
    } else {
      combined_plots <- combined_plots + wrap_elements(p1)
      combined_plots <- combined_plots + wrap_elements(p2)
    }

  }
  print(combined_plots + plot_layout(ncol = 20))
  dev.off()
}

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  combined_plots<-NULL
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  png(paste0("Plots/RNAseq/ResampledShuffledInSampleKSNormalityTest(", name,").png"), width=1200*20, height=1200*10, res = 500)
  RawcountsT<-read.table(paste0("CountsTables/RNAseqRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # calculating mean var
  groups <- unique(meta$conditions)

  group1 <- groups[1]
  group2 <- groups[2]
  sample1<-rownames(meta)[meta$conditions == group1]
  sample2<-rownames(meta)[meta$conditions == group2]

  # Apply the function to each row of the table and combine results into a data frame
  MeanSd_table <- t(apply(RawcountsT, 1, calculateMeanSd))
  MeanSd_table <- as.data.frame(MeanSd_table)
  rownames(MeanSd_table) <- rownames(table)

  # resample
  RawcountsT<-resampleRNORM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 2, sample)

    shuffledCountsT1<-shuffledCountsT[, meta$conditions == group1, drop = F]
    shuffledCountsT2<-shuffledCountsT[, meta$conditions == group2, drop = F]

    shuffledCountsT1<-shuffledCountsT1[apply(shuffledCountsT1, 1, function(row) any(row != 0)), ]
    shuffledCountsT2<-shuffledCountsT2[apply(shuffledCountsT2, 1, function(row) any(row != 0)), ]

    pvals1<-data.frame(apply(shuffledCountsT1, 1, NormalityKS))
    pvals2<-data.frame(apply(shuffledCountsT2, 1, NormalityKS))

    colnames(pvals1)<-"pval"
    colnames(pvals2)<-"pval"

    p1<-PvalHistogram(pvals1, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 1)")
    p2<-PvalHistogram(pvals2, "#df676b", title = "Resampled Shuffled In Sample KS Test \n(condition 2)")

    if (is.null(combined_plots)) {
      combined_plots <- wrap_elements(p1)
    } else {
      combined_plots <- combined_plots + wrap_elements(p1)
      combined_plots <- combined_plots + wrap_elements(p2)
    }

  }
  print(combined_plots + plot_layout(ncol = 20))
  dev.off()
}
