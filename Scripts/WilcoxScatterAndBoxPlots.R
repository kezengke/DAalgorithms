# Running Assal dataset through resamping and shuffle counts in sample pipeline
rm(list = ls())
library(coin)
library(DESeq2)
library(edgeR)
library(ggplot2)

set.seed(26)

# Function for normalization
normFun <- function(table) {
  n <- colSums(table)
  sumx <- sum(table)
  for (j in 1:ncol(table)) {
    table[, j] <- table[, j]/n[j]
  }
  table <- log10(table * (sumx/ncol(table)) + 1)
  table <- data.frame(table, check.names = F)
  return(table)
}

# Function for calculating T-test results
calcTtest <- function(table, meta) {
  t_stats <- apply(table, 1, function(x) {
    t.test(unlist(x) ~ meta$conditions)$stat
  })
  t_test_p <- apply(table, 1, function(x) {
    t.test(unlist(x) ~ meta$conditions)$p.value
  })
  t_results <- cbind(t_stats, t_test_p)
  rownames(t_results) <- rownames(table)
  colnames(t_results) <- c("stats", "pval")
  t_results <- data.frame(t_results, check.names = F)
  return(t_results)
}

# Function for calculating Wilcoxon test results
calcWilcox <- function(table, meta) {
  wilcox_stats <- apply(table, 1, function(x) {
    statistic(wilcox_test(unlist(x) ~ factor(meta$conditions)))
  })
  wilcox_p <- apply(table, 1, function(x) {
    pvalue(wilcox_test(unlist(x) ~ factor(meta$conditions)))
  })
  wilcox_results <- cbind(wilcox_stats, wilcox_p)
  rownames(wilcox_results) <- rownames(table)
  colnames(wilcox_results) <- c("stats", "pval")
  wilcox_results <- data.frame(wilcox_results, check.names = F)
  return(wilcox_results)
}

# Function to resample each taxon for the entire counts table
resampleRNORM <- function(table, meta, multiple)
{
  if (nrow(table) == 0 || nrow(meta) == 0) {
    stop("Input table or meta data frame is empty.")
  }
  groups <- unique(meta$conditions)
  group1 <- groups[1]
  group2 <- groups[2]
  # function for calculating mean and std for each group of one taxon
  calculateMeanSd <- function(z) {
    g1 <- unlist(z[meta$conditions == group1])
    g2 <- unlist(z[meta$conditions == group2])
    c(mean1 = mean(g1), mean2 = mean(g2), sd1 = sd(g1), sd2 = sd(g2))
  }
  # Calculate mean and std for each group for the entire counts table
  MeanSd_table <- t(apply(table, 1, calculateMeanSd))
  MeanSd_table <- as.data.frame(MeanSd_table)
  rownames(MeanSd_table) <- rownames(table)
  # Function for resampling the a taxon for each group using the corresponding meand and std
  resample_counts <- function(row_index) {
    z <- table[row_index, ]
    g1 <- unlist(z[meta$conditions == group1])
    g2 <- unlist(z[meta$conditions == group2])
    # Multiple is for how many times to multiply the original std for std/variance increament experiemnts
    ng1 <- rnorm(n = length(g1), mean = MeanSd_table[row_index,
                                                     "mean1"], sd = sqrt(multiple * (MeanSd_table[row_index,
                                                                                                  "sd1"])^2))
    ng2 <- rnorm(n = length(g2), mean = MeanSd_table[row_index,
                                                     "mean2"], sd = sqrt(multiple * (MeanSd_table[row_index,
                                                                                                  "sd2"])^2))
    c(ng1, ng2)
  }
  # Resampling
  newT <- t(sapply(seq_len(nrow(table)), resample_counts))
  colnames(newT) <- c(rownames(meta)[meta$conditions == group1],
                      rownames(meta)[meta$conditions == group2])
  newT <- newT[, colnames(table)]
  rownames(newT) <- rownames(table)
  # Make all the counts >= 0
  newT <- newT + abs(min(newT))
  newT <- data.frame(round(newT, digits = 0), check.names = F)
  return(newT)
}

# Function for calculating fraction of significant counts
calcFraction <- function(pvalT) {
  fractions<-apply(pvalT, 2,
                   function(column){
                     sum(column < 0.05)/length(column)
                   })
  return(fractions)
}

###############################################################
countsT<-read.table("CountsTables/RDPRaw/assal.txt", sep = "\t", header = T, row.names = 1, check.names = F)
meta<-read.table("MetaData/metadata_assal.txt",
                 header = T, row.names = 1)

# Filter out low counts taxa
if (all(which(rowMeans(countsT)<2) == 0)) {
  countsT<-countsT
} else {
  countsT<-countsT[-c(which(rowMeans(countsT)<2)), , drop = F]
}

countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

rownames(meta)<-colnames(countsT)
colnames(meta)<-"conditions"


# countsT_resampled<-resampleRNORM(countsT, meta, 1)

# shuffle counts in each sample for 100 times and calculate results of each shuffled counts table
ttestpvals<-c()
wilcoxpvals<-c()
for (i in 1:100) {
  # Resampling counts table
  countsT_resampled<-resampleRNORM(countsT, meta, 1)
  # Shuffling counts in each sample (column wise)
  shuffledCountsT<-apply(countsT_resampled, 2, sample)
  # Normalizing shuffled counts table
  NormCountsT<-normFun(shuffledCountsT)
  write.table(NormCountsT, paste0("ResampledShuffledCountsTDump/", i, ".txt"), sep = "\t")

  ttestresults<-calcTtest(NormCountsT, meta)
  ttestpvals<-cbind(ttestpvals, ttestresults$pval)
  wilcoxresults<-calcWilcox(NormCountsT, meta)
  wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
}

sigW<-calcFraction(wilcoxpvals)
moreSig<-which(sigW>0.05)
pdf("Plots/CountsScatter.pdf", height = 21, width = 14)
par(mfrow = c(3, 2))
for (i in 1:length(sigW)) {
  NcountsT<-read.table(paste0("ResampledShuffledCountsTDump/", i, ".txt"),
             header = T, row.names = 1, sep = "\t")
  rownames(NcountsT)<-rownames(countsT)
  NcountsT<-data.frame(NcountsT)
  if(sigW[i] >= 0.05)
  {
    myCol<-"red"
  }
  else
  {
    myCol<-"blue"
  }
  stripchart(NcountsT,
             col = myCol, method = "jitter", vertical = F, pch = 16,
             xlab = "Counts", ylab = "Samples",
             main = paste0("Fraction of significant results: ", sigW[i]))
}

dev.off()

for(i in 1:length(moreSig)){
  NcountsT<-read.table(paste0("ResampledShuffledCountsTDump/", moreSig[i], ".txt"),
                       header = T, row.names = 1, sep = "\t")
  rownames(NcountsT)<-rownames(countsT)
  NcountsT<-data.frame(NcountsT)

  pvals<-wilcoxpvals[, moreSig[i], drop = F]
  plotOrder<-order(pvals)

  NcountsT<-NcountsT[plotOrder, , drop = F]
  pdf(paste0("Plots/WilcoxBoxPlots/", i, "-MoreSigBoxPlots.pdf"), height = 21, width = 14)
  par(mfrow = c(3, 2))
  for (j in 1:nrow(NcountsT)) {
    boxplot(unlist(NcountsT[j, ])~meta$conditions,
            main = paste0(rownames(NcountsT)[j], "\n P-val:", pvals[plotOrder[j]]),
            xlab = "Conditions", ylab = "Counts(Normalized)",
            outline = F)
  }
  dev.off()

}

