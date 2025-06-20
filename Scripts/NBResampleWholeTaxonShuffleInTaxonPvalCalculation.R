#shuffle counts in every taxon 100x and test all pval for resampled whole taxon of counts tables (negative binomial distribution)
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(coin)
library(DESeq2)
library(edgeR)

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/RDPRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  # Filter out low counts taxa
  if (all(which(rowMeans(RawcountsT)<2) == 0)) {
    RawcountsT<-RawcountsT
  } else {
    RawcountsT<-RawcountsT[-c(which(rowMeans(RawcountsT)<2)), , drop = F]
  }
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # resample
  RawcountsT<-resampleWholeTaxonRNBINOM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    colnames(shuffledCountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(shuffledCountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(shuffledCountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(shuffledCountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ttestpvals[is.na(ttestpvals)]<-1
  wilcoxpvals[is.na(wilcoxpvals)]<-1

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/dada2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  # Filter out low counts taxa
  if (all(which(rowMeans(RawcountsT)<2) == 0)) {
    RawcountsT<-RawcountsT
  } else {
    RawcountsT<-RawcountsT[-c(which(rowMeans(RawcountsT)<2)), , drop = F]
  }
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # resample
  RawcountsT<-resampleWholeTaxonRNBINOM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    colnames(shuffledCountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(shuffledCountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(shuffledCountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(shuffledCountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ttestpvals[is.na(ttestpvals)]<-1
  wilcoxpvals[is.na(wilcoxpvals)]<-1

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/WGSRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  # Filter out low counts taxa
  if (all(which(rowMeans(RawcountsT)<2) == 0)) {
    RawcountsT<-RawcountsT
  } else {
    RawcountsT<-RawcountsT[-c(which(rowMeans(RawcountsT)<2)), , drop = F]
  }
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # resample
  RawcountsT<-resampleWholeTaxonRNBINOM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    colnames(shuffledCountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(shuffledCountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(shuffledCountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(shuffledCountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ttestpvals[is.na(ttestpvals)]<-1
  wilcoxpvals[is.na(wilcoxpvals)]<-1

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}

# RNAseq
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  RawcountsT<-read.table(paste0("CountsTables/RNAseqRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
  # Filter out low counts taxa
  if (all(which(rowMeans(RawcountsT)<2) == 0)) {
    RawcountsT<-RawcountsT
  } else {
    RawcountsT<-RawcountsT[-c(which(rowMeans(RawcountsT)<2)), , drop = F]
  }
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # resample
  RawcountsT<-resampleWholeTaxonRNBINOM(RawcountsT, meta, 1)

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    colnames(shuffledCountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(shuffledCountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(shuffledCountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(shuffledCountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ttestpvals[is.na(ttestpvals)]<-1
  wilcoxpvals[is.na(wilcoxpvals)]<-1

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"NBResampleWholeTaxonShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}
