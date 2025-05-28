#shuffle counts in every taxon 100x and test all pval for resampled counts tables
rm(list = ls())
set.seed(123)
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

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    # resample
    ResampledcountsT<-resampleRNORM(shuffledCountsT, meta, 1)
    colnames(ResampledcountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(ResampledcountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(ResampledcountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(ResampledcountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  save_dir<-"ResampleShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/RDP/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/RDP/"
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

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    # resample
    ResampledcountsT<-resampleRNORM(shuffledCountsT, meta, 1)
    colnames(ResampledcountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(ResampledcountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(ResampledcountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(ResampledcountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  save_dir<-"ResampleShuffleInTaxonDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/dada2/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/dada2/"
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

  RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    # resample
    ResampledcountsT<-resampleRNORM(shuffledCountsT, meta, 1)
    colnames(ResampledcountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(ResampledcountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(ResampledcountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(ResampledcountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  save_dir<-"ResampleShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/WGS/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/WGS/"
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

  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)
  meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(RawcountsT)
  colnames(meta)<-"conditions"

  # Filter out low counts taxa
  if (all(which(rowMeans(RawcountsT)<2) == 0)) {
    RawcountsT<-RawcountsT
  } else {
    RawcountsT<-RawcountsT[-c(which(rowMeans(RawcountsT)<2)), , drop = F]
  }

  # Remove taxa or gene with 0 variance
  meanStd<-calcMeanStd(RawcountsT, meta)

  index_with_0std1<-which(meanStd$sd1 == 0)
  index_with_0std2<-which(meanStd$sd2 == 0)
  all_index<-unique(c(index_with_0std1, index_with_0std2))

  if (length(all_index) > 0) {
    RawcountsT <- RawcountsT[-all_index, , drop = FALSE]
  }

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    shuffledCountsT<-apply(RawcountsT, 1, sample)
    shuffledCountsT<-data.frame(t(shuffledCountsT))
    # resample
    ResampledcountsT<-resampleRNORM(shuffledCountsT, meta, 1)
    colnames(ResampledcountsT)<-colnames(RawcountsT)
    NormCountsT<-normFun(ResampledcountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(ResampledcountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(ResampledcountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(ResampledcountsT)
  rownames(wilcoxpvals)<-rownames(ResampledcountsT)
  rownames(deseq2pvals)<-rownames(ResampledcountsT)
  rownames(edgerpvals)<-rownames(ResampledcountsT)

  #set NA pvals to 1
  ttestpvals[is.na(ttestpvals)]<-1
  wilcoxpvals[is.na(wilcoxpvals)]<-1

  save_dir<-"ResampleShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleShuffleInTaxonDump/RNAseq/"
  if (file.exists(file.path(save_dir))){
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}
