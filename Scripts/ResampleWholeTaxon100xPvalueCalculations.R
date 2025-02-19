#resample counts in every taxon ignoring grouping 100x and test all pval
rm(list = ls())
set.seed(9527)
library(MetagenomeTools)
library(coin)
library(DESeq2)
library(edgeR)

# # RDP
# all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)
#
# for (file in all_files) {
#   name <- gsub(basename(file), pattern=".txt$", replacement="")
#   RawcountsT<-read.table(paste0("CountsTables/RDPRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
#   meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
#                    header = T, row.names = 1)
#
#   RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
#
#   meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]
#
#   rownames(meta)<-colnames(RawcountsT)
#   colnames(meta)<-"conditions"
#
#   ttestpvals<-c()
#   wilcoxpvals<-c()
#   deseq2pvals<-c()
#   edgerpvals<-c()
#   for (i in 1:100) {
#     newCountsT<-resampleWholeTaxonRNORM(RawcountsT, meta, 1)
#
#     NormCountsT<-normFun(newCountsT)
#
#     ttestresults<-calcTtest(NormCountsT, meta)
#     ttestpvals<-cbind(ttestpvals, ttestresults$pval)
#     wilcoxresults<-calcWilcox(NormCountsT, meta)
#     wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
#     deseq2results<-calcDESeq2(newCountsT, meta)
#     deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
#     edgerresults<-calcEdgeR(newCountsT, meta)
#     edgerpvals<-cbind(edgerpvals, edgerresults$pval)
#   }
#
#   rownames(ttestpvals)<-rownames(RawcountsT)
#   rownames(wilcoxpvals)<-rownames(RawcountsT)
#   rownames(deseq2pvals)<-rownames(RawcountsT)
#   rownames(edgerpvals)<-rownames(RawcountsT)
#
#   save_dir<-"ResampleDump/RDP/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/RDP/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/RDP/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/RDP/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
#   }
# }
#
# # dada2
# all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)
#
# for (file in all_files) {
#   name <- gsub(basename(file), pattern=".txt$", replacement="")
#   RawcountsT<-read.table(paste0("CountsTables/dada2Raw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
#   meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
#                    header = T, row.names = 1)
#
#   RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
#
#   meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]
#
#   rownames(meta)<-colnames(RawcountsT)
#   colnames(meta)<-"conditions"
#
#   ttestpvals<-c()
#   wilcoxpvals<-c()
#   deseq2pvals<-c()
#   edgerpvals<-c()
#   for (i in 1:100) {
#     newCountsT<-resampleWholeTaxonRNORM(RawcountsT, meta, 1)
#
#     NormCountsT<-normFun(newCountsT)
#
#     ttestresults<-calcTtest(NormCountsT, meta)
#     ttestpvals<-cbind(ttestpvals, ttestresults$pval)
#     wilcoxresults<-calcWilcox(NormCountsT, meta)
#     wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
#     deseq2results<-calcDESeq2(newCountsT, meta)
#     deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
#     edgerresults<-calcEdgeR(newCountsT, meta)
#     edgerpvals<-cbind(edgerpvals, edgerresults$pval)
#   }
#
#   rownames(ttestpvals)<-rownames(RawcountsT)
#   rownames(wilcoxpvals)<-rownames(RawcountsT)
#   rownames(deseq2pvals)<-rownames(RawcountsT)
#   rownames(edgerpvals)<-rownames(RawcountsT)
#
#   save_dir<-"ResampleDump/dada2/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/dada2/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/dada2/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/dada2/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
#   }
# }
#
# # WGS
# all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)
#
# for (file in all_files) {
#   name <- gsub(basename(file), pattern=".txt$", replacement="")
#   RawcountsT<-read.table(paste0("CountsTables/WGSRaw/", name, ".txt"), sep = "\t", header = T, row.names = 1, check.names = F)
#   meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
#                    header = T, row.names = 1)
#
#   RawcountsT<-RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
#
#   meta<-meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]
#
#   rownames(meta)<-colnames(RawcountsT)
#   colnames(meta)<-"conditions"
#
#   ttestpvals<-c()
#   wilcoxpvals<-c()
#   deseq2pvals<-c()
#   edgerpvals<-c()
#   for (i in 1:100) {
#     newCountsT<-resampleWholeTaxonRNORM(RawcountsT, meta, 1)
#
#     NormCountsT<-normFun(newCountsT)
#
#     ttestresults<-calcTtest(NormCountsT, meta)
#     ttestpvals<-cbind(ttestpvals, ttestresults$pval)
#     wilcoxresults<-calcWilcox(NormCountsT, meta)
#     wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
#     deseq2results<-calcDESeq2(newCountsT, meta)
#     deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
#     edgerresults<-calcEdgeR(newCountsT, meta)
#     edgerpvals<-cbind(edgerpvals, edgerresults$pval)
#   }
#
#   rownames(ttestpvals)<-rownames(RawcountsT)
#   rownames(wilcoxpvals)<-rownames(RawcountsT)
#   rownames(deseq2pvals)<-rownames(RawcountsT)
#   rownames(edgerpvals)<-rownames(RawcountsT)
#
#   save_dir<-"ResampleDump/WGS/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/WGS/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/WGS/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
#   }
#
#   save_dir<-"ResampleDump/WGS/rnormWholeTaxon/"
#   if (file.exists(file.path(save_dir))){
#     write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
#   } else {
#     dir.create(file.path(save_dir))
#     write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
#   }
# }

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

  ttestpvals<-c()
  wilcoxpvals<-c()
  deseq2pvals<-c()
  edgerpvals<-c()
  for (i in 1:100) {
    newCountsT<-resampleWholeTaxonRNORM(RawcountsT, meta, 1)

    NormCountsT<-normFun(newCountsT)

    ttestresults<-calcTtest(NormCountsT, meta)
    ttestpvals<-cbind(ttestpvals, ttestresults$pval)
    wilcoxresults<-calcWilcox(NormCountsT, meta)
    wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results<-calcDESeq2(newCountsT, meta)
    deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)
    edgerresults<-calcEdgeR(newCountsT, meta)
    edgerpvals<-cbind(edgerpvals, edgerresults$pval)
  }

  rownames(ttestpvals)<-rownames(RawcountsT)
  rownames(wilcoxpvals)<-rownames(RawcountsT)
  rownames(deseq2pvals)<-rownames(RawcountsT)
  rownames(edgerpvals)<-rownames(RawcountsT)

  #set NA pvals to 1
  ttestpvals[is.na(ttestpvals)]<-1
  wilcoxpvals[is.na(wilcoxpvals)]<-1

  save_dir<-"ResampleDump/RNAseq/rnormWholeTaxon/"
  if (file.exists(file.path(save_dir))){
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(ttestpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleDump/RNAseq/rnormWholeTaxon/"
  if (file.exists(file.path(save_dir))){
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(wilcoxpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleDump/RNAseq/rnormWholeTaxon/"
  if (file.exists(file.path(save_dir))){
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(deseq2pvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_deseq2.txt"), sep = "\t", row.names = T)
  }

  save_dir<-"ResampleDump/RNAseq/rnormWholeTaxon/"
  if (file.exists(file.path(save_dir))){
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(save_dir))
    write.table(edgerpvals, paste0(file.path(save_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}

