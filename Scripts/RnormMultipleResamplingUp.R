rm(list = ls())
library(MetagenomeTools)
library(DESeq2)
library(dplyr)
set.seed(9527)

times<-c(0.5, 2, seq(5, 45, 5))

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  countsT<-LoadCountsT(file)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  for (m in times) {
    newT<-resampleRNORM(countsT, meta, m)
    tRES<-BackTestTtest(countsT, newT, meta)
    dRES<-BackTestDESeq2(countsT, newT, meta)
    p<-directPFun(tRES, dRES)
    colnames(p)<-c("X1", "X2")
    rownames(p)<-rownames(countsT)
    write.table(p, paste0("ResampleDump/RDP/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", row.names = T)
  }
}

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  countsT<-LoadCountsT(file)
  meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  for (m in times) {
    newT<-resampleRNORM(countsT, meta, m)
    tRES<-BackTestTtest(countsT, newT, meta)
    dRES<-BackTestDESeq2(countsT, newT, meta)
    p<-directPFun(tRES, dRES)
    colnames(p)<-c("X1", "X2")
    rownames(p)<-rownames(countsT)
    write.table(p, paste0("ResampleDump/dada2/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", row.names = T)
  }
}

# # WGS
# all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)
#
# for (file in all_files) {
#   name <- gsub(basename(file), pattern=".txt$", replacement="")
#   countsT<-LoadCountsT(file)
#   meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))
#
#   countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
#   meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]
#
#   rownames(meta)<-colnames(countsT)
#   colnames(meta)<-"conditions"
#
#   for (m in times) {
#     newT<-resampleRNORM(countsT, meta, m)
#     tRES<-BackTestTtest(countsT, newT, meta)
#     dRES<-BackTestDESeq2(countsT, newT, meta)
#     p<-directPFun(tRES, dRES)
#     colnames(p)<-c("X1", "X2")
#     rownames(p)<-rownames(countsT)
#     write.table(p, paste0("ResampleDump/WGS/rnormUp/", name, "/x", m, "Resampled-P.txt"), sep = "\t", row.names = T)
#   }
# }
