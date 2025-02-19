# WGS
rm(list = ls())
library(MetagenomeTools)
library(DESeq2)
library(dplyr)
set.seed(9527)

file<-"CountsTables/WGSRaw/stoolswab.txt"

m<-0.5

name <- gsub(basename(file), pattern=".txt$", replacement="")
countsT<-LoadCountsT(file)
meta<-LoadMeta(paste0("MetaData/metadata_", name, ".txt"))

countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

rownames(meta)<-colnames(countsT)
colnames(meta)<-"conditions"

newT<-resampleRNORM(countsT, meta, m)
tRES<-BackTestTtest(countsT, newT, meta)
dRES<-BackTestDESeq2(countsT, newT, meta)
p<-directPFun(tRES, dRES)
colnames(p)<-c("X1", "X2")
rownames(p)<-rownames(countsT)

