#shuffle sample and re-run results
rm(list = ls())
library(MetagenomeTools)

# RDP
rm(list = ls())
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  meta$conditions<-sample(meta$conditions)

  results<-calcEdgeR(countsT, meta)


  main_dir<-"PkgResults/RDP"
  sub_dir<-"shufflededgeR/"
  if (file.exists(file.path(main_dir, sub_dir))){
    write.table(results, paste0(file.path(main_dir, sub_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(main_dir, sub_dir))
    write.table(results, paste0(file.path(main_dir, sub_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}

# dada2
rm(list = ls())
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  meta$conditions<-sample(meta$conditions)

  results<-calcEdgeR(countsT, meta)

  main_dir<-"PkgResults/dada2"
  sub_dir<-"shufflededgeR/"
  if (file.exists(file.path(main_dir, sub_dir))){
    write.table(results, paste0(file.path(main_dir, sub_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(main_dir, sub_dir))
    write.table(results, paste0(file.path(main_dir, sub_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}

# WGS
rm(list = ls())
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  meta$conditions<-sample(meta$conditions)

  results<-calcEdgeR(countsT, meta)

  main_dir<-"PkgResults/WGS"
  sub_dir<-"shufflededgeR/"
  if (file.exists(file.path(main_dir, sub_dir))){
    write.table(results, paste0(file.path(main_dir, sub_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  } else {
    dir.create(file.path(main_dir, sub_dir))
    write.table(results, paste0(file.path(main_dir, sub_dir), gsub(basename(file), pattern=".txt$", replacement=""), "_edgeR.txt"), sep = "\t", row.names = T)
  }
}
