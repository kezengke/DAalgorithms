rm(list = ls())
library(MetagenomeTools)

# normalize rdp countsT
rm(list = ls())
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {

  countsT<-read.table(file, sep = "\t", header = T, row.names = 1)
  normT<-normFun(countsT)

  write.table(normT, paste0("CountsTables/RDPNorm/", basename(file)), sep = "\t", row.names = T)
}

# normalize DADA2 countsT
rm(list = ls())
all_files<-list.files("CountsTables/dadaRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {

  countsT<-read.table(file, sep = "\t", header = T, row.names = 1)
  normT<-normFun(countsT)

  write.table(normT, paste0("CountsTables/dadaNorm/", basename(file)), sep = "\t", row.names = T)
}

# normalize wgs countsT
rm(list = ls())
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {

  countsT<-read.table(file, sep = "\t", header = T, row.names = 1)
  normT<-normFun(countsT)

  write.table(normT, paste0("CountsTables/WGSNorm/", basename(file)), sep = "\t", row.names = T)
}