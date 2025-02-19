rm(list = ls())
library(ggplot2)
library(patchwork)
library(ALDEx2)
library(MetagenomeTools)

testTaxon<-"Unclassified.Burkholderiales.Order"
rawT<-read.table("CountsTables/RDPRaw/rhizo.txt", header = T, row.names = 1, sep = "\t")
meta<-read.table("MetaData/metadata_rhizo.txt", header = T, row.names = 1, sep = "\t")


rawT<-rawT[, intersect(colnames(rawT), rownames(meta)), drop = F]
meta<-meta[intersect(colnames(rawT), rownames(meta)), , drop = F]

rownames(meta)<-colnames(rawT)
colnames(meta)<-"conditions"

all_Tpvals<-c()
all_Wpvals<-c()
for (i in 1:100) {
  TtestRES<-calcALDEx2Ttest(rawT, meta)
  pval<-TtestRES[testTaxon, 2]
  all_Tpvals<-c(all_Tpvals, pval)

  WilcoxRES<-calcALDEx2Wilcoxon(rawT, meta)
  pval<-WilcoxRES[testTaxon, 2]
  all_Wpvals<-c(all_Wpvals, pval)
}

norm_Tpvals<-log10(all_Tpvals)
norm_Tpvals<-norm_Tpvals * -1

norm_Wpvals<-log10(all_Wpvals)
norm_Wpvals<-norm_Wpvals * -1

p1<-ggplot(data = data.frame(x = norm_Tpvals), aes(x = x)) +
  geom_histogram(binwidth = 0.05, fill = "pink", color = "black", alpha = 0.7) +
  labs(x = "-log10(p-value)", y = "Frequency", title = paste(testTaxon, "\n100x ALDEx2 t-test P-values")) +
  theme_minimal()

p2<-ggplot(data = data.frame(x = norm_Wpvals), aes(x = x)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(x = "-log10(p-value)", y = "Frequency", title = paste(testTaxon, "\n100x ALDEx2 Wilcoxon P-values")) +
  theme_minimal()

png("Plots/RDP/100xALDEx2PvalHistogram.png", width= 3000, height=1200, res = 300)
combined_plots <- p1 + p2
print(combined_plots + plot_layout(ncol = 2))
dev.off()
