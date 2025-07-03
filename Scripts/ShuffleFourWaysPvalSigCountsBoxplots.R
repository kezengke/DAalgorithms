#shuffled sample tags significant pvalue visulization
rm(list = ls())
library(MetagenomeTools)
library(ggExtra)
library(patchwork)
library(ggplot2)

makePlot <- function(fractT, classifier, name, shuffleType) {

  df_long <- stack(fractT)

  p<-ggplot(df_long, aes(x = ind, y = values, fill = ind)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("coral1", "lightslateblue", "olivedrab3", "goldenrod1")) +
    theme_classic() +
    labs(title = paste0("(", classifier, "-", name, ")\n", shuffleType),
         x = "DAA methods",
         y = "Fraction of significant results") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    ylim(0, 0.75)

  return(p)
}

calcFraction <- function(pvalT) {
  fractions<-apply(pvalT, 2,
                   function(column){
    sum(column < 0.05)/length(column)
  })
  return(fractions)
}

# RDP
combined_plots<-NULL
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RDP/FractionOfSignificantResults(RDP).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RDP", name, "Shuffle sample tags")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInSampleDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInSampleDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInSampleDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInSampleDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RDP", name, "Shuffle counts in sample")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInTaxonDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInTaxonDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInTaxonDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInTaxonDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RDP", name, "Shuffle counts in taxon")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleCountsTDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RDP", name, "Shuffle counts table")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

print(combined_plots + plot_layout(ncol = length(all_files)))
dev.off()

# DADA2
combined_plots<-NULL
all_files<-list.files("CountsTables/DADA2Raw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/DADA2/FractionOfSignificantResults(DADA2).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "DADA2", name, "Shuffle sample tags")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInSampleDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInSampleDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInSampleDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInSampleDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "DADA2", name, "Shuffle counts in sample")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInTaxonDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInTaxonDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInTaxonDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInTaxonDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "DADA2", name, "Shuffle counts in taxon")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleCountsTDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleCountsTDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleCountsTDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleCountsTDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "DADA2", name, "Shuffle counts table")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

print(combined_plots + plot_layout(ncol = length(all_files)))
dev.off()

# WGS
combined_plots<-NULL
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/WGS/FractionOfSignificantResults(WGS).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "WGS", name, "Shuffle sample tags")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInSampleDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInSampleDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInSampleDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInSampleDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "WGS", name, "Shuffle counts in sample")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInTaxonDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInTaxonDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInTaxonDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInTaxonDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "WGS", name, "Shuffle counts in taxon")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleCountsTDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "WGS", name, "Shuffle counts table")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

print(combined_plots + plot_layout(ncol = length(all_files)))
dev.off()

# RNAseq
combined_plots<-NULL
all_files<-list.files("CountsTables/RNAseqRaw", pattern = "*.txt",full.names = TRUE)

png(paste0("Plots/RNAseq/FractionOfSignificantResults(RNAseq).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RNAseq", name, "Shuffle sample tags")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInSampleDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInSampleDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInSampleDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInSampleDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RNAseq", name, "Shuffle counts in sample")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleInTaxonDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleInTaxonDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleInTaxonDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleInTaxonDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RNAseq", name, "Shuffle counts in taxon")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleCountsTDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleCountsTDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleCountsTDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleCountsTDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  # significant fraction
  sigT<-calcFraction(ttestpvals)
  sigW<-calcFraction(wilcoxpvals)
  sigD<-calcFraction(deseq2pvals)
  sigE<-calcFraction(edgerpvals)

  fractionT<-cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT<-data.frame(fractionT, check.names = F)

  p<-makePlot(fractionT, "RNAseq", name, "Shuffle counts table")

  if (is.null(combined_plots)) {
    combined_plots <- wrap_elements(p)
  } else {
    combined_plots <- combined_plots + wrap_elements(p)
  }

}

print(combined_plots + plot_layout(ncol = length(all_files)))
dev.off()
