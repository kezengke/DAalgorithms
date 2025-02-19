#resampled shuffled significant pvalue visulization
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
    scale_fill_manual(values = c("#2aa7de", "#25377f", "#ca0e12", "#f6bd21")) +
    theme_classic() +
    labs(title = paste0("(", classifier, "-", name, ")\n", shuffleType),
         x = "DAA methods",
         y = "Fraction of significant results") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

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

png(paste0("Plots/RDP/ResampledShuffledFractionOfSignificantResults(RDP).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ResampleShuffleTagDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleTagDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleTagDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleTagDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInTaxonDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInTaxonDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInTaxonDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInTaxonDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleCountsTDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleCountsTDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleCountsTDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleCountsTDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

png(paste0("Plots/DADA2/ResampledShuffledFractionOfSignificantResults(DADA2).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ResampleShuffleTagDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleTagDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleTagDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleTagDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInSampleDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInSampleDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInSampleDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInSampleDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInTaxonDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInTaxonDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInTaxonDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInTaxonDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleCountsTDump/DADA2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleCountsTDump/DADA2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleCountsTDump/DADA2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleCountsTDump/DADA2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

png(paste0("Plots/WGS/ResampledShuffledFractionOfSignificantResults(WGS).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ResampleShuffleTagDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleTagDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleTagDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleTagDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInSampleDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInSampleDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInSampleDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInSampleDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInTaxonDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInTaxonDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInTaxonDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInTaxonDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleCountsTDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleCountsTDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleCountsTDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleCountsTDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

png(paste0("Plots/RNAseq/ResampledShuffledFractionOfSignificantResults(RNAseq).png"), width=1200*length(all_files), height=4800, res = 300)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ResampleShuffleTagDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleTagDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleTagDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleTagDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInSampleDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInSampleDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInSampleDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInSampleDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleInTaxonDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleInTaxonDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleInTaxonDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleInTaxonDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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

  ttestpvals<-read.table(paste0("ResampleShuffleCountsTDump/RNAseq/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ResampleShuffleCountsTDump/RNAseq/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ResampleShuffleCountsTDump/RNAseq/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ResampleShuffleCountsTDump/RNAseq/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

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
