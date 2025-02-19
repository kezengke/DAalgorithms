#shuffle sample 100x and test all pval distribution
rm(list = ls())
library(MetagenomeTools)
library(patchwork)
library(ggplot2)

# RDP
all_files<-list.files("CountsTables/RDPRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleDump/RDP/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleDump/RDP/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  combined_plots<-NULL
  for (i in 1:ncol(ttestpvals)) {
    results<-data.frame(ttestpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "red", paste0("(RDP-", name, ") t-test"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/RDP/ShuffledHistograms(", name, "-ttest).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(wilcoxpvals)) {
    results<-data.frame(wilcoxpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "tan2", paste0("(RDP-", name, ") Wilcoxon"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/RDP/ShuffledHistograms(", name, "-Wilcoxon).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(deseq2pvals)) {
    results<-data.frame(deseq2pvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "purple", paste0("(RDP-", name, ") DESeq2"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/RDP/ShuffledHistograms(", name, "-DESeq2).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(edgerpvals)) {
    results<-data.frame(edgerpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "cornflowerblue", paste0("(RDP-", name, ") edgeR"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/RDP/ShuffledHistograms(", name, "-edgeR).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

}

# dada2
all_files<-list.files("CountsTables/dada2Raw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleDump/dada2/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleDump/dada2/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleDump/dada2/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleDump/dada2/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  combined_plots<-NULL
  for (i in 1:ncol(ttestpvals)) {
    results<-data.frame(ttestpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "red", paste0("(dada2-", name, ") t-test"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/dada2/ShuffledHistograms(", name, "-ttest).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(wilcoxpvals)) {
    results<-data.frame(wilcoxpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "tan2", paste0("(dada2-", name, ") Wilcoxon"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/dada2/ShuffledHistograms(", name, "-Wilcoxon).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(deseq2pvals)) {
    results<-data.frame(deseq2pvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "purple", paste0("(dada2-", name, ") DESeq2"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/dada2/ShuffledHistograms(", name, "-DESeq2).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(edgerpvals)) {
    results<-data.frame(edgerpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "cornflowerblue", paste0("(dada2-", name, ") edgeR"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/dada2/ShuffledHistograms(", name, "-edgeR).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

}

# WGS
all_files<-list.files("CountsTables/WGSRaw", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")

  ttestpvals<-read.table(paste0("ShuffleDump/WGS/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
  wilcoxpvals<-read.table(paste0("ShuffleDump/WGS/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")
  deseq2pvals<-read.table(paste0("ShuffleDump/WGS/", name, "_deseq2.txt"), header = T, row.names = 1, sep = "\t")
  edgerpvals<-read.table(paste0("ShuffleDump/WGS/", name, "_edgeR.txt"), header = T, row.names = 1, sep = "\t")

  combined_plots<-NULL
  for (i in 1:ncol(ttestpvals)) {
    results<-data.frame(ttestpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "red", paste0("(WGS-", name, ") t-test"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/WGS/ShuffledHistograms(", name, "-ttest).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(wilcoxpvals)) {
    results<-data.frame(wilcoxpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "tan2", paste0("(WGS-", name, ") Wilcoxon"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/WGS/ShuffledHistograms(", name, "-Wilcoxon).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(deseq2pvals)) {
    results<-data.frame(deseq2pvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "purple", paste0("(WGS-", name, ") DESeq2"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/WGS/ShuffledHistograms(", name, "-DESeq2).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

  combined_plots<-NULL
  for (i in 1:ncol(edgerpvals)) {
    results<-data.frame(edgerpvals[,i])
    colnames(results)<-"pval"
    p<-PvalHistogram(results, "cornflowerblue", paste0("(WGS-", name, ") edgeR"))
    if (is.null(combined_plots)) {
      combined_plots <- p
    } else {
      combined_plots <- combined_plots + p
    }
  }
  png(paste0("Plots/WGS/ShuffledHistograms(", name, "-edgeR).png"), width=7500, height=3500, res = 300)
  print(combined_plots + plot_layout(ncol = 15))
  dev.off()

}

