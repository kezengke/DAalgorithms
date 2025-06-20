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
    scale_fill_manual(values = c("#2aa7de", "#25377f")) +
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

name <- "bsurgery"

ttestpvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "_t.txt"), header = T, row.names = 1, sep = "\t")
wilcoxpvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "_wilcox.txt"), header = T, row.names = 1, sep = "\t")

# significant fraction
sigT<-calcFraction(ttestpvals)
sigW<-calcFraction(wilcoxpvals)

fractionT<-cbind(sigT, sigW)
colnames(fractionT)<-c("t-test", "Wilcoxon")
fractionT<-data.frame(fractionT, check.names = F)

p<-makePlot(fractionT, "RDP", name, "Shuffle counts in sample")

ttestpvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "Python_t.txt"), header = T, row.names = 1, sep = "\t")
wilcoxpvals<-read.table(paste0("ResampleShuffleInSampleDump/RDP/", name, "Python_wilcox.txt"), header = T, row.names = 1, sep = "\t")

# significant fraction
sigT<-calcFraction(ttestpvals)
sigW<-calcFraction(wilcoxpvals)

fractionT<-cbind(sigT, sigW)
colnames(fractionT)<-c("t-test", "Wilcoxon")
fractionT<-data.frame(fractionT, check.names = F)

p<-makePlot(fractionT, "RDP", name, "(Python)Shuffle counts in sample")
p
