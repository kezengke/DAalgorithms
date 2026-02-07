rm(list = ls())
library(phyloseq)
library(coin)
library(DESeq2)
library(edgeR)
library(ALDEx2)
library(ANCOMBC)
library(metagenomeSeq)
library(ggplot2)

normFun <- function(table) {
  n<-colSums(table)
  sumx<-sum(table)
  for (j in 1:ncol(table)) {
    table[,j]<-table[,j]/n[j]
  }
  table<-log10(table*(sumx/ncol(table))+1)
  table<-data.frame(table, check.names = F)
  return(table)
}

calcTtest <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")

  t_stats<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$stat})
  t_test_p<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$p.value})

  t_results<-cbind(t_stats, t_test_p)
  rownames(t_results)<-rownames(table)
  colnames(t_results)<-c("stats", "pval")
  t_results<-data.frame(t_results, check.names = F)
  return(t_results)
}

calcWilcox <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")

  wilcox_stats<-apply(table, 1, function(x){statistic(wilcox_test(unlist(x)~factor(meta$conditions)))})
  wilcox_p<-apply(table, 1, function(x){pvalue(wilcox_test(unlist(x)~factor(meta$conditions)))})

  wilcox_results<-cbind(wilcox_stats, wilcox_p)
  rownames(wilcox_results)<-rownames(table)
  colnames(wilcox_results)<-c("stats", "pval")
  wilcox_results<-data.frame(wilcox_results, check.names = F)
  return(wilcox_results)
}

#' Function to calculate DESeq2 results
calcDESeq2 <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")

  # Solve deseq2 all 0 issue
  table<-table+1

  meta$conditions<-factor(meta$conditions)
  dds1 <- DESeqDataSetFromMatrix(countData=table,
                                 colData=meta,
                                 design=~conditions)

  dds2 <- tryCatch({
    DESeq(dds1)
  }, error = function(e) {
    # If error occurs, return NA instead of running DESeq
    return(NULL)
  })

  # If dds2 is NULL, skip the remaining lines and create deseq_results as NA
  if (is.null(dds2)) {
    deseq_results <- matrix(NA, nrow = nrow(table), ncol = 2)
  } else {
    res <- results(dds2, cooksCutoff=FALSE, independentFiltering=FALSE)
    deseq_results<-cbind(res$stat, res$pvalue)
  }

  rownames(deseq_results)<-rownames(table)
  colnames(deseq_results)<-c("stats", "pval")
  deseq_results<-data.frame(deseq_results, check.names = F)
  return(deseq_results)
}

#' Function to calculate edgeR results
calcEdgeR <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
  group <- meta$condition
  dgList <- DGEList(counts=table, group = group)
  dgList <- edgeR::calcNormFactors(dgList, method="TMM")

  dgList <- tryCatch({
    estimateDisp(dgList)
  }, error = function(e) {
    # If error occurs, return NA instead of running estimateDisp
    return(NULL)
  })

  # If dgList is NULL, skip the remaining lines and create edger_results as NA
  if (is.null(dgList)) {
    edger_results <- matrix(NA, nrow = nrow(table), ncol = 2)
  } else {
    et <- exactTest(dgList)
    res <- et$table
    edger_results <- cbind(res$logFC, res$PValue)
  }

  rownames(edger_results) <- rownames(table)
  colnames(edger_results) <- c("stats", "pval")
  edger_results <- data.frame(edger_results, check.names = F)

  return(edger_results)
}

#' Function to calculate ALDEx2 t-test results
calcALDEx2Ttest <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
  CLR<-aldex.clr(table, meta$conditions, mc.samples = 128, denom = "all", verbose = F, gamma = 0.5)
  results<-aldex.ttest(CLR, hist.plot = F, paired.test = F, verbose = F)
  effectSizes<-aldex.effect(CLR)

  aldex2T_results<-cbind(effectSizes$diff.btw, results$we.ep)

  rownames(aldex2T_results)<-rownames(table)
  colnames(aldex2T_results)<-c("stats", "pval")
  aldex2T_results<-data.frame(aldex2T_results, check.names = F)
  return(aldex2T_results)
}

#' Function to calculate ALDEx2 wilcoxon results
calcALDEx2Wilcoxon <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
  CLR<-aldex.clr(table, meta$conditions, mc.samples = 128, denom = "all", verbose = F, gamma = 0.5)
  results<-aldex.ttest(CLR, hist.plot = F, paired.test = F, verbose = F)
  effectSizes<-aldex.effect(CLR)

  aldex2W_results<-cbind(effectSizes$diff.btw, results$wi.ep)

  rownames(aldex2W_results)<-rownames(table)
  colnames(aldex2W_results)<-c("stats", "pval")
  aldex2W_results<-data.frame(aldex2W_results, check.names = F)
  return(aldex2W_results)
}

#' Function to calculate ANCOMBC2 results
calcANCOMBC2 <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
  meta$conditions<-factor(meta$conditions, levels = unique(meta$conditions))

  ps <- phyloseq(
    otu_table(table, taxa_are_rows = TRUE),
    sample_data(meta)
  )

  out <- ancombc2(
    data         = ps,
    fix_formula  = names(meta),
    p_adj_method = "holm",
    prv_cut      = 0,
    lib_cut      = 0,
    s0_perc      = 0.05,
    alpha        = 0.05
  )

  res <- out$res
  ancombc2_results <- data.frame(res[, 7, drop = F], res[, 9, drop = F])
  colnames(ancombc2_results)<-c("stats", "pval")

  ancombc2_results$stats[is.na(ancombc2_results$stats)] <- 1
  rownames(ancombc2_results)<-rownames(table)

  return(ancombc2_results)
}

#' Function to calculate metagenomeSeq results
calcMetagenomeSeq <- function(table, meta) {
  group<-unique(meta$conditions)
  meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
  meta<-AnnotatedDataFrame(meta)
  obj<-newMRexperiment(table, phenoData = meta)
  percentile<-cumNormStatFast(obj)
  obj<-cumNorm(obj, p = percentile)
  pd<-pData(obj)
  mod<-model.matrix(~1 + conditions, data = pd)
  res<-fitFeatureModel(obj, mod)

  mtgnms_results<-cbind(res@fitZeroLogNormal$logFC, res@pvalues)
  rownames(mtgnms_results)<-rownames(table)
  colnames(mtgnms_results)<-c("stats", "pval")
  mtgnms_results<-data.frame(mtgnms_results, check.names = F)
  mtgnms_results$pval[is.na(mtgnms_results$pval)] <- 1
  return(mtgnms_results)
}

# 'Function to calculate fraction of significant results
calcFraction <- function(pvalT) {
  fractions<-apply(pvalT, 2,
                   function(column){
                     sum(column < 0.05)/length(column)
                   })
  return(fractions)
}

# 'Function to make the box plot for all the fraction of significant results
makePlot <- function(fractT, classifier, name, shuffleType) {

  df_long <- stack(fractT)

  p<-ggplot(df_long, aes(x = ind, y = values, fill = ind)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("coral1", "lightslateblue", "olivedrab3", "goldenrod1", "skyblue", "pink", "orchid4", "darkorange3")) +
    theme_classic() +
    labs(title = paste0("(", classifier, "-", name, ")\n", shuffleType),
         x = "DAA methods",
         y = "Fraction of significant results") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  return(p)
}

png(paste0("Plots/RDP/TestingScript_rhizo(R).png"), width=1200, height=1200, res = 300)
# Run script for rhizo RDP
rawT<-read.table("CountsTables/RDPRaw//rhizo.txt", sep = "\t", header = T, row.names = 1)
metaData<-read.table("MetaData/metadata_rhizo.txt", sep = "\t", header = T, row.names = 1)

rawT<-rawT[, intersect(colnames(rawT), rownames(metaData)), drop = F]
metaData<-metaData[intersect(colnames(rawT), rownames(metaData)), , drop = F]

normT<-normFun(rawT)

rownames(metaData)<-colnames(rawT)
colnames(metaData)<-"conditions"

ttestpvals<-c()
wilcoxpvals<-c()
deseq2pvals<-c()
edgerpvals<-c()
aldexttestpvals<-c()
aldexwilcoxonpvals<-c()
ancombc2pvals<-c()
metagenomeseqpvals<-c()

for (i in 1:10) {
  shuffleMeta<-metaData
  shuffleMeta$conditions<-sample(shuffleMeta$conditions)

  ttestresults<-calcTtest(normT, shuffleMeta)
  ttestpvals<-cbind(ttestpvals, ttestresults$pval)

  wilcoxresults<-calcWilcox(normT, shuffleMeta)
  wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)

  deseq2results<-calcDESeq2(rawT, shuffleMeta)
  deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)

  edgerresults<-calcEdgeR(rawT, shuffleMeta)
  edgerpvals<-cbind(edgerpvals, edgerresults$pval)

  aldexttestresults<-calcALDEx2Ttest(rawT, shuffleMeta)
  aldexttestpvals<-cbind(aldexttestpvals, aldexttestresults$pval)

  aldexwilcoxonresults<-calcALDEx2Wilcoxon(rawT, shuffleMeta)
  aldexwilcoxonpvals<-cbind(aldexwilcoxonpvals, aldexwilcoxonresults$pval)

  ancombc2results<-calcANCOMBC2(rawT, shuffleMeta)
  ancombc2pvals<-cbind(ancombc2pvals, ancombc2results$pval)

  metagenomeseqresults<-calcMetagenomeSeq(rawT, shuffleMeta)
  metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)
}

rownames(ttestpvals)<-rownames(rawT)
rownames(wilcoxpvals)<-rownames(rawT)
rownames(deseq2pvals)<-rownames(rawT)
rownames(edgerpvals)<-rownames(rawT)
rownames(aldexttestpvals)<-rownames(rawT)
rownames(aldexwilcoxonpvals)<-rownames(rawT)
rownames(ancombc2pvals)<-rownames(rawT)
rownames(metagenomeseqpvals)<-rownames(rawT)

# significant fraction
sigT<-calcFraction(ttestpvals)
sigW<-calcFraction(wilcoxpvals)
sigD<-calcFraction(deseq2pvals)
sigE<-calcFraction(edgerpvals)
sigAT<-calcFraction(aldexttestpvals)
sigAW<-calcFraction(aldexwilcoxonpvals)
sigB<-calcFraction(ancombc2pvals)
sigM<-calcFraction(metagenomeseqpvals)

fractionT<-cbind(sigT, sigW, sigD, sigE, sigAT, sigAW, sigB, sigM)
colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR", "ALDEx2t-test", "ALDEx2Wilcoxon", "ANCOMBC2", "metagenomeSeq")
fractionT<-data.frame(fractionT, check.names = F)

p<-makePlot(fractionT, "RDP", "Rhizo", "Shuffle sample tags")

print(p)
dev.off()

# Run script for tomatoPollen RNAseq
png(paste0("Plots/RNAseq/TestingScript_tomatoPollen(R).png"), width=1200, height=1200, res = 300)
rawT<-read.table("CountsTables/RNAseqRaw/tomatoPollen.txt", sep = "\t", header = T, row.names = 1)
metaData<-read.table("MetaData/metadata_tomatoPollen.txt", sep = "\t", header = T, row.names = 1)

rawT<-rawT[, intersect(colnames(rawT), rownames(metaData)), drop = F]
metaData<-metaData[intersect(colnames(rawT), rownames(metaData)), , drop = F]

normT<-normFun(rawT)

rownames(metaData)<-colnames(rawT)
colnames(metaData)<-"conditions"

ttestpvals<-c()
wilcoxpvals<-c()
deseq2pvals<-c()
edgerpvals<-c()
aldexttestpvals<-c()
aldexwilcoxonpvals<-c()
ancombc2pvals<-c()
metagenomeseqpvals<-c()

for (i in 1:10) {
  shuffleMeta<-metaData
  shuffleMeta$conditions<-sample(shuffleMeta$conditions)

  ttestresults<-calcTtest(normT, shuffleMeta)
  ttestpvals<-cbind(ttestpvals, ttestresults$pval)

  wilcoxresults<-calcWilcox(normT, shuffleMeta)
  wilcoxpvals<-cbind(wilcoxpvals, wilcoxresults$pval)

  deseq2results<-calcDESeq2(rawT, shuffleMeta)
  deseq2pvals<-cbind(deseq2pvals, deseq2results$pval)

  edgerresults<-calcEdgeR(rawT, shuffleMeta)
  edgerpvals<-cbind(edgerpvals, edgerresults$pval)

  aldexttestresults<-calcALDEx2Ttest(rawT, shuffleMeta)
  aldexttestpvals<-cbind(aldexttestpvals, aldexttestresults$pval)

  aldexwilcoxonresults<-calcALDEx2Wilcoxon(rawT, shuffleMeta)
  aldexwilcoxonpvals<-cbind(aldexwilcoxonpvals, aldexwilcoxonresults$pval)

  ancombc2results<-calcANCOMBC2(rawT, shuffleMeta)
  ancombc2pvals<-cbind(ancombc2pvals, ancombc2results$pval)

  metagenomeseqresults<-calcMetagenomeSeq(rawT, shuffleMeta)
  metagenomeseqpvals<-cbind(metagenomeseqpvals, metagenomeseqresults$pval)
}

rownames(ttestpvals)<-rownames(rawT)
rownames(wilcoxpvals)<-rownames(rawT)
rownames(deseq2pvals)<-rownames(rawT)
rownames(edgerpvals)<-rownames(rawT)
rownames(aldexttestpvals)<-rownames(rawT)
rownames(aldexwilcoxonpvals)<-rownames(rawT)
rownames(ancombc2pvals)<-rownames(rawT)
rownames(metagenomeseqpvals)<-rownames(rawT)

# significant fraction
sigT<-calcFraction(ttestpvals)
sigW<-calcFraction(wilcoxpvals)
sigD<-calcFraction(deseq2pvals)
sigE<-calcFraction(edgerpvals)
sigAT<-calcFraction(aldexttestpvals)
sigAW<-calcFraction(aldexwilcoxonpvals)
sigB<-calcFraction(ancombc2pvals)
sigM<-calcFraction(metagenomeseqpvals)

fractionT<-cbind(sigT, sigW, sigD, sigE, sigAT, sigAW, sigB, sigM)
colnames(fractionT)<-c("t-test", "Wilcoxon", "DESeq2", "edgeR", "ALDEx2t-test", "ALDEx2Wilcoxon", "ANCOMBC2", "metagenomeSeq")
fractionT<-data.frame(fractionT, check.names = F)

p<-makePlot(fractionT, "RNAseq", "tomatoPollen", "Shuffle sample tags")

print(p)
dev.off()


