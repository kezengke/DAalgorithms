#!/usr/bin/env python3
"""
Python script for the R testing script:
- Python: t-test, Wilcoxon, normalization, data loading, plotting
- R functions: DESeq2, edgeR, ALDEx2 t-test, ALDEx2 Wilcoxon, ANCOMBC2, metagenomeSeq
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu
import subprocess
import tempfile
import os
import warnings
warnings.filterwarnings('ignore')

def norm_fun(table):
    """Normalization function for t-test and Wilcoxon"""
    data = table.values.copy().astype(np.float64)
    n = np.sum(data, axis=0)
    sumx = np.sum(data)
    
    for j in range(data.shape[1]):
        if n[j] > 0:
            data[:, j] = data[:, j] / n[j]
    
    data = np.log10(data * (sumx / data.shape[1]) + 1)
    result = pd.DataFrame(data, index=table.index, columns=table.columns)
    return result

def calc_ttest(table, meta):
    """Python t-test"""
    conditions = meta['conditions'].unique()
    group_labels = np.where(meta['conditions'] == conditions[0], 'group1', 'group2')
    
    t_stats, t_pvals = [], []
    
    for taxon in table.index:
        values = table.loc[taxon].values
        group1_values = values[group_labels == 'group1']
        group2_values = values[group_labels == 'group2']
        
        group1_values = group1_values[~np.isnan(group1_values)]
        group2_values = group2_values[~np.isnan(group2_values)]
        
        if len(group1_values) > 1 and len(group2_values) > 1:
            try:
                t_stat, p_val = ttest_ind(group1_values, group2_values)
                t_stats.append(t_stat)
                t_pvals.append(p_val)
            except:
                t_stats.append(np.nan)
                t_pvals.append(1.0)
        else:
            t_stats.append(np.nan)
            t_pvals.append(1.0)
    
    results = pd.DataFrame({'stats': t_stats, 'pval': t_pvals}, index=table.index)
    return results

def calc_wilcox(table, meta):
    """Python Wilcoxon"""
    conditions = meta['conditions'].unique()
    group_labels = np.where(meta['conditions'] == conditions[0], 'group1', 'group2')
    
    wilcox_stats, wilcox_pvals = [], []
    
    for taxon in table.index:
        values = table.loc[taxon].values
        group1_values = values[group_labels == 'group1']
        group2_values = values[group_labels == 'group2']
        
        group1_values = group1_values[~np.isnan(group1_values)]
        group2_values = group2_values[~np.isnan(group2_values)]
        
        if len(group1_values) > 1 and len(group2_values) > 1:
            try:
                u_stat, p_val = mannwhitneyu(group1_values, group2_values, alternative='two-sided')
                wilcox_stats.append(u_stat)
                wilcox_pvals.append(p_val)
            except:
                wilcox_stats.append(np.nan)
                wilcox_pvals.append(1.0)
        else:
            wilcox_stats.append(np.nan)
            wilcox_pvals.append(1.0)
    
    results = pd.DataFrame({'stats': wilcox_stats, 'pval': wilcox_pvals}, index=table.index)
    return results

def call_r_function(function_name, table, meta, temp_dir, iteration, total_iterations):
    """Call specific R function via subprocess"""
    
    # Write data to temporary files
    counts_file = os.path.join(temp_dir, f"counts_{function_name}.txt")
    meta_file = os.path.join(temp_dir, f"meta_{function_name}.txt")
    results_file = os.path.join(temp_dir, f"results_{function_name}.txt")
    
    table.to_csv(counts_file, sep='\t')
    meta.to_csv(meta_file, sep='\t')
    
    # R functions
    if function_name == "DESeq2":
        r_script = f"""
        library(DESeq2)
        
        calcDESeq2 <- function(table, meta) {{
          group<-unique(meta$conditions)
          meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
          table<-table+1
          meta$conditions<-factor(meta$conditions)
          dds1 <- DESeqDataSetFromMatrix(countData=table, colData=meta, design=~conditions)
          dds2 <- tryCatch({{
            DESeq(dds1)
          }}, error = function(e) {{
            return(NULL)
          }})
          if (is.null(dds2)) {{
            deseq_results <- matrix(NA, nrow = nrow(table), ncol = 2)
          }} else {{
            res <- results(dds2, cooksCutoff=FALSE, independentFiltering=FALSE)
            deseq_results<-cbind(res$stat, res$pvalue)
          }}
          rownames(deseq_results)<-rownames(table)
          colnames(deseq_results)<-c("stats", "pval")
          deseq_results<-data.frame(deseq_results, check.names = F)
          return(deseq_results)
        }}
        
        rawT <- read.table("{counts_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        metaData <- read.table("{meta_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        rawT <- rawT[, intersect(colnames(rawT), rownames(metaData)), drop = FALSE]
        metaData <- metaData[intersect(colnames(rawT), rownames(metaData)), , drop = FALSE]
        results <- calcDESeq2(rawT, metaData)
        write.table(results, "{results_file}", sep = "\\t", quote = FALSE)
        """
        
    elif function_name == "edgeR":
        r_script = f"""
        library(edgeR)
        
        calcEdgeR <- function(table, meta) {{
          group<-unique(meta$conditions)
          meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
          group <- meta$conditions
          dgList <- DGEList(counts=table, group = group)
          dgList <- edgeR::calcNormFactors(dgList, method="TMM")
          dgList <- tryCatch({{
            estimateDisp(dgList)
          }}, error = function(e) {{
            return(NULL)
          }})
          if (is.null(dgList)) {{
            edger_results <- matrix(NA, nrow = nrow(table), ncol = 2)
          }} else {{
            et <- exactTest(dgList)
            res <- et$table
            edger_results <- cbind(res$logFC, res$PValue)
          }}
          rownames(edger_results) <- rownames(table)
          colnames(edger_results) <- c("stats", "pval")
          edger_results <- data.frame(edger_results, check.names = F)
          return(edger_results)
        }}
        
        rawT <- read.table("{counts_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        metaData <- read.table("{meta_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        rawT <- rawT[, intersect(colnames(rawT), rownames(metaData)), drop = FALSE]
        metaData <- metaData[intersect(colnames(rawT), rownames(metaData)), , drop = FALSE]
        results <- calcEdgeR(rawT, metaData)
        write.table(results, "{results_file}", sep = "\\t", quote = FALSE)
        """
        
    elif function_name == "ALDEx2Ttest":
        r_script = f"""
        library(ALDEx2)
        
        calcALDEx2Ttest <- function(table, meta) {{
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
        }}
        
        rawT <- read.table("{counts_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        metaData <- read.table("{meta_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        rawT <- rawT[, intersect(colnames(rawT), rownames(metaData)), drop = FALSE]
        metaData <- metaData[intersect(colnames(rawT), rownames(metaData)), , drop = FALSE]
        results <- calcALDEx2Ttest(rawT, metaData)
        write.table(results, "{results_file}", sep = "\\t", quote = FALSE)
        """
        
    elif function_name == "ALDEx2Wilcoxon":
        r_script = f"""
        library(ALDEx2)
        
        calcALDEx2Wilcoxon <- function(table, meta) {{
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
        }}
        
        rawT <- read.table("{counts_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        metaData <- read.table("{meta_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        rawT <- rawT[, intersect(colnames(rawT), rownames(metaData)), drop = FALSE]
        metaData <- metaData[intersect(colnames(rawT), rownames(metaData)), , drop = FALSE]
        results <- calcALDEx2Wilcoxon(rawT, metaData)
        write.table(results, "{results_file}", sep = "\\t", quote = FALSE)
        """
        
    elif function_name == "ANCOMBC2":
        r_script = f"""
        library(phyloseq)
        library(ANCOMBC)
        
        calcANCOMBC2 <- function(table, meta) {{
          group<-unique(meta$conditions)
          meta$conditions<-ifelse(meta$conditions == group[1], "group1", "group2")
          meta$conditions<-factor(meta$conditions, levels = unique(meta$conditions))
          ps <- phyloseq(otu_table(table, taxa_are_rows = TRUE), sample_data(meta))
          out <- ancombc2(data = ps, fix_formula = names(meta), p_adj_method = "holm", 
                         prv_cut = 0, lib_cut = 0, s0_perc = 0.05, alpha = 0.05)
          res <- out$res
          ancombc2_results <- data.frame(res[, 7, drop = F], res[, 9, drop = F])
          colnames(ancombc2_results)<-c("stats", "pval")
          ancombc2_results$stats[is.na(ancombc2_results$stats)] <- 1
          rownames(ancombc2_results)<-rownames(table)
          return(ancombc2_results)
        }}
        
        rawT <- read.table("{counts_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        metaData <- read.table("{meta_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        rawT <- rawT[, intersect(colnames(rawT), rownames(metaData)), drop = FALSE]
        metaData <- metaData[intersect(colnames(rawT), rownames(metaData)), , drop = FALSE]
        results <- calcANCOMBC2(rawT, metaData)
        write.table(results, "{results_file}", sep = "\\t", quote = FALSE)
        """
        
    elif function_name == "metagenomeSeq":
        r_script = f"""
        library(metagenomeSeq)
        
        calcMetagenomeSeq <- function(table, meta) {{
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
        }}
        
        rawT <- read.table("{counts_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        metaData <- read.table("{meta_file}", sep = "\\t", header = TRUE, row.names = 1, check.names = FALSE)
        rawT <- rawT[, intersect(colnames(rawT), rownames(metaData)), drop = FALSE]
        metaData <- metaData[intersect(colnames(rawT), rownames(metaData)), , drop = FALSE]
        results <- calcMetagenomeSeq(rawT, metaData)
        write.table(results, "{results_file}", sep = "\\t", quote = FALSE)
        """
    
    # Write R script to file
    r_script_file = os.path.join(temp_dir, f"script_{function_name}.R")
    with open(r_script_file, 'w') as f:
        f.write(r_script)
    
    # Run R script
    try:
        result = subprocess.run(['Rscript', r_script_file], 
                              capture_output=True, text=True, timeout=120)
        
        if result.returncode == 0:
            # Read results
            results_df = pd.read_csv(results_file, sep='\t', index_col=0)
            print(f"{function_name} completed ({iteration}/{total_iterations})")
            return results_df
        else:
            print(f"{function_name} failed ({iteration}/{total_iterations}): {result.stderr}")
            return None
            
    except subprocess.TimeoutExpired:
        print(f"{function_name} timed out ({iteration}/{total_iterations})")
        return None
    except Exception as e:
        print(f"{function_name} error ({iteration}/{total_iterations}): {e}")
        return None

def calc_fraction(pval_df):
    """Calculate fraction of significant results"""
    fractions = []
    for col in pval_df.columns:
        significant_count = (pval_df[col] < 0.05).sum()
        total_count = len(pval_df)
        fraction = significant_count / total_count
        fractions.append(fraction)
    
    return pd.Series(fractions, index=pval_df.columns)

def make_plot(fraction_df, classifier, name, shuffle_type):
    """Create boxplot"""
    plot_data = []
    for method in fraction_df.columns:
        for fraction in fraction_df[method]:
            plot_data.append({'Method': method, 'Fraction_Significant': fraction})
    
    df_plot = pd.DataFrame(plot_data)
    
    plt.figure(figsize=(12, 8))
    
    # Colors matching R script
    colors = ['coral', 'slateblue', 'olivedrab', 'gold', 'skyblue', 'pink', 'orchid', 'darkorange']
    
    sns.boxplot(data=df_plot, x='Method', y='Fraction_Significant', 
                hue='Method', palette=colors[:len(fraction_df.columns)], legend=False)
    
    plt.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, linewidth=2)
    
    plt.title(f'({classifier}-{name})\n{shuffle_type}', fontsize=14, fontweight='bold')
    plt.xlabel('DAA methods', fontsize=12)
    plt.ylabel('Fraction of significant results', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    # Fix y-axis limits to show all methods
    y_max = max(0.25, df_plot['Fraction_Significant'].max() * 1.2)
    y_min = -0.01  # Go slightly below 0 to show ALDEx2 boxes fully
    plt.ylim(y_min, y_max)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return plt.gcf()

def main():
    
    # ==================== Dataset 1: RDP - Rhizo ====================
    print("="*60)
    print("Dataset 1: RDP - Rhizo")
    print("="*60)
    
    # Load data
    rawT = pd.read_csv("CountsTables/RDPRaw/rhizo.txt", sep='\t', index_col=0)
    metaData = pd.read_csv("MetaData/metadata_rhizo.txt", sep='\t', index_col=0)
    
    rawT.columns = rawT.columns.str.strip('"')
    
    common_samples = rawT.columns.intersection(metaData.index)
    rawT = rawT[common_samples]
    metaData = metaData.loc[common_samples]
    
    
    # Normalize data
    normT = norm_fun(rawT)
    
    # Prepare metadata
    metaData.index = rawT.columns
    metaData.columns = ['conditions']
    
    # Create temporary directory for R files
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Using temporary directory: {temp_dir}")
        
        # Test with 10 iterations
        n_iterations = 10
        print(f"Running {n_iterations} iterations...")
        
        # Initialize result DataFrames
        results_storage = {
            't-test': pd.DataFrame(index=rawT.index),
            'Wilcoxon': pd.DataFrame(index=rawT.index),
            'DESeq2': pd.DataFrame(index=rawT.index),
            'edgeR': pd.DataFrame(index=rawT.index),
            'ALDEx2t-test': pd.DataFrame(index=rawT.index),
            'ALDEx2Wilcoxon': pd.DataFrame(index=rawT.index),
            'ANCOMBC2': pd.DataFrame(index=rawT.index),
            'metagenomeSeq': pd.DataFrame(index=rawT.index)
        }
        
        # Run iterations
        for i in range(n_iterations):
            
            # Shuffle metadata
            shuffleMeta = metaData.copy()
            shuffleMeta['conditions'] = np.random.permutation(shuffleMeta['conditions'])
            
            # Python methods
            ttestresults = calc_ttest(normT, shuffleMeta)
            results_storage['t-test'][f'iter_{i+1}'] = ttestresults['pval']
            
            wilcoxresults = calc_wilcox(normT, shuffleMeta)
            results_storage['Wilcoxon'][f'iter_{i+1}'] = wilcoxresults['pval']
            
            # R methods
            for r_method in ['DESeq2', 'edgeR', 'ALDEx2Ttest', 'ALDEx2Wilcoxon', 'ANCOMBC2', 'metagenomeSeq']:
                r_results = call_r_function(r_method, rawT, shuffleMeta, temp_dir, i+1, n_iterations)
                if r_results is not None:
                    # Map R method names to storage keys
                    storage_key = r_method
                    if r_method == 'ALDEx2Ttest':
                        storage_key = 'ALDEx2t-test'
                    elif r_method == 'ALDEx2Wilcoxon':
                        storage_key = 'ALDEx2Wilcoxon'
                    
                    results_storage[storage_key][f'iter_{i+1}'] = r_results['pval']
                else:
                    # Fallback to t-test if R method fails
                    results_storage[storage_key][f'iter_{i+1}'] = ttestresults['pval']
        
        
        # Calculate fractions
        fractions = {}
        for method, pval_df in results_storage.items():
            fractions[method] = calc_fraction(pval_df)
        
        fractionT = pd.DataFrame(fractions)
        
        # Create plot
        p = make_plot(fractionT, "RDP", "Rhizo", "Shuffle sample tags")
        
        # Save plot
        plt.savefig("Plots/RDP/TestingScript_rhizo(Python).png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # ==================== Dataset 2: RNAseq - tomatoPollen ====================
    print("="*60)
    print("Dataset 2: RNAseq - tomatoPollen")
    print("="*60)
    
    # Load data
    rawT = pd.read_csv("CountsTables/RNAseqRaw/tomatoPollen.txt", sep='\t', index_col=0)
    metaData = pd.read_csv("MetaData/metadata_tomatoPollen.txt", sep='\t', index_col=0)
    
    rawT.columns = rawT.columns.str.strip('"')
    
    common_samples = rawT.columns.intersection(metaData.index)
    rawT = rawT[common_samples]
    metaData = metaData.loc[common_samples]
    
    # Normalize data
    normT = norm_fun(rawT)
    
    # Prepare metadata
    metaData.index = rawT.columns
    metaData.columns = ['conditions']
    
    # Create temporary directory for R files
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Using temporary directory: {temp_dir}")
        
        # Test with 10 iterations
        n_iterations = 10
        print(f"Running {n_iterations} iterations...")
        
        # Initialize result DataFrames
        results_storage = {
            't-test': pd.DataFrame(index=rawT.index),
            'Wilcoxon': pd.DataFrame(index=rawT.index),
            'DESeq2': pd.DataFrame(index=rawT.index),
            'edgeR': pd.DataFrame(index=rawT.index),
            'ALDEx2t-test': pd.DataFrame(index=rawT.index),
            'ALDEx2Wilcoxon': pd.DataFrame(index=rawT.index),
            'ANCOMBC2': pd.DataFrame(index=rawT.index),
            'metagenomeSeq': pd.DataFrame(index=rawT.index)
        }
        
        # Run iterations
        for i in range(n_iterations):
            
            # Shuffle metadata
            shuffleMeta = metaData.copy()
            shuffleMeta['conditions'] = np.random.permutation(shuffleMeta['conditions'])
            
            # Python methods
            ttestresults = calc_ttest(normT, shuffleMeta)
            results_storage['t-test'][f'iter_{i+1}'] = ttestresults['pval']
            
            wilcoxresults = calc_wilcox(normT, shuffleMeta)
            results_storage['Wilcoxon'][f'iter_{i+1}'] = wilcoxresults['pval']
            
            # R methods
            for r_method in ['DESeq2', 'edgeR', 'ALDEx2Ttest', 'ALDEx2Wilcoxon', 'ANCOMBC2', 'metagenomeSeq']:
                r_results = call_r_function(r_method, rawT, shuffleMeta, temp_dir, i+1, n_iterations)
                if r_results is not None:
                    # Map R method names to storage keys
                    storage_key = r_method
                    if r_method == 'ALDEx2Ttest':
                        storage_key = 'ALDEx2t-test'
                    elif r_method == 'ALDEx2Wilcoxon':
                        storage_key = 'ALDEx2Wilcoxon'
                    
                    results_storage[storage_key][f'iter_{i+1}'] = r_results['pval']
                else:
                    # Fallback to t-test if R method fails
                    results_storage[storage_key][f'iter_{i+1}'] = ttestresults['pval']
        
        
        # Calculate fractions
        fractions = {}
        for method, pval_df in results_storage.items():
            fractions[method] = calc_fraction(pval_df)
        
        fractionT = pd.DataFrame(fractions)
        
        # Create plot
        p = make_plot(fractionT, "RNAseq", "tomatoPollen", "Shuffle sample tags")
        
        # Save plot
        plt.savefig("Plots/RNAseq/TestingScript_tomatoPollen(Python).png", dpi=300, bbox_inches='tight')
        plt.close()
        
if __name__ == "__main__":
    main()
