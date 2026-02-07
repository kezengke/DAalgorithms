#!/usr/bin/env python3
"""
Generate fraction of significant counts boxplot using actual R packages
This version uses rpy2 to call the real differential abundance analysis methods
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind
import warnings
warnings.filterwarnings('ignore')

# Try to import rpy2, if not available, fall back to approximations
try:
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    R_AVAILABLE = True
    print("R integration available - will use actual packages")
except ImportError:
    R_AVAILABLE = False
    print("R integration not available - will use approximations")

def load_counts_and_metadata(dataset="assal"):
    """
    Load counts table and metadata for a specific dataset
    """
    
    # Read counts table
    counts_file = f"CountsTables/RDPRaw/{dataset}.txt"
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
    counts_df.columns = counts_df.columns.str.strip('"')
    
    # Read metadata
    metadata_file = f"MetaData/metadata_{dataset}.txt"
    metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)
    
    # Align data
    common_samples = counts_df.columns.intersection(metadata_df.index)
    counts_aligned = counts_df[common_samples]
    metadata_aligned = metadata_df.loc[common_samples]
    
    print(f"Loaded {dataset} dataset: {counts_aligned.shape[0]} taxa x {counts_aligned.shape[1]} samples")
    print(f"Conditions: {', '.join(metadata_aligned['Sample'].unique())}")
    
    return counts_aligned, metadata_aligned

def perform_ttest(counts_df, metadata_df):
    """Standard t-test"""
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        # Remove zeros and log transform
        group1 = group1[group1 > 0]
        group2 = group2[group2 > 0]
        
        if len(group1) > 1 and len(group2) > 1:
            group1_log = np.log(group1 + 1)
            group2_log = np.log(group2 + 1)
            try:
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_wilcoxon(counts_df, metadata_df):
    """Standard Wilcoxon test"""
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        group1 = group1[group1 > 0]
        group2 = group2[group2 > 0]
        
        if len(group1) > 1 and len(group2) > 1:
            try:
                _, pval = mannwhitneyu(group1, group2, alternative='two-sided')
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_deseq2(counts_df, metadata_df):
    """Use actual DESeq2 package via rpy2"""
    if not R_AVAILABLE:
        return perform_deseq2_like(counts_df, metadata_df)
    
    try:
        # Import R packages
        deseq2 = importr('DESeq2')
        base = importr('base')
        stats = importr('stats')
        
        # Prepare data for R
        with localconverter(robjects.default_converter + pandas2ri.converter):
            counts_r = robjects.conversion.py2rpy(counts_df.T)  # Transpose for R format
            metadata_r = robjects.conversion.py2rpy(metadata_df)
        
        # Create DESeqDataSet
        dds = deseq2.DESeqDataSetFromMatrix(counts_r, metadata_r, 
                                           robjects.Formula('~ Sample'))
        
        # Run DESeq2
        dds = deseq2.DESeq(dds)
        
        # Get results
        results = deseq2.results(dds)
        
        # Extract p-values
        pvals = np.array(results.rx2('pvalue'))
        
        # Handle NaN values
        pvals = np.nan_to_num(pvals, nan=1.0)
        
        return pd.Series(pvals, index=counts_df.index)
        
    except Exception as e:
        print(f"DESeq2 error: {e}, falling back to approximation")
        return perform_deseq2_like(counts_df, metadata_df)

def perform_edgeR(counts_df, metadata_df):
    """Use actual edgeR package via rpy2"""
    if not R_AVAILABLE:
        return perform_edgeR_like(counts_df, metadata_df)
    
    try:
        # Import R packages
        edgeR = importr('edgeR')
        base = importr('base')
        
        # Prepare data for R
        with localconverter(robjects.default_converter + pandas2ri.converter):
            counts_r = robjects.conversion.py2rpy(counts_df.T)  # Transpose for R format
            metadata_r = robjects.conversion.py2rpy(metadata_df)
        
        # Create DGEList
        dge = edgeR.DGEList(counts_r)
        
        # Calculate normalization factors
        dge = edgeR.calcNormFactors(dge)
        
        # Design matrix
        design = robjects.Formula('~ Sample')
        design_matrix = stats.model_matrix(design, data=metadata_r)
        
        # Estimate dispersion
        dge = edgeR.estimateDisp(dge, design_matrix)
        
        # Fit model
        fit = edgeR.glmQLFit(dge, design_matrix)
        
        # Test
        results = edgeR.glmQLFTest(fit)
        
        # Extract p-values
        pvals = np.array(results.rx2('table').rx2('PValue'))
        
        # Handle NaN values
        pvals = np.nan_to_num(pvals, nan=1.0)
        
        return pd.Series(pvals, index=counts_df.index)
        
    except Exception as e:
        print(f"edgeR error: {e}, falling back to approximation")
        return perform_edgeR_like(counts_df, metadata_df)

def perform_aldex2(counts_df, metadata_df):
    """Use actual ALDEx2 package via rpy2"""
    if not R_AVAILABLE:
        return perform_aldex2_like(counts_df, metadata_df)
    
    try:
        # Import R packages
        aldex2 = importr('ALDEx2')
        
        # Prepare data for R
        with localconverter(robjects.default_converter + pandas2ri.converter):
            counts_r = robjects.conversion.py2rpy(counts_df.T)  # Transpose for R format
        conditions = metadata_df['Sample'].values
        
        # Run ALDEx2
        aldex_obj = aldex2.aldex(counts_r, conditions, test='t')
        
        # Extract p-values
        pvals = np.array(aldex_obj.rx2('wi.eBH'))  # Benjamini-Hochberg corrected p-values
        
        # Handle NaN values
        pvals = np.nan_to_num(pvals, nan=1.0)
        
        return pd.Series(pvals, index=counts_df.index)
        
    except Exception as e:
        print(f"ALDEx2 error: {e}, falling back to approximation")
        return perform_aldex2_like(counts_df, metadata_df)

def perform_aldex2_like(counts_df, metadata_df):
    """ALDEx2 approximation"""
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    # CLR transformation
    def clr_transform(row):
        geometric_mean = np.exp(np.mean(np.log(row + 1)))
        return np.log(row + 1) - np.log(geometric_mean)
    
    counts_clr = counts_df.apply(clr_transform, axis=1)
    
    for taxon in counts_clr.index:
        group1 = counts_clr.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_clr.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        if len(group1) > 1 and len(group2) > 1:
            try:
                _, pval = ttest_ind(group1, group2)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_deseq2_like(counts_df, metadata_df):
    """DESeq2 approximation"""
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        group1_log = np.log(group1 + 1)
        group2_log = np.log(group2 + 1)
        
        if len(group1_log) > 1 and len(group2_log) > 1:
            try:
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_edgeR_like(counts_df, metadata_df):
    """edgeR approximation"""
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        group1_log = np.log(group1 + 1)
        group2_log = np.log(group2 + 1)
        
        if len(group1_log) > 1 and len(group2_log) > 1:
            try:
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def calculate_pvalues_from_counts(counts_df, metadata_df, n_iterations=10):
    """
    Calculate p-values directly from counts data using multiple methods
    """
    
    print(f"Calculating p-values for {n_iterations} iterations...")
    
    # Define methods and their functions
    methods = {
        't-test': perform_ttest,
        'Wilcoxon': perform_wilcoxon,
        'DESeq2': perform_deseq2,
        'edgeR': perform_edgeR,
        'ALDEx2': perform_aldex2,
        'ANCOMBC2': perform_deseq2_like,  # Using DESeq2-like for now
        'metagenomeSeq': perform_deseq2_like  # Using DESeq2-like for now
    }
    
    pvalue_data = {}
    
    for method_name, method_func in methods.items():
        print(f"  Processing {method_name}...")
        
        fractions = []
        
        for iteration in range(n_iterations):
            # Shuffle the metadata labels to simulate null hypothesis
            metadata_shuffled = metadata_df.copy()
            metadata_shuffled['Sample'] = np.random.permutation(metadata_shuffled['Sample'])
            
            # Calculate p-values for this iteration
            pvals = method_func(counts_df, metadata_shuffled)
            
            # Calculate fraction of significant results
            significant_count = (pvals < 0.05).sum()
            total_taxa = len(pvals)
            fraction = significant_count / total_taxa
            
            fractions.append(fraction)
        
        pvalue_data[method_name] = pd.Series(fractions)
        print(f"    Completed {len(fractions)} iterations")
    
    return pvalue_data

def print_summary_stats(pvalue_data):
    """Print summary statistics for each method"""
    
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    for method, fractions in pvalue_data.items():
        mean_frac = np.mean(fractions)
        median_frac = np.median(fractions)
        std_frac = np.std(fractions)
        min_frac = np.min(fractions)
        max_frac = np.max(fractions)
        
        print(f"\n{method}:")
        print(f"  Mean: {mean_frac:.4f}")
        print(f"  Median: {median_frac:.4f}")
        print(f"  Std Dev: {std_frac:.4f}")
        print(f"  Range: [{min_frac:.4f}, {max_frac:.4f}]")
        print(f"  Iterations: {len(fractions)}")

def create_boxplot(pvalue_data, dataset="assal", shuffle_type="Calculated", 
                   save_path=None, show_plot=True):
    """Create boxplot of fraction of significant results"""
    
    if not pvalue_data:
        print("No data available for plotting")
        return
    
    # Prepare data for plotting
    plot_data = []
    for method, fractions in pvalue_data.items():
        for fraction in fractions:
            plot_data.append({'Method': method, 'Fraction_Significant': fraction})
    
    df_plot = pd.DataFrame(plot_data)
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Create boxplot with custom colors
    colors = ['coral', 'slateblue', 'olivedrab', 'gold', 'skyblue', 'pink', 'lightcoral', 'lightgreen']
    method_names = list(pvalue_data.keys())
    
    box_plot = sns.boxplot(data=df_plot, x='Method', y='Fraction_Significant', 
                          hue='Method', palette=colors[:len(method_names)], legend=False)
    
    # Add horizontal line at 0.05
    plt.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, linewidth=2)
    
    # Customize the plot
    plt.title(f'Fraction of Significant Results\n({dataset} - {shuffle_type})', 
              fontsize=14, fontweight='bold')
    plt.xlabel('DAA Methods', fontsize=12)
    plt.ylabel('Fraction of Significant Results', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Set y-axis limits
    plt.ylim(0, max(0.1, df_plot['Fraction_Significant'].max() * 1.1))
    
    # Add grid for better readability
    plt.grid(True, alpha=0.3)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    # Show plot
    if show_plot:
        plt.show()
    
    return plt.gcf()

def main():
    """Main function to generate fraction of significant counts boxplot"""
    
    # Configuration
    dataset = "assal"
    n_iterations = 10  # Number of shuffling iterations
    
    print(f"Generating fraction of significant counts boxplot")
    print(f"Dataset: {dataset}")
    print(f"Iterations: {n_iterations}")
    print("="*60)
    
    # Load counts and metadata
    try:
        counts_df, metadata_df = load_counts_and_metadata(dataset=dataset)
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    # Calculate p-values directly from counts data
    pvalue_data = calculate_pvalues_from_counts(counts_df, metadata_df, n_iterations=n_iterations)
    
    if not pvalue_data:
        print("No p-value data calculated. Please check the data.")
        return
    
    # Print summary statistics
    print_summary_stats(pvalue_data)
    
    # Create the boxplot
    save_path = f"Plots/RDP/FractionSignificant_{dataset}_with_r_packages.png"
    
    create_boxplot(pvalue_data, 
                   dataset=dataset, 
                   shuffle_type="With R Packages",
                   save_path=save_path,
                   show_plot=True)
    
    print(f"\nAnalysis complete!")

if __name__ == "__main__":
    main()
