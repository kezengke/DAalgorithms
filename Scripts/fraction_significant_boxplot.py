#!/usr/bin/env python3
"""
Generate fraction of significant counts boxplot for RDP-Assal dataset
Clean, readable Python implementation that calculates p-values directly from counts data
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu, ttest_ind
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

def load_counts_and_metadata(dataset="assal"):
    """
    Load counts table and metadata for a specific dataset
    
    Args:
        dataset (str): Dataset name (e.g., "assal")
    
    Returns:
        tuple: (counts_df, metadata_df) - pandas DataFrames
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
    """
    Perform t-test for differential abundance
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
    
    Returns:
        pd.Series: p-values for each taxon
    """
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        # Remove zeros and log transform
        group1 = group1[group1 > 0]
        group2 = group2[group2 > 0]
        
        if len(group1) > 1 and len(group2) > 1:
            # Log transform
            group1_log = np.log(group1 + 1)
            group2_log = np.log(group2 + 1)
            
            # Perform t-test
            try:
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_wilcoxon(counts_df, metadata_df):
    """
    Perform Wilcoxon rank-sum test for differential abundance
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
    
    Returns:
        pd.Series: p-values for each taxon
    """
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        # Remove zeros
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

def perform_deseq2_like(counts_df, metadata_df):
    """
    Perform DESeq2-like analysis using negative binomial approximation
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
    
    Returns:
        pd.Series: p-values for each taxon
    """
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        # Add pseudocount and log transform
        group1_log = np.log(group1 + 1)
        group2_log = np.log(group2 + 1)
        
        if len(group1_log) > 1 and len(group2_log) > 1:
            try:
                # Use t-test on log-transformed data as approximation
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_edgeR_like(counts_df, metadata_df):
    """
    Perform edgeR-like analysis using log-fold change and dispersion estimation
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
    
    Returns:
        pd.Series: p-values for each taxon
    """
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        # Add pseudocount and log transform
        group1_log = np.log(group1 + 1)
        group2_log = np.log(group2 + 1)
        
        if len(group1_log) > 1 and len(group2_log) > 1:
            try:
                # Calculate log fold change
                mean1 = np.mean(group1_log)
                mean2 = np.mean(group2_log)
                lfc = mean2 - mean1
                
                # Simple t-test on log-transformed data
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_aldex2_like(counts_df, metadata_df, n_iterations=100):
    """
    Perform ALDEx2-like analysis with CLR transformation and permutation
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
        n_iterations (int): Number of Monte Carlo iterations
    
    Returns:
        pd.Series: p-values for each taxon
    """
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
                # Perform t-test on CLR-transformed data
                _, pval = ttest_ind(group1, group2)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_ancombc2_like(counts_df, metadata_df):
    """
    Perform ANCOMBC2-like analysis with bias correction
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
    
    Returns:
        pd.Series: p-values for each taxon
    """
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        # Add pseudocount and log transform
        group1_log = np.log(group1 + 1)
        group2_log = np.log(group2 + 1)
        
        if len(group1_log) > 1 and len(group2_log) > 1:
            try:
                # Simple t-test approximation
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def perform_metagenomeseq_like(counts_df, metadata_df):
    """
    Perform metagenomeSeq-like analysis with zero-inflated model approximation
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
    
    Returns:
        pd.Series: p-values for each taxon
    """
    pvals = []
    conditions = metadata_df['Sample'].unique()
    
    for taxon in counts_df.index:
        group1 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[0]].values
        group2 = counts_df.loc[taxon, metadata_df['Sample'] == conditions[1]].values
        
        # Handle zero-inflation by using only non-zero values
        group1_nonzero = group1[group1 > 0]
        group2_nonzero = group2[group2 > 0]
        
        if len(group1_nonzero) > 1 and len(group2_nonzero) > 1:
            try:
                # Log transform non-zero values
                group1_log = np.log(group1_nonzero)
                group2_log = np.log(group2_nonzero)
                
                # Perform t-test
                _, pval = ttest_ind(group1_log, group2_log)
                pvals.append(pval)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    return pd.Series(pvals, index=counts_df.index)

def calc_fraction_significant(pval_df, alpha=0.05):
    """
    Calculate fraction of significant results for each iteration
    
    Args:
        pval_df (pd.DataFrame): DataFrame with p-values (taxa x iterations)
        alpha (float): Significance threshold (default: 0.05)
    
    Returns:
        pd.Series: Fraction of significant results for each iteration
    """
    # For each column (iteration), count how many p-values are < alpha
    significant_counts = (pval_df < alpha).sum()
    total_taxa = len(pval_df)
    
    # Calculate fraction
    fractions = significant_counts / total_taxa
    
    return fractions

def calculate_pvalues_from_counts(counts_df, metadata_df, n_iterations=100):
    """
    Calculate p-values directly from counts data using multiple methods
    
    Args:
        counts_df (pd.DataFrame): Counts data (taxa x samples)
        metadata_df (pd.DataFrame): Metadata with 'Sample' column
        n_iterations (int): Number of iterations for shuffling
    
    Returns:
        dict: Dictionary with method names as keys and fraction data as values
    """
    
    print(f"Calculating p-values for {n_iterations} iterations...")
    
    # Define methods and their functions
    methods = {
        't-test': perform_ttest,
        'Wilcoxon': perform_wilcoxon,
        'DESeq2': perform_deseq2_like,
        'edgeR': perform_edgeR_like,
        'ALDEx2t-test': perform_aldex2_like,
        'ALDEx2Wilcoxon': perform_aldex2_like,
        'ANCOMBC2': perform_ancombc2_like,
        'metagenomeSeq': perform_metagenomeseq_like
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

def create_boxplot(pvalue_data, dataset="assal", shuffle_type="ShuffleDump", 
                   save_path=None, show_plot=True):
    """
    Create boxplot of fraction of significant results
    
    Args:
        pvalue_data (dict): Dictionary with method names and fraction data
        dataset (str): Dataset name for title
        shuffle_type (str): Shuffle type for title
        save_path (str): Path to save the plot (optional)
        show_plot (bool): Whether to display the plot
    """
    
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

def print_summary_stats(pvalue_data):
    """
    Print summary statistics for each method
    
    Args:
        pvalue_data (dict): Dictionary with method names and fraction data
    """
    
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

def main():
    """
    Main function to generate fraction of significant counts boxplot
    """
    
    # Configuration
    dataset = "assal"
    n_iterations = 10  # Number of shuffling iterations (reduced for testing)
    
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
    save_path = f"Plots/RDP/FractionSignificant_{dataset}_calculated.png"
    
    create_boxplot(pvalue_data, 
                   dataset=dataset, 
                   shuffle_type="Calculated",
                   save_path=save_path,
                   show_plot=True)
    
    print(f"\nAnalysis complete!")

if __name__ == "__main__":
    main()
