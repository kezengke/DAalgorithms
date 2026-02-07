#!/usr/bin/env python3
"""
Script to read in RDP assal counts table and corresponding metadata
"""

import pandas as pd
import numpy as np

def read_assal_data():
    """
    Read in the RDP assal counts table and metadata
    
    """
    
    # Read in the counts table
    counts_file = "CountsTables/RDPRaw/assal.txt"
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
    
    # Remove quotes from column names if they exist
    counts_df.columns = counts_df.columns.str.strip('"')
    
    # Read in the metadata
    metadata_file = "MetaData/metadata_assal.txt"
    metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)
    
    # Align the data (keep only samples that exist in both)
    common_samples = counts_df.columns.intersection(metadata_df.index)
    
    counts_aligned = counts_df[common_samples]
    metadata_aligned = metadata_df.loc[common_samples]
    
    print(f"Common samples: {len(common_samples)}")
    print(f"Aligned counts shape: {counts_aligned.shape}")
    print(f"Aligned metadata shape: {metadata_aligned.shape}")
    
    # Check for any missing values
    print(f"\nMissing values in counts: {counts_aligned.isnull().sum().sum()}")
    print(f"Missing values in metadata: {metadata_aligned.isnull().sum().sum()}")
    
    return counts_aligned, metadata_aligned

def explore_data(counts_df, metadata_df):
    """
    Explore the loaded data
    
    Args:
        counts_df (pd.DataFrame): Counts data
        metadata_df (pd.DataFrame): Metadata
    """
    
    print("\n" + "="*50)
    print("DATA EXPLORATION")
    print("="*50)
    
    # Basic statistics
    print(f"\nCounts summary:")
    print(f"  Total features (taxa): {counts_df.shape[0]}")
    print(f"  Total samples: {counts_df.shape[1]}")
    print(f"  Mean counts per feature: {counts_df.mean().mean():.3f}")
    print(f"  Median counts per feature: {counts_df.median().median():.3f}")
    
    # Condition distribution
    condition_counts = metadata_df['Sample'].value_counts()
    print(f"\nCondition distribution:")
    for condition, count in condition_counts.items():
        print(f"  {condition}: {count} samples")
    
    # Zero counts
    print(f"\nFeatures with zero counts:")
    print(f"  Features with all zeros: {(counts_df == 0).all(axis=1).sum()}")
    print(f"  Mean percentage of zeros per feature: {(counts_df == 0).mean(axis=1).mean()*100:.1f}%")
    
    # Top features by mean abundance
    print(f"\nTop 10 features by mean abundance:")
    top_features = counts_df.mean(axis=1).sort_values(ascending=False).head(10)
    for feature, mean_abundance in top_features.items():
        print(f"  {feature}: {mean_abundance:.3f}")

if __name__ == "__main__":
    # Read the data
    counts_df, metadata_df = read_assal_data()
    
    # Explore the data
    explore_data(counts_df, metadata_df)
    
    print("\n" + "="*50)
    print("DATA LOADED SUCCESSFULLY!")
    print("="*50)
    print(f"Counts DataFrame: {counts_df.shape[0]} features x {counts_df.shape[1]} samples")
    print(f"Metadata DataFrame: {metadata_df.shape[0]} samples x {metadata_df.shape[1]} variables")
    print(f"Conditions: {', '.join(metadata_df['Sample'].unique())}")
    
    # Example of how to access the data
    print(f"\nExample data access:")
    first_taxon = counts_df.index[0]
    first_sample = counts_df.columns[0]
    print(f"  {first_taxon} in {first_sample}: {counts_df.loc[first_taxon, first_sample]:.3f}")
    print(f"  {first_sample} condition: {metadata_df.loc[first_sample, 'Sample']}")
    
    # You can now use counts_df and metadata_df for further analysis
