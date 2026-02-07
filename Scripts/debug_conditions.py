#!/usr/bin/env python3
"""
Debug script to check condition values and group sizes
"""

import pandas as pd
import numpy as np

# Load data
rawT = pd.read_csv("CountsTables/RDPRaw/assal.txt", sep='\t', index_col=0)
metaData = pd.read_csv("MetaData/metadata_assal.txt", sep='\t', index_col=0)

# Clean column names
rawT.columns = rawT.columns.str.strip('"')

# Align data
common_samples = rawT.columns.intersection(metaData.index)
rawT = rawT[common_samples]
metaData = metaData.loc[common_samples]

# Prepare metadata
metaData.index = rawT.columns
metaData.columns = ['conditions']

print("Original conditions:")
print(metaData['conditions'].value_counts())
print(f"\nTotal samples: {len(metaData)}")
print(f"Unique conditions: {metaData['conditions'].unique()}")

# Test shuffling
print("\n" + "="*50)
print("After shuffling:")
shuffled_meta = metaData.copy()
shuffled_meta['conditions'] = np.random.permutation(shuffled_meta['conditions'])
print(shuffled_meta['conditions'].value_counts())

# Test group creation
conditions = shuffled_meta['conditions'].unique()
print(f"\nUnique conditions after shuffle: {conditions}")
group_labels = np.where(shuffled_meta['conditions'] == conditions[0], 'group1', 'group2')
print(f"Group labels: {np.unique(group_labels, return_counts=True)}")

# Test with a few taxa
print("\n" + "="*50)
print("Testing t-test with first few taxa:")
from scipy.stats import ttest_ind

for i, taxon in enumerate(rawT.index[:3]):
    values = rawT.loc[taxon].values
    group1_values = values[group_labels == 'group1']
    group2_values = values[group_labels == 'group2']
    
    print(f"\nTaxon {i+1}: {taxon}")
    print(f"  Group1 values: {group1_values[:5]}... (n={len(group1_values)})")
    print(f"  Group2 values: {group2_values[:5]}... (n={len(group2_values)})")
    
    if len(group1_values) > 1 and len(group2_values) > 1:
        try:
            t_stat, p_val = ttest_ind(group1_values, group2_values)
            print(f"  t-test: stat={t_stat:.4f}, pval={p_val:.4f}")
        except Exception as e:
            print(f"  t-test error: {e}")
    else:
        print(f"  Insufficient groups: group1={len(group1_values)}, group2={len(group2_values)}")

