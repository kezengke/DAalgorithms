#!/usr/bin/env python3
"""
Debug script to check fraction calculations
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

# Load data
rawT = pd.read_csv("CountsTables/RDPRaw/assal.txt", sep='\t', index_col=0)
metaData = pd.read_csv("MetaData/metadata_assal.txt", sep='\t', index_col=0)

# Clean and align data
rawT.columns = rawT.columns.str.strip('"')
common_samples = rawT.columns.intersection(metaData.index)
rawT = rawT[common_samples]
metaData = metaData.loc[common_samples]
metaData.index = rawT.columns
metaData.columns = ['conditions']

# Test multiple iterations
print("Testing fraction calculations across iterations:")
print("="*60)

all_pvals = []
for i in range(5):
    # Shuffle metadata
    shuffled_meta = metaData.copy()
    shuffled_meta['conditions'] = np.random.permutation(shuffled_meta['conditions'])
    
    # Calculate t-test p-values
    conditions = shuffled_meta['conditions'].unique()
    group_labels = np.where(shuffled_meta['conditions'] == conditions[0], 'group1', 'group2')
    
    pvals = []
    for taxon in rawT.index:
        values = rawT.loc[taxon].values
        group1_values = values[group_labels == 'group1']
        group2_values = values[group_labels == 'group2']
        
        if len(group1_values) > 1 and len(group2_values) > 1:
            try:
                _, p_val = ttest_ind(group1_values, group2_values)
                pvals.append(p_val)
            except:
                pvals.append(1.0)
        else:
            pvals.append(1.0)
    
    # Calculate fraction significant
    significant_count = sum(1 for p in pvals if p < 0.05)
    fraction = significant_count / len(pvals)
    
    print(f"Iteration {i+1}: {significant_count}/{len(pvals)} significant, fraction = {fraction:.4f}")
    print(f"  P-value range: {min(pvals):.4f} - {max(pvals):.4f}")
    print(f"  P-values < 0.05: {[f'{p:.4f}' for p in pvals if p < 0.05][:5]}...")
    
    all_pvals.extend(pvals)

print(f"\nOverall: {sum(1 for p in all_pvals if p < 0.05)}/{len(all_pvals)} significant")
print(f"Overall fraction: {sum(1 for p in all_pvals if p < 0.05) / len(all_pvals):.4f}")
print(f"Expected under null hypothesis: ~0.05")

