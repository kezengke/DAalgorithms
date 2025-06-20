import os
import pandas as pd
import numpy as np
from glob import glob
from scipy.stats import ttest_ind, ranksums

def resample_rnorm(table, meta, multiple):
    if table.shape[0] == 0 or meta.shape[0] == 0:
        raise ValueError("Input table or meta data frame is empty.")

    # Get the two sample groups
    groups = meta['conditions'].unique()
    group1, group2 = groups[0], groups[1]

    # Calculate group-wise mean and SD for one row
    def calculate_mean_sd(z):
        g1 = z[meta['conditions'] == group1]
        g2 = z[meta['conditions'] == group2]
        return {
            'mean1': np.mean(g1),
            'mean2': np.mean(g2),
            'sd1': np.std(g1),
            'sd2': np.std(g2)
        }

    # Calculate Mean/SD table
    mean_sd_list = []
    for _, row in table.iterrows():
        stats = calculate_mean_sd(row)
        mean_sd_list.append(stats)
    MeanSd_table = pd.DataFrame(mean_sd_list, index=table.index)

    # Resample counts for a single taxon
    def resample_counts(row_index):
        z = table.iloc[row_index]
        g1 = z[meta['conditions'] == group1]
        g2 = z[meta['conditions'] == group2]

        stats = MeanSd_table.iloc[row_index]
        mean1, sd1 = stats['mean1'], stats['sd1']
        mean2, sd2 = stats['mean2'], stats['sd2']

        ng1 = np.random.normal(loc=mean1, scale=np.sqrt(multiple * sd1 ** 2), size=len(g1))
        ng2 = np.random.normal(loc=mean2, scale=np.sqrt(multiple * sd2 ** 2), size=len(g2))

        nz = np.concatenate([ng1, ng2])
        nz = nz + abs(min(nz))
        return nz

    # Apply resampling across all taxa (rows)
    resampled_data = np.array([
        resample_counts(i) for i in range(table.shape[0])
    ])

    # Create new DataFrame
    sample_names = list(meta.index[meta['conditions'] == group1]) + \
                   list(meta.index[meta['conditions'] == group2])
    newT = pd.DataFrame(resampled_data, index=table.index, columns=sample_names)

    # Reorder columns to match original table
    newT = newT[table.columns]

    # Round to integers
    newT = newT.round(0).astype(int)

    return newT


def norm_fun(table):
    # Copy to avoid modifying original
    table = table.copy()

    # Column sums (total per sample)
    n = table.sum(axis=0)

    # Total sum of all values
    sumx = table.values.sum()

    # Normalize each column by its own sum
    table = table.div(n, axis=1)

    # Multiply by average sample total and add 1, then log10 transform
    table = np.log10(table * (sumx / table.shape[1]) + 1)

    # Return as DataFrame with original column names preserved
    return pd.DataFrame(table, index=table.index, columns=table.columns)


def calc_ttest(norm_counts_df, meta_df):
    groups = meta_df['conditions'].unique()
    pvals = []
    for taxon in norm_counts_df.index:
        group1 = norm_counts_df.loc[taxon, meta_df[meta_df['conditions'] == groups[0]].index]
        group2 = norm_counts_df.loc[taxon, meta_df[meta_df['conditions'] == groups[1]].index]
        stat, p = ttest_ind(group1, group2, equal_var=False)
        pvals.append(p)
    return pd.Series(pvals, index=norm_counts_df.index)


def calc_wilcox(norm_counts_df, meta_df):
    groups = meta_df['conditions'].unique()
    pvals = []
    for taxon in norm_counts_df.index:
        group1 = norm_counts_df.loc[taxon, meta_df[meta_df['conditions'] == groups[0]].index]
        group2 = norm_counts_df.loc[taxon, meta_df[meta_df['conditions'] == groups[1]].index]
        stat, p = ranksums(group1, group2)
        pvals.append(p)
    return pd.Series(pvals, index=norm_counts_df.index)


if __name__ == '__main__':
    # Load all files
    all_files = glob("CountsTables/RDPRaw/*.txt")
    print("Found files:", all_files)
    print(os.getcwd())

    for file in all_files:
        name = os.path.basename(file).replace(".txt", "")
        raw_counts = pd.read_csv(file, sep='\t', index_col=0)

        # Filter low-count taxa
        raw_counts = raw_counts[raw_counts.mean(axis=1) >= 2]

        meta = pd.read_csv(f"MetaData/metadata_{name}.txt", sep='\t', index_col=0)
        common_samples = raw_counts.columns.intersection(meta.index)
        meta = meta.loc[common_samples]
        raw_counts = raw_counts[common_samples]
        meta = meta.rename(columns={meta.columns[0]: 'conditions'})  # rename first column

        ttest_pvals = pd.DataFrame(index=raw_counts.index)
        wilcox_pvals = pd.DataFrame(index=raw_counts.index)

        for i in range(100):
            resampled = resample_rnorm(raw_counts, meta, 1)
            shuffled = pd.DataFrame({col: np.random.permutation(resampled[col].values) for col in resampled.columns},
                                    index=resampled.index)

            norm_counts = norm_fun(shuffled)

            ttest_pvals[i] = calc_ttest(norm_counts, meta)
            wilcox_pvals[i] = calc_wilcox(norm_counts, meta)

        save_dir = "ResampleShuffleInSampleDump/RDP/"
        os.makedirs(save_dir, exist_ok=True)
        ttest_pvals.to_csv(os.path.join(save_dir, f"{name}Python_t.txt"), sep='\t')
        wilcox_pvals.to_csv(os.path.join(save_dir, f"{name}Python_wilcox.txt"), sep='\t')

    print("Done.")
