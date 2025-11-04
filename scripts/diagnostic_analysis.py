"""
Diagnostic Analysis - Investigating Set 1 vs Set 2 Differences
==============================================================
This script investigates why there's low overlap between old and fresh samples.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
sys.path.insert(0, project_dir)

from scripts.data_loader import ProteomicsDataLoader

def analyze_filtering_stage(df, loader, group_cols):
    """Analyze what happens at the filtering stage"""

    print("\n" + "="*80)
    print("STEP 2/6 ANALYSIS - PROTEIN FILTERING STAGE")
    print("="*80)

    abundance_cols = loader.abundance_cols

    # Calculate detection rates for each group
    set1_sarc_det = (~df[group_cols['set1_sarc']].isna()).sum(axis=1) / len(group_cols['set1_sarc'])
    set1_ctrl_det = (~df[group_cols['set1_ctrl']].isna()).sum(axis=1) / len(group_cols['set1_ctrl'])
    set2_sarc_det = (~df[group_cols['set2_sarc']].isna()).sum(axis=1) / len(group_cols['set2_sarc'])
    set2_ctrl_det = (~df[group_cols['set2_ctrl']].isna()).sum(axis=1) / len(group_cols['set2_ctrl'])

    # Current filtering (requires BOTH sets to pass)
    set1_keep = (set1_sarc_det >= 0.4) | (set1_ctrl_det >= 0.4)
    set2_keep = (set2_sarc_det >= 0.4) | (set2_ctrl_det >= 0.4)
    both_keep = set1_keep & set2_keep

    print(f"\n1. PROTEIN DETECTION WITH 40% THRESHOLD:")
    print(f"   Proteins passing in Set 1 (Old): {set1_keep.sum()}")
    print(f"   Proteins passing in Set 2 (Fresh): {set2_keep.sum()}")
    print(f"   Proteins passing in BOTH sets: {both_keep.sum()}")
    print(f"   Proteins ONLY in Set 1: {(set1_keep & ~set2_keep).sum()}")
    print(f"   Proteins ONLY in Set 2: {(set2_keep & ~set1_keep).sum()}")

    # Calculate mean detection rates
    print(f"\n2. MEAN DETECTION RATES:")
    print(f"   Set 1 Sarcoidosis: {set1_sarc_det.mean():.1%}")
    print(f"   Set 1 Control: {set1_ctrl_det.mean():.1%}")
    print(f"   Set 2 Sarcoidosis: {set2_sarc_det.mean():.1%}")
    print(f"   Set 2 Control: {set2_ctrl_det.mean():.1%}")

    return {
        'set1_keep': set1_keep,
        'set2_keep': set2_keep,
        'both_keep': both_keep,
        'set1_sarc_det': set1_sarc_det,
        'set1_ctrl_det': set1_ctrl_det,
        'set2_sarc_det': set2_sarc_det,
        'set2_ctrl_det': set2_ctrl_det
    }


def analyze_overall_correlation(df, loader, filter_info):
    """Calculate correlation of all proteins between two sets"""

    print("\n" + "="*80)
    print("OVERALL CORRELATION ANALYSIS")
    print("="*80)

    group_cols = loader.get_group_columns()

    # Get proteins that pass in both sets
    df_both = df[filter_info['both_keep']].copy()

    # Calculate mean abundance for each set (log2 scale)
    set1_cols = group_cols['set1_sarc'] + group_cols['set1_ctrl']
    set2_cols = group_cols['set2_sarc'] + group_cols['set2_ctrl']

    # Replace 0 with NaN and log2 transform
    df_both_set1 = df_both[set1_cols].replace(0, np.nan)
    df_both_set2 = df_both[set2_cols].replace(0, np.nan)

    df_both_set1_log = np.log2(df_both_set1)
    df_both_set2_log = np.log2(df_both_set2)

    # Calculate mean for each protein in each set
    mean_set1 = df_both_set1_log.mean(axis=1)
    mean_set2 = df_both_set2_log.mean(axis=1)

    # Remove proteins with all missing values in either set
    valid_mask = mean_set1.notna() & mean_set2.notna()
    mean_set1_valid = mean_set1[valid_mask]
    mean_set2_valid = mean_set2[valid_mask]

    pearson_r = mean_set1_valid.corr(mean_set2_valid, method='pearson')
    spearman_r = mean_set1_valid.corr(mean_set2_valid, method='spearman')

    print(f"\n1. OVERALL PROTEIN ABUNDANCE CORRELATION (n={len(mean_set1_valid)} proteins):")
    print(f"   Pearson r: {pearson_r:.3f}")
    print(f"   Spearman rho: {spearman_r:.3f}")

    # Plot correlation
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(mean_set1_valid, mean_set2_valid, alpha=0.3, s=10)
    ax.set_xlabel('Set 1 (Old) - Mean Log2 Abundance', fontsize=12)
    ax.set_ylabel('Set 2 (Fresh) - Mean Log2 Abundance', fontsize=12)
    ax.set_title(f'Overall Protein Abundance Correlation\n(Pearson r = {pearson_r:.3f}, n = {len(mean_set1_valid)})',
                 fontsize=13, fontweight='bold')

    # Add diagonal line
    min_val = min(mean_set1_valid.min(), mean_set2_valid.min())
    max_val = max(mean_set1_valid.max(), mean_set2_valid.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, alpha=0.5, label='y=x')
    ax.legend()
    ax.grid(True, alpha=0.3)

    figures_dir = os.path.join(project_dir, 'figures')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'overall_abundance_correlation.png'), dpi=300, bbox_inches='tight')
    print(f"\n   Saved: figures/overall_abundance_correlation.png")
    plt.close()

    return {
        'mean_set1': mean_set1_valid,
        'mean_set2': mean_set2_valid,
        'pearson': pearson_r,
        'spearman': spearman_r
    }


def compare_top_proteins(df, loader, filter_info):
    """Compare top abundant proteins between sets"""

    print("\n" + "="*80)
    print("TOP PROTEIN COMPARISON")
    print("="*80)

    group_cols = loader.get_group_columns()
    df_both = df[filter_info['both_keep']].copy()

    set1_cols = group_cols['set1_sarc'] + group_cols['set1_ctrl']
    set2_cols = group_cols['set2_sarc'] + group_cols['set2_ctrl']

    # Calculate mean abundance (raw scale)
    df_both['Mean_Set1'] = df_both[set1_cols].replace(0, np.nan).mean(axis=1)
    df_both['Mean_Set2'] = df_both[set2_cols].replace(0, np.nan).mean(axis=1)

    # Get top 20 proteins by abundance in each set
    top_set1 = df_both.nlargest(20, 'Mean_Set1')[['Gene Symbol', 'Description', 'Mean_Set1', 'Mean_Set2']]
    top_set2 = df_both.nlargest(20, 'Mean_Set2')[['Gene Symbol', 'Description', 'Mean_Set1', 'Mean_Set2']]

    print("\n1. TOP 20 PROTEINS IN SET 1 (OLD):")
    print(top_set1.to_string(index=False))

    print("\n\n2. TOP 20 PROTEINS IN SET 2 (FRESH):")
    print(top_set2.to_string(index=False))

    # Check overlap
    top_genes_set1 = set(top_set1['Gene Symbol'].tolist())
    top_genes_set2 = set(top_set2['Gene Symbol'].tolist())
    overlap = top_genes_set1 & top_genes_set2

    print(f"\n3. OVERLAP IN TOP 20:")
    print(f"   Common proteins: {len(overlap)}")
    print(f"   Genes: {', '.join(sorted(overlap))}")

    return {
        'top_set1': top_set1,
        'top_set2': top_set2,
        'overlap': overlap
    }


def analyze_fold_change_distribution(df, loader, filter_info):
    """Analyze fold change distribution for Sarc vs Control in each set"""

    print("\n" + "="*80)
    print("FOLD CHANGE DISTRIBUTION ANALYSIS")
    print("="*80)

    group_cols = loader.get_group_columns()
    df_both = df[filter_info['both_keep']].copy()

    # Log2 transform
    set1_sarc = df_both[group_cols['set1_sarc']].replace(0, np.nan)
    set1_ctrl = df_both[group_cols['set1_ctrl']].replace(0, np.nan)
    set2_sarc = df_both[group_cols['set2_sarc']].replace(0, np.nan)
    set2_ctrl = df_both[group_cols['set2_ctrl']].replace(0, np.nan)

    set1_sarc_log = np.log2(set1_sarc)
    set1_ctrl_log = np.log2(set1_ctrl)
    set2_sarc_log = np.log2(set2_sarc)
    set2_ctrl_log = np.log2(set2_ctrl)

    # Calculate fold changes (Sarc - Control)
    fc_set1 = set1_sarc_log.mean(axis=1) - set1_ctrl_log.mean(axis=1)
    fc_set2 = set2_sarc_log.mean(axis=1) - set2_ctrl_log.mean(axis=1)

    # Remove NaN
    fc_set1_valid = fc_set1.dropna()
    fc_set2_valid = fc_set2.dropna()

    print(f"\n1. FOLD CHANGE STATISTICS:")
    print(f"   Set 1 (Old) - Mean: {fc_set1_valid.mean():.3f}, Median: {fc_set1_valid.median():.3f}, SD: {fc_set1_valid.std():.3f}")
    print(f"   Set 2 (Fresh) - Mean: {fc_set2_valid.mean():.3f}, Median: {fc_set2_valid.median():.3f}, SD: {fc_set2_valid.std():.3f}")

    print(f"\n2. PROTEINS WITH |Log2FC| > 0.5:")
    print(f"   Set 1: {(abs(fc_set1_valid) > 0.5).sum()} proteins")
    print(f"   Set 2: {(abs(fc_set2_valid) > 0.5).sum()} proteins")

    print(f"\n3. PROTEINS WITH |Log2FC| > 1.0:")
    print(f"   Set 1: {(abs(fc_set1_valid) > 1.0).sum()} proteins")
    print(f"   Set 2: {(abs(fc_set2_valid) > 1.0).sum()} proteins")

    # Plot distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].hist(fc_set1_valid, bins=50, alpha=0.7, edgecolor='black')
    axes[0].axvline(x=0, color='red', linestyle='--', linewidth=2)
    axes[0].set_xlabel('Log2 Fold Change (Sarc vs Control)', fontsize=11)
    axes[0].set_ylabel('Frequency', fontsize=11)
    axes[0].set_title(f'Set 1 (Old) - FC Distribution\nMean={fc_set1_valid.mean():.3f}, SD={fc_set1_valid.std():.3f}',
                     fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3)

    axes[1].hist(fc_set2_valid, bins=50, alpha=0.7, edgecolor='black', color='orange')
    axes[1].axvline(x=0, color='red', linestyle='--', linewidth=2)
    axes[1].set_xlabel('Log2 Fold Change (Sarc vs Control)', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_title(f'Set 2 (Fresh) - FC Distribution\nMean={fc_set2_valid.mean():.3f}, SD={fc_set2_valid.std():.3f}',
                     fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3)

    figures_dir = os.path.join(project_dir, 'figures')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'fold_change_distribution.png'), dpi=300, bbox_inches='tight')
    print(f"\n   Saved: figures/fold_change_distribution.png")
    plt.close()


def main():
    print("="*80)
    print("DIAGNOSTIC ANALYSIS - SET 1 VS SET 2 COMPARISON")
    print("="*80)

    # Load data
    data_path = os.path.join(project_dir, 'proteome proteins.xlsx')
    loader = ProteomicsDataLoader(data_path)
    df = loader.load_data()
    loader.identify_sample_columns()
    group_cols = loader.get_group_columns()

    summary = loader.get_summary()
    print(f"\nTotal proteins: {summary['total_proteins']}")
    print(f"Total samples: {summary['total_samples']}")

    # Analysis 1: Filtering stage
    filter_info = analyze_filtering_stage(df, loader, group_cols)

    # Analysis 2: Overall correlation
    corr_info = analyze_overall_correlation(df, loader, filter_info)

    # Analysis 3: Top proteins
    top_info = compare_top_proteins(df, loader, filter_info)

    # Analysis 4: Fold change distribution
    analyze_fold_change_distribution(df, loader, filter_info)

    print("\n" + "="*80)
    print("DIAGNOSTIC ANALYSIS COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()
