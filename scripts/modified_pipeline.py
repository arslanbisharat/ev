"""
Modified Analysis Pipeline - Separate Filtering Per Set
========================================================
This pipeline analyzes Set 1 and Set 2 independently, allowing proteins
unique to each set to be retained.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
sys.path.insert(0, project_dir)

from scripts.data_loader import ProteomicsDataLoader


def filter_proteins_per_set(df, group_cols, detection_threshold=0.4):
    """Filter proteins separately for each set"""

    # Set 1 filtering
    set1_sarc_det = (~df[group_cols['set1_sarc']].isna()).sum(axis=1) / len(group_cols['set1_sarc'])
    set1_ctrl_det = (~df[group_cols['set1_ctrl']].isna()).sum(axis=1) / len(group_cols['set1_ctrl'])
    set1_keep = (set1_sarc_det >= detection_threshold) | (set1_ctrl_det >= detection_threshold)

    # Set 2 filtering
    set2_sarc_det = (~df[group_cols['set2_sarc']].isna()).sum(axis=1) / len(group_cols['set2_sarc'])
    set2_ctrl_det = (~df[group_cols['set2_ctrl']].isna()).sum(axis=1) / len(group_cols['set2_ctrl'])
    set2_keep = (set2_sarc_det >= detection_threshold) | (set2_ctrl_det >= detection_threshold)

    return {
        'df_set1': df[set1_keep].copy(),
        'df_set2': df[set2_keep].copy(),
        'n_set1': set1_keep.sum(),
        'n_set2': set2_keep.sum(),
        'n_both': (set1_keep & set2_keep).sum(),
        'n_set1_only': (set1_keep & ~set2_keep).sum(),
        'n_set2_only': (set2_keep & ~set1_keep).sum()
    }


def preprocess_data(df, abundance_cols):
    """Log2 transform and impute missing values"""

    # Log2 transform
    data = df[abundance_cols].replace(0, np.nan)
    data_log2 = np.log2(data)

    # Impute missing values (minimum value - 1)
    data_imputed = data_log2.copy()
    for i in range(len(data_imputed)):
        protein_values = data_log2.iloc[i]
        if protein_values.notna().sum() > 0:
            min_val = protein_values.min()
            data_imputed.iloc[i] = protein_values.fillna(min_val - 1)

    return data_imputed


def perform_differential_analysis(data, df_info, group1_cols, group2_cols):
    """Perform t-test for differential expression"""

    results = []

    for i in range(len(data)):
        g1 = data.iloc[i][group1_cols].values
        g2 = data.iloc[i][group2_cols].values

        t_stat, p_val = stats.ttest_ind(g1, g2)

        pooled_std = np.sqrt(((len(g1)-1)*np.var(g1, ddof=1) + (len(g2)-1)*np.var(g2, ddof=1)) / (len(g1)+len(g2)-2))
        cohens_d = (np.mean(g1) - np.mean(g2)) / pooled_std if pooled_std > 0 else 0

        log2fc = np.mean(g1) - np.mean(g2)

        results.append({
            'Gene': df_info.iloc[i]['Gene Symbol'] if pd.notna(df_info.iloc[i]['Gene Symbol']) else 'Unknown',
            'Accession': df_info.iloc[i]['Accession'],
            'Description': df_info.iloc[i]['Description'],
            'Mean_Sarc': np.mean(g1),
            'Mean_Control': np.mean(g2),
            'Log2FC': log2fc,
            'P_value': p_val,
            'T_statistic': t_stat,
            'Cohens_d': cohens_d
        })

    df_results = pd.DataFrame(results)
    df_results['FDR'] = multipletests(df_results['P_value'], method='fdr_bh')[1]

    return df_results


def plot_volcano_comparison(results_set1, results_set2):
    """Create side-by-side volcano plots"""

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for idx, (results, title, ax) in enumerate([
        (results_set1, 'Set 1 (Old Samples)', axes[0]),
        (results_set2, 'Set 2 (Fresh Samples)', axes[1])
    ]):
        # Prepare data
        results = results.copy()
        results['-log10(FDR)'] = -np.log10(results['FDR'])
        results['Significant'] = (results['FDR'] < 0.05)

        # Separate significant from non-significant
        sig = results[results['Significant']]
        non_sig = results[~results['Significant']]

        # Plot
        ax.scatter(non_sig['Log2FC'], non_sig['-log10(FDR)'],
                  alpha=0.3, s=10, color='gray', label='Not significant')
        ax.scatter(sig['Log2FC'], sig['-log10(FDR)'],
                  alpha=0.7, s=20, color='red', label=f'Significant (n={len(sig)})')

        # Add threshold lines
        ax.axhline(y=-np.log10(0.05), color='blue', linestyle='--', linewidth=1.5, alpha=0.7, label='FDR = 0.05')
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1, alpha=0.3)

        # Labels and title
        ax.set_xlabel('Log2 Fold Change (Sarcoidosis vs Control)', fontsize=12, fontweight='bold')
        ax.set_ylabel('-Log10(FDR)', fontsize=12, fontweight='bold')
        ax.set_title(f'{title}\n{len(sig)} significant proteins', fontsize=13, fontweight='bold')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    figures_dir = os.path.join(project_dir, 'figures')
    plt.savefig(os.path.join(figures_dir, 'modified_volcano_comparison.png'), dpi=300, bbox_inches='tight')
    print(f"\n   Saved: figures/modified_volcano_comparison.png")
    plt.close()


def analyze_overlap(results_set1, results_set2, fdr_threshold=0.05):
    """Analyze overlap between significant proteins"""

    sig1 = results_set1[results_set1['FDR'] < fdr_threshold]
    sig2 = results_set2[results_set2['FDR'] < fdr_threshold]

    genes1 = set(sig1['Gene'].tolist())
    genes2 = set(sig2['Gene'].tolist())

    overlap = genes1 & genes2
    only_set1 = genes1 - genes2
    only_set2 = genes2 - genes1

    jaccard = len(overlap) / len(genes1 | genes2) if len(genes1 | genes2) > 0 else 0

    # Create Venn diagram
    fig, ax = plt.subplots(figsize=(10, 8))

    from matplotlib.patches import Circle

    # Draw circles
    circle1 = Circle((0.35, 0.5), 0.3, alpha=0.5, color='blue', label='Set 1 (Old)')
    circle2 = Circle((0.65, 0.5), 0.3, alpha=0.5, color='orange', label='Set 2 (Fresh)')
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # Add text
    ax.text(0.25, 0.5, f'{len(only_set1)}', fontsize=24, fontweight='bold', ha='center', va='center')
    ax.text(0.5, 0.5, f'{len(overlap)}', fontsize=24, fontweight='bold', ha='center', va='center')
    ax.text(0.75, 0.5, f'{len(only_set2)}', fontsize=24, fontweight='bold', ha='center', va='center')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    ax.text(0.5, 0.85, 'Overlap of Significant Proteins\n(Modified Analysis)',
            fontsize=16, fontweight='bold', ha='center')
    ax.text(0.5, 0.15, f'Jaccard Index: {jaccard:.3f}\nSet 1: {len(genes1)} | Set 2: {len(genes2)} | Overlap: {len(overlap)}',
            fontsize=13, ha='center')

    ax.legend([circle1, circle2], ['Set 1 (Old)', 'Set 2 (Fresh)'], loc='upper left', fontsize=12)

    plt.tight_layout()
    figures_dir = os.path.join(project_dir, 'figures')
    plt.savefig(os.path.join(figures_dir, 'modified_venn_diagram.png'), dpi=300, bbox_inches='tight')
    print(f"   Saved: figures/modified_venn_diagram.png")
    plt.close()

    return {
        'overlap': overlap,
        'only_set1': only_set1,
        'only_set2': only_set2,
        'jaccard': jaccard,
        'n_overlap': len(overlap),
        'n_only_set1': len(only_set1),
        'n_only_set2': len(only_set2),
        'sig1': sig1,
        'sig2': sig2
    }


def calculate_fc_correlation(results_set1, results_set2):
    """Calculate fold-change correlation between sets"""

    merged = results_set1[['Gene', 'Log2FC']].merge(
        results_set2[['Gene', 'Log2FC']],
        on='Gene',
        suffixes=('_Set1', '_Set2')
    )

    pearson_r = merged['Log2FC_Set1'].corr(merged['Log2FC_Set2'], method='pearson')
    spearman_r = merged['Log2FC_Set1'].corr(merged['Log2FC_Set2'], method='spearman')

    # Plot correlation
    fig, ax = plt.subplots(figsize=(9, 9))

    ax.scatter(merged['Log2FC_Set1'], merged['Log2FC_Set2'], alpha=0.4, s=15)
    ax.set_xlabel('Set 1 (Old) - Log2FC (Sarc vs Control)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Set 2 (Fresh) - Log2FC (Sarc vs Control)', fontsize=12, fontweight='bold')
    ax.set_title(f'Fold-Change Correlation Between Sets\n(Pearson r = {pearson_r:.3f}, n = {len(merged)})',
                 fontsize=13, fontweight='bold')

    # Add diagonal line
    min_val = min(merged['Log2FC_Set1'].min(), merged['Log2FC_Set2'].min())
    max_val = max(merged['Log2FC_Set1'].max(), merged['Log2FC_Set2'].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, alpha=0.5, label='y=x')

    # Add horizontal and vertical lines at 0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1, alpha=0.3)

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    figures_dir = os.path.join(project_dir, 'figures')
    plt.savefig(os.path.join(figures_dir, 'modified_fc_correlation.png'), dpi=300, bbox_inches='tight')
    print(f"   Saved: figures/modified_fc_correlation.png")
    plt.close()

    return {
        'pearson': pearson_r,
        'spearman': spearman_r,
        'n_shared': len(merged)
    }


def main():
    print("="*80)
    print("MODIFIED ANALYSIS - SEPARATE FILTERING PER SET")
    print("="*80)

    # Load data
    print("\n[1/6] Loading data")
    data_path = os.path.join(project_dir, 'proteome proteins.xlsx')
    loader = ProteomicsDataLoader(data_path)
    df = loader.load_data()
    loader.identify_sample_columns()
    group_cols = loader.get_group_columns()

    summary = loader.get_summary()
    print(f"  Proteins: {summary['total_proteins']}")
    print(f"  Set 1 (Old): {summary['set1_samples']} samples")
    print(f"  Set 2 (Fresh): {summary['set2_samples']} samples")

    # Filter proteins separately per set
    print("\n[2/6] Preprocessing (separate filtering per set)")
    filter_results = filter_proteins_per_set(df, group_cols, detection_threshold=0.4)

    print(f"  Set 1 kept: {filter_results['n_set1']} proteins")
    print(f"  Set 2 kept: {filter_results['n_set2']} proteins")
    print(f"  Both sets: {filter_results['n_both']} proteins")
    print(f"  Set 1 only: {filter_results['n_set1_only']} proteins")
    print(f"  Set 2 only: {filter_results['n_set2_only']} proteins")

    # Preprocess each set
    abundance_cols = loader.abundance_cols
    set1_cols = group_cols['set1_sarc'] + group_cols['set1_ctrl']
    set2_cols = group_cols['set2_sarc'] + group_cols['set2_ctrl']

    data_set1 = preprocess_data(filter_results['df_set1'], set1_cols)
    data_set2 = preprocess_data(filter_results['df_set2'], set2_cols)

    # Differential analysis
    print("\n[3/6] Differential analysis")
    results_set1 = perform_differential_analysis(
        data_set1, filter_results['df_set1'],
        group_cols['set1_sarc'], group_cols['set1_ctrl']
    )
    results_set2 = perform_differential_analysis(
        data_set2, filter_results['df_set2'],
        group_cols['set2_sarc'], group_cols['set2_ctrl']
    )

    sig1_count = (results_set1['FDR'] < 0.05).sum()
    sig2_count = (results_set2['FDR'] < 0.05).sum()

    print(f"  Set 1: {sig1_count} significant proteins (FDR < 0.05)")
    print(f"  Set 2: {sig2_count} significant proteins (FDR < 0.05)")

    # Overlap analysis
    print("\n[4/6] Overlap analysis")
    overlap_info = analyze_overlap(results_set1, results_set2)
    print(f"  Overlap: {overlap_info['n_overlap']} proteins")
    print(f"  Only Set 1: {overlap_info['n_only_set1']} proteins")
    print(f"  Only Set 2: {overlap_info['n_only_set2']} proteins")
    print(f"  Jaccard: {overlap_info['jaccard']:.3f}")

    # Fold-change correlation
    print("\n[5/6] Fold-change correlation")
    corr_info = calculate_fc_correlation(results_set1, results_set2)
    print(f"  Pearson r: {corr_info['pearson']:.3f} (n={corr_info['n_shared']} shared proteins)")

    # Visualizations
    print("\n[6/6] Creating visualizations")
    plot_volcano_comparison(results_set1, results_set2)

    # Save results
    print("\n[7/6] Saving results")
    results_dir = os.path.join(project_dir, 'results')
    results_set1.to_csv(os.path.join(results_dir, 'modified_results_set1_old.csv'), index=False)
    results_set2.to_csv(os.path.join(results_dir, 'modified_results_set2_fresh.csv'), index=False)

    # Save storage-stable proteins (significant in both)
    storage_stable = overlap_info['sig1'][overlap_info['sig1']['Gene'].isin(overlap_info['overlap'])]
    storage_stable.to_csv(os.path.join(results_dir, 'modified_storage_stable_proteins.csv'), index=False)

    print("\n" + "="*80)
    print("MODIFIED ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nKey Findings:")
    print(f"  - Set 1 (Old): {sig1_count} significant proteins")
    print(f"  - Set 2 (Fresh): {sig2_count} significant proteins")
    print(f"  - Overlap: {overlap_info['n_overlap']} proteins")
    print(f"  - Fold-change correlation: {corr_info['pearson']:.3f}")
    print(f"  - Jaccard index: {overlap_info['jaccard']:.3f}")


if __name__ == "__main__":
    main()
