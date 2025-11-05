"""
Recommended Analysis Pipeline
Using 30% detection threshold + OR logic (literature-based best practices)
This demonstrates the improved approach without modifying original code.
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def load_data(data_path):
    """Load the proteomics dataset"""
    df = pd.read_excel(data_path)
    return df


def identify_sample_columns(df):
    """Identify and categorize sample columns"""
    abundance_cols = [col for col in df.columns if 'Abundances (Normalized)' in col]

    # Set 1 markers (old samples)
    set1_markers = ['302-', '310-', '311-', '327-', '338-', '343-', '346-', '364-', '368-', '371-',
                    '375-', '432-', '435-', '436-', '443-', '465-', '467-', '481-', '485-', '486-',
                    '500-', '510-', '533-', '537-', '545-', '561-', '567-', '577-', '584-', '613-']

    set1_cols = [col for col in abundance_cols if any(x in col for x in set1_markers)]
    set2_cols = [col for col in abundance_cols if ('Sarc' in col or ', C' in col) and col not in set1_cols]
    sarc_cols = [col for col in abundance_cols if 'Sarcoidosis' in col]
    control_cols = [col for col in abundance_cols if 'Control' in col]

    set1_sarc = [col for col in set1_cols if col in sarc_cols]
    set1_ctrl = [col for col in set1_cols if col in control_cols]
    set2_sarc = [col for col in set2_cols if col in sarc_cols]
    set2_ctrl = [col for col in set2_cols if col in control_cols]

    return {
        'abundance_cols': abundance_cols,
        'set1_sarc': set1_sarc,
        'set1_ctrl': set1_ctrl,
        'set2_sarc': set2_sarc,
        'set2_ctrl': set2_ctrl
    }


def filter_proteins_recommended(df, group_cols, detection_threshold=0.3):
    """
    RECOMMENDED APPROACH: Filter using OR logic with 30% threshold
    Keep proteins detected in ≥30% of at least ONE group
    """
    abundance_cols = group_cols['abundance_cols']

    # Calculate detection rates for each group
    set1_sarc_det = (~df[group_cols['set1_sarc']].isna()).sum(axis=1) / len(group_cols['set1_sarc'])
    set1_ctrl_det = (~df[group_cols['set1_ctrl']].isna()).sum(axis=1) / len(group_cols['set1_ctrl'])
    set2_sarc_det = (~df[group_cols['set2_sarc']].isna()).sum(axis=1) / len(group_cols['set2_sarc'])
    set2_ctrl_det = (~df[group_cols['set2_ctrl']].isna()).sum(axis=1) / len(group_cols['set2_ctrl'])

    # OR logic: Keep if ANY group passes threshold
    proteins_keep = ((set1_sarc_det >= detection_threshold) |
                     (set1_ctrl_det >= detection_threshold) |
                     (set2_sarc_det >= detection_threshold) |
                     (set2_ctrl_det >= detection_threshold))

    df_filtered = df[proteins_keep].copy()

    print(f"  Detection rates per group:")
    print(f"    Set1 Sarc: {set1_sarc_det[proteins_keep].mean():.1%} average")
    print(f"    Set1 Ctrl: {set1_ctrl_det[proteins_keep].mean():.1%} average")
    print(f"    Set2 Sarc: {set2_sarc_det[proteins_keep].mean():.1%} average")
    print(f"    Set2 Ctrl: {set2_ctrl_det[proteins_keep].mean():.1%} average")

    return df_filtered, proteins_keep.sum()


def preprocess_data(df_filtered, abundance_cols):
    """Log2 transform and impute missing values"""
    # Log2 transform
    abundance_data = df_filtered[abundance_cols].replace(0, np.nan)
    data_log2 = np.log2(abundance_data)

    # Impute with min-1 method
    data_imputed = data_log2.copy()
    for i in range(len(data_imputed)):
        protein_values = data_log2.iloc[i]
        if protein_values.notna().sum() > 0:
            min_val = protein_values.min()
            data_imputed.iloc[i] = protein_values.fillna(min_val - 1)

    return data_imputed


def perform_differential_analysis(data, protein_info, group1_cols, group2_cols):
    """Perform t-tests with FDR correction"""
    results = []

    for i in range(len(data)):
        g1 = data.iloc[i][group1_cols].values
        g2 = data.iloc[i][group2_cols].values

        t_stat, p_val = stats.ttest_ind(g1, g2)

        # Calculate effect size
        pooled_std = np.sqrt(((len(g1)-1)*np.var(g1, ddof=1) + (len(g2)-1)*np.var(g2, ddof=1)) / (len(g1)+len(g2)-2))
        cohens_d = (np.mean(g1) - np.mean(g2)) / pooled_std if pooled_std > 0 else 0

        log2fc = np.mean(g1) - np.mean(g2)

        results.append({
            'Gene': protein_info.iloc[i]['Gene Symbol'] if pd.notna(protein_info.iloc[i]['Gene Symbol']) else 'Unknown',
            'Accession': protein_info.iloc[i]['Accession'],
            'Description': protein_info.iloc[i]['Description'],
            'Mean_Group1': np.mean(g1),
            'Mean_Group2': np.mean(g2),
            'Log2FC': log2fc,
            'P_value': p_val,
            'T_statistic': t_stat,
            'Cohens_d': cohens_d
        })

    df_results = pd.DataFrame(results)
    df_results['FDR'] = multipletests(df_results['P_value'], method='fdr_bh')[1]

    return df_results


def analyze_thresholds(results, fdr_thresholds, fc_thresholds):
    """Analyze different FDR and FC threshold combinations"""
    summary = []

    for fdr in fdr_thresholds:
        for fc in fc_thresholds:
            sig = results[(results['FDR'] < fdr) & (abs(results['Log2FC']) > fc)]
            up = sig[sig['Log2FC'] > 0]
            down = sig[sig['Log2FC'] < 0]

            summary.append({
                'FDR_cutoff': fdr,
                'FC_cutoff': fc,
                'Total': len(sig),
                'Upregulated': len(up),
                'Downregulated': len(down)
            })

    return pd.DataFrame(summary)


def main():
    print("="*80)
    print("RECOMMENDED ANALYSIS: 30% THRESHOLD + OR LOGIC")
    print("="*80)

    # Setup paths
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(project_dir, 'proteome proteins.xlsx')
    output_dir = os.path.join(project_dir, 'results', 'recommended_approach')
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    print("\n[1/7] Loading data...")
    df = load_data(data_path)
    print(f"  Total proteins: {len(df)}")

    # Identify samples
    print("\n[2/7] Identifying sample groups...")
    group_cols = identify_sample_columns(df)
    print(f"  Set1 Sarc: {len(group_cols['set1_sarc'])} samples")
    print(f"  Set1 Ctrl: {len(group_cols['set1_ctrl'])} samples")
    print(f"  Set2 Sarc: {len(group_cols['set2_sarc'])} samples")
    print(f"  Set2 Ctrl: {len(group_cols['set2_ctrl'])} samples")

    # Filter with recommended approach
    print("\n[3/7] Filtering proteins (30% OR logic)...")
    df_filtered, n_kept = filter_proteins_recommended(df, group_cols, detection_threshold=0.3)
    print(f"  Proteins kept: {n_kept} / {len(df)} ({n_kept/len(df)*100:.1f}%)")
    print(f"  Proteins removed: {len(df) - n_kept} ({(len(df)-n_kept)/len(df)*100:.1f}%)")

    # Preprocess
    print("\n[4/7] Preprocessing (log2 + imputation)...")
    data_processed = preprocess_data(df_filtered, group_cols['abundance_cols'])
    print(f"  Log2 transformation: Done")
    print(f"  Missing value imputation: Done")

    # Differential analysis
    print("\n[5/7] Differential analysis...")
    print("  Running t-tests for Set 1 (Old)...")
    results_set1 = perform_differential_analysis(
        data_processed, df_filtered,
        group_cols['set1_sarc'], group_cols['set1_ctrl']
    )

    print("  Running t-tests for Set 2 (Fresh)...")
    results_set2 = perform_differential_analysis(
        data_processed, df_filtered,
        group_cols['set2_sarc'], group_cols['set2_ctrl']
    )

    # Threshold analysis
    print("\n[6/7] Testing FDR + FC threshold combinations...")
    fdr_thresholds = [0.001, 0.01, 0.05, 0.1, 0.2]
    fc_thresholds = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]

    summary_set1 = analyze_thresholds(results_set1, fdr_thresholds, fc_thresholds)
    summary_set2 = analyze_thresholds(results_set2, fdr_thresholds, fc_thresholds)

    # Key results
    print("\n[7/7] Saving results...")

    # Save full results
    results_set1.to_csv(os.path.join(output_dir, 'results_set1_recommended.csv'), index=False)
    results_set2.to_csv(os.path.join(output_dir, 'results_set2_recommended.csv'), index=False)

    # Save filtered results (FDR < 0.05, FC > 0.5)
    set1_sig = results_set1[(results_set1['FDR'] < 0.05) & (abs(results_set1['Log2FC']) > 0.5)]
    set2_sig = results_set2[(results_set2['FDR'] < 0.05) & (abs(results_set2['Log2FC']) > 0.5)]

    set1_sig.to_csv(os.path.join(output_dir, 'set1_significant_FDR0.05_FC0.5.csv'), index=False)
    set2_sig.to_csv(os.path.join(output_dir, 'set2_significant_FDR0.05_FC0.5.csv'), index=False)

    # Save overlap
    set1_genes = set(set1_sig['Gene'].tolist())
    set2_genes = set(set2_sig['Gene'].tolist())
    overlap_genes = set1_genes & set2_genes

    overlap_proteins = set1_sig[set1_sig['Gene'].isin(overlap_genes)]
    overlap_proteins.to_csv(os.path.join(output_dir, 'storage_stable_recommended.csv'), index=False)

    # Save threshold analysis
    summary_set1.to_csv(os.path.join(output_dir, 'threshold_analysis_set1.csv'), index=False)
    summary_set2.to_csv(os.path.join(output_dir, 'threshold_analysis_set2.csv'), index=False)

    # Print summary
    print("\n" + "="*80)
    print("RESULTS SUMMARY (30% OR LOGIC)")
    print("="*80)

    print(f"\nTotal proteins analyzed: {n_kept}")

    print(f"\nWith FDR < 0.05 + |Log2FC| > 0.5:")
    print(f"  Set 1 (Old):    {len(set1_sig)} significant proteins")
    print(f"    ├─ Up:        {len(set1_sig[set1_sig['Log2FC'] > 0])} proteins")
    print(f"    └─ Down:      {len(set1_sig[set1_sig['Log2FC'] < 0])} proteins")
    print(f"\n  Set 2 (Fresh):  {len(set2_sig)} significant proteins")
    print(f"    ├─ Up:        {len(set2_sig[set2_sig['Log2FC'] > 0])} proteins")
    print(f"    └─ Down:      {len(set2_sig[set2_sig['Log2FC'] < 0])} proteins")
    print(f"\n  Overlap:        {len(overlap_proteins)} storage-stable proteins")

    print(f"\nTop 5 proteins by |Log2FC| in Set 1:")
    top5_set1 = set1_sig.nlargest(5, lambda x: abs(x), columns='Log2FC')
    for idx, row in top5_set1.iterrows():
        direction = "↑" if row['Log2FC'] > 0 else "↓"
        print(f"  {direction} {row['Gene']:<12} Log2FC: {row['Log2FC']:>6.2f}, FDR: {row['FDR']:.2e}")

    print(f"\nTop 5 proteins by |Log2FC| in Set 2:")
    top5_set2 = set2_sig.nlargest(5, lambda x: abs(x), columns='Log2FC')
    for idx, row in top5_set2.iterrows():
        direction = "↑" if row['Log2FC'] > 0 else "↓"
        print(f"  {direction} {row['Gene']:<12} Log2FC: {row['Log2FC']:>6.2f}, FDR: {row['FDR']:.2e}")

    print(f"\nFiles saved to: {output_dir}/")
    print("  ✓ results_set1_recommended.csv")
    print("  ✓ results_set2_recommended.csv")
    print("  ✓ set1_significant_FDR0.05_FC0.5.csv")
    print("  ✓ set2_significant_FDR0.05_FC0.5.csv")
    print("  ✓ storage_stable_recommended.csv")
    print("  ✓ threshold_analysis_set1.csv")
    print("  ✓ threshold_analysis_set2.csv")

    # Print threshold matrix
    print("\n" + "="*80)
    print("THRESHOLD COMBINATION MATRIX")
    print("="*80)

    print("\nSet 1 (Old) - Significant Protein Counts:")
    print(f"{'FDR':<8} " + "  ".join(f"FC>{fc:<4.1f}" for fc in fc_thresholds))
    print("-" * 70)
    for fdr in fdr_thresholds:
        counts = summary_set1[summary_set1['FDR_cutoff'] == fdr]['Total'].tolist()
        print(f"{fdr:<8} " + "  ".join(f"{c:>6}" for c in counts))

    print("\nSet 2 (Fresh) - Significant Protein Counts:")
    print(f"{'FDR':<8} " + "  ".join(f"FC>{fc:<4.1f}" for fc in fc_thresholds))
    print("-" * 70)
    for fdr in fdr_thresholds:
        counts = summary_set2[summary_set2['FDR_cutoff'] == fdr]['Total'].tolist()
        print(f"{fdr:<8} " + "  ".join(f"{c:>6}" for c in counts))

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()
