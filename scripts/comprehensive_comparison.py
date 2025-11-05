"""
Comprehensive Comparison Between Set 1 (Old) and Set 2 (Fresh)
================================================================
This script answers:
1. How many proteins detected in Controls vs Cases in each set?
2. Overall correlation between sets
3. Similarities/differences in protein abundance (Controls vs Controls, Cases vs Cases)
4. Missing % differences for individual proteins between sets
5. Step-by-step preprocessing documentation
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
sys.path.insert(0, project_dir)

from scripts.data_loader import ProteomicsDataLoader


def analyze_detection_by_group(df, group_cols):
    """
    Question 1: How many proteins detected in Controls vs Cases in each set?
    """
    print("\n" + "="*80)
    print("QUESTION 1: PROTEIN DETECTION BY GROUP")
    print("="*80)

    # Calculate detection (non-missing values) for each protein in each group
    set1_sarc = df[group_cols['set1_sarc']]
    set1_ctrl = df[group_cols['set1_ctrl']]
    set2_sarc = df[group_cols['set2_sarc']]
    set2_ctrl = df[group_cols['set2_ctrl']]

    # Detection = at least one non-missing value
    n_proteins = len(df)

    # Count how many proteins are detected in each group
    set1_sarc_detected = (~set1_sarc.isna()).any(axis=1).sum()
    set1_ctrl_detected = (~set1_ctrl.isna()).any(axis=1).sum()
    set2_sarc_detected = (~set2_sarc.isna()).any(axis=1).sum()
    set2_ctrl_detected = (~set2_ctrl.isna()).any(axis=1).sum()

    # Average detection rate per protein
    set1_sarc_det_rate = (~set1_sarc.isna()).sum(axis=1).mean() / len(group_cols['set1_sarc'])
    set1_ctrl_det_rate = (~set1_ctrl.isna()).sum(axis=1).mean() / len(group_cols['set1_ctrl'])
    set2_sarc_det_rate = (~set2_sarc.isna()).sum(axis=1).mean() / len(group_cols['set2_sarc'])
    set2_ctrl_det_rate = (~set2_ctrl.isna()).sum(axis=1).mean() / len(group_cols['set2_ctrl'])

    print(f"\n1. PROTEINS DETECTED (at least 1 sample):")
    print(f"   Total proteins in dataset: {n_proteins}")
    print(f"")
    print(f"   Set 1 (Old) - Sarcoidosis:  {set1_sarc_detected} ({set1_sarc_detected/n_proteins*100:.1f}%)")
    print(f"   Set 1 (Old) - Control:      {set1_ctrl_detected} ({set1_ctrl_detected/n_proteins*100:.1f}%)")
    print(f"   Set 2 (Fresh) - Sarcoidosis: {set2_sarc_detected} ({set2_sarc_detected/n_proteins*100:.1f}%)")
    print(f"   Set 2 (Fresh) - Control:     {set2_ctrl_detected} ({set2_ctrl_detected/n_proteins*100:.1f}%)")

    print(f"\n2. AVERAGE DETECTION RATE PER PROTEIN:")
    print(f"   Set 1 (Old) - Sarcoidosis:  {set1_sarc_det_rate*100:.1f}%")
    print(f"   Set 1 (Old) - Control:      {set1_ctrl_det_rate*100:.1f}%")
    print(f"   Set 2 (Fresh) - Sarcoidosis: {set2_sarc_det_rate*100:.1f}%")
    print(f"   Set 2 (Fresh) - Control:     {set2_ctrl_det_rate*100:.1f}%")

    # Comparison
    print(f"\n3. DIFFERENCES BETWEEN SETS:")
    print(f"   Old vs Fresh (Sarcoidosis): {set1_sarc_detected - set2_sarc_detected:+d} proteins")
    print(f"   Old vs Fresh (Control):     {set1_ctrl_detected - set2_ctrl_detected:+d} proteins")

    return {
        'set1_sarc_detected': set1_sarc_detected,
        'set1_ctrl_detected': set1_ctrl_detected,
        'set2_sarc_detected': set2_sarc_detected,
        'set2_ctrl_detected': set2_ctrl_detected,
        'set1_sarc_det_rate': set1_sarc_det_rate,
        'set1_ctrl_det_rate': set1_ctrl_det_rate,
        'set2_sarc_det_rate': set2_sarc_det_rate,
        'set2_ctrl_det_rate': set2_ctrl_det_rate
    }


def analyze_overall_correlation(df, group_cols):
    """
    Question 2: Overall correlation between two sets
    """
    print("\n" + "="*80)
    print("QUESTION 2: OVERALL CORRELATION BETWEEN SETS")
    print("="*80)

    # Get columns for each set
    set1_cols = group_cols['set1_sarc'] + group_cols['set1_ctrl']
    set2_cols = group_cols['set2_sarc'] + group_cols['set2_ctrl']

    # Calculate mean abundance for each protein in each set (log2 scale)
    df_set1 = df[set1_cols].replace(0, np.nan)
    df_set2 = df[set2_cols].replace(0, np.nan)

    df_set1_log = np.log2(df_set1)
    df_set2_log = np.log2(df_set2)

    mean_set1 = df_set1_log.mean(axis=1)
    mean_set2 = df_set2_log.mean(axis=1)

    # Remove proteins with all missing values
    valid_mask = mean_set1.notna() & mean_set2.notna()
    mean_set1_valid = mean_set1[valid_mask]
    mean_set2_valid = mean_set2[valid_mask]

    pearson_r, pearson_p = pearsonr(mean_set1_valid, mean_set2_valid)
    spearman_r, spearman_p = spearmanr(mean_set1_valid, mean_set2_valid)

    print(f"\n1. OVERALL PROTEIN ABUNDANCE CORRELATION:")
    print(f"   Number of proteins: {len(mean_set1_valid)}")
    print(f"   Pearson r:  {pearson_r:.4f} (p={pearson_p:.2e})")
    print(f"   Spearman ρ: {spearman_r:.4f} (p={spearman_p:.2e})")
    print(f"\n   Interpretation: {'VERY STRONG' if pearson_r > 0.9 else 'STRONG' if pearson_r > 0.7 else 'MODERATE'} correlation")

    return {
        'mean_set1': mean_set1_valid,
        'mean_set2': mean_set2_valid,
        'pearson_r': pearson_r,
        'spearman_r': spearman_r,
        'n_proteins': len(mean_set1_valid)
    }


def analyze_group_specific_correlations(df, group_cols):
    """
    Question 3: Controls vs Controls, Cases vs Cases comparisons
    """
    print("\n" + "="*80)
    print("QUESTION 3: GROUP-SPECIFIC ABUNDANCE COMPARISONS")
    print("="*80)

    # Calculate mean abundance for each group (log2 scale)
    set1_sarc_log = np.log2(df[group_cols['set1_sarc']].replace(0, np.nan))
    set1_ctrl_log = np.log2(df[group_cols['set1_ctrl']].replace(0, np.nan))
    set2_sarc_log = np.log2(df[group_cols['set2_sarc']].replace(0, np.nan))
    set2_ctrl_log = np.log2(df[group_cols['set2_ctrl']].replace(0, np.nan))

    mean_set1_sarc = set1_sarc_log.mean(axis=1)
    mean_set1_ctrl = set1_ctrl_log.mean(axis=1)
    mean_set2_sarc = set2_sarc_log.mean(axis=1)
    mean_set2_ctrl = set2_ctrl_log.mean(axis=1)

    # Cases vs Cases (Sarcoidosis Set1 vs Set2)
    valid_sarc = mean_set1_sarc.notna() & mean_set2_sarc.notna()
    sarc_r, sarc_p = pearsonr(mean_set1_sarc[valid_sarc], mean_set2_sarc[valid_sarc])

    # Controls vs Controls (Control Set1 vs Set2)
    valid_ctrl = mean_set1_ctrl.notna() & mean_set2_ctrl.notna()
    ctrl_r, ctrl_p = pearsonr(mean_set1_ctrl[valid_ctrl], mean_set2_ctrl[valid_ctrl])

    print(f"\n1. SARCOIDOSIS (Cases) - Old vs Fresh:")
    print(f"   Number of proteins: {valid_sarc.sum()}")
    print(f"   Pearson r: {sarc_r:.4f} (p={sarc_p:.2e})")

    print(f"\n2. CONTROLS - Old vs Fresh:")
    print(f"   Number of proteins: {valid_ctrl.sum()}")
    print(f"   Pearson r: {ctrl_r:.4f} (p={ctrl_p:.2e})")

    print(f"\n3. COMPARISON:")
    print(f"   Cases correlation:    {sarc_r:.4f}")
    print(f"   Controls correlation: {ctrl_r:.4f}")
    print(f"   Difference: {abs(sarc_r - ctrl_r):.4f}")

    if abs(sarc_r - ctrl_r) < 0.05:
        print(f"   → Similar correlations for both groups")
    else:
        better = "Cases" if sarc_r > ctrl_r else "Controls"
        print(f"   → {better} show higher correlation between sets")

    return {
        'sarc_r': sarc_r,
        'ctrl_r': ctrl_r,
        'mean_set1_sarc': mean_set1_sarc[valid_sarc],
        'mean_set2_sarc': mean_set2_sarc[valid_sarc],
        'mean_set1_ctrl': mean_set1_ctrl[valid_ctrl],
        'mean_set2_ctrl': mean_set2_ctrl[valid_ctrl]
    }


def analyze_missing_data(df, group_cols):
    """
    Question 4: Missing % differences between sets
    """
    print("\n" + "="*80)
    print("QUESTION 4: MISSING DATA ANALYSIS")
    print("="*80)

    set1_cols = group_cols['set1_sarc'] + group_cols['set1_ctrl']
    set2_cols = group_cols['set2_sarc'] + group_cols['set2_ctrl']

    # Calculate missing % for each protein in each set
    set1_missing_pct = df[set1_cols].isna().sum(axis=1) / len(set1_cols) * 100
    set2_missing_pct = df[set2_cols].isna().sum(axis=1) / len(set2_cols) * 100

    # Overall statistics
    print(f"\n1. OVERALL MISSING DATA:")
    print(f"   Set 1 (Old):  Mean={set1_missing_pct.mean():.1f}%, Median={set1_missing_pct.median():.1f}%")
    print(f"   Set 2 (Fresh): Mean={set2_missing_pct.mean():.1f}%, Median={set2_missing_pct.median():.1f}%")
    print(f"   Difference: {set1_missing_pct.mean() - set2_missing_pct.mean():.1f}% (Old - Fresh)")

    # Proteins with large differences
    missing_diff = set1_missing_pct - set2_missing_pct

    # More missing in Set 1 (Old)
    more_missing_set1 = df[missing_diff > 30].copy()
    more_missing_set1['Missing_Set1'] = set1_missing_pct[missing_diff > 30]
    more_missing_set1['Missing_Set2'] = set2_missing_pct[missing_diff > 30]
    more_missing_set1['Difference'] = missing_diff[missing_diff > 30]

    # More missing in Set 2 (Fresh)
    more_missing_set2 = df[missing_diff < -30].copy()
    more_missing_set2['Missing_Set1'] = set1_missing_pct[missing_diff < -30]
    more_missing_set2['Missing_Set2'] = set2_missing_pct[missing_diff < -30]
    more_missing_set2['Difference'] = missing_diff[missing_diff < -30]

    print(f"\n2. PROTEINS WITH LARGE MISSING DATA DIFFERENCES (>30%):")
    print(f"   More missing in Set 1 (Old): {len(more_missing_set1)} proteins")
    print(f"   More missing in Set 2 (Fresh): {len(more_missing_set2)} proteins")

    # Proteins completely missing in one set but present in another
    only_set1 = (set1_missing_pct < 100) & (set2_missing_pct == 100)
    only_set2 = (set1_missing_pct == 100) & (set2_missing_pct < 100)

    print(f"\n3. PROTEINS UNIQUE TO ONE SET:")
    print(f"   Only in Set 1 (Old): {only_set1.sum()} proteins")
    print(f"   Only in Set 2 (Fresh): {only_set2.sum()} proteins")

    # Save detailed results
    results_dir = os.path.join(project_dir, 'results')

    missing_summary = pd.DataFrame({
        'Gene': df['Gene Symbol'],
        'Accession': df['Accession'],
        'Missing_Set1_pct': set1_missing_pct,
        'Missing_Set2_pct': set2_missing_pct,
        'Difference_pct': missing_diff
    })
    missing_summary = missing_summary.sort_values('Difference_pct', key=abs, ascending=False)
    missing_summary.to_csv(os.path.join(results_dir, 'missing_data_comparison.csv'), index=False)
    print(f"\n   Saved: results/missing_data_comparison.csv")

    # Show top examples
    print(f"\n4. TOP 10 PROTEINS MORE MISSING IN SET 1 (OLD):")
    if len(more_missing_set1) > 0:
        top_set1 = more_missing_set1.nlargest(10, 'Difference')[['Gene Symbol', 'Missing_Set1', 'Missing_Set2', 'Difference']]
        print(top_set1.to_string(index=False))

    print(f"\n5. TOP 10 PROTEINS MORE MISSING IN SET 2 (FRESH):")
    if len(more_missing_set2) > 0:
        top_set2 = more_missing_set2.nsmallest(10, 'Difference')[['Gene Symbol', 'Missing_Set1', 'Missing_Set2', 'Difference']]
        print(top_set2.to_string(index=False))

    return {
        'set1_missing_pct': set1_missing_pct,
        'set2_missing_pct': set2_missing_pct,
        'more_missing_set1': len(more_missing_set1),
        'more_missing_set2': len(more_missing_set2),
        'only_set1': only_set1.sum(),
        'only_set2': only_set2.sum()
    }


def create_comprehensive_plots(overall_corr, group_corr, missing_info):
    """Create comprehensive visualization"""

    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # 1. Overall correlation
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(overall_corr['mean_set1'], overall_corr['mean_set2'], alpha=0.3, s=10)
    ax1.set_xlabel('Set 1 (Old) - Mean Log2 Abundance', fontsize=10)
    ax1.set_ylabel('Set 2 (Fresh) - Mean Log2 Abundance', fontsize=10)
    ax1.set_title(f'Overall Correlation\nr = {overall_corr["pearson_r"]:.3f}', fontsize=11, fontweight='bold')
    min_val = min(overall_corr['mean_set1'].min(), overall_corr['mean_set2'].min())
    max_val = max(overall_corr['mean_set1'].max(), overall_corr['mean_set2'].max())
    ax1.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5)
    ax1.grid(True, alpha=0.3)

    # 2. Cases correlation
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(group_corr['mean_set1_sarc'], group_corr['mean_set2_sarc'], alpha=0.3, s=10, color='red')
    ax2.set_xlabel('Set 1 (Old) - Sarcoidosis Mean', fontsize=10)
    ax2.set_ylabel('Set 2 (Fresh) - Sarcoidosis Mean', fontsize=10)
    ax2.set_title(f'Cases (Sarcoidosis)\nr = {group_corr["sarc_r"]:.3f}', fontsize=11, fontweight='bold')
    min_val = min(group_corr['mean_set1_sarc'].min(), group_corr['mean_set2_sarc'].min())
    max_val = max(group_corr['mean_set1_sarc'].max(), group_corr['mean_set2_sarc'].max())
    ax2.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
    ax2.grid(True, alpha=0.3)

    # 3. Controls correlation
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.scatter(group_corr['mean_set1_ctrl'], group_corr['mean_set2_ctrl'], alpha=0.3, s=10, color='blue')
    ax3.set_xlabel('Set 1 (Old) - Control Mean', fontsize=10)
    ax3.set_ylabel('Set 2 (Fresh) - Control Mean', fontsize=10)
    ax3.set_title(f'Controls\nr = {group_corr["ctrl_r"]:.3f}', fontsize=11, fontweight='bold')
    min_val = min(group_corr['mean_set1_ctrl'].min(), group_corr['mean_set2_ctrl'].min())
    max_val = max(group_corr['mean_set1_ctrl'].max(), group_corr['mean_set2_ctrl'].max())
    ax3.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
    ax3.grid(True, alpha=0.3)

    # 4. Missing data histogram
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.hist(missing_info['set1_missing_pct'], bins=50, alpha=0.5, label='Set 1 (Old)', color='orange')
    ax4.hist(missing_info['set2_missing_pct'], bins=50, alpha=0.5, label='Set 2 (Fresh)', color='cyan')
    ax4.set_xlabel('Missing Data %', fontsize=10)
    ax4.set_ylabel('Number of Proteins', fontsize=10)
    ax4.set_title('Missing Data Distribution', fontsize=11, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Missing data difference
    ax5 = fig.add_subplot(gs[1, 1])
    missing_diff = missing_info['set1_missing_pct'] - missing_info['set2_missing_pct']
    ax5.hist(missing_diff, bins=50, edgecolor='black', alpha=0.7)
    ax5.axvline(x=0, color='red', linestyle='--', linewidth=2)
    ax5.set_xlabel('Missing % Difference (Set1 - Set2)', fontsize=10)
    ax5.set_ylabel('Number of Proteins', fontsize=10)
    ax5.set_title('Missing Data Difference\n(Positive = More missing in Old)', fontsize=11, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # 6. Summary statistics
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')

    summary_text = f"""
SUMMARY STATISTICS

Overall Correlation: {overall_corr['pearson_r']:.3f}
  Cases correlation: {group_corr['sarc_r']:.3f}
  Controls correlation: {group_corr['ctrl_r']:.3f}

Missing Data (Mean %):
  Set 1 (Old): {missing_info['set1_missing_pct'].mean():.1f}%
  Set 2 (Fresh): {missing_info['set2_missing_pct'].mean():.1f}%

Proteins with >30% difference:
  More in Set 1: {missing_info['more_missing_set1']}
  More in Set 2: {missing_info['more_missing_set2']}

Unique proteins:
  Only Set 1: {missing_info['only_set1']}
  Only Set 2: {missing_info['only_set2']}
"""
    ax6.text(0.1, 0.5, summary_text, fontsize=10, family='monospace',
             verticalalignment='center')

    # Save figure
    figures_dir = os.path.join(project_dir, 'figures')
    plt.savefig(os.path.join(figures_dir, 'comprehensive_comparison.png'), dpi=300, bbox_inches='tight')
    print(f"\n   Saved: figures/comprehensive_comparison.png")
    plt.close()


def document_preprocessing_methods():
    """
    Question 5: Step-by-step preprocessing methods documentation
    """
    print("\n" + "="*80)
    print("QUESTION 5: PREPROCESSING PIPELINE DOCUMENTATION")
    print("="*80)

    methods = """
STEP-BY-STEP PREPROCESSING METHODS
===================================

1. DATA LOADING
   - Input: Excel file with protein abundance data (2,800 proteins × 120 samples)
   - Sample identification: Automatic detection of sample columns based on naming patterns
   - Group assignment: Set 1/2, Sarcoidosis/Control based on sample metadata

2. QUALITY CONTROL (BEFORE FILTERING)
   Purpose: Identify low-quality samples and proteins

   Sample QC metrics:
   - Missing data % per sample (flag if >50%)
   - Total intensity per sample (flag outliers)
   - PCA outlier detection (Mahalanobis distance)

   Protein QC metrics:
   - Detection rate across samples
   - Coefficient of variation (CV)
   - Missing value patterns

   Decision: Review flagged samples/proteins before proceeding

3. PROTEIN FILTERING
   Purpose: Remove proteins with too much missing data

   Criteria: Keep proteins detected in ≥40% of samples in AT LEAST ONE CONDITION
   - For Set 1: Keep if ≥40% in Sarc OR ≥40% in Control
   - For Set 2: Keep if ≥40% in Sarc OR ≥40% in Control

   Rationale:
   - 40% threshold balances data quality vs. coverage
   - "OR" logic preserves disease-specific proteins
   - Applied per-set to avoid bias

   Alternative (original): Require presence in BOTH sets
   - More conservative, but loses set-specific proteins

   Result: ~1,600-1,800 proteins retained per set

4. LOG2 TRANSFORMATION
   Purpose: Stabilize variance and normalize distribution

   Method:
   - Replace 0 values with NA (true missing)
   - Apply log2 transformation: log2(abundance)

   Criteria for use:
   - Proteomics data is typically log-normally distributed
   - Fold-changes become symmetric (2x up = -2x down)
   - Stabilizes variance across abundance range

   When to skip:
   - Data already log-transformed
   - Using rank-based methods that don't assume normality

5. IMPUTATION
   Purpose: Handle missing values for statistical testing

   Method: Minimum value - 1 approach
   - For each protein, find minimum observed value
   - Impute missing as min - 1 (on log2 scale)

   Rationale:
   - Missing values in proteomics often = below detection limit
   - This assumes missing = low abundance
   - Conservative approach that doesn't inflate significance

   Alternatives:
   - KNN imputation (use similar proteins)
   - Half-minimum (min/2 on linear scale)
   - Random draws from low end of distribution

   When to skip:
   - Using methods that handle missing data directly
   - Too much missing data (>70% may be unreliable)

6. BATCH EFFECT CORRECTION
   Status: NOT APPLIED in current analysis

   Why not:
   - Set 1 and Set 2 are analyzed SEPARATELY
   - No attempt to merge or directly compare samples
   - Comparison is done at summary statistic level (fold-changes)

   When to apply:
   - If combining sets for joint analysis
   - If technical batch effects detected within a set

   Method (if needed):
   - ComBat (from sva package in R)
   - Remove batch variation while preserving biological variation
   - Requires batch labels and ideally balanced design

7. NORMALIZATION
   Status: NOT EXPLICITLY APPLIED

   Current approach:
   - Log2 transformation provides implicit normalization
   - Within-sample normalization assumed from data provider

   When to apply:
   - If samples have different total protein amounts
   - If systematic differences in abundance distributions

   Methods:
   - Median centering (shift each sample to same median)
   - Quantile normalization (make distributions identical)
   - Total intensity normalization (divide by total per sample)

   Check needed:
   - Compare total intensity distributions
   - If similar, normalization may not be needed

8. STATISTICAL TESTING
   Method: Student's t-test with FDR correction

   Test: Welch's t-test (unequal variance)
   - Compare Sarcoidosis vs Control within each set
   - Assumes: independence, approximate normality (CLT applies with n=30)

   Multiple testing correction:
   - Benjamini-Hochberg FDR (False Discovery Rate)
   - Threshold: FDR < 0.05

   Effect size:
   - Log2 fold-change (Sarc - Control)
   - Cohen's d for standardized effect size

DECISION TREE FOR PREPROCESSING
================================

1. Start with raw data
   ↓
2. QC → Remove bad samples? (if >50% missing or clear outlier)
   ↓
3. Filter proteins → Use 30-40% detection threshold
   ↓
4. Check distribution → Skewed? → YES → Log2 transform
   ↓                              ↓ NO
5. Missing data <70%? → YES → Impute (min-1 or KNN)
   ↓                     ↓ NO → Consider removing protein
6. Multiple sets/batches? → YES → Check for batch effects → Correct if needed
   ↓                         ↓ NO
7. Check total intensity → Different? → YES → Normalize
   ↓                         ↓ NO
8. Statistical testing (t-test, FDR correction)

"""
    print(methods)

    # Save to file
    with open(os.path.join(project_dir, 'PREPROCESSING_METHODS.md'), 'w') as f:
        f.write("# Preprocessing Methods for EV Proteomics Analysis\n\n")
        f.write(methods)

    print("\n   Saved: PREPROCESSING_METHODS.md")


def main():
    print("="*80)
    print("COMPREHENSIVE COMPARISON: SET 1 (OLD) vs SET 2 (FRESH)")
    print("="*80)

    # Load data
    data_path = os.path.join(project_dir, 'proteome proteins.xlsx')
    loader = ProteomicsDataLoader(data_path)
    df = loader.load_data()
    loader.identify_sample_columns()
    group_cols = loader.get_group_columns()

    # Question 1: Detection by group
    detection_info = analyze_detection_by_group(df, group_cols)

    # Question 2: Overall correlation
    overall_corr = analyze_overall_correlation(df, group_cols)

    # Question 3: Group-specific correlations
    group_corr = analyze_group_specific_correlations(df, group_cols)

    # Question 4: Missing data analysis
    missing_info = analyze_missing_data(df, group_cols)

    # Create comprehensive plots
    print("\n" + "="*80)
    print("CREATING COMPREHENSIVE VISUALIZATIONS")
    print("="*80)
    create_comprehensive_plots(overall_corr, group_corr, missing_info)

    # Question 5: Methods documentation
    document_preprocessing_methods()

    print("\n" + "="*80)
    print("COMPREHENSIVE ANALYSIS COMPLETE")
    print("="*80)
    print("\nGenerated files:")
    print("  - results/missing_data_comparison.csv")
    print("  - figures/comprehensive_comparison.png")
    print("  - PREPROCESSING_METHODS.md")


if __name__ == "__main__":
    main()
