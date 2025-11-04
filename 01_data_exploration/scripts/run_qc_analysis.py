import os
import sys
import pandas as pd
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, script_dir)

from data_loader import DataLoader
from quality_control import ProteomicsQC
from descriptive_stats import DescriptiveStatistics
from qc_visualization import QCVisualizer


def run_comprehensive_qc():

    print("="*80)
    print("DATA EXPLORATION & QUALITY CONTROL ANALYSIS")
    print("="*80)
    print()

    data_file = os.path.join(project_dir, 'proteome proteins.xlsx')
    figures_dir = os.path.join(script_dir, '..', 'figures')
    reports_dir = os.path.join(script_dir, '..', 'reports')

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    print("[STEP 1/8] Loading Data")
    print("-" * 80)

    loader = DataLoader(data_file)
    df = loader.load_excel_data()
    abundance_cols = loader.identify_abundance_columns()

    set1_markers = ['302-', '310-', '311-', '327-', '338-', '343-', '346-', '364-', '368-', '371-',
                    '375-', '432-', '435-', '436-', '443-', '465-', '467-', '481-', '485-', '486-',
                    '500-', '510-', '533-', '537-', '545-', '561-', '567-', '577-', '584-', '613-']

    group_definitions = {
        'Set1_Old': set1_markers,
        'Set2_Fresh': 'Sarc',
        'Sarcoidosis': 'Sarcoidosis',
        'Control': 'Control'
    }

    groups = loader.identify_sample_groups(group_definitions)

    set1_cols = groups['Set1_Old']
    set2_cols = [col for col in abundance_cols if col not in set1_cols]

    set1_sarc = [col for col in set1_cols if 'Sarcoidosis' in col]
    set1_ctrl = [col for col in set1_cols if 'Control' in col]
    set2_sarc = [col for col in set2_cols if 'Sarcoidosis' in col]
    set2_ctrl = [col for col in set2_cols if 'Control' in col]

    group_info = {
        'Set1_Old': set1_cols,
        'Set2_Fresh': set2_cols,
        'Set1_Sarcoidosis': set1_sarc,
        'Set1_Control': set1_ctrl,
        'Set2_Sarcoidosis': set2_sarc,
        'Set2_Control': set2_ctrl
    }

    print(f"\nSample Distribution:")
    for group_name, group_cols in group_info.items():
        print(f"  {group_name}: {len(group_cols)} samples")

    validation = loader.validate_data_integrity()
    print(f"\nData Validation:")
    print(f"  Total proteins: {validation['total_rows']}")
    print(f"  Total columns: {validation['total_columns']}")
    print(f"  Duplicate rows: {validation['duplicate_rows']}")
    print(f"  Empty rows: {validation['all_nan_rows']}")

    print("\n[STEP 2/8] Sample Quality Control")
    print("-" * 80)

    qc = ProteomicsQC(df, abundance_cols)

    sample_qc = qc.check_sample_completeness()
    print(f"\nSample QC Results:")
    print(f"  Mean missing per sample: {sample_qc['Missing_Percent'].mean():.2f}%")
    print(f"  Median proteins detected: {sample_qc['Detected_Proteins'].median():.0f}")
    print(f"  Samples flagged: {sample_qc['QC_Flag'].sum()}")
    print(f"    - High missing (>50%): {sample_qc['High_Missing'].sum()}")
    print(f"    - Low intensity: {sample_qc['Low_Intensity'].sum()}")
    print(f"    - High intensity: {sample_qc['High_Intensity'].sum()}")

    if sample_qc['QC_Flag'].sum() > 0:
        print("\n  Flagged samples:")
        flagged_samples = sample_qc[sample_qc['QC_Flag']][['Sample', 'Missing_Percent', 'Detected_Proteins']]
        for idx, row in flagged_samples.iterrows():
            print(f"    - {row['Sample'][:50]}: {row['Missing_Percent']:.1f}% missing, {row['Detected_Proteins']} detected")

    sample_qc.to_csv(os.path.join(reports_dir, 'sample_qc_report.csv'), index=False)

    print("\n[STEP 3/8] Protein Quality Control")
    print("-" * 80)

    protein_qc = qc.check_protein_completeness()
    print(f"\nProtein QC Results:")
    print(f"  Mean detection rate: {protein_qc['Detection_Rate'].mean():.2f}%")
    print(f"  Proteins detected in all samples: {(protein_qc['Detection_Rate'] == 100).sum()}")
    print(f"  Proteins detected in >75% samples: {(protein_qc['Detection_Rate'] >= 75).sum()}")
    print(f"  Proteins detected in >50% samples: {(protein_qc['Detection_Rate'] >= 50).sum()}")
    print(f"  Proteins detected in >25% samples: {(protein_qc['Detection_Rate'] >= 25).sum()}")
    print(f"  Median CV: {protein_qc['CV_Percent'].median():.2f}%")

    protein_qc.to_csv(os.path.join(reports_dir, 'protein_qc_report.csv'), index=False)

    high_cv = qc.identify_high_cv_proteins(cv_threshold=100)
    print(f"\n  High CV proteins (>100%): {len(high_cv)}")

    print("\n[STEP 4/8] Sample Correlation Analysis")
    print("-" * 80)

    corr_matrix = qc.calculate_sample_correlations()
    mean_corr = corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].mean()
    print(f"\nSample Correlations:")
    print(f"  Mean pairwise correlation: {mean_corr:.3f}")

    print("\n[STEP 5/8] PCA & Outlier Detection")
    print("-" * 80)

    pca_results, variance_explained = qc.perform_pca(min_detection_rate=0.3)
    print(f"\nPCA Results:")
    print(f"  PC1 variance: {variance_explained[0]*100:.2f}%")
    print(f"  PC2 variance: {variance_explained[1]*100:.2f}%")
    print(f"  Cumulative (PC1+PC2): {(variance_explained[0] + variance_explained[1])*100:.2f}%")
    print(f"  PCA outliers detected: {pca_results['PCA_Outlier'].sum()}")

    if pca_results['PCA_Outlier'].sum() > 0:
        outlier_samples = pca_results[pca_results['PCA_Outlier']]['Sample']
        print(f"\n  Outlier samples:")
        for sample in outlier_samples:
            print(f"    - {sample[:70]}")

    pca_results.to_csv(os.path.join(reports_dir, 'pca_results.csv'), index=False)

    print("\n[STEP 6/8] Normality Testing")
    print("-" * 80)

    normality = qc.check_normality()
    print(f"\nNormality Test Results (Shapiro-Wilk on 10 samples):")
    print(f"  Normally distributed: {normality['Is_Normal'].sum()}/10")
    print(f"  Non-normal: {(~normality['Is_Normal']).sum()}/10")

    normality.to_csv(os.path.join(reports_dir, 'normality_test.csv'), index=False)

    qc_summary = qc.generate_qc_summary()
    print(f"\n[STEP 7/8] Descriptive Statistics")
    print("-" * 80)

    desc_stats = DescriptiveStatistics(df, abundance_cols, group_info)

    group_stats = desc_stats.calculate_group_statistics()
    summary_table = desc_stats.create_summary_table()
    print("\nGroup Summary:")
    print(summary_table.to_string(index=False))

    summary_table.to_csv(os.path.join(reports_dir, 'group_summary.csv'), index=False)

    missing_summary = desc_stats.calculate_missing_value_summary()
    print(f"\nMissing Value Summary:")
    print(f"  Total missing: {missing_summary['missing_percent']:.2f}%")
    print(f"  Proteins with no missing: {missing_summary['proteins_complete']}")
    print(f"  Proteins >50% missing: {missing_summary['proteins_50pct_missing']}")
    print(f"  Proteins >75% missing: {missing_summary['proteins_75pct_missing']}")

    filtering_analysis = desc_stats.analyze_filtering_thresholds()
    print(f"\nFiltering Threshold Analysis:")
    print(f"  At 30% threshold:")
    for group in ['Set1_Sarcoidosis', 'Set1_Control', 'Set2_Sarcoidosis', 'Set2_Control']:
        group_data = filtering_analysis[(filtering_analysis['Threshold'] == 0.3) &
                                       (filtering_analysis['Group'] == group)]
        if len(group_data) > 0:
            print(f"    {group}: {group_data['Proteins_Kept'].values[0]} proteins")

    filtering_analysis.to_csv(os.path.join(reports_dir, 'filtering_threshold_analysis.csv'), index=False)

    percentiles = desc_stats.calculate_abundance_percentiles()
    print(f"\nAbundance Percentiles:")
    for _, row in percentiles.iterrows():
        print(f"  {row['Percentile']}th: {row['Value']:.2e}")

    print("\n[STEP 8/8] Generating Visualizations")
    print("-" * 80)

    visualizer = QCVisualizer(figures_dir)

    abundance_data = df[abundance_cols]

    print("  Creating missing values heatmap...")
    visualizer.plot_missing_values_heatmap(abundance_data)

    print("  Creating sample completeness plots...")
    visualizer.plot_sample_completeness(sample_qc)

    print("  Creating protein detection plots...")
    visualizer.plot_protein_detection_rates(protein_qc)

    print("  Creating PCA plot...")
    group_labels = []
    for sample in abundance_cols:
        if sample in set1_sarc:
            group_labels.append('Set1_Sarc')
        elif sample in set1_ctrl:
            group_labels.append('Set1_Ctrl')
        elif sample in set2_sarc:
            group_labels.append('Set2_Sarc')
        elif sample in set2_ctrl:
            group_labels.append('Set2_Ctrl')
        else:
            group_labels.append('Unknown')

    visualizer.plot_pca(pca_results, variance_explained, group_labels)

    print("  Creating correlation heatmap...")
    visualizer.plot_sample_correlation_heatmap(corr_matrix, group_labels)

    print("  Creating intensity distribution plots...")
    visualizer.plot_intensity_distributions(abundance_data, group_info)

    print("  Creating filtering threshold analysis plot...")
    visualizer.plot_filtering_threshold_analysis(filtering_analysis)

    print("  Creating clustered heatmap...")
    visualizer.plot_clustered_heatmap(abundance_data, n_proteins=100)

    summary_report = pd.DataFrame([qc_summary])
    summary_report.to_csv(os.path.join(reports_dir, 'qc_summary_report.csv'), index=False)

    print("\n" + "="*80)
    print("QC ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nKey Findings:")
    print(f"  Total proteins: {qc_summary['total_proteins']}")
    print(f"  Total samples: {qc_summary['total_samples']}")
    print(f"  Mean detection rate: {qc_summary['mean_detection_rate']:.2f}%")
    print(f"  Samples flagged: {qc_summary['samples_with_flags']}")
    print(f"  PCA outliers: {qc_summary['pca_outliers']}")
    print(f"  PC1 explains: {qc_summary['variance_pc1']*100:.2f}% variance")
    print(f"  Non-normal samples: {qc_summary['samples_non_normal']}/10")
    print(f"\nOutputs saved to:")
    print(f"  Figures: {figures_dir}")
    print(f"  Reports: {reports_dir}")


if __name__ == "__main__":
    run_comprehensive_qc()
