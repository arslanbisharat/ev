import os
import sys
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
sys.path.insert(0, project_dir)

from scripts.data_loader import ProteomicsDataLoader
from scripts.preprocessing import ProteomicsPreprocessor
from scripts.differential_analysis import DifferentialAnalyzer
from scripts.visualization import ProteomicsVisualizer

def run_analysis(data_path):

    print("="*80)
    print("EV PROTEOMICS STORAGE EFFECT ANALYSIS")
    print("="*80)

    print("[1/6] Loading data")
    loader = ProteomicsDataLoader(data_path)
    df = loader.load_data()
    loader.identify_sample_columns()
    group_cols = loader.get_group_columns()
    summary = loader.get_summary()

    print(f"  Proteins: {summary['total_proteins']}")
    print(f"  Set 1 (Old): {summary['set1_samples']} samples")
    print(f"  Set 2 (Fresh): {summary['set2_samples']} samples")

    print("[2/6] Preprocessing")
    preprocessor = ProteomicsPreprocessor(df, loader.abundance_cols)
    df_filt, n_kept = preprocessor.filter_proteins(group_cols)
    print(f"  Kept: {n_kept} / {summary['total_proteins']} proteins")

    data_log2 = preprocessor.log2_transform()
    data_processed = preprocessor.impute_missing()

    print("[3/6] Differential analysis")
    analyzer = DifferentialAnalyzer(data_processed, df_filt)

    results_set1 = analyzer.perform_ttest(group_cols['set1_sarc'], group_cols['set1_ctrl'])
    results_set2 = analyzer.perform_ttest(group_cols['set2_sarc'], group_cols['set2_ctrl'])

    sig1 = analyzer.get_significant_proteins(results_set1)
    sig2 = analyzer.get_significant_proteins(results_set2)

    print(f"  Set 1: {len(sig1)} significant proteins")
    print(f"  Set 2: {len(sig2)} significant proteins")

    print("[4/6] Overlap analysis")
    overlap_data = analyzer.calculate_overlap(results_set1, results_set2)
    print(f"  Overlap: {overlap_data['n_overlap']} proteins")
    print(f"  Jaccard: {overlap_data['jaccard']:.3f}")

    corr_data = analyzer.calculate_correlation(results_set1, results_set2)
    print(f"  Pearson r: {corr_data['pearson']:.3f}")

    print("[5/6] Creating visualizations")
    visualizer = ProteomicsVisualizer()

    abundance_before = df_filt[loader.abundance_cols].replace(0, np.nan)
    visualizer.plot_log_transformation(abundance_before, data_log2)

    visualizer.plot_volcano(results_set1, 'Set 1 (Old Samples)', 'volcano_set1.png')
    visualizer.plot_volcano(results_set2, 'Set 2 (Fresh Samples)', 'volcano_set2.png')

    visualizer.plot_venn_diagram(overlap_data, len(sig1), len(sig2))
    visualizer.plot_correlation(corr_data['merged_data'], corr_data['pearson'], overlap_data)

    print("[6/6] Saving results")
    results_dir = os.path.join(project_dir, 'results')
    results_set1.to_csv(os.path.join(results_dir, 'results_set1_old.csv'), index=False)
    results_set2.to_csv(os.path.join(results_dir, 'results_set2_fresh.csv'), index=False)

    storage_stable = sig1[sig1['Gene'].isin(overlap_data['overlap'])]
    storage_stable.to_csv(os.path.join(results_dir, 'storage_stable_proteins.csv'), index=False)

    print("" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"Key Findings:")
    print(f"  - {len(sig2)/len(sig1):.1f}x more proteins in fresh samples")
    print(f"  - Correlation: {corr_data['pearson']:.3f}")
    print(f"  - {overlap_data['n_overlap']} storage-stable biomarkers identified")

if __name__ == "__main__":
    data_file = os.path.join(project_dir, 'proteome proteins.xlsx')
    run_analysis(data_file)
