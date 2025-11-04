import os
import sys
import pandas as pd
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, script_dir)

from data_loader import DataLoader
from transformation import DataTransformation
from transformation_visualization import TransformationVisualizer


def run_log2_transformation():

    print("="*80)
    print("LOG2 TRANSFORMATION & NORMALIZATION VERIFICATION")
    print("="*80)
    print()

    data_file = os.path.join(project_dir, 'proteome proteins.xlsx')
    figures_dir = os.path.join(script_dir, '..', 'figures')
    reports_dir = os.path.join(script_dir, '..', 'reports')
    processed_dir = os.path.join(script_dir, '..', 'processed_data')

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(processed_dir, exist_ok=True)

    print("[STEP 1/6] Loading Data")
    print("-" * 80)

    loader = DataLoader(data_file)
    df = loader.load_excel_data()
    abundance_cols = loader.identify_abundance_columns()

    abundance_data = df[abundance_cols]

    print(f"Data shape: {abundance_data.shape}")
    print(f"Total values: {abundance_data.size}")
    print(f"Missing values: {abundance_data.isna().sum().sum()} ({abundance_data.isna().sum().sum() / abundance_data.size * 100:.2f}%)")
    print(f"Zero values: {(abundance_data == 0).sum().sum()}")

    print("\n[STEP 2/6] Applying Log2 Transformation")
    print("-" * 80)

    transformer = DataTransformation(abundance_data)

    transformed_data = transformer.log2_transform(pseudocount=1)

    print(f"\nTransformation complete:")
    print(f"  Pseudocount used: 1")
    print(f"  Zeros replaced with NaN: {transformer.transformation_stats['n_zeros_replaced']}")
    print(f"  Negative values found: {transformer.transformation_stats['n_negative_values']}")
    print(f"  Missing in transformed: {transformed_data.isna().sum().sum()}")

    print("\n[STEP 3/6] Calculating Distribution Metrics")
    print("-" * 80)

    dist_metrics = transformer.calculate_distribution_metrics()

    print("\nRaw Data Metrics:")
    print(f"  Mean: {dist_metrics['raw']['mean']:.2e}")
    print(f"  Median: {dist_metrics['raw']['median']:.2e}")
    print(f"  SD: {dist_metrics['raw']['std']:.2e}")
    print(f"  Skewness: {dist_metrics['raw']['skewness']:.3f}")
    print(f"  Kurtosis: {dist_metrics['raw']['kurtosis']:.3f}")

    print("\nLog2 Transformed Metrics:")
    print(f"  Mean: {dist_metrics['transformed']['mean']:.2f}")
    print(f"  Median: {dist_metrics['transformed']['median']:.2f}")
    print(f"  SD: {dist_metrics['transformed']['std']:.2f}")
    print(f"  Skewness: {dist_metrics['transformed']['skewness']:.3f}")
    print(f"  Kurtosis: {dist_metrics['transformed']['kurtosis']:.3f}")

    print(f"\nSkewness reduction: {abs(dist_metrics['raw']['skewness']) - abs(dist_metrics['transformed']['skewness']):.3f}")
    print(f"Kurtosis reduction: {abs(dist_metrics['raw']['kurtosis']) - abs(dist_metrics['transformed']['kurtosis']):.3f}")

    print("\n[STEP 4/6] Verifying Normality Improvement")
    print("-" * 80)

    normality_results = transformer.verify_normality_improvement(n_samples=10)

    print(f"\nShapiro-Wilk test results (10 samples):")
    print(f"  Normal before transformation: {normality_results['Raw_Is_Normal'].sum()}/10")
    print(f"  Normal after transformation: {normality_results['Transformed_Is_Normal'].sum()}/10")
    print(f"  Samples showing improvement: {(normality_results['Improvement'] == 'Yes').sum()}/10")

    normality_results.to_csv(os.path.join(reports_dir, 'transformation_normality_tests.csv'), index=False)

    print("\n[STEP 5/6] Checking Sample Balance")
    print("-" * 80)

    balance_stats = transformer.check_sample_distribution_balance()

    print("\nSample Median Statistics:")
    print(f"  Raw data:")
    print(f"    Median of sample medians: {balance_stats['raw']['median_of_medians']:.2e}")
    print(f"    SD of sample medians: {balance_stats['raw']['std_of_medians']:.2e}")
    print(f"    CV of sample medians: {balance_stats['raw']['cv_of_medians']:.2f}%")

    print(f"\n  Transformed data:")
    print(f"    Median of sample medians: {balance_stats['transformed']['median_of_medians']:.2f}")
    print(f"    SD of sample medians: {balance_stats['transformed']['std_of_medians']:.2f}")
    print(f"    CV of sample medians: {balance_stats['transformed']['cv_of_medians']:.2f}%")

    cv_comparison = transformer.compare_cv_before_after()

    print("\nCoefficient of Variation Analysis:")
    print(f"  Raw data median CV: {cv_comparison['raw_median_cv']:.2f}%")
    print(f"  Transformed data median CV: {cv_comparison['transformed_median_cv']:.2f}%")
    print(f"  CV reduction: {cv_comparison['cv_reduction']:.2f}%")

    outlier_stats = transformer.detect_transformation_outliers()

    print("\nOutlier Detection (IQR method):")
    print(f"  Raw data - mean outliers per sample: {outlier_stats['raw_mean_outliers_per_sample']:.1f}")
    print(f"  Transformed data - mean outliers per sample: {outlier_stats['transformed_mean_outliers_per_sample']:.1f}")

    summary = transformer.get_transformation_summary()
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(os.path.join(reports_dir, 'transformation_summary.csv'), index=False)

    print("\n[STEP 6/6] Generating Visualizations")
    print("-" * 80)

    visualizer = TransformationVisualizer(figures_dir)

    print("  Creating transformation comparison plots...")
    visualizer.plot_transformation_comparison(abundance_data, transformed_data)

    print("  Creating sample medians comparison...")
    visualizer.plot_sample_medians_comparison(abundance_data, transformed_data)

    print("  Creating CV comparison...")
    visualizer.plot_cv_comparison(abundance_data, transformed_data)

    print("  Creating variance stabilization plots...")
    visualizer.plot_variance_stabilization(abundance_data, transformed_data)

    print("  Creating normality test comparison...")
    visualizer.plot_normality_test_comparison(normality_results)

    print("  Creating skewness/kurtosis comparison...")
    visualizer.plot_skewness_kurtosis(dist_metrics)

    print("\nSaving transformed data...")
    transformed_with_info = pd.concat([df[['Accession', 'Gene Symbol', 'Description']], transformed_data], axis=1)
    transformed_with_info.to_csv(os.path.join(processed_dir, 'log2_transformed_data.csv'), index=False)

    abundance_data_with_info = pd.concat([df[['Accession', 'Gene Symbol', 'Description']], abundance_data], axis=1)
    abundance_data_with_info.to_csv(os.path.join(processed_dir, 'raw_abundance_data.csv'), index=False)

    print("\n" + "="*80)
    print("LOG2 TRANSFORMATION COMPLETE")
    print("="*80)
    print(f"\nKey Results:")
    print(f"  Transformation method: Log2 (pseudocount = 1)")
    print(f"  Skewness reduced by: {summary.get('skewness_reduction', 0):.3f}")
    print(f"  Kurtosis reduced by: {summary.get('kurtosis_reduction', 0):.3f}")
    print(f"  CV median before: {summary.get('median_cv_before', 0):.2f}%")
    print(f"  CV median after: {summary.get('median_cv_after', 0):.2f}%")
    print(f"  Normal samples before: {summary.get('normal_samples_before', 0)}/10")
    print(f"  Normal samples after: {summary.get('normal_samples_after', 0)}/10")

    print(f"\nOutputs saved to:")
    print(f"  Figures: {figures_dir}")
    print(f"  Reports: {reports_dir}")
    print(f"  Processed data: {processed_dir}")
    print(f"\nTransformed data file: log2_transformed_data.csv")
    print(f"Raw data file: raw_abundance_data.csv")


if __name__ == "__main__":
    run_log2_transformation()
