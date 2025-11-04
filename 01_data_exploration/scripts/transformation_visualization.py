import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats


class TransformationVisualizer:

    def __init__(self, output_dir):
        self.output_dir = output_dir
        sns.set_style("whitegrid")
        plt.rcParams['figure.dpi'] = 300

    def plot_transformation_comparison(self, raw_data, transformed_data, filename='transformation_comparison.png'):
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))

        raw_flat = raw_data.values.flatten()
        raw_flat = raw_flat[~np.isnan(raw_flat)]

        trans_flat = transformed_data.values.flatten()
        trans_flat = trans_flat[~np.isnan(trans_flat)]

        axes[0, 0].hist(raw_flat, bins=100, color='steelblue', edgecolor='black', alpha=0.7)
        axes[0, 0].set_xlabel('Raw Abundance')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Raw Data Distribution')
        axes[0, 0].set_yscale('log')

        axes[1, 0].hist(trans_flat, bins=100, color='green', edgecolor='black', alpha=0.7)
        axes[1, 0].set_xlabel('Log2 Abundance')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Log2 Transformed Distribution')

        stats.probplot(raw_flat[::100], dist="norm", plot=axes[0, 1])
        axes[0, 1].set_title('Q-Q Plot: Raw Data')
        axes[0, 1].get_lines()[0].set_markerfacecolor('steelblue')
        axes[0, 1].get_lines()[0].set_markersize(3)

        stats.probplot(trans_flat[::100], dist="norm", plot=axes[1, 1])
        axes[1, 1].set_title('Q-Q Plot: Log2 Transformed')
        axes[1, 1].get_lines()[0].set_markerfacecolor('green')
        axes[1, 1].get_lines()[0].set_markersize(3)

        boxplot_data_raw = [raw_data[col].dropna().values for col in raw_data.columns[:20]]
        bp_raw = axes[0, 2].boxplot(boxplot_data_raw, patch_artist=True, showfliers=False)
        for patch in bp_raw['boxes']:
            patch.set_facecolor('lightblue')
        axes[0, 2].set_xlabel('Sample (first 20)')
        axes[0, 2].set_ylabel('Abundance')
        axes[0, 2].set_title('Raw Data by Sample')
        axes[0, 2].set_xticklabels([])

        boxplot_data_trans = [transformed_data[col].dropna().values for col in transformed_data.columns[:20]]
        bp_trans = axes[1, 2].boxplot(boxplot_data_trans, patch_artist=True, showfliers=False)
        for patch in bp_trans['boxes']:
            patch.set_facecolor('lightgreen')
        axes[1, 2].set_xlabel('Sample (first 20)')
        axes[1, 2].set_ylabel('Log2 Abundance')
        axes[1, 2].set_title('Log2 Transformed by Sample')
        axes[1, 2].set_xticklabels([])

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_sample_medians_comparison(self, raw_data, transformed_data, filename='sample_medians_comparison.png'):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        raw_medians = raw_data.median().values
        trans_medians = transformed_data.median().values

        axes[0].bar(range(len(raw_medians)), raw_medians, color='steelblue', alpha=0.7)
        axes[0].axhline(np.median(raw_medians), color='red', linestyle='--',
                       label=f'Median: {np.median(raw_medians):.2e}')
        axes[0].set_xlabel('Sample Index')
        axes[0].set_ylabel('Median Abundance')
        axes[0].set_title('Raw Data: Sample Medians')
        axes[0].legend()

        axes[1].bar(range(len(trans_medians)), trans_medians, color='green', alpha=0.7)
        axes[1].axhline(np.median(trans_medians), color='red', linestyle='--',
                       label=f'Median: {np.median(trans_medians):.2f}')
        axes[1].set_xlabel('Sample Index')
        axes[1].set_ylabel('Median Log2 Abundance')
        axes[1].set_title('Log2 Transformed: Sample Medians')
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_cv_comparison(self, raw_data, transformed_data, filename='cv_comparison.png'):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        cv_raw = []
        cv_transformed = []

        for i in range(len(raw_data)):
            raw_values = raw_data.iloc[i].dropna()
            trans_values = transformed_data.iloc[i].dropna()

            if len(raw_values) > 0 and raw_values.mean() > 0:
                cv_raw.append((raw_values.std() / raw_values.mean()) * 100)

            if len(trans_values) > 0 and trans_values.mean() != 0:
                cv_transformed.append((trans_values.std() / trans_values.mean()) * 100)

        cv_raw = np.array(cv_raw)
        cv_transformed = np.array(cv_transformed)

        cv_raw_clean = cv_raw[~np.isnan(cv_raw)]
        cv_trans_clean = cv_transformed[~np.isnan(cv_transformed)]

        axes[0].hist(cv_raw_clean, bins=50, color='steelblue', edgecolor='black', alpha=0.7, range=(0, 200))
        axes[0].axvline(np.median(cv_raw_clean), color='red', linestyle='--',
                       label=f'Median: {np.median(cv_raw_clean):.1f}%')
        axes[0].set_xlabel('Coefficient of Variation (%)')
        axes[0].set_ylabel('Number of Proteins')
        axes[0].set_title('CV Distribution: Raw Data')
        axes[0].legend()

        axes[1].hist(cv_trans_clean, bins=50, color='green', edgecolor='black', alpha=0.7, range=(0, 200))
        axes[1].axvline(np.median(cv_trans_clean), color='red', linestyle='--',
                       label=f'Median: {np.median(cv_trans_clean):.1f}%')
        axes[1].set_xlabel('Coefficient of Variation (%)')
        axes[1].set_ylabel('Number of Proteins')
        axes[1].set_title('CV Distribution: Log2 Transformed')
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_variance_stabilization(self, raw_data, transformed_data, filename='variance_stabilization.png'):
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        raw_means = raw_data.mean(axis=1).values
        raw_vars = raw_data.var(axis=1).values

        trans_means = transformed_data.mean(axis=1).values
        trans_vars = transformed_data.var(axis=1).values

        raw_means_clean = raw_means[~np.isnan(raw_means) & ~np.isnan(raw_vars)]
        raw_vars_clean = raw_vars[~np.isnan(raw_means) & ~np.isnan(raw_vars)]

        trans_means_clean = trans_means[~np.isnan(trans_means) & ~np.isnan(trans_vars)]
        trans_vars_clean = trans_vars[~np.isnan(trans_means) & ~np.isnan(trans_vars)]

        axes[0].scatter(np.log10(raw_means_clean + 1), np.log10(raw_vars_clean + 1),
                       alpha=0.3, s=10, color='steelblue')
        axes[0].set_xlabel('Log10 Mean Abundance')
        axes[0].set_ylabel('Log10 Variance')
        axes[0].set_title('Mean-Variance Relationship: Raw Data')
        axes[0].grid(alpha=0.3)

        from scipy.stats import pearsonr
        if len(raw_means_clean) > 0:
            corr_raw, _ = pearsonr(np.log10(raw_means_clean + 1), np.log10(raw_vars_clean + 1))
            axes[0].text(0.05, 0.95, f'Correlation: {corr_raw:.3f}',
                        transform=axes[0].transAxes, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        axes[1].scatter(trans_means_clean, trans_vars_clean,
                       alpha=0.3, s=10, color='green')
        axes[1].set_xlabel('Mean Log2 Abundance')
        axes[1].set_ylabel('Variance')
        axes[1].set_title('Mean-Variance Relationship: Log2 Transformed')
        axes[1].grid(alpha=0.3)

        if len(trans_means_clean) > 0:
            corr_trans, _ = pearsonr(trans_means_clean, trans_vars_clean)
            axes[1].text(0.05, 0.95, f'Correlation: {corr_trans:.3f}',
                        transform=axes[1].transAxes, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_normality_test_comparison(self, normality_results, filename='normality_comparison.png'):
        fig, ax = plt.subplots(figsize=(10, 6))

        x = np.arange(len(normality_results))
        width = 0.35

        ax.bar(x - width/2, -np.log10(normality_results['Raw_Shapiro_Pval']),
               width, label='Raw Data', color='steelblue', alpha=0.7)
        ax.bar(x + width/2, -np.log10(normality_results['Transformed_Shapiro_Pval']),
               width, label='Log2 Transformed', color='green', alpha=0.7)

        ax.axhline(-np.log10(0.05), color='red', linestyle='--', linewidth=2,
                   label='p = 0.05 threshold')

        ax.set_xlabel('Sample')
        ax.set_ylabel('-Log10(p-value)')
        ax.set_title('Shapiro-Wilk Normality Test: Before vs After Transformation')
        ax.set_xticks(x)
        ax.set_xticklabels([f'S{i+1}' for i in range(len(normality_results))])
        ax.legend()
        ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_skewness_kurtosis(self, distribution_metrics, filename='skewness_kurtosis.png'):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        categories = ['Raw Data', 'Log2 Transformed']
        skewness = [abs(distribution_metrics['raw']['skewness']),
                   abs(distribution_metrics['transformed']['skewness'])]
        kurtosis = [abs(distribution_metrics['raw']['kurtosis']),
                   abs(distribution_metrics['transformed']['kurtosis'])]

        axes[0].bar(categories, skewness, color=['steelblue', 'green'], alpha=0.7)
        axes[0].set_ylabel('Absolute Skewness')
        axes[0].set_title('Distribution Skewness')
        axes[0].grid(alpha=0.3)

        for i, v in enumerate(skewness):
            axes[0].text(i, v + 0.1, f'{v:.2f}', ha='center', fontweight='bold')

        axes[1].bar(categories, kurtosis, color=['steelblue', 'green'], alpha=0.7)
        axes[1].set_ylabel('Absolute Kurtosis')
        axes[1].set_title('Distribution Kurtosis')
        axes[1].grid(alpha=0.3)

        for i, v in enumerate(kurtosis):
            axes[1].text(i, v + 0.1, f'{v:.2f}', ha='center', fontweight='bold')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()
