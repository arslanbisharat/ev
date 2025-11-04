import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


class QCVisualizer:

    def __init__(self, output_dir):
        self.output_dir = output_dir
        sns.set_style("whitegrid")
        plt.rcParams['figure.dpi'] = 300

    def plot_missing_values_heatmap(self, abundance_data, filename='missing_values_heatmap.png'):
        fig, ax = plt.subplots(figsize=(16, 10))

        missing_matrix = abundance_data.isna().astype(int)

        sns.heatmap(missing_matrix, cmap=['lightblue', 'darkred'],
                    cbar_kws={'label': 'Missing (1) vs Present (0)'},
                    xticklabels=False, yticklabels=False, ax=ax)

        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('Proteins', fontsize=12)
        ax.set_title('Missing Value Pattern Across All Samples', fontsize=14, fontweight='bold')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_sample_completeness(self, sample_qc, filename='sample_completeness.png'):
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        axes[0, 0].bar(range(len(sample_qc)), sample_qc['Detected_Proteins'], color='steelblue')
        axes[0, 0].set_xlabel('Sample Index')
        axes[0, 0].set_ylabel('Number of Proteins Detected')
        axes[0, 0].set_title('Proteins Detected per Sample')
        axes[0, 0].axhline(sample_qc['Detected_Proteins'].median(), color='red',
                          linestyle='--', label=f'Median: {sample_qc["Detected_Proteins"].median():.0f}')
        axes[0, 0].legend()

        axes[0, 1].bar(range(len(sample_qc)), sample_qc['Missing_Percent'], color='coral')
        axes[0, 1].set_xlabel('Sample Index')
        axes[0, 1].set_ylabel('Missing Data (%)')
        axes[0, 1].set_title('Missing Data per Sample')
        axes[0, 1].axhline(50, color='red', linestyle='--', label='50% threshold')
        axes[0, 1].legend()

        axes[1, 0].hist(sample_qc['Log_Total_Intensity'], bins=30, color='green', edgecolor='black')
        axes[1, 0].set_xlabel('Log10 Total Intensity')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Distribution of Total Intensity per Sample')
        mean_int = sample_qc['Log_Total_Intensity'].mean()
        std_int = sample_qc['Log_Total_Intensity'].std()
        axes[1, 0].axvline(mean_int - 2*std_int, color='red', linestyle='--', label='-2 SD')
        axes[1, 0].axvline(mean_int + 2*std_int, color='red', linestyle='--', label='+2 SD')
        axes[1, 0].legend()

        flagged = sample_qc[sample_qc['QC_Flag']]
        axes[1, 1].scatter(range(len(sample_qc)), sample_qc['Log_Total_Intensity'],
                          c=['red' if flag else 'blue' for flag in sample_qc['QC_Flag']],
                          alpha=0.6, s=50)
        axes[1, 1].set_xlabel('Sample Index')
        axes[1, 1].set_ylabel('Log10 Total Intensity')
        axes[1, 1].set_title(f'Sample QC Flags (Red = Flagged, n={len(flagged)})')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_protein_detection_rates(self, protein_qc, filename='protein_detection_rates.png'):
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        axes[0, 0].hist(protein_qc['Detection_Rate'], bins=50, color='purple', edgecolor='black')
        axes[0, 0].set_xlabel('Detection Rate (%)')
        axes[0, 0].set_ylabel('Number of Proteins')
        axes[0, 0].set_title('Distribution of Protein Detection Rates')
        axes[0, 0].axvline(50, color='red', linestyle='--', label='50% threshold')
        axes[0, 0].legend()

        detection_bins = pd.cut(protein_qc['Detection_Rate'],
                               bins=[0, 25, 50, 75, 100],
                               labels=['0-25%', '25-50%', '50-75%', '75-100%'])
        bin_counts = detection_bins.value_counts().sort_index()
        axes[0, 1].bar(range(len(bin_counts)), bin_counts.values, color='teal')
        axes[0, 1].set_xticks(range(len(bin_counts)))
        axes[0, 1].set_xticklabels(bin_counts.index)
        axes[0, 1].set_ylabel('Number of Proteins')
        axes[0, 1].set_title('Proteins by Detection Rate Categories')

        cv_data = protein_qc['CV_Percent'].dropna()
        axes[1, 0].hist(cv_data, bins=50, color='orange', edgecolor='black', range=(0, 150))
        axes[1, 0].set_xlabel('Coefficient of Variation (%)')
        axes[1, 0].set_ylabel('Number of Proteins')
        axes[1, 0].set_title('Distribution of Protein CV')
        axes[1, 0].axvline(cv_data.median(), color='red', linestyle='--',
                          label=f'Median: {cv_data.median():.1f}%')
        axes[1, 0].legend()

        mean_abundance = protein_qc['Mean_Abundance'].dropna()
        axes[1, 1].hist(np.log10(mean_abundance + 1), bins=50, color='brown', edgecolor='black')
        axes[1, 1].set_xlabel('Log10 Mean Abundance')
        axes[1, 1].set_ylabel('Number of Proteins')
        axes[1, 1].set_title('Distribution of Mean Protein Abundance')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_pca(self, pca_results, variance_explained, group_labels=None, filename='pca_plot.png'):
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        if group_labels is None:
            colors = ['blue'] * len(pca_results)
        else:
            unique_groups = list(set(group_labels))
            color_map = dict(zip(unique_groups, plt.cm.tab10(range(len(unique_groups)))))
            colors = [color_map[g] for g in group_labels]

        scatter = axes[0].scatter(pca_results['PC1'], pca_results['PC2'],
                                  c=colors, alpha=0.6, s=80, edgecolors='black')
        axes[0].set_xlabel(f'PC1 ({variance_explained[0]*100:.1f}% variance)')
        axes[0].set_ylabel(f'PC2 ({variance_explained[1]*100:.1f}% variance)')
        axes[0].set_title('PCA of Samples')
        axes[0].grid(alpha=0.3)

        outliers = pca_results[pca_results['PCA_Outlier']]
        if len(outliers) > 0:
            axes[0].scatter(outliers['PC1'], outliers['PC2'],
                          s=200, facecolors='none', edgecolors='red',
                          linewidths=2, label=f'Outliers (n={len(outliers)})')
            axes[0].legend()

        var_cumsum = np.cumsum(variance_explained[:10])
        axes[1].bar(range(1, len(var_cumsum)+1), variance_explained[:10]*100,
                   color='steelblue', edgecolor='black')
        axes[1].plot(range(1, len(var_cumsum)+1), var_cumsum*100,
                    marker='o', color='red', linewidth=2, label='Cumulative')
        axes[1].set_xlabel('Principal Component')
        axes[1].set_ylabel('Variance Explained (%)')
        axes[1].set_title('Scree Plot')
        axes[1].legend()
        axes[1].grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_sample_correlation_heatmap(self, corr_matrix, group_labels=None, filename='correlation_heatmap.png'):
        fig, ax = plt.subplots(figsize=(14, 12))

        if group_labels is not None:
            unique_groups = sorted(set(group_labels))
            colors = sns.color_palette("tab10", len(unique_groups))
            group_colors = [colors[unique_groups.index(g)] for g in group_labels]

            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=colors[i], label=g)
                             for i, g in enumerate(unique_groups)]
        else:
            group_colors = None
            legend_elements = None

        sns.heatmap(corr_matrix, cmap='RdYlBu_r', center=0.5,
                   vmin=0, vmax=1, square=True,
                   xticklabels=False, yticklabels=False,
                   cbar_kws={'label': 'Pearson Correlation'}, ax=ax)

        ax.set_title('Sample-to-Sample Correlation Matrix', fontsize=14, fontweight='bold')

        if legend_elements:
            ax.legend(handles=legend_elements, bbox_to_anchor=(1.15, 1), loc='upper left')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_intensity_distributions(self, abundance_data, group_info, filename='intensity_distributions.png'):
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        all_values = abundance_data.values.flatten()
        all_values = all_values[~np.isnan(all_values)]

        axes[0, 0].hist(all_values, bins=100, color='steelblue', edgecolor='black', alpha=0.7)
        axes[0, 0].set_xlabel('Abundance')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Raw Abundance Distribution (All Samples)')
        axes[0, 0].set_yscale('log')

        axes[0, 1].hist(np.log2(all_values + 1), bins=100, color='green', edgecolor='black', alpha=0.7)
        axes[0, 1].set_xlabel('Log2 Abundance')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Log2 Abundance Distribution (All Samples)')

        boxplot_data = []
        labels = []
        for sample in abundance_data.columns:
            sample_values = abundance_data[sample].dropna()
            if len(sample_values) > 0:
                boxplot_data.append(sample_values.values)
                labels.append(sample)

        bp = axes[1, 0].boxplot(boxplot_data, patch_artist=True, showfliers=False)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        axes[1, 0].set_xlabel('Sample Index')
        axes[1, 0].set_ylabel('Abundance')
        axes[1, 0].set_title('Abundance Distribution per Sample (Boxplot)')
        axes[1, 0].set_xticklabels([])

        if group_info:
            for group_name, group_cols in group_info.items():
                group_data = abundance_data[group_cols].values.flatten()
                group_data = group_data[~np.isnan(group_data)]
                axes[1, 1].hist(np.log2(group_data + 1), bins=50, alpha=0.5, label=group_name, edgecolor='black')

            axes[1, 1].set_xlabel('Log2 Abundance')
            axes[1, 1].set_ylabel('Frequency')
            axes[1, 1].set_title('Log2 Abundance by Group')
            axes[1, 1].legend()

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_filtering_threshold_analysis(self, filtering_df, filename='filtering_threshold_analysis.png'):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        for group in filtering_df['Group'].unique():
            group_data = filtering_df[filtering_df['Group'] == group]
            axes[0].plot(group_data['Threshold'], group_data['Proteins_Kept'],
                        marker='o', linewidth=2, label=group)

        axes[0].set_xlabel('Detection Threshold')
        axes[0].set_ylabel('Number of Proteins Kept')
        axes[0].set_title('Proteins Retained at Different Thresholds')
        axes[0].legend()
        axes[0].grid(alpha=0.3)

        for group in filtering_df['Group'].unique():
            group_data = filtering_df[filtering_df['Group'] == group]
            axes[1].plot(group_data['Threshold'], group_data['Percent_Kept'],
                        marker='o', linewidth=2, label=group)

        axes[1].set_xlabel('Detection Threshold')
        axes[1].set_ylabel('Percent of Proteins Kept (%)')
        axes[1].set_title('Percentage of Proteins Retained')
        axes[1].legend()
        axes[1].grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_clustered_heatmap(self, abundance_data, n_proteins=100, filename='clustered_heatmap.png'):
        variance_per_protein = abundance_data.var(axis=1)
        top_variable_indices = variance_per_protein.nlargest(n_proteins).index
        plot_data = abundance_data.loc[top_variable_indices]

        plot_data_filled = plot_data.T.fillna(plot_data.min().min())

        fig, ax = plt.subplots(figsize=(12, 10))

        sns.heatmap(plot_data_filled.T, cmap='RdBu_r', center=plot_data_filled.values.mean(),
                   xticklabels=False, yticklabels=False,
                   cbar_kws={'label': 'Abundance'}, ax=ax)

        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('Proteins', fontsize=12)
        ax.set_title(f'Heatmap of Top {n_proteins} Most Variable Proteins', fontsize=14, fontweight='bold')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()
