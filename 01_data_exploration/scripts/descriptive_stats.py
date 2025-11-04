import pandas as pd
import numpy as np
from scipy import stats


class DescriptiveStatistics:

    def __init__(self, df, abundance_cols, group_info):
        self.df = df
        self.abundance_cols = abundance_cols
        self.group_info = group_info
        self.stats_results = {}

    def calculate_group_statistics(self):
        abundance_data = self.df[self.abundance_cols]

        group_stats = {}

        for group_name, group_cols in self.group_info.items():
            group_data = abundance_data[group_cols]

            stats_dict = {
                'n_samples': len(group_cols),
                'mean_per_protein': group_data.mean(axis=1),
                'median_per_protein': group_data.median(axis=1),
                'std_per_protein': group_data.std(axis=1),
                'min_per_protein': group_data.min(axis=1),
                'max_per_protein': group_data.max(axis=1),
                'missing_per_protein': group_data.isna().sum(axis=1),
                'detection_rate_per_protein': (~group_data.isna()).sum(axis=1) / len(group_cols) * 100
            }

            group_stats[group_name] = stats_dict

        self.stats_results['group_statistics'] = group_stats

        return group_stats

    def create_summary_table(self):
        if 'group_statistics' not in self.stats_results:
            self.calculate_group_statistics()

        summary_rows = []

        for group_name, stats_dict in self.stats_results['group_statistics'].items():
            summary_rows.append({
                'Group': group_name,
                'N_Samples': stats_dict['n_samples'],
                'Mean_Proteins_Detected': stats_dict['detection_rate_per_protein'].mean(),
                'Median_Proteins_Detected': stats_dict['detection_rate_per_protein'].median(),
                'Mean_Abundance_Overall': stats_dict['mean_per_protein'].mean(),
                'Median_Abundance_Overall': stats_dict['median_per_protein'].median(),
                'SD_Abundance_Overall': stats_dict['std_per_protein'].mean()
            })

        summary_table = pd.DataFrame(summary_rows)

        return summary_table

    def calculate_missing_value_summary(self):
        abundance_data = self.df[self.abundance_cols]

        total_values = len(self.df) * len(self.abundance_cols)
        missing_values = abundance_data.isna().sum().sum()
        missing_percent = (missing_values / total_values) * 100

        missing_per_sample = abundance_data.isna().sum() / len(abundance_data) * 100
        missing_per_protein = abundance_data.isna().sum(axis=1) / len(self.abundance_cols) * 100

        missing_summary = {
            'total_values': total_values,
            'missing_values': missing_values,
            'missing_percent': missing_percent,
            'mean_missing_per_sample': missing_per_sample.mean(),
            'median_missing_per_sample': missing_per_sample.median(),
            'max_missing_per_sample': missing_per_sample.max(),
            'mean_missing_per_protein': missing_per_protein.mean(),
            'median_missing_per_protein': missing_per_protein.median(),
            'max_missing_per_protein': missing_per_protein.max(),
            'proteins_complete': (missing_per_protein == 0).sum(),
            'proteins_50pct_missing': (missing_per_protein >= 50).sum(),
            'proteins_75pct_missing': (missing_per_protein >= 75).sum()
        }

        self.stats_results['missing_summary'] = missing_summary

        return missing_summary

    def calculate_detection_rate_distribution(self):
        abundance_data = self.df[self.abundance_cols]

        detection_per_protein = (~abundance_data.isna()).sum(axis=1) / len(self.abundance_cols) * 100

        detection_bins = pd.cut(detection_per_protein, bins=[0, 10, 25, 50, 75, 90, 100],
                                labels=['0-10%', '10-25%', '25-50%', '50-75%', '75-90%', '90-100%'])

        detection_distribution = detection_bins.value_counts().sort_index()

        return detection_distribution

    def analyze_filtering_thresholds(self, thresholds=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]):
        abundance_data = self.df[self.abundance_cols]

        group_stats = self.stats_results.get('group_statistics', {})
        if not group_stats:
            self.calculate_group_statistics()
            group_stats = self.stats_results['group_statistics']

        filtering_results = []

        for threshold in thresholds:
            proteins_kept_overall = 0

            for group_name, stats_dict in group_stats.items():
                detection_rates = stats_dict['detection_rate_per_protein']
                proteins_above_threshold = (detection_rates >= threshold * 100).sum()

                filtering_results.append({
                    'Threshold': threshold,
                    'Group': group_name,
                    'Proteins_Kept': proteins_above_threshold,
                    'Percent_Kept': (proteins_above_threshold / len(detection_rates)) * 100
                })

        filtering_df = pd.DataFrame(filtering_results)

        self.stats_results['filtering_analysis'] = filtering_df

        return filtering_df

    def compare_groups(self, group1_name, group2_name):
        if 'group_statistics' not in self.stats_results:
            self.calculate_group_statistics()

        group_stats = self.stats_results['group_statistics']

        group1_stats = group_stats[group1_name]
        group2_stats = group_stats[group2_name]

        comparison = pd.DataFrame({
            'Mean_Group1': group1_stats['mean_per_protein'],
            'Mean_Group2': group2_stats['mean_per_protein'],
            'Median_Group1': group1_stats['median_per_protein'],
            'Median_Group2': group2_stats['median_per_protein'],
            'SD_Group1': group1_stats['std_per_protein'],
            'SD_Group2': group2_stats['std_per_protein'],
            'Detection_Group1': group1_stats['detection_rate_per_protein'],
            'Detection_Group2': group2_stats['detection_rate_per_protein']
        })

        if 'Gene Symbol' in self.df.columns:
            comparison['Gene'] = self.df['Gene Symbol'].values
        if 'Accession' in self.df.columns:
            comparison['Accession'] = self.df['Accession'].values

        return comparison

    def calculate_abundance_percentiles(self, percentiles=[5, 10, 25, 50, 75, 90, 95]):
        abundance_data = self.df[self.abundance_cols]

        all_values = abundance_data.values.flatten()
        all_values = all_values[~np.isnan(all_values)]

        percentile_values = np.percentile(all_values, percentiles)

        percentile_df = pd.DataFrame({
            'Percentile': percentiles,
            'Value': percentile_values
        })

        return percentile_df

    def identify_consistently_detected_proteins(self, min_detection=0.8):
        if 'group_statistics' not in self.stats_results:
            self.calculate_group_statistics()

        group_stats = self.stats_results['group_statistics']

        consistent_proteins = None

        for group_name, stats_dict in group_stats.items():
            detection_rates = stats_dict['detection_rate_per_protein']
            group_consistent = detection_rates >= (min_detection * 100)

            if consistent_proteins is None:
                consistent_proteins = group_consistent
            else:
                consistent_proteins = consistent_proteins & group_consistent

        consistent_protein_indices = np.where(consistent_proteins)[0]

        consistent_df = self.df.iloc[consistent_protein_indices][['Gene Symbol', 'Accession', 'Description']].copy()
        consistent_df['N_Proteins'] = len(consistent_protein_indices)

        return consistent_df
