import pandas as pd
import numpy as np
from scipy import stats


class DataTransformation:

    def __init__(self, abundance_data):
        self.raw_data = abundance_data.copy()
        self.transformed_data = None
        self.transformation_stats = {}

    def log2_transform(self, pseudocount=1):
        transformed = self.raw_data.copy()

        transformed = transformed.replace(0, np.nan)

        transformed = np.log2(transformed + pseudocount)

        self.transformed_data = transformed

        self.transformation_stats['method'] = 'log2'
        self.transformation_stats['pseudocount'] = pseudocount
        self.transformation_stats['n_zeros_replaced'] = (self.raw_data == 0).sum().sum()
        self.transformation_stats['n_negative_values'] = (self.raw_data < 0).sum().sum()

        return transformed

    def verify_normality_improvement(self, n_samples=10):
        if self.transformed_data is None:
            raise ValueError("No transformation applied yet")

        results = []

        sample_cols = self.raw_data.columns[:n_samples]

        for col in sample_cols:
            raw_values = self.raw_data[col].dropna()
            transformed_values = self.transformed_data[col].dropna()

            if len(raw_values) >= 3:
                raw_stat, raw_pval = stats.shapiro(raw_values)
                transformed_stat, transformed_pval = stats.shapiro(transformed_values)

                results.append({
                    'Sample': col,
                    'Raw_Shapiro_Stat': raw_stat,
                    'Raw_Shapiro_Pval': raw_pval,
                    'Raw_Is_Normal': raw_pval > 0.05,
                    'Transformed_Shapiro_Stat': transformed_stat,
                    'Transformed_Shapiro_Pval': transformed_pval,
                    'Transformed_Is_Normal': transformed_pval > 0.05,
                    'Improvement': 'Yes' if transformed_pval > raw_pval else 'No'
                })

        normality_df = pd.DataFrame(results)
        self.transformation_stats['normality_tests'] = normality_df

        return normality_df

    def calculate_distribution_metrics(self):
        if self.transformed_data is None:
            raise ValueError("No transformation applied yet")

        raw_flat = self.raw_data.values.flatten()
        raw_flat = raw_flat[~np.isnan(raw_flat)]

        trans_flat = self.transformed_data.values.flatten()
        trans_flat = trans_flat[~np.isnan(trans_flat)]

        metrics = {
            'raw': {
                'mean': np.mean(raw_flat),
                'median': np.median(raw_flat),
                'std': np.std(raw_flat),
                'skewness': stats.skew(raw_flat),
                'kurtosis': stats.kurtosis(raw_flat),
                'min': np.min(raw_flat),
                'max': np.max(raw_flat),
                'q25': np.percentile(raw_flat, 25),
                'q75': np.percentile(raw_flat, 75)
            },
            'transformed': {
                'mean': np.mean(trans_flat),
                'median': np.median(trans_flat),
                'std': np.std(trans_flat),
                'skewness': stats.skew(trans_flat),
                'kurtosis': stats.kurtosis(trans_flat),
                'min': np.min(trans_flat),
                'max': np.max(trans_flat),
                'q25': np.percentile(trans_flat, 25),
                'q75': np.percentile(trans_flat, 75)
            }
        }

        self.transformation_stats['distribution_metrics'] = metrics

        return metrics

    def check_sample_distribution_balance(self):
        if self.transformed_data is None:
            raise ValueError("No transformation applied yet")

        sample_medians_raw = self.raw_data.median()
        sample_medians_transformed = self.transformed_data.median()

        balance_stats = {
            'raw': {
                'median_of_medians': sample_medians_raw.median(),
                'std_of_medians': sample_medians_raw.std(),
                'cv_of_medians': (sample_medians_raw.std() / sample_medians_raw.mean()) * 100,
                'range_of_medians': sample_medians_raw.max() - sample_medians_raw.min()
            },
            'transformed': {
                'median_of_medians': sample_medians_transformed.median(),
                'std_of_medians': sample_medians_transformed.std(),
                'cv_of_medians': (sample_medians_transformed.std() / sample_medians_transformed.mean()) * 100,
                'range_of_medians': sample_medians_transformed.max() - sample_medians_transformed.min()
            }
        }

        self.transformation_stats['sample_balance'] = balance_stats

        return balance_stats

    def compare_cv_before_after(self):
        if self.transformed_data is None:
            raise ValueError("No transformation applied yet")

        cv_raw = []
        cv_transformed = []

        for i in range(len(self.raw_data)):
            raw_values = self.raw_data.iloc[i].dropna()
            trans_values = self.transformed_data.iloc[i].dropna()

            if len(raw_values) > 0 and raw_values.mean() > 0:
                cv_raw.append((raw_values.std() / raw_values.mean()) * 100)
            else:
                cv_raw.append(np.nan)

            if len(trans_values) > 0 and trans_values.mean() != 0:
                cv_transformed.append((trans_values.std() / trans_values.mean()) * 100)
            else:
                cv_transformed.append(np.nan)

        cv_comparison = {
            'raw_median_cv': np.nanmedian(cv_raw),
            'raw_mean_cv': np.nanmean(cv_raw),
            'transformed_median_cv': np.nanmedian(cv_transformed),
            'transformed_mean_cv': np.nanmean(cv_transformed),
            'cv_reduction': np.nanmedian(cv_raw) - np.nanmedian(cv_transformed)
        }

        self.transformation_stats['cv_comparison'] = cv_comparison

        return cv_comparison

    def detect_transformation_outliers(self):
        if self.transformed_data is None:
            raise ValueError("No transformation applied yet")

        outlier_counts_raw = []
        outlier_counts_transformed = []

        for col in self.raw_data.columns:
            raw_values = self.raw_data[col].dropna()
            trans_values = self.transformed_data[col].dropna()

            if len(raw_values) > 0:
                q1_raw = raw_values.quantile(0.25)
                q3_raw = raw_values.quantile(0.75)
                iqr_raw = q3_raw - q1_raw
                outliers_raw = ((raw_values < (q1_raw - 1.5*iqr_raw)) |
                               (raw_values > (q3_raw + 1.5*iqr_raw))).sum()
                outlier_counts_raw.append(outliers_raw)

            if len(trans_values) > 0:
                q1_trans = trans_values.quantile(0.25)
                q3_trans = trans_values.quantile(0.75)
                iqr_trans = q3_trans - q1_trans
                outliers_trans = ((trans_values < (q1_trans - 1.5*iqr_trans)) |
                                 (trans_values > (q3_trans + 1.5*iqr_trans))).sum()
                outlier_counts_transformed.append(outliers_trans)

        outlier_stats = {
            'raw_mean_outliers_per_sample': np.mean(outlier_counts_raw),
            'raw_median_outliers_per_sample': np.median(outlier_counts_raw),
            'transformed_mean_outliers_per_sample': np.mean(outlier_counts_transformed),
            'transformed_median_outliers_per_sample': np.median(outlier_counts_transformed)
        }

        self.transformation_stats['outlier_analysis'] = outlier_stats

        return outlier_stats

    def get_transformation_summary(self):
        if self.transformed_data is None:
            raise ValueError("No transformation applied yet")

        summary = {
            'transformation_method': self.transformation_stats.get('method', 'Unknown'),
            'pseudocount': self.transformation_stats.get('pseudocount', 0),
            'zeros_replaced': self.transformation_stats.get('n_zeros_replaced', 0),
            'negative_values': self.transformation_stats.get('n_negative_values', 0)
        }

        if 'distribution_metrics' in self.transformation_stats:
            metrics = self.transformation_stats['distribution_metrics']
            summary['skewness_reduction'] = abs(metrics['raw']['skewness']) - abs(metrics['transformed']['skewness'])
            summary['kurtosis_reduction'] = abs(metrics['raw']['kurtosis']) - abs(metrics['transformed']['kurtosis'])

        if 'cv_comparison' in self.transformation_stats:
            cv_comp = self.transformation_stats['cv_comparison']
            summary['median_cv_before'] = cv_comp['raw_median_cv']
            summary['median_cv_after'] = cv_comp['transformed_median_cv']

        if 'normality_tests' in self.transformation_stats:
            norm_tests = self.transformation_stats['normality_tests']
            summary['normal_samples_before'] = norm_tests['Raw_Is_Normal'].sum()
            summary['normal_samples_after'] = norm_tests['Transformed_Is_Normal'].sum()

        return summary

    def get_transformed_data(self):
        if self.transformed_data is None:
            raise ValueError("No transformation applied yet")

        return self.transformed_data
