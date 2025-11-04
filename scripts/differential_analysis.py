import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

class DifferentialAnalyzer:

    def __init__(self, data, protein_info):
        self.data = data
        self.protein_info = protein_info

    def perform_ttest(self, group1_cols, group2_cols):
        results = []

        for i in range(len(self.data)):
            g1 = self.data.iloc[i][group1_cols].values
            g2 = self.data.iloc[i][group2_cols].values

            t_stat, p_val = stats.ttest_ind(g1, g2)

            pooled_std = np.sqrt(((len(g1)-1)*np.var(g1, ddof=1) + (len(g2)-1)*np.var(g2, ddof=1)) / (len(g1)+len(g2)-2))
            cohens_d = (np.mean(g1) - np.mean(g2)) / pooled_std if pooled_std > 0 else 0

            log2fc = np.mean(g1) - np.mean(g2)

            results.append({
                'Gene': self.protein_info.iloc[i]['Gene Symbol'] if pd.notna(self.protein_info.iloc[i]['Gene Symbol']) else 'Unknown',
                'Accession': self.protein_info.iloc[i]['Accession'],
                'Description': self.protein_info.iloc[i]['Description'],
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

    def get_significant_proteins(self, results_df, fdr_threshold=0.05):
        return results_df[results_df['FDR'] < fdr_threshold]

    def calculate_overlap(self, results1, results2, fdr_threshold=0.05):
        sig1 = results1[results1['FDR'] < fdr_threshold]
        sig2 = results2[results2['FDR'] < fdr_threshold]

        genes1 = set(sig1['Gene'].tolist())
        genes2 = set(sig2['Gene'].tolist())

        overlap = genes1 & genes2
        only_set1 = genes1 - genes2
        only_set2 = genes2 - genes1

        jaccard = len(overlap) / len(genes1 | genes2) if len(genes1 | genes2) > 0 else 0

        return {
            'overlap': overlap,
            'only_set1': only_set1,
            'only_set2': only_set2,
            'jaccard': jaccard,
            'n_overlap': len(overlap),
            'n_only_set1': len(only_set1),
            'n_only_set2': len(only_set2)
        }

    def calculate_correlation(self, results1, results2):
        merged = results1[['Gene', 'Log2FC']].merge(
            results2[['Gene', 'Log2FC']],
            on='Gene',
            suffixes=('_Set1', '_Set2')
        )

        pearson_r = merged['Log2FC_Set1'].corr(merged['Log2FC_Set2'], method='pearson')
        spearman_r = merged['Log2FC_Set1'].corr(merged['Log2FC_Set2'], method='spearman')

        return {
            'pearson': pearson_r,
            'spearman': spearman_r,
            'merged_data': merged
        }
