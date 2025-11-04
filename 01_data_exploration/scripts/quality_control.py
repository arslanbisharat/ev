import pandas as pd
import numpy as np
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')


class ProteomicsQC:

    def __init__(self, df, abundance_cols, sample_info=None):
        self.df = df
        self.abundance_cols = abundance_cols
        self.sample_info = sample_info
        self.qc_results = {}

    def check_sample_completeness(self):
        abundance_data = self.df[self.abundance_cols]

        missing_per_sample = abundance_data.isna().sum() / len(abundance_data) * 100
        detected_per_sample = (~abundance_data.isna()).sum()
        total_intensity_per_sample = abundance_data.sum()

        sample_qc = pd.DataFrame({
            'Sample': self.abundance_cols,
            'Missing_Percent': missing_per_sample.values,
            'Detected_Proteins': detected_per_sample.values,
            'Total_Intensity': total_intensity_per_sample.values,
            'Log_Total_Intensity': np.log10(total_intensity_per_sample.values + 1)
        })

        mean_intensity = sample_qc['Log_Total_Intensity'].mean()
        std_intensity = sample_qc['Log_Total_Intensity'].std()

        sample_qc['High_Missing'] = sample_qc['Missing_Percent'] > 50
        sample_qc['Low_Intensity'] = sample_qc['Log_Total_Intensity'] < (mean_intensity - 2*std_intensity)
        sample_qc['High_Intensity'] = sample_qc['Log_Total_Intensity'] > (mean_intensity + 2*std_intensity)
        sample_qc['QC_Flag'] = sample_qc['High_Missing'] | sample_qc['Low_Intensity'] | sample_qc['High_Intensity']

        self.qc_results['sample_completeness'] = sample_qc

        return sample_qc

    def check_protein_completeness(self):
        abundance_data = self.df[self.abundance_cols]

        detection_per_protein = (~abundance_data.isna()).sum(axis=1) / len(self.abundance_cols) * 100

        protein_stats = []
        for i in range(len(abundance_data)):
            values = abundance_data.iloc[i].dropna()
            if len(values) > 0:
                protein_stats.append({
                    'Detection_Rate': detection_per_protein.iloc[i],
                    'N_Detected': len(values),
                    'Mean_Abundance': values.mean(),
                    'Median_Abundance': values.median(),
                    'SD_Abundance': values.std(),
                    'CV_Percent': (values.std() / values.mean() * 100) if values.mean() > 0 else np.nan,
                    'Min_Abundance': values.min(),
                    'Max_Abundance': values.max()
                })
            else:
                protein_stats.append({
                    'Detection_Rate': 0,
                    'N_Detected': 0,
                    'Mean_Abundance': np.nan,
                    'Median_Abundance': np.nan,
                    'SD_Abundance': np.nan,
                    'CV_Percent': np.nan,
                    'Min_Abundance': np.nan,
                    'Max_Abundance': np.nan
                })

        protein_qc = pd.DataFrame(protein_stats)

        if 'Gene Symbol' in self.df.columns:
            protein_qc['Gene'] = self.df['Gene Symbol'].values
        if 'Accession' in self.df.columns:
            protein_qc['Accession'] = self.df['Accession'].values
        if 'Description' in self.df.columns:
            protein_qc['Description'] = self.df['Description'].values

        self.qc_results['protein_completeness'] = protein_qc

        return protein_qc

    def calculate_sample_correlations(self):
        abundance_data = self.df[self.abundance_cols]
        corr_matrix = abundance_data.T.corr(method='pearson')
        self.qc_results['sample_correlation'] = corr_matrix
        return corr_matrix

    def perform_pca(self, min_detection_rate=0.3):
        abundance_data = self.df[self.abundance_cols]

        detection_rates = (~abundance_data.isna()).sum(axis=1) / len(self.abundance_cols)
        filtered_data = abundance_data[detection_rates >= min_detection_rate]

        print(f"  Using {len(filtered_data)} proteins with >{min_detection_rate*100}% detection for PCA")

        pca_data = filtered_data.T.copy()
        for col in pca_data.columns:
            if pca_data[col].isna().any():
                min_val = pca_data[col].min()
                if pd.notna(min_val):
                    pca_data[col].fillna(min_val - 1, inplace=True)
                else:
                    pca_data[col].fillna(0, inplace=True)

        scaler = StandardScaler()
        pca_data_scaled = scaler.fit_transform(pca_data)

        pca = PCA()
        pca_coords = pca.fit_transform(pca_data_scaled)

        pca_results = pd.DataFrame({
            'Sample': self.abundance_cols,
            'PC1': pca_coords[:, 0],
            'PC2': pca_coords[:, 1],
            'PC3': pca_coords[:, 2] if pca_coords.shape[1] > 2 else 0
        })

        from scipy.spatial.distance import mahalanobis
        coords_3d = pca_coords[:, :min(3, pca_coords.shape[1])]
        centroid = coords_3d.mean(axis=0)

        try:
            cov_matrix = np.cov(coords_3d.T)
            inv_cov = np.linalg.inv(cov_matrix)
            mahal_dist = [mahalanobis(sample, centroid, inv_cov) for sample in coords_3d]

            mean_dist = np.mean(mahal_dist)
            std_dist = np.std(mahal_dist)
            pca_results['Mahalanobis_Distance'] = mahal_dist
            pca_results['PCA_Outlier'] = mahal_dist > (mean_dist + 3*std_dist)
        except:
            pca_results['Mahalanobis_Distance'] = 0
            pca_results['PCA_Outlier'] = False

        self.qc_results['pca'] = {
            'coordinates': pca_results,
            'variance_explained': pca.explained_variance_ratio_,
            'n_proteins_used': len(filtered_data)
        }

        return pca_results, pca.explained_variance_ratio_

    def check_normality(self, sample_size=1000):
        abundance_data = self.df[self.abundance_cols]

        normality_results = []

        proteins_with_data = abundance_data.dropna(thresh=len(self.abundance_cols)*0.5).index

        if len(proteins_with_data) > sample_size:
            sample_proteins = np.random.choice(proteins_with_data, sample_size, replace=False)
        else:
            sample_proteins = proteins_with_data

        for col in self.abundance_cols[:10]:
            sample_data = abundance_data.loc[sample_proteins, col].dropna()

            if len(sample_data) >= 3:
                stat, pval = stats.shapiro(sample_data)
                normality_results.append({
                    'Sample': col,
                    'Shapiro_Stat': stat,
                    'P_value': pval,
                    'Is_Normal': pval > 0.05
                })

        normality_df = pd.DataFrame(normality_results)
        self.qc_results['normality'] = normality_df

        return normality_df

    def identify_high_cv_proteins(self, cv_threshold=50):
        if 'protein_completeness' not in self.qc_results:
            self.check_protein_completeness()

        protein_qc = self.qc_results['protein_completeness']
        high_cv = protein_qc[protein_qc['CV_Percent'] > cv_threshold].copy()
        high_cv = high_cv.sort_values('CV_Percent', ascending=False)

        return high_cv

    def generate_qc_summary(self):
        summary = {
            'total_proteins': len(self.df),
            'total_samples': len(self.abundance_cols),
        }

        if 'sample_completeness' in self.qc_results:
            sample_qc = self.qc_results['sample_completeness']
            summary['samples_with_flags'] = sample_qc['QC_Flag'].sum()
            summary['samples_high_missing'] = sample_qc['High_Missing'].sum()
            summary['samples_low_intensity'] = sample_qc['Low_Intensity'].sum()
            summary['samples_high_intensity'] = sample_qc['High_Intensity'].sum()
            summary['mean_missing_per_sample'] = sample_qc['Missing_Percent'].mean()
            summary['median_proteins_per_sample'] = sample_qc['Detected_Proteins'].median()

        if 'protein_completeness' in self.qc_results:
            protein_qc = self.qc_results['protein_completeness']
            summary['proteins_detected_all_samples'] = (protein_qc['Detection_Rate'] == 100).sum()
            summary['proteins_detected_50pct'] = (protein_qc['Detection_Rate'] >= 50).sum()
            summary['proteins_detected_25pct'] = (protein_qc['Detection_Rate'] >= 25).sum()
            summary['mean_detection_rate'] = protein_qc['Detection_Rate'].mean()
            summary['median_cv'] = protein_qc['CV_Percent'].median()

        if 'pca' in self.qc_results:
            pca_results = self.qc_results['pca']['coordinates']
            summary['pca_outliers'] = pca_results['PCA_Outlier'].sum()
            summary['variance_pc1'] = self.qc_results['pca']['variance_explained'][0]
            summary['variance_pc2'] = self.qc_results['pca']['variance_explained'][1]

        if 'normality' in self.qc_results:
            normality = self.qc_results['normality']
            summary['samples_normally_distributed'] = normality['Is_Normal'].sum()
            summary['samples_non_normal'] = (~normality['Is_Normal']).sum()

        self.qc_results['summary'] = summary

        return summary
