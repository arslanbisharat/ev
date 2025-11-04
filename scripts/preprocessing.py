import pandas as pd
import numpy as np

class ProteomicsPreprocessor:

    def __init__(self, df, abundance_cols):
        self.df = df
        self.abundance_cols = abundance_cols
        self.df_filtered = None
        self.data_log2 = None
        self.data_imputed = None

    def filter_proteins(self, group_cols_dict, detection_threshold=0.4):
        set1_sarc_det = (~self.df[group_cols_dict['set1_sarc']].isna()).sum(axis=1) / len(group_cols_dict['set1_sarc'])
        set1_ctrl_det = (~self.df[group_cols_dict['set1_ctrl']].isna()).sum(axis=1) / len(group_cols_dict['set1_ctrl'])
        set1_keep = (set1_sarc_det >= detection_threshold) | (set1_ctrl_det >= detection_threshold)

        set2_sarc_det = (~self.df[group_cols_dict['set2_sarc']].isna()).sum(axis=1) / len(group_cols_dict['set2_sarc'])
        set2_ctrl_det = (~self.df[group_cols_dict['set2_ctrl']].isna()).sum(axis=1) / len(group_cols_dict['set2_ctrl'])
        set2_keep = (set2_sarc_det >= detection_threshold) | (set2_ctrl_det >= detection_threshold)

        proteins_keep = set1_keep & set2_keep
        self.df_filtered = self.df[proteins_keep].copy()

        return self.df_filtered, proteins_keep.sum()

    def log2_transform(self):
        abundance_data = self.df_filtered[self.abundance_cols].replace(0, np.nan)
        self.data_log2 = np.log2(abundance_data)
        return self.data_log2

    def impute_missing(self):
        self.data_imputed = self.data_log2.copy()

        for i in range(len(self.data_imputed)):
            protein_values = self.data_log2.iloc[i]
            if protein_values.notna().sum() > 0:
                min_val = protein_values.min()
                self.data_imputed.iloc[i] = protein_values.fillna(min_val - 1)

        return self.data_imputed

    def get_processed_data(self):
        return self.data_imputed
