import pandas as pd
import numpy as np
import os


class DataLoader:

    def __init__(self, data_path):
        self.data_path = data_path
        self.df = None
        self.abundance_cols = None
        self.metadata = {}

    def load_excel_data(self):
        if not os.path.exists(self.data_path):
            raise FileNotFoundError(f"Data file not found: {self.data_path}")

        try:
            self.df = pd.read_excel(self.data_path)
            print(f"Loaded data: {len(self.df)} proteins x {len(self.df.columns)} columns")
            return self.df
        except Exception as e:
            raise Exception(f"Error loading data: {str(e)}")

    def identify_abundance_columns(self, pattern='Abundances (Normalized)'):
        if self.df is None:
            raise ValueError("Data not loaded. Call load_excel_data() first.")

        self.abundance_cols = [col for col in self.df.columns if pattern in col]

        if len(self.abundance_cols) == 0:
            raise ValueError(f"No abundance columns found matching pattern: {pattern}")

        print(f"Found {len(self.abundance_cols)} abundance columns")
        return self.abundance_cols

    def identify_sample_groups(self, group_definitions):
        if self.abundance_cols is None:
            raise ValueError("Abundance columns not identified. Call identify_abundance_columns() first.")

        identified_groups = {}

        for group_name, criteria in group_definitions.items():
            if isinstance(criteria, list):
                matching_cols = [col for col in self.abundance_cols
                               if any(marker in col for marker in criteria)]
            elif isinstance(criteria, str):
                matching_cols = [col for col in self.abundance_cols if criteria in col]
            else:
                matching_cols = []

            identified_groups[group_name] = matching_cols

        self.metadata['groups'] = identified_groups

        return identified_groups

    def get_protein_info_columns(self):
        if self.df is None:
            raise ValueError("Data not loaded. Call load_excel_data() first.")

        protein_info_cols = []

        standard_cols = ['Accession', 'Gene Symbol', 'Description', 'Protein',
                        'Gene', 'UniProt', 'Entry']

        for col in standard_cols:
            if col in self.df.columns:
                protein_info_cols.append(col)

        return protein_info_cols

    def validate_data_integrity(self):
        if self.df is None:
            raise ValueError("Data not loaded.")

        validation_report = {
            'total_rows': len(self.df),
            'total_columns': len(self.df.columns),
            'abundance_columns': len(self.abundance_cols) if self.abundance_cols else 0,
            'protein_info_columns': len(self.get_protein_info_columns()),
            'duplicate_rows': self.df.duplicated().sum(),
            'all_nan_rows': self.df.isna().all(axis=1).sum(),
            'all_nan_columns': self.df.isna().all(axis=0).sum()
        }

        return validation_report

    def create_sample_metadata(self, group_mappings):
        sample_metadata = []

        for sample_col in self.abundance_cols:
            metadata_row = {'Sample': sample_col}

            for attr_name, group_dict in group_mappings.items():
                for group_value, group_cols in group_dict.items():
                    if sample_col in group_cols:
                        metadata_row[attr_name] = group_value
                        break

            sample_metadata.append(metadata_row)

        sample_metadata_df = pd.DataFrame(sample_metadata)

        return sample_metadata_df

    def get_abundance_matrix(self):
        if self.abundance_cols is None:
            raise ValueError("Abundance columns not identified.")

        return self.df[self.abundance_cols]

    def get_protein_annotations(self):
        if self.df is None:
            raise ValueError("Data not loaded.")

        annotation_cols = self.get_protein_info_columns()

        if len(annotation_cols) == 0:
            raise ValueError("No protein annotation columns found.")

        return self.df[annotation_cols]

    def export_processed_data(self, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        abundance_matrix = self.get_abundance_matrix()
        abundance_matrix.to_csv(os.path.join(output_dir, 'abundance_matrix.csv'), index=False)

        protein_annotations = self.get_protein_annotations()
        protein_annotations.to_csv(os.path.join(output_dir, 'protein_annotations.csv'), index=False)

        if 'groups' in self.metadata:
            groups_df = pd.DataFrame({
                'Group': list(self.metadata['groups'].keys()),
                'N_Samples': [len(cols) for cols in self.metadata['groups'].values()]
            })
            groups_df.to_csv(os.path.join(output_dir, 'group_definitions.csv'), index=False)

        print(f"Exported processed data to: {output_dir}")
