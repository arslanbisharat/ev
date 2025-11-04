import pandas as pd
import numpy as np

class ProteomicsDataLoader:

    def __init__(self, data_path):
        self.data_path = data_path
        self.df = None
        self.abundance_cols = None
        self.set1_cols = None
        self.set2_cols = None
        self.sarc_cols = None
        self.control_cols = None

    def load_data(self):
        self.df = pd.read_excel(self.data_path)
        return self.df

    def identify_sample_columns(self):
        self.abundance_cols = [col for col in self.df.columns if 'Abundances (Normalized)' in col]

        set1_markers = ['302-', '310-', '311-', '327-', '338-', '343-', '346-', '364-', '368-', '371-',
                        '375-', '432-', '435-', '436-', '443-', '465-', '467-', '481-', '485-', '486-',
                        '500-', '510-', '533-', '537-', '545-', '561-', '567-', '577-', '584-', '613-']

        self.set1_cols = [col for col in self.abundance_cols if any(x in col for x in set1_markers)]
        self.set2_cols = [col for col in self.abundance_cols if ('Sarc' in col or ', C' in col) and col not in self.set1_cols]
        self.sarc_cols = [col for col in self.abundance_cols if 'Sarcoidosis' in col]
        self.control_cols = [col for col in self.abundance_cols if 'Control' in col]

        return {
            'set1': self.set1_cols,
            'set2': self.set2_cols,
            'sarcoidosis': self.sarc_cols,
            'control': self.control_cols
        }

    def get_group_columns(self):
        set1_sarc = [col for col in self.set1_cols if col in self.sarc_cols]
        set1_ctrl = [col for col in self.set1_cols if col in self.control_cols]
        set2_sarc = [col for col in self.set2_cols if col in self.sarc_cols]
        set2_ctrl = [col for col in self.set2_cols if col in self.control_cols]

        return {
            'set1_sarc': set1_sarc,
            'set1_ctrl': set1_ctrl,
            'set2_sarc': set2_sarc,
            'set2_ctrl': set2_ctrl
        }

    def get_summary(self):
        groups = self.get_group_columns()
        return {
            'total_proteins': len(self.df),
            'total_samples': len(self.abundance_cols),
            'set1_samples': len(self.set1_cols),
            'set2_samples': len(self.set2_cols),
            'set1_sarc': len(groups['set1_sarc']),
            'set1_ctrl': len(groups['set1_ctrl']),
            'set2_sarc': len(groups['set2_sarc']),
            'set2_ctrl': len(groups['set2_ctrl'])
        }
