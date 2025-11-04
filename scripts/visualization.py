import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import patheffects
import numpy as np
import os

class ProteomicsVisualizer:

    def __init__(self, output_dir=None):
        if output_dir is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            project_dir = os.path.dirname(script_dir)
            output_dir = os.path.join(project_dir, 'figures')
        self.output_dir = output_dir

    def plot_log_transformation(self, data_before, data_after, filename='log_transformation.png'):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

        vals_before = data_before.values.flatten()[~np.isnan(data_before.values.flatten())]
        vals_after = data_after.values.flatten()[~np.isnan(data_after.values.flatten())]

        ax1.hist(vals_before, bins=50, edgecolor='black', alpha=0.7)
        ax1.set_xlabel('Raw Abundance')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Before Log2 Transformation')
        ax1.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

        ax2.hist(vals_after, bins=50, edgecolor='black', alpha=0.7, color='green')
        ax2.set_xlabel('Log2 Abundance')
        ax2.set_ylabel('Frequency')
        ax2.set_title('After Log2 Transformation')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_volcano(self, results_df, title, filename, n_label=10):
        fig, ax = plt.subplots(figsize=(10, 8))

        ax.scatter(results_df['Log2FC'], -np.log10(results_df['P_value']),
                   c='gray', alpha=0.5, s=20, label='Not significant')

        sig_mask = results_df['FDR'] < 0.05
        ax.scatter(results_df.loc[sig_mask, 'Log2FC'],
                   -np.log10(results_df.loc[sig_mask, 'P_value']),
                   c='red', alpha=0.7, s=30, label=f'Significant (n={sig_mask.sum()})')

        top_proteins = results_df.nsmallest(n_label, 'FDR')
        for _, row in top_proteins.iterrows():
            ax.annotate(row['Gene'],
                        xy=(row['Log2FC'], -np.log10(row['P_value'])),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, alpha=0.7)

        ax.axhline(-np.log10(0.05), color='blue', linestyle='--', linewidth=1, alpha=0.5)
        ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-Log10 P-value')
        ax.set_title(title)
        ax.legend()
        ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    def plot_venn_diagram(self, overlap_data, n_sig1, n_sig2, filename='venn_diagram.png'):
        fig, ax = plt.subplots(figsize=(12, 9), facecolor='white')

        circle1 = Circle((-1.5, 0), 2.5, color='#E74C3C', alpha=0.4, zorder=1)
        circle2 = Circle((1.5, 0), 2.5, color='#3498DB', alpha=0.4, zorder=1)
        ax.add_patch(circle1)
        ax.add_patch(circle2)

        circle1_border = Circle((-1.5, 0), 2.5, fill=False, edgecolor='#C0392B', linewidth=3, zorder=2)
        circle2_border = Circle((1.5, 0), 2.5, fill=False, edgecolor='#2980B9', linewidth=3, zorder=2)
        ax.add_patch(circle1_border)
        ax.add_patch(circle2_border)

        def add_text_with_outline(ax, x, y, text, fontsize, color):
            txt = ax.text(x, y, text, fontsize=fontsize, color=color,
                          ha='center', va='center', weight='bold', zorder=5)
            txt.set_path_effects([patheffects.Stroke(linewidth=3, foreground='white'),
                                  patheffects.Normal()])

        add_text_with_outline(ax, -2.5, 0, f'{overlap_data["n_only_set1"]}', 42, '#8B0000')
        add_text_with_outline(ax, 2.5, 0, f'{overlap_data["n_only_set2"]}', 42, '#00008B')
        add_text_with_outline(ax, 0, 0, f'{overlap_data["n_overlap"]}', 46, '#006400')

        ax.text(-1.5, 3.5, f'Set 1 (Old){n_sig1} proteins',
                fontsize=15, ha='center', weight='bold',
                bbox=dict(boxstyle='round,pad=0.6', facecolor='#FFE5E5',
                          edgecolor='#E74C3C', linewidth=2))

        ax.text(1.5, 3.5, f'Set 2 (Fresh){n_sig2} proteins',
                fontsize=15, ha='center', weight='bold',
                bbox=dict(boxstyle='round,pad=0.6', facecolor='#E5F2FF',
                          edgecolor='#3498DB', linewidth=2))

        ax.text(0, 4.8, 'Differentially Expressed ProteinsSarcoidosis vs Control (FDR < 0.05)',
                fontsize=17, ha='center', weight='bold',
                bbox=dict(boxstyle='round,pad=0.7', facecolor='white',
                          edgecolor='#2C3E50', linewidth=2.5))

        ax.text(0, 4.2, f'Jaccard Similarity = {overlap_data["jaccard"]:.3f}',
                fontsize=13, ha='center', style='italic',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#F8F9FA',
                          edgecolor='gray', linewidth=1.5))

        ax.set_xlim(-6, 6)
        ax.set_ylim(-3.5, 6)
        ax.set_aspect('equal')
        ax.axis('off')

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()

    def plot_correlation(self, merged_data, pearson_r, overlap_info, filename='correlation.png'):
        fig, ax = plt.subplots(figsize=(10, 10))

        colors = []
        for _, row in merged_data.iterrows():
            gene = row['Gene']
            if gene in overlap_info['overlap']:
                colors.append('red')
            elif gene in overlap_info['only_set1']:
                colors.append('orange')
            elif gene in overlap_info['only_set2']:
                colors.append('blue')
            else:
                colors.append('gray')

        ax.scatter(merged_data['Log2FC_Set1'], merged_data['Log2FC_Set2'],
                   c=colors, alpha=0.6, s=25)

        ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
        ax.axvline(0, color='black', linestyle='--', linewidth=0.5)

        from scipy.stats import linregress
        slope, intercept, r_value, _, _ = linregress(merged_data['Log2FC_Set1'], merged_data['Log2FC_Set2'])
        x_line = np.array([merged_data['Log2FC_Set1'].min(), merged_data['Log2FC_Set1'].max()])
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, 'k-', linewidth=2, label=f'R² = {r_value**2:.3f}')

        ax.set_xlabel('Log2 FC - Set 1 (Old)', fontsize=12)
        ax.set_ylabel('Log2 FC - Set 2 (Fresh)', fontsize=12)
        ax.set_title(f'Effect Size Correlation (r={pearson_r:.3f})', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()
