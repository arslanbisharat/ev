# EV Proteomics Storage Effect Analysis

## Research Question
Does long-term storage of plasma samples affect the integrity of extracellular vesicle (EV) proteomics data?

## Study Design
- **Set 1 (ACCESS trial)**: Old/stored samples - 30 Sarcoidosis + 30 Control
- **Set 2 (LUC trial)**: Fresh samples - 30 Sarcoidosis + 30 Control
- **Total**: 2,800 proteins measured across 120 samples

## Project Structure
```
EV_Proteomics_Analysis/
├── scripts/
│   ├── data_loader.py            # Data loading and sample identification
│   ├── preprocessing.py          # Filtering, log2 transform, imputation
│   ├── differential_analysis.py  # T-tests, overlap, correlation
│   ├── visualization.py          # All plotting functions
│   └── main_pipeline.py          # Main analysis pipeline
├── data/                         # Input data files
├── results/                      # Output CSV files
├── figures/                      # Output figures
└── README.md
```

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
cd EV_Proteomics_Analysis
python scripts/main_pipeline.py
```

## Workflow

1. **Data Loading**: Identify sample groups (Set 1/2, Sarc/Control)
2. **Preprocessing**: Filter (≥40% detection), log2 transform, impute missing values
3. **Differential Analysis**: T-tests with FDR correction for each set
4. **Overlap Analysis**: Compare significant proteins between sets
5. **Visualization**: Volcano plots, Venn diagram, correlation plot
6. **Export Results**: CSV files with all statistics

## Key Outputs

### Results Files
- `results_set1_old.csv` - Differential analysis for old samples
- `results_set2_fresh.csv` - Differential analysis for fresh samples
- `storage_stable_proteins.csv` - Proteins significant in both sets

### Figures
- `log_transformation.png` - Before/after log2 transformation
- `volcano_set1.png` - Volcano plot for old samples
- `volcano_set2.png` - Volcano plot for fresh samples
- `venn_diagram.png` - Overlap of significant proteins
- `correlation.png` - Fold-change correlation between sets

## Main Findings

- **4x** more proteins detected in fresh vs old samples (347 vs 85)
- **r = 0.451** correlation of fold-changes (moderate preservation)
- **41 storage-stable proteins** identified (robust biomarkers)
- **Jaccard = 0.106** (10.6% overlap)

## Interpretation

Storage has a **moderate** impact on EV proteomics:
- Reduces detection sensitivity (~75% reduction)
- Preserves core disease biology (moderate correlation)
- 41 proteins are storage-resistant (high-confidence biomarkers)

## Authors
Proteomics Analysis Pipeline
