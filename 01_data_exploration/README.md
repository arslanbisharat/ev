# Data Exploration & Quality Control Analysis

## Overview

This directory contains comprehensive quality control and exploratory analysis performed on the EV proteomics dataset **BEFORE** any statistical testing or differential analysis.

## Analysis Performed

### 1. Sample Quality Control
- Detection of samples with high missing data (>50%)
- Identification of samples with abnormal total intensity
- Outlier detection using PCA and Mahalanobis distance

### 2. Protein Quality Control
- Assessment of protein detection rates across samples
- Calculation of coefficient of variation (CV) for each protein
- Identification of consistently detected vs. sporadically detected proteins

### 3. Missing Value Analysis
- Overall missing data percentage calculation
- Missing value patterns across samples and proteins
- Visualization of missing data structure

### 4. Sample Correlation Analysis
- Pairwise correlation matrix of all samples
- Assessment of within-group vs. between-group correlations

### 5. Principal Component Analysis
- PCA for dimensionality reduction and sample clustering
- Variance explained by top principal components
- Detection of outlier samples

### 6. Normality Testing
- Shapiro-Wilk tests on sample distributions
- Assessment of whether parametric tests are appropriate

### 7. Descriptive Statistics
- Group-wise summary statistics (mean, median, SD)
- Detection rate analysis by group
- Abundance percentile calculations

### 8. Filtering Threshold Analysis
- Evaluation of different detection thresholds (10-70%)
- Assessment of how many proteins would be retained at each threshold

## Key Findings

### Data Structure
- **Total Proteins**: 2,800
- **Total Samples**: 120
  - Set 1 (Old samples): 60
  - Set 2 (Fresh samples): 60
  - Sarcoidosis cases: 60 (30 per set)
  - Controls: 60 (30 per set)

### Data Quality
- **Mean Detection Rate**: 57.06%
- **Missing Data**: 42.94% overall
- **Proteins Detected in All Samples**: 972 (35%)
- **Proteins Detected in >50% Samples**: 1,561 (56%)
- **Median CV**: 45.48%

### Quality Flags
- **Samples Flagged**: 5 out of 120 (4.2%)
  - High missing data: 2 samples
  - Low total intensity: 2 samples
  - High total intensity: 1 sample

- **PCA Outliers**: 2 samples
  - F1: Sample, Sarcoidosis, 302-1
  - F101: Sample, Sarcoidosis, sarc9

### Distribution Characteristics
- **All samples fail normality tests** (Shapiro-Wilk p < 0.05)
- This suggests need for:
  - Log transformation before analysis
  - Non-parametric tests OR
  - Large sample size for CLT to apply

### PCA Results
- **PC1 explains**: 19.44% of variance
- **PC2 explains**: 11.60% of variance
- **Cumulative (PC1+PC2)**: 31.05%

Low variance explained suggests:
- High biological heterogeneity
- Multiple sources of variation (storage time, case/control, technical batch)
- Need for variance partitioning analysis

### Group Comparisons
Group statistics show:
- Set2 (Fresh) has slightly higher mean abundance than Set1 (Old)
- Detection rates similar between Sarcoidosis and Control within each set
- Slightly lower detection in Set2 compared to Set1

### Filtering Recommendations
At 30% detection threshold:
- Set1 Sarcoidosis: 1,828 proteins (65%)
- Set1 Control: 1,797 proteins (64%)
- Set2 Sarcoidosis: 1,642 proteins (59%)
- Set2 Control: 1,644 proteins (59%)

## Generated Outputs

### Figures (01_data_exploration/figures/)
1. `missing_values_heatmap.png` - Heatmap showing missing value patterns
2. `sample_completeness.png` - 4-panel plot of sample QC metrics
3. `protein_detection_rates.png` - 4-panel plot of protein detection statistics
4. `pca_plot.png` - PCA scatter plot with outlier identification
5. `correlation_heatmap.png` - Sample-to-sample correlation matrix
6. `intensity_distributions.png` - Distribution of abundance values
7. `filtering_threshold_analysis.png` - Effect of different filtering thresholds
8. `clustered_heatmap.png` - Heatmap of top 100 most variable proteins

### Reports (01_data_exploration/reports/)
1. `sample_qc_report.csv` - Detailed QC metrics for each sample
2. `protein_qc_report.csv` - Detection rates and CV for each protein
3. `pca_results.csv` - PCA coordinates and outlier flags
4. `normality_test.csv` - Shapiro-Wilk test results
5. `group_summary.csv` - Summary statistics by group
6. `filtering_threshold_analysis.csv` - Proteins retained at each threshold
7. `qc_summary_report.csv` - Overall QC summary statistics

## Recommendations for Next Steps

### 1. Sample Exclusion Decisions
Review flagged samples and PCA outliers:
- F50 (Control, 545): 53.5% missing - consider excluding
- F111 (Sarcoidosis): 51.6% missing - consider excluding
- F1 and F101 are PCA outliers - investigate further

### 2. Filtering Strategy
Recommend 30-40% detection threshold in at least one condition:
- Balances protein coverage vs. data quality
- Retains ~1,600-1,800 proteins

### 3. Statistical Approach
Due to non-normal distributions:
- Apply log2 transformation
- Use linear models robust to non-normality OR
- Use non-parametric tests (Wilcoxon) for small subgroups

### 4. Confounders to Address
- Account for batch effects (Set1 vs Set2 processed at different times)
- Include covariates in models (age, sex if available)
- Variance partitioning to quantify storage vs. biological effects

### 5. Imputation Strategy
- 43% missing data requires careful imputation
- Recommend: KNN imputation OR minimum value imputation
- Sensitivity analysis with different methods

## Scripts

All analysis code is in `01_data_exploration/scripts/`:
- `data_loader.py` - Data loading and validation
- `quality_control.py` - QC metrics and outlier detection
- `descriptive_stats.py` - Summary statistics and filtering analysis
- `qc_visualization.py` - All plotting functions
- `run_qc_analysis.py` - Main pipeline script

## Running the Analysis

```bash
cd 01_data_exploration/scripts
python run_qc_analysis.py
```

This will regenerate all figures and reports.
