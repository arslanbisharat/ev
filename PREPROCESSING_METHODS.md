# Preprocessing Methods for EV Proteomics Analysis


STEP-BY-STEP PREPROCESSING METHODS
===================================

1. DATA LOADING
   - Input: Excel file with protein abundance data (2,800 proteins × 120 samples)
   - Sample identification: Automatic detection of sample columns based on naming patterns
   - Group assignment: Set 1/2, Sarcoidosis/Control based on sample metadata

2. QUALITY CONTROL (BEFORE FILTERING)
   Purpose: Identify low-quality samples and proteins

   Sample QC metrics:
   - Missing data % per sample (flag if >50%)
   - Total intensity per sample (flag outliers)
   - PCA outlier detection (Mahalanobis distance)

   Protein QC metrics:
   - Detection rate across samples
   - Coefficient of variation (CV)
   - Missing value patterns

   Decision: Review flagged samples/proteins before proceeding

3. PROTEIN FILTERING
   Purpose: Remove proteins with too much missing data

   Criteria: Keep proteins detected in ≥40% of samples in AT LEAST ONE CONDITION
   - For Set 1: Keep if ≥40% in Sarc OR ≥40% in Control
   - For Set 2: Keep if ≥40% in Sarc OR ≥40% in Control

   Rationale:
   - 40% threshold balances data quality vs. coverage
   - "OR" logic preserves disease-specific proteins
   - Applied per-set to avoid bias

   Alternative (original): Require presence in BOTH sets
   - More conservative, but loses set-specific proteins

   Result: ~1,600-1,800 proteins retained per set

4. LOG2 TRANSFORMATION
   Purpose: Stabilize variance and normalize distribution

   Method:
   - Replace 0 values with NA (true missing)
   - Apply log2 transformation: log2(abundance)

   Criteria for use:
   - Proteomics data is typically log-normally distributed
   - Fold-changes become symmetric (2x up = -2x down)
   - Stabilizes variance across abundance range

   When to skip:
   - Data already log-transformed
   - Using rank-based methods that don't assume normality

5. IMPUTATION
   Purpose: Handle missing values for statistical testing

   Method: Minimum value - 1 approach
   - For each protein, find minimum observed value
   - Impute missing as min - 1 (on log2 scale)

   Rationale:
   - Missing values in proteomics often = below detection limit
   - This assumes missing = low abundance
   - Conservative approach that doesn't inflate significance

   Alternatives:
   - KNN imputation (use similar proteins)
   - Half-minimum (min/2 on linear scale)
   - Random draws from low end of distribution

   When to skip:
   - Using methods that handle missing data directly
   - Too much missing data (>70% may be unreliable)

6. BATCH EFFECT CORRECTION
   Status: NOT APPLIED in current analysis

   Why not:
   - Set 1 and Set 2 are analyzed SEPARATELY
   - No attempt to merge or directly compare samples
   - Comparison is done at summary statistic level (fold-changes)

   When to apply:
   - If combining sets for joint analysis
   - If technical batch effects detected within a set

   Method (if needed):
   - ComBat (from sva package in R)
   - Remove batch variation while preserving biological variation
   - Requires batch labels and ideally balanced design

7. NORMALIZATION
   Status: NOT EXPLICITLY APPLIED

   Current approach:
   - Log2 transformation provides implicit normalization
   - Within-sample normalization assumed from data provider

   When to apply:
   - If samples have different total protein amounts
   - If systematic differences in abundance distributions

   Methods:
   - Median centering (shift each sample to same median)
   - Quantile normalization (make distributions identical)
   - Total intensity normalization (divide by total per sample)

   Check needed:
   - Compare total intensity distributions
   - If similar, normalization may not be needed

8. STATISTICAL TESTING
   Method: Student's t-test with FDR correction

   Test: Welch's t-test (unequal variance)
   - Compare Sarcoidosis vs Control within each set
   - Assumes: independence, approximate normality (CLT applies with n=30)

   Multiple testing correction:
   - Benjamini-Hochberg FDR (False Discovery Rate)
   - Threshold: FDR < 0.05

   Effect size:
   - Log2 fold-change (Sarc - Control)
   - Cohen's d for standardized effect size

DECISION TREE FOR PREPROCESSING
================================

1. Start with raw data
   ↓
2. QC → Remove bad samples? (if >50% missing or clear outlier)
   ↓
3. Filter proteins → Use 30-40% detection threshold
   ↓
4. Check distribution → Skewed? → YES → Log2 transform
   ↓                              ↓ NO
5. Missing data <70%? → YES → Impute (min-1 or KNN)
   ↓                     ↓ NO → Consider removing protein
6. Multiple sets/batches? → YES → Check for batch effects → Correct if needed
   ↓                         ↓ NO
7. Check total intensity → Different? → YES → Normalize
   ↓                         ↓ NO
8. Statistical testing (t-test, FDR correction)

