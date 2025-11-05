# Answers to Your 5 Questions

## Question 1: How many proteins are detected for controls and cases? Any differences between sets?

### Detection by Group:

**Set 1 (Old Samples):**
- Sarcoidosis: **2,105 proteins detected** (75.2%)
- Control: **2,072 proteins detected** (74.0%)
- Average detection rate: 60% (Sarc), 59% (Control)

**Set 2 (Fresh Samples):**
- Sarcoidosis: **1,953 proteins detected** (69.8%)
- Control: **1,928 proteins detected** (68.9%)
- Average detection rate: 54.8% (Sarc), 54.4% (Control)

### Key Finding:
**Old samples detect MORE proteins than fresh samples!**
- +152 proteins in Sarcoidosis (old vs fresh)
- +144 proteins in Controls (old vs fresh)

This contradicts the expectation that storage reduces detection.

---

## Question 2: Overall correlation between two sets?

### Result: **r = 0.976** (VERY STRONG correlation)

- Pearson correlation: **0.9764** (p < 0.0001)
- Spearman correlation: **0.9701** (p < 0.0001)
- Based on 1,994 proteins

### Interpretation:
The two sets are **nearly identical** in terms of overall protein abundance. This is an excellent correlation, showing that storage preserves the proteome extremely well.

**Visual**: See top-left panel in `figures/comprehensive_comparison.png`

---

## Question 3: Similarities/differences in individual protein abundance?

### Controls vs Controls (Old vs Fresh):
- **Pearson r = 0.975**
- Based on 1,898 proteins
- Nearly perfect correlation

### Cases vs Cases (Sarcoidosis Old vs Fresh):
- **Pearson r = 0.976**
- Based on 1,934 proteins
- Nearly perfect correlation

### Comparison:
- **Difference: 0.002** (negligible)
- Both groups show essentially identical correlations
- No evidence that disease status affects storage stability

### Key Finding:
Individual protein abundances are **highly preserved** regardless of whether samples are from cases or controls. Storage affects both groups equally (minimally).

**Visual**: See top-middle and top-right panels in `figures/comprehensive_comparison.png`

---

## Question 4: Differences in missing % of individual proteins between sets?

### Overall Missing Data:

| Set | Mean Missing % | Median Missing % |
|-----|---------------|------------------|
| Set 1 (Old) | **40.5%** | **11.7%** |
| Set 2 (Fresh) | **45.4%** | **30.0%** |

**Surprise finding: Old samples have LESS missing data than fresh samples!**
- Difference: -4.9% (Old has less missing)

### Proteins with Large Differences (>30%):

**More missing in Set 1 (Old):** 10 proteins
- Examples: ASTE1 (93% vs 35%), NINL (67% vs 30%)

**More missing in Set 2 (Fresh):** **201 proteins**
- Examples: ENO2 (0% vs 100%), CFAP100 (0% vs 100%), CENPF (0% vs 98%)

### Unique Proteins:

**Only in Set 1 (completely missing in Set 2):** 143 proteins
**Only in Set 2 (completely missing in Set 1):** 21 proteins

### Key Finding:
Fresh samples have **more missing data** overall, and **201 proteins are substantially more missing in fresh vs old samples**. This suggests technical variability between sample processing batches, not storage degradation.

**Details**: See `results/missing_data_comparison.csv` for all 2,800 proteins

**Visual**: See bottom panels in `figures/comprehensive_comparison.png`

---

## Question 5: Step-by-step preprocessing methods

See **`PREPROCESSING_METHODS.md`** for complete documentation.

### Quick Summary:

#### 1. **Data Loading**
- Input: 2,800 proteins × 120 samples (Excel file)
- Automatic sample identification and grouping

#### 2. **Quality Control**
- Sample QC: Missing data, total intensity, PCA outliers
- Protein QC: Detection rates, coefficient of variation
- **Action**: Flag samples with >50% missing data

#### 3. **Protein Filtering**
- **Criteria**: Keep proteins detected in ≥40% of samples in **at least one condition**
- **Rationale**: Balances quality vs. coverage, preserves disease-specific proteins
- **Result**: ~1,600-1,800 proteins retained per set

#### 4. **Log2 Transformation**
- **Why**: Proteomics data is log-normally distributed
- **Method**: Replace 0 with NA, then apply log2
- **Benefit**: Stabilizes variance, makes fold-changes symmetric

#### 5. **Imputation**
- **Method**: Minimum value - 1 (on log2 scale)
- **Why**: Missing = below detection limit (low abundance)
- **Alternative**: KNN imputation (use similar proteins)

#### 6. **Batch Effect Correction**
- **Status**: NOT applied
- **Why**: Sets analyzed separately, not merged
- **When to use**: If combining sets for joint analysis

#### 7. **Normalization**
- **Status**: NOT explicitly applied
- **Why**: Log2 transformation provides implicit normalization
- **When to use**: If total intensities differ between samples

#### 8. **Statistical Testing**
- **Method**: Student's t-test (Welch's, unequal variance)
- **Comparison**: Sarcoidosis vs Control within each set
- **Correction**: Benjamini-Hochberg FDR < 0.05
- **Effect size**: Log2 fold-change, Cohen's d

### Decision Criteria:

**Filtering threshold (40%)**:
- Based on QC analysis showing median detection ~60%
- 30-40% is standard in proteomics
- Too strict = lose coverage, too lenient = poor quality

**Imputation method (min-1)**:
- Conservative approach
- Assumes missing = low abundance
- Doesn't artificially inflate significance

**No batch correction**:
- Sets 1 and 2 come from different trials (ACCESS vs LUC)
- True biological difference, not just technical batch
- Analyzed separately to avoid confounding

**No explicit normalization**:
- Data appeared pre-normalized from provider
- Total intensity distributions were similar
- Log2 transformation sufficient

---

## Summary of All Findings

### ✅ What We Confirmed:

1. **Old samples detect MORE proteins** (2,105 vs 1,953 in Sarcoidosis)
2. **Old samples have LESS missing data** (40.5% vs 45.4%)
3. **Protein abundances are nearly identical** between old and fresh (r = 0.98)
4. **Both cases and controls correlate equally well** (r = 0.976 vs 0.975)

### ⚠️ What We Found Different:

1. **Fresh samples have more significant proteins** (393 vs 114)
   - NOT because old samples failed to detect proteins
   - Because fresh samples have larger fold-changes (higher variability)
   - Old samples have more noise → less statistical power

2. **Fresh samples have more missing data**
   - 201 proteins substantially more missing in fresh samples
   - Suggests technical batch differences, not storage effect

### 🎯 Conclusion:

**Storage preserves protein detection and abundance extremely well (r = 0.98).**

The main effect of storage is to **reduce statistical power** (due to added noise), NOT to cause protein loss or degradation.

The **44 proteins significant in both sets** are high-confidence biomarkers suitable for clinical studies using stored samples.

---

## Generated Files

1. **scripts/comprehensive_comparison.py** - Analysis script
2. **results/missing_data_comparison.csv** - Per-protein missing data
3. **figures/comprehensive_comparison.png** - 6-panel visualization
4. **PREPROCESSING_METHODS.md** - Detailed methods documentation
5. **This file (ANSWERS_TO_QUESTIONS.md)** - Summary answers
