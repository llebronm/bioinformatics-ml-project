# End-to-End BRCA RNA-Seq Classification Pipeline

This document outlines the leak-free bioinformatics and machine learning workflow for classifying Breast Invasive Carcinoma (`Cancer` vs. `Healthy`) using bulk RNA-seq count data.

---

## Visual Architecture Overview

```text
[Raw Merged Counts: 232 Samples x 59,427 Genes]
                     │
                     ▼
   ┌───────────────────────────────────┐
   │  PHASE 1: Unsupervised QC &       │
   │  Low-Expression Filtering         │
   │  (Safe before splitting)          │
   └─────────────────┬─────────────────┘
                     │
       [232 Samples x 27,249 Genes]
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼ (Branch A: Biology)     ▼ (Branch B: Machine Learning)
┌───────────────────┐     ┌───────────────────────────────────────┐
│  Full-Cohort DEA  │     │  PHASE 2: Stratified Train/Test Split │
│  (All 232 samples)│     └───────────┬───────────────────────────┘
│  • Volcano Plots  │                 │
│  • GO / Pathways  │        ┌────────┴────────┐
│  (DO NOT USE FOR  │        ▼                 ▼
│   ML TRAINING!)   │  [Train: 80%]     [Locked Test: 20%]
└───────────────────┘  (~186 samples)    (~46 samples)
                             │                 │
                             ▼                 │ (Locked in vault
                   ┌───────────────────┐       │  until Phase 4)
                   │  PHASE 3: 5-Fold  │       │
                   │  Cross-Validation │       │
                   │  (Repeat 5 times) │       │
                   └─────────┬─────────┘       │
                             │                 │
            ┌────────────────┴────────────────┐│
            ▼                                 ▼│
  [4 Training Folds: ~149]          [1 Val Fold: ~37]
  1. Fit PyDESeq2 DEA               1. Subset to top N genes
  2. Select top N DEGs              2. Apply Train VST/Scaler
  3. Fit VST / Scaler               3. Predict & Score
  4. Train ML Classifier            (Average 5 fold scores)
            │                                  │
            └────────────────┬─────────────────┘
                             │
                             ▼
                   ┌───────────────────┐
                   │  PHASE 4: Final   │
                   │  Model & Testing  │
                   └─────────┬─────────┘
                             │
     1. Fit PyDESeq2 on ALL 80% Train (~186 samples)
     2. Select final top N genes & train best model
     3. Unlock 20% Test (~46 samples) ────────┘
     4. Evaluate ONCE for final metrics
```

---

## Step-by-Step Execution Plan

### Phase 1: Data Preprocessing & Quality Control (`02_data_qc.ipynb`)
*Note: Unsupervised operations that do not look at `Cancer` vs. `Healthy` target labels can be safely performed on the full dataset.*

1. **Sample & Metadata Alignment:** Verify 1-to-1 matching between sample IDs in the count matrix and clinical metadata (`Condition`).
2. **Sample-Level QC:** 
   * Inspect library sizes (total read counts per sample) to catch sequencing failures.
   * Inspect zero-count sparsity per sample (flagging extreme outliers with abnormally low sparsity, e.g., <15%, for potential removal).
3. **Gene-Level Unsupervised Filtering:**
   * Drop genes with 0 counts across all samples.
   * Retain genes with $\ge 10$ counts in at least 10% of the cohort ($\ge 23$ samples), reducing the feature space from `59,427` to `27,249` genes.
4. **Export Clean Raw Matrix:** Save `raw_merged_filtered.csv` (integer counts are required as input for `PyDESeq2`).

---

### Branch A: Biological Discovery (`03_dea_full_cohort.ipynb`)
*Purpose: Purely for biological interpretation, biomarker exploration, and visualization. **Never pass these selected genes into your ML evaluation loop.***

1. Run `PyDESeq2` (`Cancer` vs. `Healthy`) on all 232 samples.
2. Filter significant Differentially Expressed Genes (DEGs) using standard thresholds (e.g., $|\log_2\text{FoldChange}| \ge 1.0$ and $p_{\text{adj}} \le 0.01$).
3. Generate publication figures: Volcano plots, MA plots, and hierarchical clustering heatmaps.
4. Perform Gene Ontology (GO) and KEGG pathway enrichment analysis on the resulting DEG list.

---

### Branch B: Leak-Free Machine Learning Pipeline (`04_ml_pipeline.ipynb`)

#### Phase 2: The Holdout Split
1. **Stratified Split:** Divide `raw_merged_filtered.csv` into **80% Training (~186 samples)** and **20% Holdout Test (~46 samples)** using `StratifiedShuffleSplit` or `train_test_split(..., stratify=y)` to preserve the exact `Cancer`/`Healthy` ratio in both sets.
2. **Lock the Test Set:** Do not touch the 20% test set during feature selection, normalization fitting, or hyperparameter tuning.

#### Phase 3: 5-Fold Cross-Validation (Model & Feature Selection)
Split the **80% Training Set (~186 samples)** into 5 stratified folds (~37 samples per fold). Loop through each fold $k \in \{1, 2, 3, 4, 5\}$:

1. **Partition:** Assign 4 folds (~149 samples) as `Train_Fold` and 1 fold (~37 samples) as `Val_Fold`.
2. **In-Fold DEA (Feature Selection):**
   * Initialize and fit `PyDESeq2` **exclusively on `Train_Fold` raw counts**.
   * Rank genes by statistical significance ($p_{\text{adj}}$) and effect size ($|\log_2\text{FoldChange}|$).
   * Select the top $N$ genes (e.g., $N \in \{50, 100, 250, 500\}$).
3. **In-Fold Transformation & Scaling:**
   * Apply Variance Stabilizing Transformation (VST) or $\log_2(\text{CPM} + 1)$ to `Train_Fold` and `Val_Fold` independently (or using dispersion trends learned from `Train_Fold`).
   * Subset both matrices to the exact $N$ genes selected in Step 2.
   * Fit a `StandardScaler` on `Train_Fold` and transform `Val_Fold`.
4. **Train & Validate Classifier:**
   * Train candidate models (e.g., SVM-RBF, Random Forest, Logistic Regression/LASSO, MLP) on `Train_Fold`.
   * Predict on `Val_Fold` and record **Accuracy, Sensitivity, Specificity, F1-Score, and ROC-AUC**.
5. **Aggregate Results:** Average the metrics across all 5 folds to identify the best classifier architecture, hyperparameters, and optimal gene count $N$.

#### Phase 4: Final Model Training & Unseen Test Evaluation
Once the winning model configuration is locked in from Phase 3:

1. **Final Feature Selection:** Run `PyDESeq2` one final time on the **entire 80% Training Set (~186 samples)** to extract the final consensus list of top $N$ genes.
2. **Final Training:** Transform and scale the 80% Training Set (subsetted to the final $N$ genes) and fit the winning model architecture.
3. **Final Evaluation:** 
   * Subset the **Locked 20% Test Set (~46 samples)** to those exact $N$ genes.
   * Apply the transformation and scaler learned from the 80% Training Set.
   * Generate predictions **once** and report final metrics + Confusion Matrix + ROC Curve.