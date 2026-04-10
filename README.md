# Pain AUC & Opioid Persistence in Head and Neck Cancer

**ASTRO 2026 Abstract: Computational Companion Repository**

> *An analysis of acute-on-chronic pain and opioid persistence in a prospective head and neck cancer cohort treated with curative radiation therapy*

Md Mahin, Oliver Nkuku, Renji He, Clifton Dave Fuller, Amy Moreno, Saba Javed  
MD Anderson Cancer Center

---

## Overview

This repository contains the statistical analysis code for evaluating whether **acute pain burden** (measured as normalized Pain AUC during radiation therapy) predicts **opioid persistence** at subacute (FU1: 6–15 weeks) and chronic (FU2: 21–30 weeks) follow-up in head and neck cancer patients.

### Key Methods

| Component | Description |
|-----------|-------------|
| **Primary predictor** | Normalized Pain AUC (%), scaled per 10% increase |
| **Primary outcome** | Active opioid prescription at FU1 / FU2 (binary) |
| **FU1 model** | Full multivariable logistic regression |
| **FU2 model** | Prespecified parsimonious model with Firth bias-reduced penalized logistic regression |
| **Threshold analysis** | Grid search + optional Bayesian optimization for clinical cutpoints |
| **Sensitivity** | L2-regularized logistic regression; imputation vs. complete-case |

---

## Repository Structure

```
├── README.md
├── R/
│   ├── 01_data_preparation.R        # Cohort filtering, variable collapsing, AUC computation
│   ├── 02_firth_primary_models.R     # Firth penalized LR (logistf) — primary inference
│   ├── 03_l2_sensitivity.R           # L2 regularized LR (glmnet) — sensitivity analysis
│   ├── 04_threshold_search.R         # Grid search for AUC cutpoints
│   └── 05_bayesian_threshold.R       # Bayesian optimization for threshold (optional)
├── python/
│   ├── firth_logistic.py             # Pure NumPy Firth implementation
│   ├── run_models.py                 # Full pipeline: data prep → Firth → L2 → threshold
│   ├── bayesian_threshold.py         # Bayesian optimization with scikit-optimize
│   └── utils.py                      # AUC computation, variable collapsing helpers
├── matlab/
│   ├── firth_logistic.m              # Firth penalized LR (Newton/IRLS)
│   └── demo_firth.m                  # Example usage with simulated data
└── docs/
    └── statistical_notes.md          # Methodological rationale & reviewer-ready language
```

---

## Quick Start

### R (recommended for primary analysis)

```r
# Install dependencies
install.packages(c("logistf", "glmnet", "dplyr"))

# Run primary Firth models
source("R/01_data_preparation.R")
source("R/02_firth_primary_models.R")
```

### Python

```bash
pip install numpy scipy scikit-learn scikit-optimize
python python/run_models.py --data path/to/fu2_cohort.csv
```

### MATLAB

```matlab
% Add matlab/ to path, then:
demo_firth
```

---

## Methodological Notes

### Why Firth Penalized Logistic Regression?

The FU2 cohort (N≈106, 47 events) with multiple categorical covariates produces **quasi-complete separation** under standard MLE — coefficients diverge, ORs explode to millions, and standard errors become meaningless. Firth regression adds a Jeffreys-prior penalty:

```
l*(β) = l(β) + ½ log|I(β)|
```

This guarantees finite estimates, corrects small-sample bias, and preserves likelihood-based inference. See `docs/statistical_notes.md` for full rationale.

### Modeling Hierarchy

1. ✅ **Firth penalized LR** — primary inference (FU2)
2. ✅ **Parsimonious reduced model** — collapsed categories
3. ⚠️ **L2 regularized LR** — sensitivity/prediction only
4. ❌ **Standard MLE with separation** — not reportable

---

## Citation

If you use this code, please cite the ASTRO 2026 abstract:

> Mahin M, Nkuku O, He R, Fuller CD, Moreno A, Javed S. An analysis of acute-on-chronic pain and opioid persistence in a prospective head and neck cancer cohort treated with curative radiation therapy. ASTRO 2026.

---

## License

MIT License. See individual files for details.
