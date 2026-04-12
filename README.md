# Pain AUC & Opioid Persistence in Head and Neck Cancer

**ASTRO 2026 — Computational Companion Repository**

> *An analysis of acute-on-chronic pain and opioid persistence in a prospective
> head and neck cancer cohort treated with curative radiation therapy*

Md Mahin, Oliver Nkuku, Renji He, Clifton Dave Fuller, Amy Moreno, Saba Javed  
Department of Radiation Oncology, The University of Texas MD Anderson Cancer Center

---

## Scientific Context

Head and neck cancer (HNC) patients undergoing curative radiation therapy (RT)
frequently develop severe acute mucositis, dysphagia, and pain that require
opioid analgesics (Elting et al., 2007; Trotti et al., 2003). A substantial
proportion of patients initiated on opioids during the acute treatment phase
continue to fill prescriptions well beyond expected toxicity resolution — a
phenomenon termed **opioid persistence** (Patel et al., 2021). This is clinically
distinct from opioid use disorder as defined by the DSM-5 (APA, 2013) but
carries significant risks including dose escalation, physiologic dependence, and
adverse health consequences (Dowell et al., 2022).

This repository implements the complete statistical analysis pipeline for a study
evaluating whether **total acute pain burden** — quantified as normalized Pain
Area Under the Curve (AUC) via weekly MDASI-HN assessments (Rosenthal et al.,
2007; Cleeland et al., 2000) — independently predicts opioid persistence at
subacute (FU1: 6–15 weeks post-RT) and chronic (FU2: 21–30 weeks post-RT)
follow-up. The biological rationale rests on central sensitization: prolonged
nociceptive input during tissue injury may induce maladaptive neuroplastic
changes that persist beyond peripheral injury resolution (Woolf, 2011;
Latremoliere & Woolf, 2009).

---

## Study Design Summary

| Component | Detail |
|-----------|--------|
| **Population** | 1,118 MOSAIC-eligible HNC patients, MD Anderson, 02/2021–06/2024 |
| **Primary predictor** | Normalized Pain AUC (%), trapezoidal rule, scaled per 10% increase |
| **Primary outcome** | Active opioid prescription within ±7 days of FU pain assessment |
| **FU1 cohort** | N = 397 (35.5%), 287 events, median 66 days post-RT |
| **FU2 cohort** | N = 106 (9.5%), 47 events, median 180 days post-RT |
| **FU1 model** | Full multivariable logistic regression (26 parameters) |
| **FU2 model** | Prespecified parsimonious Firth penalized LR (4–6 parameters) |
| **Threshold analysis** | Grid search + optional Bayesian optimization |
| **Sensitivity** | L2-regularized LR; imputation vs. complete-case |

### Key Results

- **Pain AUC** predicts opioid persistence at both timepoints: OR ≈ 1.63 per 10% increase (~63% higher odds)
- **Smoking** doubles subacute opioid odds (FU1 OR 2.15); shows temporal amplification at FU2
- **Age ≥60** is protective at chronic phase (FU2 OR 0.28)
- **Threshold divergence**: FU1 optimal cutpoint 53.2% (OR 7.8); FU2 optimal 18.3% (OR 7.0) — biologically coherent

---

## Repository Structure

```
├── README.md
├── R/
│   ├── 01_data_preparation.R          # Cohort filtering, AUC computation, variable collapsing
│   ├── 02_firth_primary_models.R       # Firth penalized LR — primary inference
│   ├── 03_l2_sensitivity.R             # L2 ridge (glmnet) — sensitivity analysis
│   ├── 04_threshold_search.R           # Grid search with bootstrap validation
│   └── 05_bayesian_threshold.R         # GP-based Bayesian optimization (optional)
├── python/
│   ├── firth_logistic.py               # Pure NumPy Firth (Jeffreys prior, Newton-Raphson)
│   ├── run_models.py                   # Full pipeline: Firth → L2 → threshold → interpretation
│   ├── bayesian_threshold.py           # Bayesian threshold optimization (scikit-optimize)
│   ├── utils.py                        # AUC trapezoidal rule, imputation, collapsing, design matrices
│   └── requirements.txt
├── matlab/
│   ├── firth_logistic.m                # Self-contained Firth (modified Newton/IRLS)
│   └── demo_firth.m                    # Full demo with separation detection
└── docs/
    ├── technical_report.md             # Complete technical report (study design through limitations)
    ├── statistical_notes.md            # Separation, penalization, EPV, scaling — with references
    ├── clinical_interpretation.md      # Translational framing & reviewer-ready language
    └── references.md                   # Full bibliography (46 references)
```

---

## Quick Start

### R (recommended for primary analysis)

```r
install.packages(c("logistf", "glmnet", "dplyr"))
source("R/01_data_preparation.R")
source("R/02_firth_primary_models.R")
source("R/04_threshold_search.R")
```

### Python

```bash
pip install -r python/requirements.txt
python python/run_models.py --demo
```

### MATLAB

```matlab
addpath('matlab/'); demo_firth
```

---

## Methodological Framework

### Pain AUC Computation

Pain AUC is computed via the trapezoidal rule on weekly MDASI-HN pain scores (0–10 NRS) and normalized to a percentage of maximum possible burden. This captures both intensity and duration (Kahn et al., 2012). Retained as a continuous predictor to preserve information (Royston et al., 2006). Reported per 10% increase for interpretability (Vittinghoff et al., 2012).

### Why Firth Penalized Logistic Regression?

The FU2 cohort produces quasi-complete separation under standard MLE. Firth regression adds a Jeffreys invariant prior (Firth, 1993; Heinze & Schemper, 2002), guaranteeing finite estimates without a tuning hyperparameter — the standard for rare-event clinical studies (King & Zeng, 2001).

### Outcome: Opioid Persistence (Not Dependency)

Active opioid prescription at follow-up — not DSM-5 OUD (APA, 2013; Volkow & McLellan, 2016).

See [`docs/`](docs/) for full technical report, statistical notes, and bibliography.

---

## Citation

> Mahin M, Nkuku O, He R, Fuller CD, Moreno A, Javed S. An analysis of acute-on-chronic pain and opioid persistence in a prospective head and neck cancer cohort treated with curative radiation therapy. *ASTRO Annual Meeting*; 2026.

## License

MIT License. See [LICENSE](LICENSE).
