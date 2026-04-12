# Statistical Notes: Methodological Rationale

**ASTRO 2026 — Pain AUC & Opioid Persistence in Head and Neck Cancer**

---

## 1. Terminology: "Persistence" not "Dependency"

The endpoint is **active opioid prescription at follow-up** (binary Y/N). This measures *opioid persistence* — ongoing prescribing behavior — not DSM-5 opioid use disorder, which requires evidence of impaired control, social impairment, risky use, and pharmacological indicators assessed over 12 months (American Psychiatric Association, 2013). Interventions to reduce persistence (earlier non-opioid transition) differ from OUD interventions (buprenorphine, behavioral therapy) (Dowell et al., 2022; Volkow & McLellan, 2016). Using "opioid dependency" without OUD diagnostic data would constitute overclaiming and would be flagged immediately by ASTRO reviewers.

---

## 2. Pain AUC: Why Continuous, Why Per 10%

### 2.1 Continuous > Binary

Dichotomizing continuous predictors discards information, reduces power, and introduces threshold-dependent bias (Royston et al., 2006; Altman & Royston, 2006). The dose–response relationship between pain burden and opioid persistence is more accurately captured by continuous specification. Threshold analysis (Section 6) is secondary and exploratory.

### 2.2 Scaling Per 10% vs. Per 1%

An OR of 1.05 per 1% sounds trivially small. Per 10%: OR = (1.05)¹⁰ ≈ 1.63, meaning 63% higher odds — clinically meaningful. The rescaling is linear:

```
β₁₀ = 10 × β₁
SE₁₀ = 10 × SE₁
Z = β₁₀ / SE₁₀ = β₁ / SE₁  (unchanged)
```

The model fit, p-values, and significance are identical (Vittinghoff et al., 2012). Apply consistently across FU1 and FU2.

### 2.3 Trapezoidal AUC Computation

The trapezoidal rule integrates both intensity and duration, adapted from pharmacokinetics (Rowland & Tozer, 2011) for symptom trajectories (Kahn et al., 2012). Normalization to percent of maximum possible AUC (pain = 10 sustained for entire window) produces an interpretable metric: 50% means average half-maximal burden.

---

## 3. Quasi-Complete Separation

### 3.1 Definition

Separation occurs when predictor values perfectly predict the outcome. In logistic regression, the MLE tries to assign P(Y=1) = 1 by pushing X^Tβ → +∞, requiring β → ±∞. The likelihood has a supremum but **no finite maximum** — the MLE does not exist in a mathematical sense (Albert & Anderson, 1984).

### 3.2 Signatures in Our Output

- ORs of 47,924 and 7,532,188
- Standard errors > 2,000
- Confidence intervals: 0 to 10²¹
- MLE optimizer failures (RuntimeError, LinAlgError, BFGS rescue)

### 3.3 Sources in This Study

| Sparse category | N in FU2 | Problem |
|----------------|----------|---------|
| Black race | 1 | Zero events or non-events in cell |
| Eye & Orbit | 1 | Perfect prediction |
| ICRT + Surgery Yes | 1 | Perfect prediction |
| Asian race | 1 | Zero events or non-events |

Additionally: FU1 fits 26 predictors with 110 non-events; FU2 fits 23 predictors with 59 non-events — violating the ~10 EPV rule (Peduzzi et al., 1996).

### 3.4 Why Pain AUC Is Unaffected

Separation is caused by sparse *categorical* variables with zero variation. Pain AUC is:

- **Continuous** — no discrete cells to separate
- **Present in all patients** — no zero-count categories
- **Broadly distributed** (0–65%+) — overlap between outcome groups
- **No perfect prediction** — both opioid users and non-users span the range

Therefore, separation-related instability does not undermine the primary finding.

---

## 4. Firth Penalized Logistic Regression

### 4.1 The Penalty

Firth (1993) modifies the log-likelihood by adding a Jeffreys invariant prior:

```
l*(β) = l(β) + ½ log|I(β)|
```

where I(β) is the Fisher information matrix. This is equivalent to modifying the score function to include a bias-correction term involving the hat-matrix diagonal.

### 4.2 Benefits

1. **Finite estimates** under complete or quasi-complete separation
2. **O(n⁻¹) bias reduction** — corrects small-sample bias of MLE coefficients
3. **Likelihood-based inference** — valid CIs and profile likelihood ratio tests preserved
4. **No tuning hyperparameter** — unlike L2/ridge, which depends on λ

### 4.3 Why Not L2 (Ridge)?

L2 regularization maximizes `l(β) − λΣβ²`. It stabilizes estimates but:

- Depends on the tuning parameter λ (chosen by cross-validation)
- Is designed for prediction, not classical hypothesis testing
- Wald p-values under penalty have uncertain theoretical grounding
- Effect sizes depend on penalty strength

Firth is specifically designed for clinical inference in sparse-data settings and is widely accepted in oncology journals (Heinze & Schemper, 2002; King & Zeng, 2001).

### 4.4 Implementation

- **R**: `logistf` package (Heinze et al., 2024) — `logistf(Y ~ X1 + X2 + ..., data = df)`
- **Python**: Custom Newton-Raphson with Jeffreys-prior hat-matrix adjustment (see `python/firth_logistic.py`)
- **MATLAB**: Custom modified Newton/IRLS (see `matlab/firth_logistic.m`)

### 4.5 Hierarchy of Approaches

| Rank | Method | Appropriate for | Issues |
|------|--------|----------------|--------|
| 1 ✅ | Firth penalized LR | Primary inference when separation exists | None (preferred) |
| 2 ✅ | Parsimonious reduced model | When full model exceeds EPV budget | Alters covariate structure |
| 3 ⚠️ | L2 regularized LR | Sensitivity analysis, prediction | Depends on λ; not classical inference |
| 4 ❌ | Standard MLE with separation | Never | Infinite coefficients, invalid p-values |

---

## 5. Events-Per-Variable (EPV) Budget

The ~10 EPV rule (Peduzzi et al., 1996) provides a practical constraint:

| Model | Events | Parameters | EPV | Status |
|-------|--------|-----------|-----|--------|
| FU1 full | 287 | 26 | 11.0 | ✅ Acceptable |
| FU2 full (original) | 47 | 23 | 2.0 | ❌ Severely underpowered |
| FU2 parsimonious | 47 | 5 | 9.4 | ✅ Near-optimal |

The FU2 parsimonious model includes: Pain AUC (continuous, per 10%), Smoking (ever/never), Age (≥60/<60), Treatment intensity (2–3 levels). Race, detailed site, and detailed regimen are excluded from FU2 inference (remain in descriptive tables).

Vittinghoff & McCulloch (2007) showed that the strict 10 EPV rule can be relaxed to 5–9 in some scenarios, but this assumes good predictor balance and no sparse cells — conditions not met in our FU2 data.

---

## 6. Threshold Search Methodology

### 6.1 Grid Search

Evenly spaced thresholds between the 5th–95th percentiles of the AUC distribution. At each cutpoint, dichotomize → fit multivariable logistic model → record adjusted OR and p-value.

### 6.2 Criticism and Mitigations

Grid search across multiple thresholds inflates type I error (Lausen & Schumacher, 1992). Mitigations:

1. **Frame as exploratory** risk-stratification, not hypothesis-testing
2. **Bootstrap internal validation** to assess threshold stability
3. **Bayesian optimization** with GP surrogate for principled selection (Snoek et al., 2012)
4. **Primary inference** uses the continuous model — thresholds are secondary

### 6.3 Binary Thresholding Is Information-Lossy

Dichotomizing throws away information and reduces power (Royston et al., 2006). Use thresholds only when clinical protocols require a categorical decision rule.

---

## 7. Imputation

### 7.1 Method

Linear interpolation between nearest non-missing neighbors. Boundary cases use half of nearest neighbor.

### 7.2 Volume

- 80.9% of patients had ≥1 imputed value
- 51.3% had Week 7 imputed
- FU1 was never imputed (required anchor)

### 7.3 Concerns

The trapezoidal AUC is partly synthetic for a majority of patients. Simple averaging may attenuate variability and bias toward the mean (Little & Rubin, 2019). Recommended: sensitivity analysis without imputation, and consideration of MICE for future analyses (van Buuren & Groothuis-Oudshoorn, 2011; Sterne et al., 2009).

---

## 8. Attrition Bias

The sequential filtering from 1,118 → 397 → 106 represents:

- **64% drop** before FU1
- **90% drop** before FU2

Patients completing extended follow-up may differ systematically — better function, greater engagement, fewer comorbidities, or conversely more persistent symptoms (Mercieca-Bebber et al., 2017; Bottomley et al., 2016). Low power simultaneously increases false-negative and false-positive rates (Button et al., 2013).

**ASTRO-ready sentence**: "Substantial attrition from the initial cohort resulted in a smaller FU2 sample size, reducing statistical power and leading to wider confidence intervals and less stable multivariable estimates compared to FU1."

---

*See [`references.md`](references.md) for the complete bibliography.*
