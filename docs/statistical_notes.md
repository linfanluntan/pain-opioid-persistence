# Statistical Notes: Methodological Rationale

**ASTRO 2026 — Pain AUC & Opioid Persistence in Head and Neck Cancer**

---

## 1. Terminology: "Persistence" not "Dependency"

The endpoint is **active opioid prescription at follow-up** (binary: Y/N). This measures **opioid persistence**, not DSM-5 opioid use disorder or clinical dependency. The distinction is important for ASTRO reviewers — overclaiming "dependency" without OUD diagnostic data invites immediate criticism.

## 2. Why Firth Penalized Logistic Regression

### The Problem: Separation

The FU2 cohort (N≈106, 47 events) with multiple categorical covariates produces **quasi-complete separation** under standard MLE:

- ORs of 47,924 and 7,532,188
- Standard errors in the thousands
- Confidence intervals spanning 0 to 10²¹
- MLE optimizer failures requiring BFGS rescue

This occurs because sparse categories (e.g., Eye & Orbit n=1, Black race n=1 in FU2) perfectly predict the outcome. The likelihood has a supremum but **no finite maximum** — mathematically, the MLE does not exist.

### The Solution: Firth's Penalty

Firth regression adds a Jeffreys-prior penalty to the log-likelihood:

```
l*(β) = l(β) + ½ log|I(β)|
```

This guarantees:
1. **Finite estimates** even under complete separation
2. **Reduced small-sample bias** (O(n⁻¹) bias correction)
3. **Valid likelihood-based inference** with interpretable ORs and CIs

### Why Not L2 Regularization?

L2 (ridge) regression also prevents divergence, but:
- Depends on a tuning hyperparameter (λ)
- Is designed for prediction, not classical inference
- Standard p-values are not well-defined under penalization
- Effect sizes depend on penalty strength

Firth is specifically designed for clinical inference in sparse-data settings and is widely accepted in oncology journals.

### Hierarchy of Approaches

| Rank | Method | Use |
|------|--------|-----|
| 1 ✅ | Firth penalized LR | Primary inference (FU2) |
| 2 ✅ | Parsimonious reduced model | Acceptable alternative |
| 3 ⚠️ | L2 regularized LR | Sensitivity/prediction only |
| 4 ❌ | Standard MLE with separation | Not reportable |

## 3. Why Pain AUC Survives Instability

Sparse categorical variables create separation because they have few patients and zero variation in one outcome group. Pain AUC is fundamentally different:

- **Continuous** — no discrete cells to separate
- **Present in all patients** — no zero-count categories
- **Broadly distributed** (0–65%+) — overlap between outcome groups
- **No perfect prediction** — both opioid users and non-users span the AUC range

Therefore, separation-related instability does not affect the primary Pain AUC finding.

## 4. FU2 Parsimonious Model Design

With 47 events, the ~10 events-per-variable rule allows 4–5 parameters:

| Variable | Coding | Parameters |
|----------|--------|------------|
| Pain AUC (acute) | Continuous, per 10% | 1 |
| Smoking | Ever vs. Never | 1 |
| Age | ≥60 vs. <60 (or continuous) | 1 |
| Treatment intensity | 2–3 levels | 1–2 |
| **Total** | | **4–5** |

Variables **not included** in FU2 inference (due to sparse strata):
- Detailed race categories
- Detailed primary site (8–10 levels)
- Detailed treatment regimen combinations

These remain in descriptive tables.

## 5. Scaling: Per 10% vs. Per 1%

Reporting OR per 1% AUC increase gives OR ≈ 1.05, which sounds small. Rescaling to **per 10% increase** yields OR ≈ 1.63 (same model, same p-value). The rescaling is a linear transformation that does not change inference:

```
β₁₀ = 10 × β₁
SE₁₀ = 10 × SE₁
Z = β₁₀ / SE₁₀ = β₁ / SE₁  (unchanged)
```

## 6. Threshold Search Framing

Grid search is exploratory risk-stratification, not primary inference. The continuous AUC model is statistically stronger. Thresholds serve clinical decision-making ("above X%, escalate pain management"). Bootstrap validation or Bayesian optimization with CV AUC can strengthen the threshold's defensibility.

## 7. Attrition & Imputation

- **Attrition**: 1,118 → 397 (FU1) → 106 (FU2) represents 64% and 90% drop-off
- **Imputation**: 80.9% had some imputed values; 51% had W7 imputed
- Both should be acknowledged as limitations
- Sensitivity analysis comparing imputed vs. complete-case results recommended

---

*Document version: April 2026*
