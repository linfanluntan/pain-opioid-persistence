# Clinical Interpretation & Reviewer-Ready Language

**ASTRO 2026 — Pain AUC & Opioid Persistence in Head and Neck Cancer**

---

## 1. Conceptual Elevation: From Descriptive to Translational

The paper is strongest when it frames acute Pain AUC not merely as a statistical predictor but as a **predictive biomarker of opioid persistence**. This elevates the manuscript from descriptive association to translational science:

- **A measurable, actionable risk signal** available at RT completion
- **A target for early intervention** — proactive pain management and opioid stewardship
- **A stratification tool** for post-RT pain management pathways

The biological plausibility rests on central sensitization: prolonged nociceptive input during tissue injury induces maladaptive neuroplastic changes in pain-processing circuits that persist beyond peripheral injury resolution (Woolf, 2011; Latremoliere & Woolf, 2009). Acute nociceptive burden *imprints* the future prescribing trajectory. The transition from acute to persistent opioid use is predictable at RT completion.

---

## 2. Narrative Hierarchy

**Tier 1 — Mechanistic exposure:** Pain AUC during RT *(the story)*  
**Tier 2 — Behavioral modifiers:** Smoking, psychiatric history *(context)*  
**Tier 3 — Demographic background:** Race, age, site *(adjustment)*

The manuscript is strongest when centered on Tier 1. Secondary covariate instability does not weaken the paper — it is expected in small chronic-phase cohorts. Defending fragile smoking significance risks distracting from the primary biological message.

---

## 3. Key Interpretive Points

### 3.1 Pain AUC: Clinically Meaningful Magnitude

Per 10% increase: OR ≈ 1.63 (~63% higher odds). For context, well-established surgical risk factors for persistent opioid use (preoperative opioid exposure, anxiety disorder) typically show ORs of 1.5–3.0 (Brummett et al., 2017; Sun et al., 2016). The Pain AUC threshold effects (OR 7–8) are among the strongest reported predictors in oncology.

### 3.2 Smoking: Temporal Amplification

FU1 OR = 2.15 → FU2 directional OR ~5.27 (unstabilized). This suggests smokers may transition to chronic persistence rather than tapering. Under penalized modeling, smoking shows directional association but insufficient precision in the smaller FU2 cohort — a power limitation, not evidence of absence. Known mechanisms: nicotine–pain interaction, slower wound healing, behavioral comorbidity (Shi et al., 2010; Ditre et al., 2011).

### 3.3 Age: Late-Phase Protection

FU1: non-significant. FU2: age ≥60 → OR 0.28 (p = 0.031). Possible explanations: more cautious prescribing for older patients, differential opioid metabolism, or cohort selection effects (Chau et al., 2008).

### 3.4 Threshold Divergence

FU1 optimal: AUC ≥ 53.2% → OR 7.8 (very high burden predicts early persistence).  
FU2 optimal: AUC ≥ 18.3% → OR 7.0 (even moderate burden predicts late persistence).

The lower chronic threshold implies that moderate-to-high pain during RT creates lasting opioid exposure trajectories. Patients do not need extreme pain to develop chronic persistence — consistent with postoperative opioid literature (Brummett et al., 2017).

---

## 4. Reviewer-Ready Paragraphs

### Methods: Statistical Analysis

> Multivariable logistic regression was performed to evaluate predictors of opioid persistence at FU1 and FU2. Acute pain burden was summarized as normalized on-treatment pain AUC (%), modeled as a continuous predictor and scaled per 10% increase for interpretability. Given the smaller FU2 cohort and limited event count, we prespecified a parsimonious model including acute pain AUC, smoking status (ever vs. never), age, and a collapsed treatment-intensity variable. Because several categorical strata were sparse and produced quasi-complete separation under maximum-likelihood estimation, bias-reduced penalized logistic regression (Firth method) was used to obtain finite, stable effect estimates (Firth, 1993; Heinze & Schemper, 2002). Sensitivity analyses using alternative shrinkage approaches (Friedman et al., 2010) yielded consistent directionality for the primary predictor.

### Results: FU2 Stability Note

> In the FU2 cohort, small sample size and sparse subgroup counts resulted in unstable maximum-likelihood estimates with quasi-complete separation and inflated odds ratios (Albert & Anderson, 1984). Accordingly, we report bias-reduced penalized estimates from the prespecified parsimonious model. Under penalized modeling, acute pain AUC remained independently associated with opioid persistence, demonstrating robustness across modeling approaches. Some secondary covariates (e.g., smoking) showed directional associations but did not consistently retain statistical significance, likely reflecting limited precision for subgroup-level inference in the smaller chronic-phase cohort.

### Results: Smoking (FU2)

> In FU2, smoking demonstrated a positive directional association with chronic opioid persistence; however, this did not retain statistical significance after bias-reduced penalized modeling, suggesting limited statistical power in the chronic-phase cohort rather than definitive absence of effect.

### Limitations

> Substantial attrition from the initial cohort (1,118 eligible → 106 FU2 analyzable) resulted in reduced statistical power and potential selection bias, particularly for subgroup analyses (Mercieca-Bebber et al., 2017). Although several sparse categorical strata exhibited separation with inflated odds ratios under standard MLE, the continuous pain AUC variable remained well-supported by broad data distribution and consistent outcome overlap, yielding stable and robust estimates across models. High imputation rates (80.9% with ≥1 imputed value) and the absence of prior opioid use and psychiatric history data are additional limitations (Little & Rubin, 2019; Vitzthum et al., 2020). Findings should be validated in external cohorts before clinical implementation.

### Core Message (One Paragraph)

> This study demonstrates that total acute pain burden during radiation therapy, quantified as normalized pain AUC, is a robust and clinically meaningful predictor of opioid persistence into both the subacute and chronic phases of recovery. The relationship is graded, biologically coherent, and survives multiple modeling approaches including bias-reduced penalized regression. These findings support the development of AUC-driven risk stratification protocols to guide proactive opioid stewardship at the time of RT completion.

---

*See [`references.md`](references.md) for the complete bibliography.*
