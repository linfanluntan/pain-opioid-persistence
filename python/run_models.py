"""
run_models.py
Full modeling pipeline: Data prep → Firth → L2 → Threshold search.
ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed
MD Anderson Cancer Center

PIPELINE (mirrors docs/technical_report.md structure):
    1. Data preparation with variable collapsing (§2, §4)
    2. Sparse factor check — identify separation risks (§11.2)
    3. Model B: Firth penalized LR — primary inference (§6, §11.4)
       Prespecified parsimonious FU2 model: Pain AUC + Smoking + Age + TxIntensity
    4. L2 sensitivity — directional consistency check (§11.5)
       [Friedman, Hastie & Tibshirani, J Stat Softw 2010]
    5. Threshold search — exploratory risk stratification (§7)
       Grid search with adjusted ORs at each cutpoint
    6. Clinical interpretation — OR per 10% scaling (§6.3)
       [Vittinghoff et al., 2012]

MODELING HIERARCHY:
    1 ✅ Firth penalized LR   [Firth 1993; Heinze & Schemper 2002]
    2 ✅ Parsimonious model    [Peduzzi et al. 1996; Vittinghoff & McCulloch 2007]
    3 ⚠️ L2 regularized LR    [Friedman et al. 2010] — sensitivity only
    4 ❌ Standard MLE          [Albert & Anderson 1984] — not reportable under separation

Usage:
    python run_models.py --data path/to/cohort.csv
    python run_models.py --demo   # run with simulated data
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

from firth_logistic import FirthLogisticRegression
from utils import (
    collapse_smoking, collapse_treatment, build_design_matrix,
    check_sparse_factors, compute_pain_auc,
)


def run_firth_model(X, y, term_names, model_name="Firth Model"):
    """Run Firth logistic regression and print results."""
    print(f"\n{'=' * 60}")
    print(f"  {model_name}")
    print(f"{'=' * 60}")
    print(f"N = {len(y)}, Events = {int(y.sum())}, Non-events = {int((1-y).sum())}")

    model = FirthLogisticRegression()
    model.fit(X, y, term_names=term_names)
    results = model.summary()

    # Classification accuracy
    preds = model.predict(X)
    acc = (preds == y).mean() * 100
    print(f"Classification accuracy: {acc:.1f}%")

    return model, results


def run_l2_model(X, y, term_names, model_name="L2 Sensitivity"):
    """Run L2 regularized logistic regression."""
    from sklearn.linear_model import LogisticRegression

    print(f"\n{'=' * 60}")
    print(f"  {model_name}")
    print(f"{'=' * 60}")

    # Remove intercept column if present (sklearn adds its own)
    if term_names[0] == "Intercept":
        X_sk = X[:, 1:]
        names_sk = term_names[1:]
    else:
        X_sk = X
        names_sk = term_names

    clf = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=5000)
    clf.fit(X_sk, y)

    print(f"{'Term':<25} {'Coef':>8} {'OR':>8}")
    print("-" * 45)
    print(f"{'Intercept':<25} {clf.intercept_[0]:>8.4f} {np.exp(clf.intercept_[0]):>8.3f}")
    for name, coef in zip(names_sk, clf.coef_[0]):
        print(f"{name:<25} {coef:>8.4f} {np.exp(coef):>8.3f}")

    acc = clf.score(X_sk, y) * 100
    print(f"\nClassification accuracy: {acc:.1f}%")

    return clf


def run_threshold_search(df, outcome_col, auc_col, adjust_X, adjust_names,
                          n_thresholds=40, model_name="Threshold Search"):
    """Grid search over AUC thresholds."""
    print(f"\n{'=' * 60}")
    print(f"  {model_name}")
    print(f"{'=' * 60}")

    auc = df[auc_col].values
    y = df[outcome_col].values.astype(float)

    lo, hi = np.percentile(auc, [5, 95])
    thresholds = np.linspace(lo, hi, n_thresholds)

    results = []
    for t in thresholds:
        above = (auc >= t).astype(float)
        n_low = int((above == 0).sum())
        n_high = int((above == 1).sum())

        if n_low < 10 or n_high < 10:
            continue

        X_t = np.column_stack([np.ones(len(y)), above, adjust_X])
        names_t = ["Intercept", "AUC_above"] + adjust_names

        try:
            model = FirthLogisticRegression(max_iter=50)
            model.fit(X_t, y, term_names=names_t)
            res = model.get_results()

            idx = 1  # AUC_above is second term
            results.append({
                "threshold": round(t, 1),
                "n_low": n_low,
                "n_high": n_high,
                "beta": round(res["beta"][idx], 3),
                "OR": round(res["OR"][idx], 2),
                "p_value": round(res["p_value"][idx], 4),
            })
        except Exception:
            continue

    results_df = pd.DataFrame(results)
    sig = results_df[results_df["p_value"] < 0.05]

    print(f"\nTested {len(results_df)} thresholds")
    print(f"Significant: {len(sig)}")

    if len(sig) > 0:
        lowest = sig.iloc[0]
        best = sig.loc[sig["OR"].abs().idxmax()]
        print(f"\nLowest significant: AUC >= {lowest['threshold']}% → OR = {lowest['OR']} (p = {lowest['p_value']})")
        print(f"Strongest: AUC >= {best['threshold']}% → OR = {best['OR']} (p = {best['p_value']})")

    return results_df


def run_demo():
    """Run full pipeline with simulated data."""
    print("\n" + "#" * 60)
    print("  DEMO: Simulated FU2-like Cohort")
    print("#" * 60)

    np.random.seed(42)
    n = 150

    # Simulate FU2-like data
    auc_pct = np.clip(np.random.normal(30, 15, n), 0, 70)
    smoking = np.random.binomial(1, 0.4, n)
    age60 = np.random.binomial(1, 0.5, n)
    tx = np.random.choice(["RT_only", "PORT", "Def_ChemoRT"], n, p=[0.2, 0.3, 0.5])

    # Outcome with true AUC effect
    logit = -1.5 + 0.04 * auc_pct + 0.5 * smoking - 0.3 * age60
    prob = 1 / (1 + np.exp(-logit))
    opioid = np.random.binomial(1, prob)

    df = pd.DataFrame({
        "Opioid_FU2": opioid,
        "AUCpain_pct": auc_pct,
        "AUC10": auc_pct / 10.0,
        "Smk_Ever": pd.Categorical(np.where(smoking, "Ever", "Never")),
        "Age60plus": pd.Categorical(np.where(age60, ">=60", "<60")),
        "TxIntensity": pd.Categorical(tx, categories=["RT_only", "PORT", "Def_ChemoRT"]),
    })

    print(f"\nCohort: N={n}, Events={opioid.sum()}, Non-events={(1-opioid).sum()}")

    # Sparse factor check
    print("\n--- Sparse Factor Check ---")
    warnings = check_sparse_factors(df, "Opioid_FU2", ["Smk_Ever", "Age60plus", "TxIntensity"])
    if len(warnings) > 0:
        print(warnings.to_string(index=False))
    else:
        print("No sparse factors detected.")

    # Build design matrix
    X, names = build_design_matrix(
        df,
        continuous_cols=["AUC10"],
        categorical_cols=["Smk_Ever", "Age60plus", "TxIntensity"],
        add_intercept=True,
    )
    y = df["Opioid_FU2"].values.astype(float)

    # ── Model B: Firth Primary ──
    firth_model, firth_res = run_firth_model(X, y, names, "Model B: FU2 Primary (Firth)")

    # ── L2 Sensitivity ──
    run_l2_model(X, y, names, "L2 Sensitivity Check")

    # ── Threshold Search ──
    # Get adjustment covariates (everything except intercept and AUC)
    adjust_idx = [i for i, n in enumerate(names) if n not in ("Intercept", "AUC10")]
    adjust_X = X[:, adjust_idx]
    adjust_names = [names[i] for i in adjust_idx]

    threshold_results = run_threshold_search(
        df, "Opioid_FU2", "AUCpain_pct", adjust_X, adjust_names,
        n_thresholds=25, model_name="AUC Threshold Search"
    )

    # ── Clinical Interpretation ──
    print(f"\n{'=' * 60}")
    print("  Clinical Interpretation")
    print(f"{'=' * 60}")
    auc_beta = firth_res["beta"][names.index("AUC10")]
    or_10pct = np.exp(auc_beta)
    print(f"OR per 10% increase in Pain AUC: {or_10pct:.3f}")
    print(f"Meaning: For every 10% increase in acute pain burden,")
    print(f"  odds of opioid persistence increase by {(or_10pct - 1) * 100:.1f}%.")
    print(f"\nOR per 1% increase: {np.exp(auc_beta / 10):.4f}")
    print(f"(Same model — just different scaling.)")


def main():
    parser = argparse.ArgumentParser(description="Pain AUC & Opioid Persistence Models")
    parser.add_argument("--data", type=str, help="Path to analytic CSV file")
    parser.add_argument("--demo", action="store_true", help="Run demo with simulated data")
    args = parser.parse_args()

    if args.demo or args.data is None:
        run_demo()
    else:
        print(f"Loading data from: {args.data}")
        print("Adapt the run_demo() pipeline to your column names.")
        print("The core functions (FirthLogisticRegression, utils) are ready to use.")

    print("\n✓ Pipeline complete.")


if __name__ == "__main__":
    main()
