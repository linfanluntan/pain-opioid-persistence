"""
bayesian_threshold.py
Bayesian Optimization for Pain AUC threshold selection.
ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed

Uses Gaussian Process surrogate to maximize cross-validated predictive
performance over the threshold space.

Usage:
    python bayesian_threshold.py --demo

REQUIRES: pip install scikit-optimize
"""

import argparse
import numpy as np
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression

try:
    from skopt import gp_minimize
    from skopt.space import Real
    HAS_SKOPT = True
except ImportError:
    HAS_SKOPT = False


def bayesian_threshold_search(
    X_auc: np.ndarray,
    X_other: np.ndarray,
    y: np.ndarray,
    auc_range: tuple = None,
    n_calls: int = 30,
    n_folds: int = 5,
    random_state: int = 42,
):
    """
    Find optimal AUC threshold using Bayesian Optimization with GP surrogate.

    Parameters
    ----------
    X_auc : ndarray of shape (n,)
        Continuous AUC values (%).
    X_other : ndarray of shape (n, k)
        Other covariates (already dummy-coded).
    y : ndarray of shape (n,)
        Binary outcome (0/1).
    auc_range : tuple (low, high), optional
        Range to search. Default: 5th–95th percentile.
    n_calls : int
        Total evaluations (init + optimization steps).
    n_folds : int
        Cross-validation folds.
    random_state : int
        Random seed for reproducibility.

    Returns
    -------
    dict with: optimal_threshold, best_cv_auc, history
    """
    if not HAS_SKOPT:
        print("scikit-optimize not installed. Falling back to manual search.")
        return manual_search(X_auc, X_other, y, auc_range, n_folds)

    if auc_range is None:
        auc_range = (np.percentile(X_auc, 5), np.percentile(X_auc, 95))

    print(f"Bayesian optimization over [{auc_range[0]:.1f}, {auc_range[1]:.1f}]")
    print(f"Using {n_calls} evaluations with {n_folds}-fold CV\n")

    def objective(params):
        threshold = params[0]

        auc_binary = (X_auc >= threshold).astype(float).reshape(-1, 1)
        X_model = np.hstack([auc_binary, X_other])

        model = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=1000)
        scores = cross_val_score(model, X_model, y, cv=n_folds, scoring="roc_auc")

        return -scores.mean()  # minimize negative AUC

    space = [Real(auc_range[0], auc_range[1], name="threshold")]

    result = gp_minimize(
        objective,
        space,
        n_calls=n_calls,
        random_state=random_state,
        verbose=False,
    )

    optimal = result.x[0]
    best_auc = -result.fun

    print(f"\n{'=' * 50}")
    print(f"  Bayesian Optimization Results")
    print(f"{'=' * 50}")
    print(f"Optimal threshold: {optimal:.1f}%")
    print(f"Best CV AUC-ROC:   {best_auc:.4f}")

    # Posterior uncertainty from evaluations
    xs = np.array([x[0] for x in result.x_iters])
    ys = -np.array(result.func_vals)

    top_mask = ys >= np.percentile(ys, 80)
    top_thresholds = xs[top_mask]

    print(f"\nTop-performing threshold range: "
          f"{top_thresholds.min():.1f} – {top_thresholds.max():.1f}")
    print(f"Mean of top 20%: {top_thresholds.mean():.1f}")

    return {
        "optimal_threshold": optimal,
        "best_cv_auc": best_auc,
        "all_thresholds": xs,
        "all_scores": ys,
        "skopt_result": result,
    }


def manual_search(X_auc, X_other, y, auc_range=None, n_folds=5, n_eval=40):
    """Fallback: stratified random search without scikit-optimize."""
    if auc_range is None:
        auc_range = (np.percentile(X_auc, 5), np.percentile(X_auc, 95))

    print(f"Manual stratified search over [{auc_range[0]:.1f}, {auc_range[1]:.1f}]")
    print(f"Evaluating {n_eval} thresholds with {n_folds}-fold CV\n")

    thresholds = np.linspace(auc_range[0], auc_range[1], n_eval)
    scores = []

    for t in thresholds:
        auc_binary = (X_auc >= t).astype(float).reshape(-1, 1)
        n_low = (auc_binary == 0).sum()
        n_high = (auc_binary == 1).sum()

        if n_low < 10 or n_high < 10:
            scores.append(0.5)
            continue

        X_model = np.hstack([auc_binary, X_other])

        try:
            model = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=1000)
            cv = cross_val_score(model, X_model, y, cv=n_folds, scoring="roc_auc")
            scores.append(cv.mean())
        except Exception:
            scores.append(0.5)

    scores = np.array(scores)
    best_idx = np.argmax(scores)

    print(f"{'=' * 50}")
    print(f"Best threshold: {thresholds[best_idx]:.1f}%")
    print(f"Best CV AUC-ROC: {scores[best_idx]:.4f}")

    top_mask = scores >= np.percentile(scores, 80)
    print(f"Top range: {thresholds[top_mask].min():.1f} – {thresholds[top_mask].max():.1f}")

    return {
        "optimal_threshold": thresholds[best_idx],
        "best_cv_auc": scores[best_idx],
        "all_thresholds": thresholds,
        "all_scores": scores,
    }


def demo():
    """Demo with simulated data."""
    print("\n" + "#" * 50)
    print("  Bayesian Threshold Optimization — Demo")
    print("#" * 50 + "\n")

    np.random.seed(42)
    n = 300

    auc = np.clip(np.random.normal(30, 15, n), 0, 70)
    smoking = np.random.binomial(1, 0.4, n).reshape(-1, 1).astype(float)
    age60 = np.random.binomial(1, 0.5, n).reshape(-1, 1).astype(float)
    X_other = np.hstack([smoking, age60])

    logit = -1.5 + 0.04 * auc + 0.5 * smoking.ravel() - 0.3 * age60.ravel()
    prob = 1 / (1 + np.exp(-logit))
    y = np.random.binomial(1, prob)

    result = bayesian_threshold_search(
        X_auc=auc,
        X_other=X_other,
        y=y,
        n_calls=20,
    )

    print(f"\n✓ Demo complete. Optimal threshold: {result['optimal_threshold']:.1f}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()

    if args.demo:
        demo()
    else:
        print("Use --demo to run with simulated data.")
        print("Import bayesian_threshold_search() for use with your data.")
