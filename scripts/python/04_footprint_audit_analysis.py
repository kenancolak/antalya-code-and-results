#!/usr/bin/env python3
"""Reproduce footprint-audit checks for the final 40-cell morphology population."""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "derived" / "morphology_footprint_audit_scl_p95_40.csv"
AUDIT47 = ROOT / "data" / "derived" / "manual_footprint_audit_47.csv"
OUT = ROOT / "results" / "reproduced_footprint_audit"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 42
BOOTSTRAPS = 10_000
OUTCOME = "delta_LST_C"


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=1)


def correlation(frame: pd.DataFrame, predictor: str, label: str) -> dict:
    clean = frame[[predictor, OUTCOME]].dropna()
    x = clean[predictor].to_numpy(float)
    y = clean[OUTCOME].to_numpy(float)
    observed = spearmanr(x, y)
    rng = np.random.default_rng(SEED)
    boot = np.full(BOOTSTRAPS, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for index in range(BOOTSTRAPS):
            sample = rng.integers(0, len(clean), len(clean))
            boot[index] = spearmanr(x[sample], y[sample]).statistic
    valid = boot[np.isfinite(boot)]
    low, high = np.quantile(valid, [0.025, 0.975])
    return {
        "measure": label,
        "n": len(clean),
        "rho": observed.statistic,
        "p_value": observed.pvalue,
        "bootstrap_ci_low": low,
        "bootstrap_ci_high": high,
        "valid_bootstrap_repetitions": len(valid),
    }


def fit_hc3(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    predictors = [
        "BCR_baseline_proxy_pct",
        "delta_BCR_pp",
        "p75_storeys",
        "IVS_median_2017_2025",
    ]
    data = frame[[OUTCOME, *predictors]].dropna().copy()
    x = np.column_stack([np.ones(len(data)), data[predictors].apply(zscore).to_numpy(float)])
    y = zscore(data[OUTCOME]).to_numpy(float)
    names = ["constant", *predictors]
    beta = np.linalg.pinv(x) @ y
    residuals = y - x @ beta
    n, k = x.shape
    df = n - k
    xtx_inv = np.linalg.pinv(x.T @ x)
    leverage = np.einsum("ij,jk,ik->i", x, xtx_inv, x)
    adjusted_sq = residuals**2 / np.maximum((1 - leverage) ** 2, 1e-15)
    covariance = xtx_inv @ (x.T @ (x * adjusted_sq[:, None])) @ xtx_inv
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    t_values = beta / se
    p_values = 2 * student_t.sf(np.abs(t_values), df)
    critical = student_t.ppf(0.975, df)
    parameters = pd.DataFrame(
        {
            "term": names,
            "standardised_beta": beta,
            "hc3_se": se,
            "hc3_t": t_values,
            "hc3_p": p_values,
            "hc3_ci95_low": beta - critical * se,
            "hc3_ci95_high": beta + critical * se,
        }
    )
    summary = pd.DataFrame({"n": [n], "predictors": [" + ".join(predictors)]})
    return parameters, summary


def main() -> None:
    data = pd.read_csv(INPUT, dtype={"GridID": str})
    audit47 = pd.read_csv(AUDIT47)
    if len(data) != 40 or len(audit47) != 47:
        raise RuntimeError("Expected 40 morphology cells and 47 audited UU cells")

    rows = [
        correlation(data, "BCR_2025_pct", "2025 BCR"),
        correlation(data, "BCR_baseline_proxy_pct", "Baseline-proxy BCR"),
        correlation(data, "delta_BCR_pp", "Net BCR change"),
        correlation(data, "turnover_pp", "Footprint turnover"),
    ]
    unchanged = data.loc[~data["any_footprint_change"].astype(bool)]
    if len(unchanged) != 17:
        raise RuntimeError("Expected 17 SCL-p95 cells with no detected footprint change")
    rows.append(correlation(unchanged, "BCR_2025_pct", "2025 BCR; no detected change"))

    correlations = pd.DataFrame(rows)
    parameters, summary = fit_hc3(data)
    audit_summary = pd.DataFrame(
        {
            "sample": ["all_retained_UU_47", "SCL_p95_morphology_40"],
            "n": [47, 40],
            "added_m2": [audit47["added_m2"].sum(), data["added_m2"].sum()],
            "removed_m2": [audit47["removed_m2"].sum(), data["removed_m2"].sum()],
            "net_change_m2": [audit47["net_change_m2"].sum(), data["net_change_m2"].sum()],
            "cells_with_change": [
                int((audit47[["added_m2", "removed_m2"]].sum(axis=1) > 0).sum()),
                int(data["any_footprint_change"].astype(bool).sum()),
            ],
        }
    )

    correlations.to_csv(OUT / "footprint_audit_correlations.csv", index=False)
    parameters.to_csv(OUT / "footprint_audit_ols_hc3_parameters.csv", index=False)
    summary.to_csv(OUT / "footprint_audit_ols_summary.csv", index=False)
    audit_summary.to_csv(OUT / "footprint_audit_summary.csv", index=False)

    print(audit_summary.to_string(index=False))
    print(correlations.to_string(index=False))
    print(parameters.to_string(index=False))


if __name__ == "__main__":
    main()

