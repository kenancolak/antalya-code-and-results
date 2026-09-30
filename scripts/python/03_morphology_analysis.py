#!/usr/bin/env python3
"""Reproduce primary continuous-morphology analyses for the 40 SCL-p95 cells."""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "derived" / "morphology_scl_p95_40.csv"
OUT = ROOT / "results" / "reproduced_morphology"
OUT.mkdir(parents=True, exist_ok=True)

PREDICTORS = ["BCR_2025_pct", "p75_storeys", "IVS_median_2017_2025"]
OUTCOME = "delta_LST_C"
SEED = 42
BOOTSTRAPS = 10_000


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=1)


def bootstrap_spearman(frame: pd.DataFrame, predictor: str) -> dict:
    x = frame[predictor].to_numpy(float)
    y = frame[OUTCOME].to_numpy(float)
    observed = spearmanr(x, y)
    rng = np.random.default_rng(SEED)
    values = np.full(BOOTSTRAPS, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for index in range(BOOTSTRAPS):
            sample = rng.integers(0, len(frame), len(frame))
            values[index] = spearmanr(x[sample], y[sample]).statistic
    valid = values[np.isfinite(values)]
    low, high = np.quantile(valid, [0.025, 0.975])
    return {
        "predictor": predictor,
        "n": len(frame),
        "rho": observed.statistic,
        "p_value": observed.pvalue,
        "bootstrap_ci_low": low,
        "bootstrap_ci_high": high,
        "valid_bootstrap_repetitions": len(valid),
    }


def fit_hc3(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = frame[[OUTCOME, *PREDICTORS]].dropna().copy()
    x_values = data[PREDICTORS].apply(zscore).to_numpy(float)
    y = zscore(data[OUTCOME]).to_numpy(float)
    x = np.column_stack([np.ones(len(data)), x_values])
    names = ["constant", *PREDICTORS]

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

    rss = float(residuals @ residuals)
    tss = float((y - y.mean()) @ (y - y.mean()))
    r_squared = 1 - rss / tss

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
    summary = pd.DataFrame(
        {"n": [n], "predictors": [" + ".join(PREDICTORS)], "r_squared": [r_squared]}
    )
    return parameters, summary


def main() -> None:
    data = pd.read_csv(INPUT, dtype={"GridID": str})
    if len(data) != 40 or data["annual_stable_urban_p95"].sum() != 40:
        raise RuntimeError("Expected the final 40-cell SCL-p95 morphology population")

    correlations = pd.DataFrame([bootstrap_spearman(data, item) for item in PREDICTORS])
    parameters, summary = fit_hc3(data)
    correlations.to_csv(OUT / "morphology_spearman_bootstrap.csv", index=False)
    parameters.to_csv(OUT / "morphology_ols_hc3_parameters.csv", index=False)
    summary.to_csv(OUT / "morphology_ols_summary.csv", index=False)

    print(correlations.to_string(index=False))
    print(parameters.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

