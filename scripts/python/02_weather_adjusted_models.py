#!/usr/bin/env python3
"""Meteorology-adjusted Antalya urban-versus-forest LST models.

The primary specification uses grid and scene fixed effects. Scene effects
absorb every meteorological state common to urban and forest cores on a given
overpass. Urban-by-weather terms allow those two surface groups to respond
differently to the same weather. Inference uses two-way CR1 clustering by
scene and prespecified 1-km spatial block.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "derived"
OUT = ROOT / "results" / "reproduced_weather_models"
OUT.mkdir(parents=True, exist_ok=True)

CORE_PATH = DATA / "core_acquisition_panel_60acq_3823obs.csv"
SCENE_PATH = DATA / "landsat_scene_meteorology_audit_75.csv"
ANNUAL_PATH = DATA / "annual_scl_core_stats_64.csv"

WEATHER_PRIMARY = [
    "t2m_C_aoi",
    "vpd_kPa_aoi",
    "wind10_mps_aoi",
    "ssrd_prev3h_Wm2_aoi",
    "sm_root_0_100_ff_background",
]
WEATHER_PRECIP = [
    "t2m_C_aoi",
    "vpd_kPa_aoi",
    "wind10_mps_aoi",
    "ssrd_prev3h_Wm2_aoi",
    "precip_prev30d_mm_aoi",
]


@dataclass
class FitResult:
    names: list[str]
    beta: np.ndarray
    residuals: np.ndarray
    vcov: np.ndarray
    rank: int
    nobs: int
    scene_clusters: int
    spatial_clusters: int
    df_inference: int
    rss: float
    r2: float
    covariance_psd_adjusted: bool


def restricted_cubic_spline(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Three-column natural/restricted cubic spline with four fixed quantile knots."""
    x = np.asarray(x, dtype=float)
    knots = np.quantile(x, [0.05, 0.35, 0.65, 0.95])
    if np.unique(knots).size != 4:
        raise RuntimeError("DOY knots are not unique")
    k1, _k2, k3, k4 = knots

    def positive_cube(z: np.ndarray) -> np.ndarray:
        return np.clip(z, 0.0, None) ** 3

    cols = [x]
    for kj in knots[:2]:
        h = (
            positive_cube(x - kj)
            - positive_cube(x - k3) * (k4 - kj) / (k4 - k3)
            + positive_cube(x - k4) * (k3 - kj) / (k4 - k3)
        ) / ((k4 - k1) ** 2)
        cols.append(h)
    basis = np.column_stack(cols)
    return basis, knots


def dummy_matrix(values: pd.Series, prefix: str) -> tuple[np.ndarray, list[str]]:
    levels = sorted(values.astype(str).unique())
    reference = levels[0]
    columns = []
    names = []
    raw = values.astype(str).to_numpy()
    for level in levels[1:]:
        columns.append((raw == level).astype(float))
        names.append(f"{prefix}[{level}]_ref[{reference}]")
    if not columns:
        return np.empty((len(values), 0)), []
    return np.column_stack(columns), names


def cluster_meat(
    x: np.ndarray,
    residuals: np.ndarray,
    groups: Iterable[str],
    rank: int,
) -> tuple[np.ndarray, int]:
    labels = pd.Series(list(groups), dtype="string")
    scores = x * residuals[:, None]
    score_frame = pd.DataFrame(scores)
    score_frame.insert(0, "cluster", labels.to_numpy())
    summed = score_frame.groupby("cluster", sort=False, observed=True).sum().to_numpy()
    meat = summed.T @ summed
    g = len(summed)
    n = len(residuals)
    if g <= 1 or n <= rank:
        raise RuntimeError("Insufficient clusters or residual degrees of freedom")
    correction = (g / (g - 1.0)) * ((n - 1.0) / (n - rank))
    return correction * meat, g


def fit_ols_two_way_cluster(
    y: np.ndarray,
    x: np.ndarray,
    names: list[str],
    scene_groups: pd.Series,
    spatial_groups: pd.Series,
) -> FitResult:
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    residuals = y - fitted
    xtx_inv = np.linalg.pinv(x.T @ x, rcond=1e-12)

    meat_scene, g_scene = cluster_meat(x, residuals, scene_groups, rank)
    meat_space, g_space = cluster_meat(x, residuals, spatial_groups, rank)
    intersection = (
        scene_groups.astype(str).to_numpy()
        + "__"
        + spatial_groups.astype(str).to_numpy()
    )
    meat_intersection, _ = cluster_meat(x, residuals, intersection, rank)
    vcov = xtx_inv @ (meat_scene + meat_space - meat_intersection) @ xtx_inv
    vcov = (vcov + vcov.T) / 2.0
    # Multiway inclusion-exclusion covariance estimates can be non-PSD in
    # finite samples. Apply the standard eigenvalue truncation used by major
    # fixed-effect implementations; record that the correction was required.
    eigenvalues, eigenvectors = np.linalg.eigh(vcov)
    psd_adjusted = bool(np.any(eigenvalues < -1e-10))
    if psd_adjusted:
        vcov = eigenvectors @ np.diag(np.clip(eigenvalues, 0.0, None)) @ eigenvectors.T
        vcov = (vcov + vcov.T) / 2.0

    tss = float(np.sum((y - np.mean(y)) ** 2))
    rss = float(residuals @ residuals)
    r2 = 1.0 - rss / tss
    return FitResult(
        names=names,
        beta=beta,
        residuals=residuals,
        vcov=vcov,
        rank=int(rank),
        nobs=len(y),
        scene_clusters=g_scene,
        spatial_clusters=g_space,
        df_inference=min(g_scene - 1, g_space - 1),
        rss=rss,
        r2=r2,
        covariance_psd_adjusted=psd_adjusted,
    )


def scalar_two_way_cluster_variance(
    influence: np.ndarray,
    residuals: np.ndarray,
    scene_groups: pd.Series,
    spatial_groups: pd.Series,
    rank: int,
) -> tuple[float, int, int, int]:
    """CR1 variance for a scalar linear estimator with crossed clusters."""

    def one_way(groups: Iterable[str]) -> tuple[float, int]:
        frame = pd.DataFrame(
            {
                "cluster": pd.Series(list(groups), dtype="string"),
                "score": influence * residuals,
            }
        )
        sums = frame.groupby("cluster", sort=False, observed=True)["score"].sum()
        g = len(sums)
        n = len(residuals)
        correction = (g / (g - 1.0)) * ((n - 1.0) / (n - rank))
        return float(correction * np.sum(sums.to_numpy() ** 2)), g

    v_scene, g_scene = one_way(scene_groups)
    v_space, g_space = one_way(spatial_groups)
    intersection = (
        scene_groups.astype(str).to_numpy()
        + "__"
        + spatial_groups.astype(str).to_numpy()
    )
    v_intersection, _ = one_way(intersection)
    return v_scene + v_space - v_intersection, g_scene, g_space, min(g_scene - 1, g_space - 1)


def scalar_one_way_cluster_variance(
    influence: np.ndarray,
    residuals: np.ndarray,
    groups: pd.Series,
    rank: int,
) -> tuple[float, int]:
    frame = pd.DataFrame(
        {
            "cluster": groups.astype(str).to_numpy(),
            "score": influence * residuals,
        }
    )
    sums = frame.groupby("cluster", sort=False, observed=True)["score"].sum()
    g = len(sums)
    n = len(residuals)
    correction = (g / (g - 1.0)) * ((n - 1.0) / (n - rank))
    return float(correction * np.sum(sums.to_numpy() ** 2)), g


def prepare_data() -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    core = pd.read_csv(CORE_PATH, dtype={"GridID": str, "scene_id": str})
    scene = pd.read_csv(SCENE_PATH, dtype={"scene_id": str})
    annual = pd.read_csv(ANNUAL_PATH, dtype={"GridID": str})

    stable = annual[
        annual["core_default"].isin(["UU", "FF"])
        & pd.to_numeric(annual["unstableRatio_p95"], errors="coerce").le(0.10)
        & pd.to_numeric(annual["annual_full_support"], errors="coerce").ge(0.90)
    ]
    stable_ids = set(stable["GridID"])
    stable_counts = stable.groupby("core_default")["GridID"].nunique().to_dict()
    if stable_counts != {"FF": 14, "UU": 40}:
        raise RuntimeError(f"Unexpected annually stable counts: {stable_counts}")

    join_cols = [
        "scene_id",
        "model_scene_keep_core95",
        "core_clear_frac",
        "aoi_clear_frac",
        "t2m_C_aoi",
        "td2m_C_aoi",
        "rh2m_pct_aoi",
        "vpd_kPa_aoi",
        "u10_mps_aoi",
        "v10_mps_aoi",
        "wind10_mps_aoi",
        "ssrd_prev3h_Wm2_aoi",
        "sm_0_7_aoi",
        "sm_0_7_ff_background",
        "sm_root_0_100_aoi",
        "sm_root_0_100_ff_background",
        "sm_100_289_ff_background",
        "precip_prev7d_mm_aoi",
        "precip_prev30d_mm_aoi",
        "t2m_localday_max_C_aoi",
    ]
    data = core.merge(scene[join_cols], on="scene_id", how="left", validate="many_to_one")
    data["urban"] = data["core_default"].eq("UU").astype(int)
    data["year"] = data["year"].astype(int)
    data["doy"] = pd.to_numeric(data["doy"], errors="raise")
    data["LST_C"] = pd.to_numeric(data["LST_C"], errors="coerce")
    data["core_scene_valid_frac"] = pd.to_numeric(
        data["core_scene_valid_frac"], errors="coerce"
    )

    base = data[
        data["core_default"].isin(["UU", "FF"])
        & data["model_scene_keep_core95"].eq(1)
        & data["core_scene_valid_frac"].ge(0.90)
        & data["LST_C"].notna()
    ].copy()
    return base, scene, stable_ids


def retain_scene_coverage(
    data: pd.DataFrame,
    expected_uu: int,
    expected_ff: int,
    threshold: float = 0.90,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    coverage = (
        data.groupby(["scene_id", "core_default"], observed=True)["GridID"]
        .nunique()
        .unstack(fill_value=0)
        .reset_index()
    )
    coverage["coverage_UU"] = coverage["UU"] / expected_uu
    coverage["coverage_FF"] = coverage["FF"] / expected_ff
    coverage["keep"] = (
        coverage["coverage_UU"].ge(threshold)
        & coverage["coverage_FF"].ge(threshold)
    )
    keep = set(coverage.loc[coverage["keep"], "scene_id"])
    return data[data["scene_id"].isin(keep)].copy(), coverage


def add_scene_covariates(data: pd.DataFrame, weather: list[str]) -> tuple[pd.DataFrame, dict]:
    scene_cols = ["scene_id", "doy", "sensor", *weather]
    scene_frame = data[scene_cols].drop_duplicates("scene_id").copy()
    if scene_frame[weather].isna().any().any():
        raise RuntimeError("Model-eligible scene has missing meteorology")

    spline, knots = restricted_cubic_spline(scene_frame["doy"].to_numpy())
    for j in range(spline.shape[1]):
        values = spline[:, j]
        sd = values.std(ddof=1)
        scene_frame[f"doy_rcs{j + 1}"] = (values - values.mean()) / sd

    scene_frame["L9_centered"] = scene_frame["sensor"].eq("LANDSAT_9").astype(float)
    scene_frame["L9_centered"] -= scene_frame["L9_centered"].mean()

    scaling = {"doy_knots": knots.tolist()}
    for column in weather:
        mean = scene_frame[column].mean()
        sd = scene_frame[column].std(ddof=1)
        if not np.isfinite(sd) or sd <= 0:
            raise RuntimeError(f"Cannot standardize {column}")
        scene_frame[f"{column}_z"] = (scene_frame[column] - mean) / sd
        scaling[column] = {"mean": float(mean), "sd": float(sd)}

    added = [
        "scene_id",
        "doy_rcs1",
        "doy_rcs2",
        "doy_rcs3",
        "L9_centered",
        *[f"{x}_z" for x in weather],
    ]
    return data.merge(scene_frame[added], on="scene_id", how="left", validate="many_to_one"), scaling


def build_design(
    data: pd.DataFrame,
    weather: list[str],
    add_ff_ndmi: bool = False,
) -> tuple[np.ndarray, list[str]]:
    n = len(data)
    matrices = [np.ones((n, 1))]
    names = ["intercept"]

    grid_matrix, grid_names = dummy_matrix(data["GridID"], "GridFE")
    scene_matrix, scene_names = dummy_matrix(data["scene_id"], "SceneFE")
    matrices.extend([grid_matrix, scene_matrix])
    names.extend(grid_names + scene_names)

    urban = data["urban"].to_numpy(dtype=float)
    for year in range(2018, 2026):
        matrices.append((urban * data["year"].eq(year).to_numpy()).reshape(-1, 1))
        names.append(f"urban_x_year_{year}")

    for column in ["doy_rcs1", "doy_rcs2", "doy_rcs3", "L9_centered"]:
        matrices.append((urban * data[column].to_numpy(dtype=float)).reshape(-1, 1))
        names.append(f"urban_x_{column}")

    for column in weather:
        z = f"{column}_z"
        matrices.append((urban * data[z].to_numpy(dtype=float)).reshape(-1, 1))
        names.append(f"urban_x_{z}")

    if add_ff_ndmi:
        matrices.append((urban * data["ff_NDMI_z"].to_numpy(dtype=float)).reshape(-1, 1))
        names.append("urban_x_ff_NDMI_z")

    return np.column_stack(matrices), names


def forest_ndmi_by_scene(data: pd.DataFrame) -> pd.DataFrame:
    ff = data[data["core_default"].eq("FF")].copy()
    valid_optical = (
        pd.to_numeric(ff["core_scene_optical_valid_frac"], errors="coerce").ge(0.80)
        & ff["NDMI"].notna()
        & ff["NDVI"].notna()
    )
    ff["NDMI_valid"] = ff["NDMI"].where(valid_optical)
    ff["NDVI_valid"] = ff["NDVI"].where(valid_optical)
    result = ff.groupby("scene_id", observed=True).agg(
        ff_LST_C=("LST_C", "mean"),
        ff_NDMI=("NDMI_valid", "mean"),
        ff_NDVI=("NDVI_valid", "mean"),
        ff_n_cores=("GridID", "nunique"),
        ff_optical_n=("NDMI_valid", "count"),
    ).reset_index()
    return result


def fit_specification(
    data: pd.DataFrame,
    sample_label: str,
    model_label: str,
    weather: list[str],
    spatial_cluster: str = "block_1000",
    add_ff_ndmi: bool = False,
) -> tuple[dict, FitResult]:
    model_data, _ = add_scene_covariates(data, weather)
    if add_ff_ndmi:
        model_data = model_data[model_data["ff_NDMI_z"].notna()].copy()
    x, names = build_design(model_data, weather, add_ff_ndmi=add_ff_ndmi)
    y = model_data["LST_C"].to_numpy(dtype=float)
    fit = fit_ols_two_way_cluster(
        y,
        x,
        names,
        model_data["scene_id"],
        model_data[spatial_cluster],
    )
    target = names.index("urban_x_year_2025")
    variance = fit.vcov[target, target]
    if not np.isfinite(variance) or variance <= 0:
        raise RuntimeError(
            f"Non-positive target variance for {sample_label}/{model_label}/{spatial_cluster}: {variance}"
        )
    se = float(np.sqrt(variance))
    estimate = float(fit.beta[target])
    critical = float(stats.t.ppf(0.975, fit.df_inference))
    t_value = estimate / se
    p_value = float(2 * stats.t.sf(abs(t_value), fit.df_inference))
    result = {
        "sample": sample_label,
        "model": model_label,
        "cluster_scheme": f"scene_id + {spatial_cluster}",
        "estimate_C": estimate,
        "std_error_C": se,
        "ci95_low_C": estimate - critical * se,
        "ci95_high_C": estimate + critical * se,
        "t_value": t_value,
        "df_inference": fit.df_inference,
        "p_value": p_value,
        "n_rows": fit.nobs,
        "n_scenes": int(model_data["scene_id"].nunique()),
        "n_cores": int(model_data["GridID"].nunique()),
        "n_scene_clusters": fit.scene_clusters,
        "n_spatial_clusters": fit.spatial_clusters,
        "model_rank": fit.rank,
        "r2_with_fixed_effects": fit.r2,
        "covariance_psd_adjusted": int(fit.covariance_psd_adjusted),
    }
    return result, fit


def standardized_endpoint_scenarios(data: pd.DataFrame) -> pd.DataFrame:
    """Standardize the original S0-S3 contrasts without changing their dates.

    Differential meteorological/calendar/sensor response is learned from the
    60-scene panel while year interactions protect the secular signal.  Each
    scenario's observed equal-core contrast is then transported to the common
    centered covariate distribution.  The reported variance treats the whole
    estimator as a linear function of the core-scene outcomes.
    """
    model_data, _ = add_scene_covariates(data, WEATHER_PRIMARY)
    x, names = build_design(model_data, WEATHER_PRIMARY)
    y = model_data["LST_C"].to_numpy(dtype=float)
    fit = fit_ols_two_way_cluster(
        y,
        x,
        names,
        model_data["scene_id"],
        model_data["block_1000"],
    )
    xtx_inv = np.linalg.pinv(x.T @ x, rcond=1e-12)
    nuisance = [
        "doy_rcs1",
        "doy_rcs2",
        "doy_rcs3",
        "L9_centered",
        *[f"{column}_z" for column in WEATHER_PRIMARY],
    ]
    rows = []
    for scenario in ["s0_member", "s1_member", "s2_current_member", "s3_member"]:
        c = np.zeros(len(model_data), dtype=float)
        unique_scene = model_data.drop_duplicates("scene_id")
        selected = unique_scene[unique_scene[scenario].eq(1)]
        counts_by_year = selected.groupby("year")["scene_id"].nunique().to_dict()
        if set(counts_by_year) != {2017, 2025}:
            raise RuntimeError(f"{scenario}: missing endpoint year")

        for year, year_sign in [(2017, -1.0), (2025, 1.0)]:
            scene_ids = selected.loc[selected["year"].eq(year), "scene_id"].tolist()
            for scene_id in scene_ids:
                scene_mask = model_data["scene_id"].eq(scene_id)
                for group, group_sign in [("UU", 1.0), ("FF", -1.0)]:
                    mask = scene_mask & model_data["core_default"].eq(group)
                    n_group = int(mask.sum())
                    if n_group == 0:
                        raise RuntimeError(f"{scenario}/{scene_id}: empty {group}")
                    c[mask.to_numpy()] = year_sign * group_sign / (
                        len(scene_ids) * n_group
                    )

        q = np.zeros(len(names), dtype=float)
        deltas = {}
        for column in nuisance:
            mean_2017 = selected.loc[selected["year"].eq(2017), column].mean()
            mean_2025 = selected.loc[selected["year"].eq(2025), column].mean()
            delta = float(mean_2025 - mean_2017)
            deltas[column] = delta
            q[names.index(f"urban_x_{column}")] = delta

        raw = float(c @ y)
        nuisance_adjustment = float(q @ fit.beta)
        standardized = raw - nuisance_adjustment
        influence = c - x @ xtx_inv @ q
        variance, g_scene, g_space, df = scalar_two_way_cluster_variance(
            influence,
            fit.residuals,
            model_data["scene_id"],
            model_data["block_1000"],
            fit.rank,
        )
        if not np.isfinite(variance) or variance <= 0:
            raise RuntimeError(f"{scenario}: non-positive standardized variance {variance}")
        se = float(np.sqrt(variance))
        critical = float(stats.t.ppf(0.975, df))
        t_value = standardized / se
        spatial_variance, spatial_g = scalar_one_way_cluster_variance(
            influence,
            fit.residuals,
            model_data["block_1000"],
            fit.rank,
        )
        spatial_se = float(np.sqrt(spatial_variance))
        spatial_critical = float(stats.t.ppf(0.975, spatial_g - 1))
        rows.append(
            {
                "scenario": scenario,
                "n_2017": counts_by_year[2017],
                "n_2025": counts_by_year[2025],
                "raw_change_in_D_C": raw,
                "estimated_nuisance_contribution_C": nuisance_adjustment,
                "weather_standardized_change_in_D_C": standardized,
                "std_error_C": se,
                "ci95_low_C": standardized - critical * se,
                "ci95_high_C": standardized + critical * se,
                "t_value": t_value,
                "df_inference": df,
                "p_value": float(2 * stats.t.sf(abs(t_value), df)),
                "n_scene_clusters": g_scene,
                "n_spatial_clusters": g_space,
                "spatial_only_std_error_C": spatial_se,
                "spatial_only_ci95_low_C": standardized - spatial_critical * spatial_se,
                "spatial_only_ci95_high_C": standardized + spatial_critical * spatial_se,
                "spatial_only_p_value": float(
                    2 * stats.t.sf(abs(standardized / spatial_se), spatial_g - 1)
                ),
                **{f"delta_{key}": value for key, value in deltas.items()},
            }
        )
    return pd.DataFrame(rows)


def endpoint_and_weather_audits(scene: pd.DataFrame, full: pd.DataFrame) -> None:
    endpoint = scene[scene["year"].isin([2017, 2025])].copy()
    weather_vars = [
        "t2m_C_aoi",
        "t2m_localday_max_C_aoi",
        "vpd_kPa_aoi",
        "wind10_mps_aoi",
        "ssrd_prev3h_Wm2_aoi",
        "sm_0_7_ff_background",
        "sm_root_0_100_ff_background",
        "precip_prev7d_mm_aoi",
        "precip_prev30d_mm_aoi",
    ]
    overlap_rows = []
    for variable in weather_vars:
        a = endpoint.loc[endpoint["year"].eq(2017), variable].dropna().to_numpy()
        b = endpoint.loc[endpoint["year"].eq(2025), variable].dropna().to_numpy()
        pooled = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
        overlap_rows.append(
            {
                "variable": variable,
                "n_2017": len(a),
                "n_2025": len(b),
                "mean_2017": np.mean(a),
                "mean_2025": np.mean(b),
                "min_2017": np.min(a),
                "max_2017": np.max(a),
                "min_2025": np.min(b),
                "max_2025": np.max(b),
                "standardized_mean_difference": (np.mean(b) - np.mean(a)) / pooled,
                "marginal_ranges_overlap": int(max(np.min(a), np.min(b)) <= min(np.max(a), np.max(b))),
            }
        )
    pd.DataFrame(overlap_rows).to_csv(OUT / "01_endpoint_weather_overlap.csv", index=False)

    hot_rows = []
    for label, variable in [
        ("overpass_T2m_ge_35C", "t2m_C_aoi"),
        ("localday_Tmax_ge_35C", "t2m_localday_max_C_aoi"),
    ]:
        for year in [2017, 2025]:
            q = endpoint[
                endpoint["year"].eq(year)
                & endpoint["s2_current_member"].eq(1)
                & endpoint[variable].ge(35)
            ]
            hot_rows.append(
                {
                    "threshold_definition": label,
                    "year": year,
                    "n_scenes": len(q),
                    "n_months": q["month"].nunique(),
                    "dates": ";".join(q["date_utc"].astype(str)),
                    "year_has_minimum_support": int(
                        len(q) >= 3 and q["month"].nunique() >= 2
                    ),
                }
            )
    hot_frame = pd.DataFrame(hot_rows)
    definition_ok = hot_frame.groupby("threshold_definition").apply(
        lambda z: int(
            z["year_has_minimum_support"].eq(1).all()
            and (z["n_scenes"].max() - z["n_scenes"].min() <= 1)
        ),
        include_groups=False,
    )
    hot_frame["feasible_for_balanced_composite"] = hot_frame[
        "threshold_definition"
    ].map(definition_ok)
    hot_frame.to_csv(OUT / "02_hot35_feasibility.csv", index=False)

    # Calendar-matched endpoint weather. The first three pairs are S0; all
    # four are S1. Differences are always 2025 minus 2017.
    pair_dates = [
        (1, "2017-06-25", "2025-06-23", 1, 1),
        (2, "2017-07-11", "2025-07-17", 1, 1),
        (3, "2017-08-12", "2025-08-10", 1, 1),
        (4, "2017-08-28", "2025-08-26", 0, 1),
    ]
    scene_date = scene.set_index("date_utc")
    paired_rows = []
    for pair_id, date_2017, date_2025, in_s0, in_s1 in pair_dates:
        row = {
            "pair_id": pair_id,
            "date_2017": date_2017,
            "date_2025": date_2025,
            "in_S0": in_s0,
            "in_S1": in_s1,
        }
        for variable in weather_vars:
            value_2017 = float(scene_date.loc[date_2017, variable])
            value_2025 = float(scene_date.loc[date_2025, variable])
            row[f"{variable}_2017"] = value_2017
            row[f"{variable}_2025"] = value_2025
            row[f"{variable}_delta_2025_minus_2017"] = value_2025 - value_2017
        paired_rows.append(row)
    pd.DataFrame(paired_rows).to_csv(
        OUT / "02b_calendar_paired_weather_balance.csv", index=False
    )

    # Reproduce the manuscript's four raw endpoint scenarios from the verified master table.
    scenario_rows = []
    for scenario in ["s0_member", "s1_member", "s2_current_member", "s3_member"]:
        q = endpoint[endpoint[scenario].eq(1)].copy()
        q["D_C"] = q["lst_C_uu"] - q["lst_C_ff"]
        means = q.groupby("year")["D_C"].mean()
        scenario_rows.append(
            {
                "scenario": scenario,
                "n_2017": int(q["year"].eq(2017).sum()),
                "n_2025": int(q["year"].eq(2025).sum()),
                "mean_D_2017_C": float(means.loc[2017]),
                "mean_D_2025_C": float(means.loc[2025]),
                "change_in_D_2025_minus_2017_C": float(means.loc[2025] - means.loc[2017]),
            }
        )
    pd.DataFrame(scenario_rows).to_csv(OUT / "03_raw_scenario_reproduction.csv", index=False)

    # Equal-core scene contrast in the primary 60-scene sample.
    contrast = (
        full.groupby(["scene_id", "date_utc", "year", "core_default"], observed=True)
        .agg(mean_LST_C=("LST_C", "mean"), n_cores=("GridID", "nunique"))
        .reset_index()
    )
    means = contrast.pivot(index=["scene_id", "date_utc", "year"], columns="core_default", values="mean_LST_C").reset_index()
    counts = contrast.pivot(index=["scene_id", "date_utc", "year"], columns="core_default", values="n_cores").reset_index()
    counts = counts.rename(columns={"UU": "n_UU", "FF": "n_FF"})
    means["urban_minus_forest_C"] = means["UU"] - means["FF"]
    means = means.merge(counts[["scene_id", "n_UU", "n_FF"]], on="scene_id", how="left")
    means.to_csv(OUT / "04_equal_core_scene_contrasts.csv", index=False)


def main() -> None:
    base, scene, stable_ids = prepare_data()
    full, full_coverage = retain_scene_coverage(base, 47, 17)
    stable_raw = base[base["GridID"].isin(stable_ids)].copy()
    stable, stable_coverage = retain_scene_coverage(stable_raw, 40, 14)
    full_coverage.to_csv(OUT / "05_full_scene_coverage.csv", index=False)
    stable_coverage.to_csv(OUT / "06_stable_scene_coverage.csv", index=False)

    if full["scene_id"].nunique() != 60 or stable["scene_id"].nunique() != 60:
        raise RuntimeError("Primary and stable samples must each retain 60 scenes")

    # Forest-canopy diagnostic and mechanism variable.
    ff_scene = forest_ndmi_by_scene(full)
    ff_scene = ff_scene.merge(
        scene,
        on="scene_id",
        how="left",
        validate="one_to_one",
    )
    ff_scene.to_csv(OUT / "07_forest_moisture_diagnostic.csv", index=False)
    ff_by_date = ff_scene.set_index("date_utc")
    forest_pair_rows = []
    forest_pair_dates = [
        (1, "2017-06-25", "2025-06-23", 1, 1),
        (2, "2017-07-11", "2025-07-17", 1, 1),
        (3, "2017-08-12", "2025-08-10", 1, 1),
        (4, "2017-08-28", "2025-08-26", 0, 1),
    ]
    forest_variables = [
        "ff_LST_C",
        "ff_NDMI",
        "ff_NDVI",
        "sm_0_7_ff_background",
        "sm_root_0_100_ff_background",
        "sm_100_289_ff_background",
        "precip_prev7d_mm_aoi",
        "precip_prev30d_mm_aoi",
    ]
    for pair_id, date_2017, date_2025, in_s0, in_s1 in forest_pair_dates:
        row = {
            "pair_id": pair_id,
            "date_2017": date_2017,
            "date_2025": date_2025,
            "in_S0": in_s0,
            "in_S1": in_s1,
        }
        for variable in forest_variables:
            value_2017 = float(ff_by_date.loc[date_2017, variable])
            value_2025 = float(ff_by_date.loc[date_2025, variable])
            row[f"{variable}_2017"] = value_2017
            row[f"{variable}_2025"] = value_2025
            row[f"{variable}_delta_2025_minus_2017"] = value_2025 - value_2017
        forest_pair_rows.append(row)
    pd.DataFrame(forest_pair_rows).to_csv(
        OUT / "07b_calendar_paired_forest_moisture.csv", index=False
    )
    ff_mean = ff_scene["ff_NDMI"].mean()
    ff_sd = ff_scene["ff_NDMI"].std(ddof=1)
    ff_scene["ff_NDMI_z"] = (ff_scene["ff_NDMI"] - ff_mean) / ff_sd
    full = full.merge(ff_scene[["scene_id", "ff_NDMI_z"]], on="scene_id", how="left", validate="many_to_one")
    stable = stable.merge(ff_scene[["scene_id", "ff_NDMI_z"]], on="scene_id", how="left", validate="many_to_one")

    endpoint_and_weather_audits(scene, full)
    standardized_endpoint_scenarios(full).to_csv(
        OUT / "05_weather_standardized_original_scenarios.csv",
        index=False,
    )

    samples = [
        ("47UU_17FF_endpoint_stable", full),
        ("40UU_14FF_SCL_exact_p95_stable", stable),
    ]
    results = []
    fits = {}
    for sample_label, data in samples:
        for model_label, weather, add_ndmi in [
            ("M0_sceneFE_calendar_sensor", [], False),
            ("M1_weather_adjusted_primary", WEATHER_PRIMARY, False),
            ("M1b_replace_soilwater_with_P30", WEATHER_PRECIP, False),
            ("M2_plus_FF_NDMI_diagnostic", WEATHER_PRIMARY, True),
        ]:
            result, fit = fit_specification(
                data,
                sample_label,
                model_label,
                weather,
                spatial_cluster="block_1000",
                add_ff_ndmi=add_ndmi,
            )
            results.append(result)
            fits[(sample_label, model_label)] = fit

    # Complete 64-core scenes hold spatial composition exactly fixed.
    complete_counts = full.groupby("scene_id")["GridID"].nunique()
    complete_ids = set(complete_counts[complete_counts.eq(64)].index)
    complete = full[full["scene_id"].isin(complete_ids)].copy()
    result, fit = fit_specification(
        complete,
        "47UU_17FF_complete_64core_scenes",
        "M1_weather_adjusted_primary",
        WEATHER_PRIMARY,
        spatial_cluster="block_1000",
    )
    results.append(result)
    fits[(result["sample"], result["model"])] = fit

    # Reviewer-4 sensitivity focused specifically on possible bare/construction
    # contamination in UU cores; all canonical FF references are retained.
    urban_bare_ids = set(
        full.drop_duplicates("GridID").loc[
            full.drop_duplicates("GridID")["core_default"].eq("FF")
            | (
                full.drop_duplicates("GridID")["bareFracMax"].le(0.05)
                & full.drop_duplicates("GridID")["bareFracIncrease"].le(0.02)
            ),
            "GridID",
        ]
    )
    urban_bare = full[full["GridID"].isin(urban_bare_ids)].copy()
    urban_bare_counts = (
        urban_bare.drop_duplicates("GridID")
        .groupby("core_default")["GridID"]
        .nunique()
        .to_dict()
    )
    urban_bare, _ = retain_scene_coverage(
        urban_bare,
        urban_bare_counts.get("UU", 0),
        urban_bare_counts.get("FF", 0),
    )
    result, fit = fit_specification(
        urban_bare,
        "UU_bare_max_le5pct_increase_le2pp_all_FF",
        "M1_weather_adjusted_primary",
        WEATHER_PRIMARY,
        spatial_cluster="block_1000",
    )
    results.append(result)
    fits[(result["sample"], result["model"])] = fit

    # Direct counterpart to the visual-audit threshold: remove every UU core
    # whose mapped bare fraction exceeds 5% at either endpoint, while retaining
    # all FF references.  Unlike the stricter check below, this does not add a
    # second bare-increase rule and therefore isolates the effect of the six
    # visually audited high-bare UU cores (41 UU + 17 FF expected).
    maxbare_ids = set(
        full.drop_duplicates("GridID").loc[
            full.drop_duplicates("GridID")["core_default"].eq("FF")
            | full.drop_duplicates("GridID")["bareFracMax"].le(0.05),
            "GridID",
        ]
    )
    maxbare = full[full["GridID"].isin(maxbare_ids)].copy()
    maxbare_counts = (
        maxbare.drop_duplicates("GridID")
        .groupby("core_default")["GridID"]
        .nunique()
        .to_dict()
    )
    if maxbare_counts != {"FF": 17, "UU": 41}:
        raise RuntimeError(
            f"Expected 41 UU and 17 FF after >5% maximum-bare exclusion; "
            f"found {maxbare_counts}"
        )
    maxbare, _ = retain_scene_coverage(maxbare, 41, 17)
    result, fit = fit_specification(
        maxbare,
        "41UU_bare_max_le5pct_17FF",
        "M1_weather_adjusted_primary",
        WEATHER_PRIMARY,
        spatial_cluster="block_1000",
    )
    results.append(result)
    fits[(result["sample"], result["model"])] = fit

    # AOI-clear sensitivity aligns the multiyear model with the conventional
    # all-AOI threshold rather than the primary core-focused clear criterion.
    aoi95 = full[full["aoi_clear_frac"].ge(0.95)].copy()
    aoi95, _ = retain_scene_coverage(aoi95, 47, 17)
    result, fit = fit_specification(
        aoi95,
        "AOI_clear_ge95pct",
        "M1_weather_adjusted_primary",
        WEATHER_PRIMARY,
        spatial_cluster="block_1000",
    )
    results.append(result)
    fits[(result["sample"], result["model"])] = fit

    # Strict bare-surface sensitivity: at most 5% bare at either endpoint and
    # no more than a 2 percentage-point increase between 2017 and 2025.
    eligible_bare_ids = set(
        full.loc[
            full["bareFracMax"].le(0.05)
            & full["bareFracIncrease"].le(0.02),
            "GridID",
        ]
    )
    bare = full[full["GridID"].isin(eligible_bare_ids)].copy()
    expected = bare.drop_duplicates("GridID").groupby("core_default")["GridID"].nunique().to_dict()
    bare, _ = retain_scene_coverage(bare, expected.get("UU", 0), expected.get("FF", 0))
    result, fit = fit_specification(
        bare,
        "bare_max_le5pct_increase_le2pp",
        "M1_weather_adjusted_primary",
        WEATHER_PRIMARY,
        spatial_cluster="block_1000",
    )
    results.append(result)
    fits[(result["sample"], result["model"])] = fit

    results_frame = pd.DataFrame(results)
    results_frame.to_csv(OUT / "08_weather_adjusted_endpoint_contrast.csv", index=False)

    # Humidity and hydroclimate definitions. These are joint-adjustment
    # sensitivities; individual weather coefficients are not causal effects.
    weather_definitions = {
        "T2m_plus_VPD_primary": WEATHER_PRIMARY,
        "T2m_plus_dewpoint": [
            "t2m_C_aoi", "td2m_C_aoi", "wind10_mps_aoi",
            "ssrd_prev3h_Wm2_aoi", "sm_root_0_100_ff_background",
        ],
        "T2m_plus_RH": [
            "t2m_C_aoi", "rh2m_pct_aoi", "wind10_mps_aoi",
            "ssrd_prev3h_Wm2_aoi", "sm_root_0_100_ff_background",
        ],
        "T2m_no_humidity": [
            "t2m_C_aoi", "wind10_mps_aoi", "ssrd_prev3h_Wm2_aoi",
            "sm_root_0_100_ff_background",
        ],
        "VPD_no_T2m": [
            "vpd_kPa_aoi", "wind10_mps_aoi", "ssrd_prev3h_Wm2_aoi",
            "sm_root_0_100_ff_background",
        ],
        "T2m_VPD_no_soilwater": [
            "t2m_C_aoi", "vpd_kPa_aoi", "wind10_mps_aoi",
            "ssrd_prev3h_Wm2_aoi",
        ],
        "T2m_VPD_replace_soilwater_P30": WEATHER_PRECIP,
    }
    weather_sensitivity_rows = []
    for definition, weather in weather_definitions.items():
        result, _ = fit_specification(
            full,
            "47UU_17FF_endpoint_stable",
            definition,
            weather,
            spatial_cluster="block_1000",
        )
        weather_sensitivity_rows.append(result)
    pd.DataFrame(weather_sensitivity_rows).to_csv(
        OUT / "12_weather_definition_sensitivity.csv", index=False
    )

    # Endpoint-scene influence audit. The primary point estimate should not be
    # driven by one 2017 or 2025 acquisition.
    endpoint_scenes = (
        full.loc[full["year"].isin([2017, 2025]), ["scene_id", "year", "date_utc"]]
        .drop_duplicates()
        .sort_values(["year", "date_utc"])
    )
    loo_rows = []
    for endpoint_row in endpoint_scenes.itertuples(index=False):
        reduced = full[full["scene_id"].ne(endpoint_row.scene_id)].copy()
        result, _ = fit_specification(
            reduced,
            "leave_one_endpoint_scene_out",
            "M1_weather_adjusted_primary",
            WEATHER_PRIMARY,
            spatial_cluster="block_1000",
        )
        result["dropped_scene_id"] = endpoint_row.scene_id
        result["dropped_year"] = endpoint_row.year
        result["dropped_date"] = endpoint_row.date_utc
        loo_rows.append(result)
    pd.DataFrame(loo_rows).to_csv(
        OUT / "13_endpoint_scene_leave_one_out.csv", index=False
    )

    # Prespecified spatial-clustering sensitivity for M1.
    cluster_rows = []
    for sample_label, data in samples:
        for cluster in ["GridID", "block_600", "block_1000", "block_1500", "block_2000"]:
            result, _ = fit_specification(
                data,
                sample_label,
                "M1_weather_adjusted_primary",
                WEATHER_PRIMARY,
                spatial_cluster=cluster,
            )
            cluster_rows.append(result)
    pd.DataFrame(cluster_rows).to_csv(OUT / "09_M1_cluster_sensitivity.csv", index=False)

    # Sample sizes by year.
    sample_rows = []
    for sample_label, data in [
        *samples,
        ("47UU_17FF_complete_64core_scenes", complete),
        ("bare_max_le5pct_increase_le2pp", bare),
        ("UU_bare_max_le5pct_increase_le2pp_all_FF", urban_bare),
        ("AOI_clear_ge95pct", aoi95),
    ]:
        for year, group in data.groupby("year"):
            sample_rows.append(
                {
                    "sample": sample_label,
                    "year": int(year),
                    "n_scenes": int(group["scene_id"].nunique()),
                    "n_cores": int(group["GridID"].nunique()),
                    "n_rows": len(group),
                }
            )
    pd.DataFrame(sample_rows).to_csv(OUT / "10_model_sample_by_year.csv", index=False)

    # Weather correlation and model-target coefficient sequence for diagnosis.
    scene_unique = full.drop_duplicates("scene_id")
    scene_unique[WEATHER_PRIMARY + ["precip_prev30d_mm_aoi"]].corr().to_csv(
        OUT / "11_weather_correlation_matrix.csv"
    )

    print(results_frame.to_string(index=False))


if __name__ == "__main__":
    main()
