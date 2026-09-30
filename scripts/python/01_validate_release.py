#!/usr/bin/env python3
"""Validate the principal sample-size and result checkpoints in this release."""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "derived"
RESULTS = ROOT / "results"


def close(value: float, target: float, tolerance: float = 5e-4) -> None:
    if not np.isclose(value, target, atol=tolerance, rtol=0):
        raise AssertionError(f"Expected {target}, obtained {value}")


def main() -> None:
    validation = pd.read_csv(
        DATA / "classification_validation_confusion_matrices.csv"
    )
    assert len(validation) == 10
    for year, expected_correct, expected_oa, expected_kappa in (
        (2017, 852, 85.2, 0.815),
        (2025, 866, 86.6, 0.8325),
    ):
        rows = validation.loc[validation.year.eq(year)]
        assert len(rows) == 5 and int(rows.reference_total.sum()) == 1000
        assert int(rows.correct.sum()) == expected_correct
        close(rows.overall_accuracy_pct.iloc[0], expected_oa, 1e-9)
        close(rows.kappa.iloc[0], expected_kappa, 1e-9)

    retained = pd.read_csv(DATA / "endpoint_thermal_retained_cores_83.csv")
    counts = retained["core_default"].value_counts().to_dict()
    assert counts == {"UU": 47, "FF": 17, "BB": 15, "AA": 4}
    close(retained.loc[retained.core_default.eq("UU"), "dlst"].mean(), 2.144, 0.001)
    close(retained.loc[retained.core_default.eq("FF"), "dlst"].mean(), 1.050, 0.001)
    close(
        retained.loc[retained.core_default.eq("UU"), "dlst"].mean()
        - retained.loc[retained.core_default.eq("FF"), "dlst"].mean(),
        1.094,
        0.001,
    )

    annual = pd.read_csv(DATA / "annual_scl_core_stats_64.csv")
    assert len(annual) == 64
    assert int(annual.loc[annual.core_default.eq("UU"), "annual_stable_urban_p90"].sum()) == 13
    assert int(annual.loc[annual.core_default.eq("UU"), "annual_stable_urban_p95"].sum()) == 40
    assert int(annual.loc[annual.core_default.eq("UU"), "annual_stable_urban_p97"].sum()) == 44

    thresholds = pd.read_csv(DATA / "annual_scl_thresholds.csv")
    assert len(thresholds) == 3
    assert thresholds.anchor_valid_pixel_count.eq(74431).all()
    p95 = thresholds.loc[thresholds.anchor_percentile.eq(95)].iloc[0]
    close(p95.threshold_maxdNDBI, 0.1511067319661378, 1e-12)
    close(p95.threshold_maxdNDVI, 0.2383859455585479, 1e-12)
    close(p95.threshold_maxdTCB, 0.274189, 1e-12)

    scene = pd.read_csv(DATA / "landsat_scene_meteorology_audit_75.csv")
    panel = pd.read_csv(DATA / "core_acquisition_panel_60acq_3823obs.csv")
    assert len(scene) == 75 and int(scene.model_scene_keep_core95.sum()) == 60
    assert len(panel) == 3823
    assert panel.scene_id.nunique() == 60 and panel.GridID.nunique() == 64

    model = pd.read_csv(RESULTS / "meteorological_model_summary.csv")
    primary = model.loc[model.display_label.eq("Primary multiyear model")].iloc[0]
    close(primary.estimate_C, 0.791, 0.001)
    close(primary.ci95_low_C, 0.037, 0.001)
    close(primary.ci95_high_C, 1.546, 0.001)

    morph = pd.read_csv(DATA / "morphology_scl_p95_40.csv")
    footprint = pd.read_csv(DATA / "morphology_footprint_audit_scl_p95_40.csv")
    assert len(morph) == 40 and len(footprint) == 40
    close(morph[["BCR_2025_pct", "delta_LST_C"]].corr(method="spearman").iloc[0, 1], 0.569, 0.001)
    close(
        footprint[["BCR_baseline_proxy_pct", "delta_LST_C"]]
        .corr(method="spearman")
        .iloc[0, 1],
        0.547,
        0.001,
    )
    unchanged = footprint.loc[~footprint.any_footprint_change.astype(bool)]
    assert len(unchanged) == 17
    close(unchanged[["BCR_2025_pct", "delta_LST_C"]].corr(method="spearman").iloc[0, 1], 0.738, 0.001)

    print(
        "Release validation passed: 1,000 validation observations per endpoint; "
        "83 retained cores; 40 morphology cells; 60 acquisitions; "
        "3,823 observations."
    )


if __name__ == "__main__":
    main()
