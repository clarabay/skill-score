#!/usr/bin/env python3
"""
Convert cBase.csv baseline model forecasts to COVID-19 Forecast Hub quantile format.

REQUIRED INPUT — not included in this repository:
  cBase.csv must be downloaded from the "Baselines_Application_Data" record
      https://zenodo.org/records/16992035
  (listed there at 94.3 MB) and placed in the baseline model comparison/ directory at the repo root
  before running this script.

cBase.csv has identical wide-format structure to iBase.csv, but covers COVID-19
case forecasts for the COVID-19 Forecast Hub (2020-04 to 2023-02).

Column mapping (identical to iBase):
  x3–x13  : lower bounds of 11 prediction intervals (98% → 10%), widest first
             → quantiles 0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250,
                          0.300, 0.350, 0.400, 0.450
  x14–x24 : upper bounds of 11 prediction intervals (98% → 10%), widest first
             → quantiles 0.990, 0.975, 0.950, 0.900, 0.850, 0.800, 0.750,
                          0.700, 0.650, 0.600, 0.550
  x1       : mean  (excluded)
  x2       : median (excluded)
  x25      : observed truth (excluded — not present in cBase, but excluded for consistency)

Date conventions:
  cBase 'date' column = as_of_date (last data point available; always a Saturday)
  reference_date = as_of_date + 7 days (the following Saturday)
  forecast_date  = as_of_date + 2 days (Monday) — the COVID-19 Forecast Hub used
                   Monday submission deadlines throughout its entire run (2020–2023),
                   unlike FluSight which switched from Monday to Wednesday in 2023-24.

Horizon conventions (identical to iBase):
  cBase horizons 1–4 are relative to as_of_date.
  Since reference_date = as_of_date + 7:
    target_end_date = as_of_date + 7 * cBase_horizon
                    = reference_date + 7 * (cBase_horizon - 1)
    output_horizon  = cBase_horizon - 1   (0–3, relative to reference_date)
  For horizon 0: target_end_date == reference_date (current epiweek ending Saturday).

Target convention:
  COVID-19 Forecast Hub: "N wk ahead inc case"
"""

import pandas as pd
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE_DIR = REPO_ROOT / "baseline model comparison"
CBASE_PATH = BASELINE_DIR / "cBase.csv"
OUTPUT_DIR = BASELINE_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# iBase/cBase integer location ID → (FluSight/COVID-hub FIPS code, state name)
# cBase location 0 = national (US); locations 1–56 map to 2-digit FIPS strings
LOCATION_MAP: dict[int, tuple[str, str]] = {
    0:  ("US", "United States"),
    1:  ("01", "Alabama"),
    2:  ("02", "Alaska"),
    4:  ("04", "Arizona"),
    5:  ("05", "Arkansas"),
    6:  ("06", "California"),
    8:  ("08", "Colorado"),
    9:  ("09", "Connecticut"),
    10: ("10", "Delaware"),
    11: ("11", "District of Columbia"),
    12: ("12", "Florida"),
    13: ("13", "Georgia"),
    15: ("15", "Hawaii"),
    16: ("16", "Idaho"),
    17: ("17", "Illinois"),
    18: ("18", "Indiana"),
    19: ("19", "Iowa"),
    20: ("20", "Kansas"),
    21: ("21", "Kentucky"),
    22: ("22", "Louisiana"),
    23: ("23", "Maine"),
    24: ("24", "Maryland"),
    25: ("25", "Massachusetts"),
    26: ("26", "Michigan"),
    27: ("27", "Minnesota"),
    28: ("28", "Mississippi"),
    29: ("29", "Missouri"),
    30: ("30", "Montana"),
    31: ("31", "Nebraska"),
    32: ("32", "Nevada"),
    33: ("33", "New Hampshire"),
    34: ("34", "New Jersey"),
    35: ("35", "New Mexico"),
    36: ("36", "New York"),
    37: ("37", "North Carolina"),
    38: ("38", "North Dakota"),
    39: ("39", "Ohio"),
    40: ("40", "Oklahoma"),
    41: ("41", "Oregon"),
    42: ("42", "Pennsylvania"),
    44: ("44", "Rhode Island"),
    45: ("45", "South Carolina"),
    46: ("46", "South Dakota"),
    47: ("47", "Tennessee"),
    48: ("48", "Texas"),
    49: ("49", "Utah"),
    50: ("50", "Vermont"),
    51: ("51", "Virginia"),
    53: ("53", "Washington"),
    54: ("54", "West Virginia"),
    55: ("55", "Wisconsin"),
    56: ("56", "Wyoming"),
}

# cBase column → quantile probability (identical mapping to iBase)
QUANTILE_MAP: dict[str, float] = {
    "x3":  0.010,  # lower bound of 98% PI
    "x4":  0.025,  # lower bound of 95% PI
    "x5":  0.050,  # lower bound of 90% PI
    "x6":  0.100,  # lower bound of 80% PI
    "x7":  0.150,  # lower bound of 70% PI
    "x8":  0.200,  # lower bound of 60% PI
    "x9":  0.250,  # lower bound of 50% PI
    "x10": 0.300,  # lower bound of 40% PI
    "x11": 0.350,  # lower bound of 30% PI
    "x12": 0.400,  # lower bound of 20% PI
    "x13": 0.450,  # lower bound of 10% PI
    "x14": 0.990,  # upper bound of 98% PI
    "x15": 0.975,  # upper bound of 95% PI
    "x16": 0.950,  # upper bound of 90% PI
    "x17": 0.900,  # upper bound of 80% PI
    "x18": 0.850,  # upper bound of 70% PI
    "x19": 0.800,  # upper bound of 60% PI
    "x20": 0.750,  # upper bound of 50% PI
    "x21": 0.700,  # upper bound of 40% PI
    "x22": 0.650,  # upper bound of 30% PI
    "x23": 0.600,  # upper bound of 20% PI
    "x24": 0.550,  # upper bound of 10% PI
}

QUANTILE_COLS = list(QUANTILE_MAP.keys())
FINAL_COLS = [
    "location", "location_name", "forecast_date", "target",
    "target_end_date", "quantile", "value", "reference_date",
    "horizon", "log_value",
]


def main() -> None:
    print(f"Reading {CBASE_PATH.name} ...")
    df = pd.read_csv(CBASE_PATH)
    print(f"  {len(df):,} rows | {df['model'].nunique()} models | "
          f"{df['location'].nunique()} locations")

    # ------------------------------------------------------------------ dates
    # cBase 'date' is the as_of_date (last data point available, always a Saturday).
    # reference_date = following Saturday (as_of_date + 7).
    # forecast_date  = Monday (as_of_date + 2): the COVID-19 Forecast Hub used
    #                  Monday deadlines throughout its entire run — no season switch.
    as_of_date = pd.to_datetime(df["date"])
    df["reference_date"] = as_of_date + pd.Timedelta(weeks=1)
    df["forecast_date"]  = as_of_date + pd.Timedelta(days=2)

    # target_end_date = as_of_date + 7 * cBase_horizon
    #                 = reference_date + 7 * (cBase_horizon - 1)
    # output_horizon  = cBase_horizon - 1  (0–3, relative to reference_date)
    df["target_end_date"] = as_of_date + df["horizon"].apply(
        lambda h: pd.Timedelta(weeks=int(h))
    )
    df["horizon"] = df["horizon"] - 1
    df["target"] = df["horizon"].apply(lambda h: f"{h} wk ahead inc case")

    # --------------------------------------------------------------- locations
    df["location_name"] = df["location"].map(
        lambda x: LOCATION_MAP.get(int(x), (None, None))[1]
    )
    df["location"] = df["location"].map(
        lambda x: LOCATION_MAP.get(int(x), (None, None))[0]
    )

    unknown = df["location"].isna()
    if unknown.any():
        print(f"  Warning: dropping {unknown.sum():,} rows with unrecognised location IDs")
        df = df[~unknown].copy()

    # --------------------------------------------------------- wide → long
    id_vars = [
        "model", "location", "location_name", "forecast_date",
        "target", "target_end_date", "reference_date", "horizon",
    ]
    long = df[id_vars + QUANTILE_COLS].melt(
        id_vars=id_vars,
        value_vars=QUANTILE_COLS,
        var_name="x_col",
        value_name="value",
    )
    long["quantile"] = long["x_col"].map(QUANTILE_MAP)
    long = long.drop(columns=["x_col"])

    # log_value: log(1 + value); safe for zero-clipped values
    long["log_value"] = np.log1p(long["value"])

    # ----------------------------------------------------------- format dates
    for col in ("forecast_date", "target_end_date", "reference_date"):
        long[col] = pd.to_datetime(long[col]).dt.strftime("%Y-%m-%d")

    # ---------------------------------------------------------- write outputs
    models = sorted(long["model"].unique())
    print(f"\nFound {len(models)} models: {models}")
    print(f"Writing one CSV per model to {OUTPUT_DIR}/\n")

    for model_id in models:
        model_df = (
            long[long["model"] == model_id]
            [FINAL_COLS]
            .sort_values(["reference_date", "location", "target_end_date", "quantile"])
            .reset_index(drop=True)
        )
        out_path = OUTPUT_DIR / f"cbase_model_{model_id}.csv"
        model_df.to_csv(out_path, index=False)
        print(f"  model {model_id:2d}: {len(model_df):>8,} rows  →  {out_path.name}")

    print("\nDone.")


if __name__ == "__main__":
    main()
