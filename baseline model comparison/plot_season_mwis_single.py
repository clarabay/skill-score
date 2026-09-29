import colorsys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from datetime import date
# ── Style definitions ──────────────────────────────────────────────────────
base_families = [
    "Constant", "Marginal", "KDE", "LSD(52,5)", "OLS(5,1)",
    "IDS(3)", "ARMA(1,1)", "INARCH(1)", "ETS+ST", "STL+ST",
]

def _desaturate(hex_color, factor=0.60):
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i+2], 16) / 255 for i in (0, 2, 4))
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s * factor)
    return "#{:02x}{:02x}{:02x}".format(round(r2 * 255), round(g2 * 255), round(b2 * 255))

_family_palette = [_desaturate(c) for c in [
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # yellow-green
    "#17becf",  # cyan
    "#ffbb78",  # light orange
    "#98df8a",  # light green
]]
family_colors = {fam: _family_palette[i % len(_family_palette)] for i, fam in enumerate(base_families)}

ETS_AD_LABEL  = "Trend Baseline"
NAIVE_LABEL   = "Naive"


def get_style(team):
    if team in ("FluSight-ensemble", "COVIDhub-4_week_ensemble"):
        return dict(color="black",     linestyle="-",  marker="o", linewidth=2.5, markersize=5)
    if team in ("FluSight-baseline", "COVIDhub-baseline"):
        return dict(color="firebrick", linestyle="-",  marker="o", linewidth=2.5, markersize=5)
    if team == "Log-baseline":
        return dict(color="navy",      linestyle="-",  marker="o", linewidth=2.5, markersize=5)
    if team == ETS_AD_LABEL:
        return dict(color="#1f77b4",   linestyle="-",  marker="o", linewidth=2.0, markersize=5)
    if team == NAIVE_LABEL:
        return dict(color="#636363",   linestyle="-",  marker="o", linewidth=2.0, markersize=5)
    base = team[4:] if team.startswith("log-") else team
    col  = family_colors.get(base, "gray")
    return dict(color=col, linestyle="--", marker="s", linewidth=1.0, markersize=4)


def team_sort_key(team):
    if team in ("FluSight-ensemble", "COVIDhub-4_week_ensemble"): return (0, team)
    if team in ("FluSight-baseline", "COVIDhub-baseline"):         return (1, team)
    if team == "Log-baseline":                                      return (2, team)
    if team == NAIVE_LABEL:                                         return (3, team)
    if team == ETS_AD_LABEL:                                        return (4, team)
    return (5, team)


def format_label(label):
    label = label.replace("-", " ")
    if label.lower().startswith("log "):
        label = "Log" + label[3:]
    return label


def legend_label(team):
    if team == "COVIDhub-4_week_ensemble":
        return "Hub Ensemble"
    if team == "COVIDhub-baseline":
        return "Hub Baseline"
    if team == "Log-baseline":
        return "RW Baseline"
    return team


def build_legend_handles(teams_present, hub_labels):
    handles = []
    for label in hub_labels:
        team = next((t for t in teams_present if legend_label(t) == label), label)
        s = get_style(team)
        handles.append(mlines.Line2D([], [],
            color=s["color"], linestyle=s["linestyle"],
            marker=s["marker"], linewidth=s["linewidth"],
            markersize=s["markersize"], label=format_label(label)))
    if ETS_AD_LABEL in teams_present:
        s = get_style(ETS_AD_LABEL)
        handles.append(mlines.Line2D([], [],
            color=s["color"], linestyle=s["linestyle"],
            marker=s["marker"], linewidth=s["linewidth"],
            markersize=s["markersize"], label=ETS_AD_LABEL))
    for fam in base_families:
        for name in (f"log-{fam}", fam):
            if name in teams_present:
                s = get_style(name)
                handles.append(mlines.Line2D([], [],
                    color=s["color"], linestyle=s["linestyle"],
                    marker=s["marker"], linewidth=s["linewidth"],
                    markersize=s["markersize"], label=format_label(name)))
                break
    return handles


# ── WIS scoring helpers ────────────────────────────────────────────────────
def next_saturday(date_series):
    d = pd.to_datetime(date_series)
    days = (5 - d.dt.dayofweek) % 7
    return (d + pd.to_timedelta(days, unit="D")).dt.date


def vectorized_log_wis(fc, group_cols):
    fc = fc.copy()
    fc["log_value"] = np.log(fc["value"].clip(lower=0) + 1)
    fc["_s"] = np.maximum(
        fc["quantile"] * (fc["log_obs_value"] - fc["log_value"]),
        (1 - fc["quantile"]) * (fc["log_value"] - fc["log_obs_value"]),
    )
    scored = (2 * fc.groupby(group_cols)["_s"].mean()).reset_index()
    scored.columns = group_cols + ["log_wis"]
    return scored


def score_ets_flu():
    """Returns row-level (team, season, year, location, target_end_date, horizon, log_wis)."""
    fc = pd.read_parquet("benchmark forecasts/trend baseline/naiveETS_damped_flu_hosp.parquet")
    fc["forecast_date"]   = pd.to_datetime(fc["forecast_date"]).dt.date
    fc["target_end_date"] = pd.to_datetime(fc["target_end_date"]).dt.date
    fc["reference_date"]  = pd.to_datetime(fc["reference_date"], errors="coerce").dt.date

    fc["horizon"] = fc["target"].str.extract(r"^(\d+) wk").astype(int) - 1
    fc = fc[fc["horizon"].isin(range(4))].copy()

    mask = fc["reference_date"].isna()
    fc.loc[mask, "reference_date"] = next_saturday(
        pd.Series(fc.loc[mask, "forecast_date"].values)
    ).values

    obs = pd.concat([
        pd.read_csv("surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv"),
        pd.read_csv("surveillance data/observed data/FluSight-hosp-observed.csv")[["target_end_date", "location", "value"]],
    ]).drop_duplicates(["location", "target_end_date"])
    obs["obs_value"]     = obs["value"].where(obs["value"] >= 0)
    obs["log_obs_value"] = np.log(obs["obs_value"] + 1)
    obs["target_end_date"] = pd.to_datetime(obs["target_end_date"]).dt.date

    fc = fc.merge(obs[["location", "target_end_date", "obs_value", "log_obs_value"]],
                  on=["location", "target_end_date"], how="inner")

    group_cols = ["location", "forecast_date", "reference_date", "target_end_date", "horizon"]
    scored = vectorized_log_wis(fc, group_cols)

    def assign_season(fd, rd, ted):
        if date(2022, 1, 10) <= fd <= date(2022, 6, 20):  return "2021-2022"
        if date(2022, 10, 17) <= fd <= date(2023, 5, 17): return "2022-2023"
        if rd is not None and date(2023, 10, 1) <= rd and ted < date(2024, 5, 4): return "2023-2024"
        if rd is not None and date(2024, 11, 20) <= rd <= date(2025, 5, 31):      return "2024-2025"
        return None

    scored["season"] = [assign_season(fd, rd, ted)
        for fd, rd, ted in zip(scored["forecast_date"], scored["reference_date"], scored["target_end_date"])]
    scored = scored.dropna(subset=["season"])

    exclude = {date(2022,1,15), date(2022,1,22), date(2025,5,31), date(2023,10,7), date(2025,1,25)}
    scored = scored[~scored["reference_date"].isin(exclude)]

    scored["n_h"] = scored.groupby(["location", "target_end_date"])["horizon"].transform("count")
    scored = scored[scored["n_h"] == 4].drop(columns="n_h")

    scored = scored[~scored["season"].isin(EXCLUDE_SEASONS)]
    scored["team"] = ETS_AD_LABEL
    scored["year"] = scored["season"].str[:4].astype(int)
    return scored[["team", "season", "year", "location", "target_end_date", "horizon", "log_wis"]]


def score_ets_covid():
    """Returns row-level (team, season, year, location, target_end_date, horizon, log_wis)."""
    fc = pd.read_parquet("benchmark forecasts/trend baseline/ets_additive_damped_covid_log.parquet")
    fc = fc[(fc["outcome"] == "case") & (fc["location"].str.len() == 2)]
    fc = fc[fc["target"].isin([f"{n} wk ahead inc cases" for n in range(1, 5)])]

    fc["forecast_date"]   = pd.to_datetime(fc["forecast_date"]).dt.date
    fc["target_end_date"] = pd.to_datetime(fc["target_end_date"]).dt.date
    fc["horizon"] = fc["target"].str.extract(r"^(\d+) wk").astype(int) - 1

    obs = pd.read_parquet("surveillance data/observed data/COVID19-observed.parquet")
    obs = obs[(obs["outcome"] == "case") & (obs["location"].str.len() == 2)]
    obs["target_end_date"] = pd.to_datetime(obs["target_end_date"]).dt.date

    fc = fc.merge(obs[["location", "target_end_date", "obs_value", "log_obs_value"]],
                  on=["location", "target_end_date"], how="inner")

    group_cols = ["location", "forecast_date", "target_end_date", "horizon"]
    scored = vectorized_log_wis(fc, group_cols)

    scored["season"] = pd.to_datetime(scored["target_end_date"]).dt.year.astype(str)
    scored = scored[
        scored["season"].isin(["2020", "2021"]) &
        (scored["forecast_date"] >= date(2020, 7, 28)) &
        (scored["forecast_date"] <= date(2021, 12, 21))
    ]

    scored["n_h"] = scored.groupby(["location", "target_end_date"])["horizon"].transform("count")
    scored = scored[scored["n_h"] == 4].drop(columns="n_h")

    scored["team"] = ETS_AD_LABEL
    scored["year"] = scored["season"].str[:4].astype(int)
    return scored[["team", "season", "year", "location", "target_end_date", "horizon", "log_wis"]]


# ── Load row-level data from pre-scored parquets ───────────────────────────
EXCLUDE_SEASONS = {"2024-2025"}

def load_rows(parquet_path, hub_teams, log_only=True):
    df = pd.read_parquet(parquet_path)
    if log_only:
        mask = df["team"].str.startswith("log-") | df["team"].isin(hub_teams)
    else:
        mask = df["team"].isin(hub_teams) | (
            (~df["team"].str.startswith("log-")) & df["team"].isin(base_families)
        )
    df = df[mask]
    df = df.dropna(subset=["log_wis"])
    df = df[~df["season"].isin(EXCLUDE_SEASONS)]
    df["target_end_date"] = pd.to_datetime(df["target_end_date"]).dt.date
    df["year"] = df["season"].str[:4].astype(int)
    return df[["team", "season", "year", "location", "target_end_date", "horizon", "log_wis"]]


def load_naive_rows(parquet_path):
    """Extract naive model scores from log_wis_naive column in combined parquet."""
    df = pd.read_parquet(parquet_path)
    df = df.dropna(subset=["log_wis_naive"])
    df = df[~df["season"].isin(EXCLUDE_SEASONS)]
    df["target_end_date"] = pd.to_datetime(df["target_end_date"]).dt.date
    df["year"] = df["season"].str[:4].astype(int)
    df = (df[["season", "year", "location", "target_end_date", "horizon", "log_wis_naive"]]
          .drop_duplicates(subset=["season", "location", "target_end_date", "horizon"]))
    df = df.rename(columns={"log_wis_naive": "log_wis"})
    df["team"] = NAIVE_LABEL
    return df[["team", "season", "year", "location", "target_end_date", "horizon", "log_wis"]]


# ── Apply ensemble reference dates and aggregate ───────────────────────────
REF_COLS = ["season", "location", "target_end_date", "horizon"]

def apply_reference_and_aggregate(rows, ets_rows, ensemble_team, label):
    ref = (rows[rows["team"] == ensemble_team][REF_COLS]
           .drop_duplicates()
           .reset_index(drop=True))

    all_rows = pd.concat([rows, ets_rows], ignore_index=True)
    filtered = all_rows.merge(ref, on=REF_COLS, how="inner")

    # Report any model missing dates relative to the ensemble
    missing_any = False
    for team in sorted(filtered["team"].unique()):
        if team == ensemble_team:
            continue
        team_dates = filtered[filtered["team"] == team][REF_COLS].drop_duplicates()
        missing = ref.merge(team_dates, on=REF_COLS, how="left", indicator=True)
        missing = missing[missing["_merge"] == "left_only"].drop(columns="_merge")
        if len(missing) > 0:
            missing_any = True
            by_season = missing.groupby("season").size()
            print(f"  [{label}] {team}: missing {len(missing)} location-date-horizon combos")
            for s, n in by_season.items():
                print(f"    {s}: {n}")
    if not missing_any:
        print(f"  [{label}] All models have complete coverage of ensemble dates.")

    mwis = (
        filtered.groupby(["team", "season", "year"])["log_wis"]
        .mean()
        .reset_index()
        .rename(columns={"log_wis": "mean_log_wis"})
    )
    mwis_overall = (
        filtered.groupby("team")["log_wis"]
        .mean()
        .reset_index()
        .rename(columns={"log_wis": "mean_log_wis"})
    )
    return mwis, mwis_overall


# ── Bar chart drawing ─────────────────────────────────────────────────────
def draw_bar_panel(ax, mwis_overall, title):
    df = mwis_overall.sort_values("mean_log_wis")

    colors  = [get_style(t)["color"] for t in df["team"]]
    hatches = ["///" if t.startswith("log-") else "" for t in df["team"]]
    labels  = [format_label(legend_label(t)) for t in df["team"]]

    x = range(len(df))
    for i, (val, col, hatch) in enumerate(zip(df["mean_log_wis"], colors, hatches)):
        ax.bar(i, val, color=col, hatch=hatch, edgecolor="white" if not hatch else col,
               linewidth=0.5)

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_title(title, fontweight="bold")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.6)
    ax.set_axisbelow(True)
    ax.set_ylim(0, df["mean_log_wis"].max() * 1.1)


# ── Panel drawing ──────────────────────────────────────────────────────────
def draw_panel(ax, mwis, hub_legend_labels, title, shared_ymax, show_legend=True):
    season_years  = sorted(mwis["year"].unique())
    season_labels = [mwis.loc[mwis["year"] == y, "season"].iloc[0] for y in season_years]

    ax.set_xlim(min(season_years) - 0.3, max(season_years) + 0.3)
    ax.set_ylim(0, shared_ymax)
    ax.set_xticks(season_years)
    ax.set_xticklabels(season_labels, rotation=45, ha="right")
    ax.set_title(title, fontweight="bold")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.6)
    ax.set_axisbelow(True)

    teams_present = sorted(mwis["team"].unique().tolist(), key=team_sort_key)
    n = len(teams_present)
    offsets = np.linspace(-0.07, 0.07, n) if n > 1 else [0.0]
    team_offset = {t: offsets[i] for i, t in enumerate(teams_present)}

    for team in teams_present:
        d = mwis[mwis["team"] == team].sort_values("year")
        s = get_style(team)
        ax.plot(d["year"] + team_offset[team], d["mean_log_wis"],
                color=s["color"], linestyle=s["linestyle"],
                linewidth=s["linewidth"], marker=s["marker"],
                markersize=s["markersize"],
                zorder=4 if team == ETS_AD_LABEL else 3)

    if show_legend:
        handles = build_legend_handles(teams_present, hub_legend_labels)
        ax.legend(handles=handles, fontsize=7, loc="upper left",
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)


# ── Score ETS models ───────────────────────────────────────────────────────
print("Scoring ETS Additive Damped (flu hosp)...")
ets_flu = score_ets_flu()
print("Scoring ETS Additive Damped (COVID cases)...")
ets_covid = score_ets_covid()

# ── Load row-level pre-scored data ─────────────────────────────────────────
ibase_rows = load_rows(
    "scored forecasts/iBase_combined_WIS.parquet",
    hub_teams=["FluSight-ensemble", "FluSight-baseline", "Log-baseline"],
)
cbase_rows = load_rows(
    "scored forecasts/cBase_combined_WIS.parquet",
    hub_teams=["COVIDhub-4_week_ensemble", "COVIDhub-baseline", "Log-baseline"],
)

naive_ibase = load_naive_rows("scored forecasts/iBase_combined_WIS.parquet")
naive_cbase = load_naive_rows("scored forecasts/cBase_combined_WIS.parquet")

# ── Filter to ensemble reference dates and aggregate ──────────────────────
print("\nChecking coverage vs ensemble reference dates:")
mwis_ibase, mwis_ibase_overall = apply_reference_and_aggregate(
    ibase_rows,
    pd.concat([ets_flu, naive_ibase], ignore_index=True),
    "FluSight-ensemble", "FluSight"
)
mwis_cbase, mwis_cbase_overall = apply_reference_and_aggregate(
    cbase_rows,
    pd.concat([ets_covid, naive_cbase], ignore_index=True),
    "COVIDhub-4_week_ensemble", "COVID"
)

shared_ymax = max(mwis_ibase["mean_log_wis"].max(), mwis_cbase["mean_log_wis"].max()) * 1.05

# ── Combined 2-panel figure ────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
fig.subplots_adjust(wspace=0.08)

draw_panel(ax1, mwis_ibase, ["FluSight-ensemble", "FluSight-baseline", "RW Baseline", NAIVE_LABEL],
           "FluSight Hospitalizations", shared_ymax, show_legend=False)
draw_panel(ax2, mwis_cbase, ["Hub Ensemble", "Hub Baseline", "RW Baseline", NAIVE_LABEL],
           "COVID-19 Cases", shared_ymax, show_legend=True)

ax1.set_ylabel("Mean log WIS")
ax2.set_ylabel("")

fig.savefig("figs/figS7-baseline_models_log_wis.pdf", bbox_inches="tight")
plt.close(fig)
print("\nSaved: figS7-baseline_models_log_wis.pdf")

