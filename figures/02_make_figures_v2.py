from __future__ import annotations

import csv
import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
from PIL import Image


PKG = Path(__file__).resolve().parent.parent
TABLES = PKG / "tables" / "supplement"
OUT = PKG / "figures_v2"
FIG_MAIN = OUT / "main"
FIG_SUPP = OUT / "supplement"
PROVENANCE_FILE = OUT / "figure_data_provenance_v2.csv"
MAIN_LEGENDS_FILE = FIG_MAIN / "Figure_legends_main.txt"


INK = "#20262E"
MUTED = "#5D6670"
AXIS = "#8C959F"
GRID = "#E3E7EA"
ROW = "#F6F7F8"
NAVY = "#1F4E79"
BLUE = "#3C78A8"
TEAL = "#2A9D8F"
VERMILLION = "#D55E00"
ORANGE = "#E76F51"
GOLD = "#D6A21F"
PURPLE = "#7B6AA8"
PALE_BLUE = "#EAF1F7"
PALE_TEAL = "#E9F5F2"
PALE_GRAY = "#F2F4F5"


mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.2,
        "axes.titlesize": 9.2,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.4,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.75,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.8,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "legend.frameon": False,
        "legend.fontsize": 7.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


PROVENANCE: list[dict[str, str]] = []


def configure_paths() -> None:
    global PKG, TABLES, OUT, FIG_MAIN, FIG_SUPP, PROVENANCE_FILE, MAIN_LEGENDS_FILE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--package-root", type=Path, default=PKG,
        help="Directory containing tables/supplement (default: repository root).",
    )
    parser.add_argument(
        "--output-root", type=Path, default=None,
        help="Figure output directory (default: <package-root>/figures_v2).",
    )
    args = parser.parse_args()
    PKG = args.package_root.resolve()
    TABLES = PKG / "tables" / "supplement"
    OUT = args.output_root.resolve() if args.output_root else PKG / "figures_v2"
    FIG_MAIN = OUT / "main"
    FIG_SUPP = OUT / "supplement"
    PROVENANCE_FILE = OUT / "figure_data_provenance_v2.csv"
    MAIN_LEGENDS_FILE = FIG_MAIN / "Figure_legends_main.txt"


def ensure_dirs() -> None:
    FIG_MAIN.mkdir(parents=True, exist_ok=True)
    FIG_SUPP.mkdir(parents=True, exist_ok=True)


def source(name: str) -> Path:
    path = TABLES / name
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def record(figure: str, panel: str, file_name: str, rule: str) -> None:
    PROVENANCE.append(
        {
            "figure": figure,
            "panel": panel,
            "source_file": str(source(file_name)),
            "selection_or_transformation": rule,
        }
    )


def save(fig: plt.Figure, stem: Path) -> None:
    png = stem.with_suffix(".png")
    fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)
    with Image.open(png) as image:
        rgb = image.convert("RGB")
        rgb.save(str(stem.with_suffix(".tif")), format="TIFF", dpi=(600, 600))
        rgb.save(str(stem.with_suffix(".jpg")), quality=95, optimize=True, dpi=(600, 600))


def is_true(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def fmt_num(x: float, digits: int = 2) -> str:
    if abs(float(x)) < 0.5 * 10 ** (-digits):
        x = 0.0
    return f"{x:.{digits}f}"


def fmt_ci(est: float, low: float, high: float, digits: int = 2) -> str:
    return (
        f"{fmt_num(est, digits)} "
        f"({fmt_num(low, digits)}, {fmt_num(high, digits)})"
    )


def clean_axis(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(AXIS)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(length=3, width=0.7)
    if grid:
        ax.grid(axis="x", color=GRID, lw=0.65, zorder=0)


def panel_title(ax: plt.Axes, letter: str, title: str) -> None:
    ax.set_title(f"{letter}   {title}", loc="left", pad=7, color=INK)


def shade_rows(ax: plt.Axes, n: int) -> None:
    for i in range(n):
        if i % 2 == 1:
            ax.axhspan(i - 0.5, i + 0.5, color=ROW, zorder=-5, lw=0)


def draw_ci(
    ax: plt.Axes,
    estimate: float,
    low: float,
    high: float,
    y: float,
    color: str,
    marker: str = "o",
    filled: bool = True,
    size: float = 27,
) -> None:
    ax.plot([low, high], [y, y], color=color, lw=1.45, solid_capstyle="round", zorder=2)
    ax.scatter(
        estimate,
        y,
        s=size,
        marker=marker,
        facecolor=color if filled else "white",
        edgecolor=color,
        linewidth=0.9,
        zorder=3,
    )


def legend_marker(color: str, label: str, marker: str = "o", filled: bool = True):
    return mpl.lines.Line2D(
        [0],
        [0],
        marker=marker,
        linestyle="none",
        markerfacecolor=color if filled else "white",
        markeredgecolor=color,
        markeredgewidth=0.9,
        markersize=5.2,
        label=label,
    )


def draw_flow_node(
    ax,
    x,
    y,
    w,
    h,
    title,
    detail,
    face="white",
    edge="#A9B1B8",
    lw=0.8,
    title_size=7.7,
    detail_size=7.2,
    compact=False,
    text_pad=0.022,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.008,rounding_size=0.009",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
    )
    ax.add_patch(patch)
    title_y = 0.73 if compact else 0.62
    detail_y = 0.27 if compact else 0.30
    ax.text(
        x + text_pad,
        y + h * title_y,
        title,
        ha="left",
        va="center",
        fontsize=title_size,
        weight="bold",
        linespacing=1.0,
    )
    ax.text(
        x + text_pad,
        y + h * detail_y,
        detail,
        ha="left",
        va="center",
        fontsize=detail_size,
        color=MUTED,
        linespacing=1.0,
    )


def flow_arrow(ax, x, y1, y2, color=AXIS):
    ax.annotate(
        "",
        xy=(x, y2),
        xytext=(x, y1),
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color=color, mutation_scale=8),
    )


def figure_1_flow() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.8), gridspec_kw={"wspace": 0.17})
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    left = axes[0]
    right = axes[1]
    left.text(0.03, 0.975, "a", fontsize=11, weight="bold", va="top")
    left.text(0.09, 0.975, "Full WebQ cohort: incident CKD", fontsize=10, weight="bold", va="top")
    right.text(0.03, 0.975, "b", fontsize=11, weight="bold", va="top")
    right.text(0.09, 0.975, "MRI cohort: multimodal brain age", fontsize=10, weight="bold", va="top")

    x, w, h = 0.055, 0.655, 0.105
    side_x, side_w = 0.755, 0.235
    left_nodes = [
        (0.82, "Valid Oxford WebQ", "n = 210,962; CKD events = 6,179", "white"),
        (0.61, "Eligible follow-up", "n = 210,053; events = 6,066", "white"),
        (0.40, "Complete 8-component LE8", "n = 156,107; events = 4,483", PALE_BLUE),
        (0.19, "Primary complete-case Cox model", "n = 155,910; events = 4,477", PALE_TEAL),
    ]
    exclusions = [
        (0.692, "Excluded", "Prevalent CKD: 889\nNo positive follow-up: 20"),
        (0.482, "Excluded", "Incomplete LE8: 53,946"),
        (0.272, "Excluded", "Incomplete analysis data: 197\nCKD events: 6"),
    ]
    for y, title, detail, face in left_nodes:
        draw_flow_node(left, x, y, w, h, title, detail, face=face, edge=NAVY if face != "white" else "#A9B1B8")
    for i in range(len(left_nodes) - 1):
        flow_arrow(left, x + w / 2, left_nodes[i][0], left_nodes[i + 1][0] + h)
    for y, title, detail in exclusions:
        draw_flow_node(
            left, side_x, y, side_w, 0.100, title, detail,
            face=PALE_GRAY, edge="#C1C7CC", lw=0.65,
            title_size=6.3, detail_size=5.7, compact=True, text_pad=0.012,
        )
        left.annotate(
            "", xy=(side_x, y + 0.046), xytext=(x + w, y + 0.046),
            arrowprops=dict(arrowstyle="-", color="#B2B9C0", lw=0.65),
        )

    right_nodes = [
        (0.84, "Brain-age prediction source", "n = 40,069", "white"),
        (0.705, "Matched valid pre-MRI WebQ", "n = 25,859", "white"),
        (0.57, "Eligible MRI follow-up", "n = 25,643", "white"),
        (0.435, "Complete LE8 and available BAG", "n = 19,326", PALE_GRAY),
        (0.30, "Out-of-sample combined BAG", "n = 16,536", PALE_BLUE),
        (0.165, "Primary complete-case MRI model", "n = 16,509", PALE_TEAL),
    ]
    for y, title, detail, face in right_nodes:
        draw_flow_node(right, x, y, w, 0.087, title, detail, face=face, edge=NAVY if face != "white" else "#A9B1B8")
    for i in range(len(right_nodes) - 1):
        flow_arrow(right, x + w / 2, right_nodes[i][0], right_nodes[i + 1][0] + 0.087)
    right_excl = [
        (0.737, "Unmatched WebQ", "n = 14,210"),
        (0.602, "Excluded", "Pre-MRI CKD: 214\nNo follow-up: 2"),
        (0.467, "Incomplete LE8", "n = 6,317"),
        (0.197, "Excluded", "Training leakage/\nmissing covariates: 27"),
    ]
    for y, title, detail in right_excl:
        draw_flow_node(
            right, side_x, y, side_w, 0.094, title, detail,
            face=PALE_GRAY, edge="#C1C7CC", lw=0.65,
            title_size=6.3, detail_size=5.7, compact=True, text_pad=0.012,
        )
        right.annotate(
            "", xy=(side_x, y + 0.042), xytext=(x + w, y + 0.042),
            arrowprops=dict(arrowstyle="-", color="#B2B9C0", lw=0.65),
        )

    save(fig, FIG_MAIN / "Figure_1_study_flow")
    record("Figure 1", "A-B", "cohort_flow_detailed.csv", "Cohort counts and exclusions shown without statistical transformation; final complete-case counts added from frozen main-model outputs.")


COMPONENT_LABELS = {
    "Diet": "Diet",
    "Physical_activity": "Physical activity",
    "Nicotine_exposure": "Nicotine exposure",
    "Sleep": "Sleep",
    "BMI": "Body mass index",
    "Non_HDL_cholesterol": "Non-HDL cholesterol",
    "Blood_glucose": "Blood glucose",
    "Blood_pressure": "Blood pressure",
}


def figure_2_components() -> None:
    file_name = "Table_3_LE8_component_profiles.csv"
    d = pd.read_csv(source(file_name))
    d = d[d["model"].eq("Mutually_adjusted_8_components")].copy()
    order = list(COMPONENT_LABELS)
    bag = d[d["outcome"].eq("BAG")].set_index("component").loc[order]
    ckd = d[d["outcome"].eq("Incident_CKD")].set_index("component").loc[order]
    y = np.arange(len(order))[::-1]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.65), gridspec_kw={"wspace": 0.08})
    for ax in axes:
        shade_rows(ax, len(order))

    ax = axes[0]
    for yy, key in zip(y, order):
        r = bag.loc[key]
        draw_ci(ax, r["estimate"], r["CI_low"], r["CI_high"], yy, BLUE, filled=float(r["p_FDR"]) < 0.05)
        ax.text(0.075, yy, fmt_ci(r["estimate"], r["CI_low"], r["CI_high"], 3), va="center", fontsize=6.8, color=MUTED)
    ax.axvline(0, color=AXIS, lw=0.85)
    ax.set_xlim(-0.225, 0.16)
    ax.set_yticks(y)
    ax.set_yticklabels([COMPONENT_LABELS[k] for k in order])
    ax.set_xlabel("Difference in combined BAG, years\nper 10-point higher component score")
    panel_title(ax, "a", "Combined multimodal BAG")
    ax.text(0.075, y.max() + 0.52, "Estimate (95% CI)", fontsize=6.8, color=MUTED)
    clean_axis(ax)

    ax = axes[1]
    for yy, key in zip(y, order):
        r = ckd.loc[key]
        draw_ci(ax, r["HR"], r["CI_low"], r["CI_high"], yy, VERMILLION, marker="s", filled=float(r["p_FDR"]) < 0.05)
        ax.text(1.072, yy, fmt_ci(r["HR"], r["CI_low"], r["CI_high"], 3), va="center", fontsize=6.8, color=MUTED)
    ax.axvline(1, color=AXIS, lw=0.85)
    ax.set_xlim(0.875, 1.135)
    ax.set_xticks([0.90, 0.95, 1.00, 1.05, 1.10])
    ax.set_xticklabels(["0.90", "0.95", "1.00", "1.05", "1.10"])
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("Hazard ratio for incident CKD\nper 10-point higher component score")
    panel_title(ax, "b", "Incident CKD")
    ax.text(1.072, y.max() + 0.52, "HR (95% CI)", fontsize=6.8, color=MUTED)
    clean_axis(ax)

    fig.subplots_adjust(bottom=0.22)
    fig.text(0.01, 0.012, "Filled symbols indicate FDR q < 0.05. All eight component scores were mutually adjusted.", fontsize=6.9, color=MUTED)
    save(fig, FIG_MAIN / "Figure_3_component_profiles")
    record("Figure 3", "A-B", file_name, "All eight mutually adjusted component estimates and 95% CIs; marker fill reflects the reported outcome-specific FDR value.")


def figure_3_imaging() -> None:
    p9 = "Table_P9_standardized_modality_specific_slopes.csv"
    p11 = "Table_P11_diffusion_metric_composites.csv"
    p12 = "Table_P12_sex_specific_slopes.csv"
    ds = pd.read_csv(source(p9))
    ds = ds[
        ds["component"].isin(["Blood glucose", "Blood pressure"])
        & ds["model"].str.contains("Mutually", case=False, na=False)
        & (~is_true(ds["selection_IPW"]))
    ].copy()
    dm = pd.read_csv(source(p11))
    dm = dm[
        dm["analysis"].isin(["SBP full sample", "SBP no BP medication"])
        & (~is_true(dm["selection_IPW"]))
        & (~dm["outcome"].str.contains("overall", case=False, na=False))
    ].copy()
    dx = pd.read_csv(source(p12))
    dx = dx[
        dx["analysis"].eq("Total LE8 per 10 points")
        & dx["modality"].isin(["combined", "t1", "diffusion"])
        & (~is_true(dx["selection_IPW"]))
    ].copy()

    fig, axes = plt.subplots(1, 3, figsize=(11.7, 3.55), gridspec_kw={"width_ratios": [0.92, 1.08, 1.0], "wspace": 0.42})

    ax = axes[0]
    comp_order = ["Blood glucose", "Blood pressure"]
    y = np.arange(2)[::-1]
    for mod, off, color, marker in [("T1", 0.12, BLUE, "o"), ("Diffusion", -0.12, TEAL, "s")]:
        sub = ds[ds["modality"].eq(mod)].set_index("component")
        for yy, comp in zip(y, comp_order):
            r = sub.loc[comp]
            draw_ci(ax, r["beta_BAG_SD"], r["CI_low"], r["CI_high"], yy + off, color, marker=marker)
    ax.axvline(0, color=AXIS, lw=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(comp_order)
    ax.set_xlabel("Standardized BAG difference\nper 10-point higher score")
    panel_title(ax, "a", "Modality-specific\ncomponent associations")
    ax.legend(handles=[legend_marker(BLUE, "T1-BAG", "o"), legend_marker(TEAL, "Diffusion-BAG", "s")], loc="upper right", ncol=1, handletextpad=0.4)
    clean_axis(ax)

    ax = axes[1]
    preferred = ["FA lower-direction", "MD higher-direction", "ICVF lower-direction", "ISOVF higher-direction", "ISOV higher-direction", "OD original direction"]
    outcomes = [x for x in preferred if x in set(dm["outcome"])]
    y = np.arange(len(outcomes))[::-1]
    for analysis, off, color, marker in [("SBP full sample", 0.12, NAVY, "o"), ("SBP no BP medication", -0.12, ORANGE, "s")]:
        sub = dm[dm["analysis"].eq(analysis)].set_index("outcome")
        for yy, outcome in zip(y, outcomes):
            if outcome not in sub.index:
                continue
            r = sub.loc[outcome]
            draw_ci(ax, r["beta"], r["CI_low"], r["CI_high"], yy + off, color, marker=marker)
    ax.axvline(0, color=AXIS, lw=0.85)
    ax.set_yticks(y)
    metric_labels = [x.split(" ")[0] for x in outcomes]
    ax.set_yticklabels(["ISOVF" if x == "ISOV" else x for x in metric_labels])
    ax.set_xlabel("Adverse metric difference, SD\nper 10-mm Hg higher SBP")
    panel_title(ax, "b", "Systolic blood pressure and\ndiffusion MRI metric composites")
    ax.legend(handles=[legend_marker(NAVY, "Full sample", "o"), legend_marker(ORANGE, "No BP medication", "s")], loc="upper left", ncol=1, handletextpad=0.4)
    clean_axis(ax)

    ax = axes[2]
    mod_order = ["combined", "t1", "diffusion"]
    y = np.arange(3)[::-1]
    for sex, off, color, marker in [("Female", 0.12, PURPLE, "o"), ("Male", -0.12, TEAL, "s")]:
        sub = dx[dx["sex"].eq(sex)].set_index("modality")
        for yy, mod in zip(y, mod_order):
            r = sub.loc[mod]
            draw_ci(ax, r["beta"], r["CI_low"], r["CI_high"], yy + off, color, marker=marker)
    ax.axvline(0, color=AXIS, lw=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(["Combined", "T1", "Diffusion"])
    ax.set_xlabel("BAG difference, years\nper 10-point higher LE8")
    panel_title(ax, "c", "Sex-specific LE8–BAG\nassociations")
    ax.legend(
        handles=[legend_marker(PURPLE, "Female", "o"), legend_marker(TEAL, "Male", "s")],
        loc="lower right",
        ncol=1,
        handletextpad=0.4,
    )
    clean_axis(ax)

    save(fig, FIG_MAIN / "Figure_4_brain_imaging_deep_phenotyping")
    record("Figure 4", "A", p9, "Mutually adjusted unweighted blood-glucose and blood-pressure component slopes.")
    record("Figure 4", "B", p11, "Unweighted SBP metric-composite models in the full and untreated samples.")
    record("Figure 4", "C", p12, "Unweighted sex-specific total-LE8 slopes after sex-specific brain-age calibration.")


def figure_s1_performance() -> None:
    file_name = "modality_model_performance.csv"
    d = pd.read_csv(source(file_name))
    d = d[d["correction"].eq("Raw") & d["split"].isin(["Healthy validation", "Healthy train", "Non-healthy test"])].copy()
    splits = ["Healthy train", "Healthy validation", "Non-healthy test"]
    split_labels = ["Training", "Internal validation", "Non-healthy OOS test"]
    y = np.arange(3)[::-1]
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    offsets = {"T1": 0.18, "Diffusion": 0.0, "Combined": -0.18}
    colors = {"T1": BLUE, "Diffusion": ORANGE, "Combined": TEAL}
    markers = {"T1": "o", "Diffusion": "^", "Combined": "s"}
    for mod in ["T1", "Diffusion", "Combined"]:
        sub = d[d["modality"].eq(mod)].set_index("split")
        for yy, split in zip(y, splits):
            value = float(sub.loc[split, "MAE"])
            ax.scatter(value, yy + offsets[mod], s=32, marker=markers[mod], color=colors[mod], zorder=3)
            ax.text(value + 0.035, yy + offsets[mod], f"{value:.2f}", va="center", fontsize=6.7, color=colors[mod])
    shade_rows(ax, 3)
    ax.set_yticks(y)
    ax.set_yticklabels(split_labels)
    ax.set_xlabel("Mean absolute error, years\n(lower is better)")
    ax.set_xlim(2.45, 4.05)
    ax.set_title("Brain-age model performance", loc="left", pad=25)
    ax.legend(handles=[legend_marker(BLUE, "T1", "o"), legend_marker(ORANGE, "Diffusion", "^"), legend_marker(TEAL, "Combined", "s")], loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3)
    clean_axis(ax)
    save(fig, FIG_SUPP / "Figure_S1_brain_age_performance")
    record("Figure S1", "all", file_name, "Raw MAE displayed for each predefined split and modality; no axis truncation that changes ordering.")


def figure_s2_absolute_risk() -> None:
    file_name = "Table_P4_CKD_adjusted_FG_CIF.csv"
    d = pd.read_csv(source(file_name))
    colors = {"Low (<50)": VERMILLION, "Moderate (50-79)": GOLD, "High (>=80)": TEAL}
    labels = {"Low (<50)": "Low LE8", "Moderate (50-79)": "Moderate LE8", "High (>=80)": "High LE8"}
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    for group in ["Low (<50)", "Moderate (50-79)", "High (>=80)"]:
        sub = d[d["LE8_category"].eq(group)].sort_values("horizon_years")
        x = sub["horizon_years"].to_numpy()
        y = sub["adjusted_risk_per_1000"].to_numpy()
        ax.plot(x, y, color=colors[group], lw=1.8, marker="o", ms=4.6)
        for xx, yy in zip(x, y):
            ax.text(xx, yy + 1.0, f"{yy:.1f}", ha="center", fontsize=6.8, color=colors[group])
        ax.text(10.12, y[-1], labels[group], va="center", fontsize=7.0, color=colors[group], weight="bold")
    ax.set_xticks([5, 10])
    ax.set_xticklabels(["5 years", "10 years"])
    ax.set_xlim(4.65, 10.90)
    ax.set_ylim(0, max(d["adjusted_risk_per_1000"]) * 1.16)
    ax.set_ylabel("Adjusted CKD risk per 1,000")
    ax.set_title("Adjusted CKD risk by LE8 category", loc="left")
    clean_axis(ax, grid=False)
    ax.grid(axis="y", color=GRID, lw=0.65)
    save(fig, FIG_SUPP / "Figure_S2_adjusted_CKD_absolute_risk")
    record("Figure S2", "all", file_name, "Model-standardized cumulative-incidence point estimates only.")


def figure_s3_ckd_sensitivity() -> None:
    file_name = "Table_P2_total_LE8_CKD_sensitivities.csv"
    d = pd.read_csv(source(file_name))
    order = ["Primary", "Plus_baseline_eGFR", "Plus_eGFR_UACR_complete_case", "Two_year_landmark"]
    d = d.set_index("model").loc[order].reset_index()
    labels = ["Primary Cox", "+ baseline eGFR", "+ eGFR and UACR", "2-year landmark"]
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    shade_rows(ax, len(d))
    for yy, (_, r) in zip(y, d.iterrows()):
        draw_ci(ax, r["HR"], r["CI_low"], r["CI_high"], yy, NAVY, filled=True)
        ax.text(0.925, yy, fmt_ci(r["HR"], r["CI_low"], r["CI_high"], 2), va="center", fontsize=7.0, color=MUTED)
    ax.axvline(1, color=AXIS, lw=0.85)
    ax.set_xlim(0.65, 1.02)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Hazard ratio\nper 10-point higher LE8")
    ax.set_title("CKD sensitivity analyses", loc="left")
    ax.text(0.925, y.max() + 0.48, "HR (95% CI)", fontsize=6.8, color=MUTED)
    clean_axis(ax)
    save(fig, FIG_SUPP / "Figure_S3_CKD_sensitivity_forest")
    record("Figure S3", "all", file_name, "All prespecified total-LE8 Cox sensitivity models, displayed in protocol order.")


def figure_s4_components_by_modality() -> None:
    file_name = "Table_P9_standardized_modality_specific_slopes.csv"
    d = pd.read_csv(source(file_name))
    d = d[d["model"].str.contains("Mutually", case=False, na=False) & (~is_true(d["selection_IPW"]))].copy()
    order = ["Diet", "Physical activity", "Nicotine exposure", "Sleep", "BMI", "Non-HDL cholesterol", "Blood glucose", "Blood pressure"]
    y = np.arange(len(order))[::-1]
    fig, ax = plt.subplots(figsize=(5.9, 4.6))
    shade_rows(ax, len(order))
    for mod, off, color, marker in [("T1", 0.13, BLUE, "o"), ("Diffusion", -0.13, TEAL, "s")]:
        sub = d[d["modality"].eq(mod)].set_index("component")
        for yy, comp in zip(y, order):
            r = sub.loc[comp]
            draw_ci(ax, r["beta_BAG_SD"], r["CI_low"], r["CI_high"], yy + off, color, marker=marker, size=23)
    ax.axvline(0, color=AXIS, lw=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.set_xlabel("Standardized BAG difference\nper 10-point higher component score")
    ax.set_title("LE8 components by MRI modality", loc="left", pad=27)
    ax.legend(handles=[legend_marker(BLUE, "T1-BAG", "o"), legend_marker(TEAL, "Diffusion-BAG", "s")], loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    clean_axis(ax)
    save(fig, FIG_SUPP / "Figure_S4_all_components_by_modality")
    record("Figure S4", "all", file_name, "All eight mutually adjusted, unweighted modality-specific slopes.")


def figure_s5_continuous_exposures() -> None:
    file_name = "Table_P10_continuous_BP_glycaemia_modality_slopes.csv"
    d = pd.read_csv(source(file_name))
    d = d[~is_true(d["selection_IPW"])].copy()
    d["row_key"] = d["component"].astype(str) + " | " + d["model"].astype(str)
    order = list(dict.fromkeys(d["row_key"].tolist()))
    labels = []
    for key in order:
        label = key.replace(" | Full sample adjusted for treatment/disease", "\nFull sample")
        label = label.replace(" | No BP medication", "\nNo BP medication")
        label = label.replace(" | No diabetes/insulin", "\nNo diabetes/insulin")
        label = label.replace(" | No self-reported diabetes/insulin", "\nNo diabetes/insulin")
        label = label.replace("Systolic BP", "SBP").replace("Diastolic BP", "DBP")
        labels.append(label)
    y = np.arange(len(order))[::-1]
    fig, ax = plt.subplots(figsize=(6.0, 5.3))
    shade_rows(ax, len(order))
    for mod, off, color, marker in [("T1", 0.13, BLUE, "o"), ("Diffusion", -0.13, TEAL, "s")]:
        sub = d[d["modality"].eq(mod)].set_index("row_key")
        for yy, key in zip(y, order):
            r = sub.loc[key]
            draw_ci(ax, r["beta_BAG_SD"], r["CI_low"], r["CI_high"], yy + off, color, marker=marker, size=23)
    ax.axvline(0, color=AXIS, lw=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Standardized BAG difference")
    ax.set_title("Continuous BP and glycaemia by MRI modality", loc="left", pad=27)
    ax.legend(handles=[legend_marker(BLUE, "T1-BAG", "o"), legend_marker(TEAL, "Diffusion-BAG", "s")], loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    clean_axis(ax)
    save(fig, FIG_SUPP / "Figure_S5_continuous_BP_glycaemia_modalities")
    record("Figure S5", "all", file_name, "All unweighted continuous blood-pressure and glycaemia slopes, with prespecified untreated/disease-free subsets retained.")


def tract_name(description: str) -> str:
    name = re.sub(r"^Weighted-mean\s+\S+\s+in tract\s+", "", str(description), flags=re.I)
    name = name.replace(" (left)", " (L)").replace(" (right)", " (R)")
    name = name.replace("part of ", "")
    return name


def figure_s6_tract_heatmap() -> None:
    file_name = "Table_P11_significant_tracts_FDR134.csv"
    d = pd.read_csv(source(file_name))
    d = d[d["analysis"].eq("SBP no BP medication") & (~is_true(d["selection_IPW"]))].copy()
    d["metric"] = d["metric"].replace({"ISOV": "ISOVF"})
    d["tract"] = d["description"].map(tract_name)
    metrics = [m for m in ["FA", "MD", "ICVF", "ISOVF", "OD"] if m in set(d["metric"])]
    d["plot_value"] = d["adverse_beta"]
    d.loc[d["metric"].eq("OD"), "plot_value"] = d.loc[d["metric"].eq("OD"), "beta"]
    pivot = d.pivot_table(index="tract", columns="metric", values="plot_value", aggfunc="first")
    pivot = pivot.reindex(columns=metrics)
    pivot = pivot.assign(_mean=pivot.mean(axis=1)).sort_values("_mean", ascending=False).drop(columns="_mean")

    fig_h = max(6.4, 0.255 * len(pivot) + 1.5)
    fig, ax = plt.subplots(figsize=(7.3, fig_h))
    cmap = LinearSegmentedColormap.from_list("signed", [PURPLE, "#F6F7F8", NAVY])
    max_abs = float(np.nanmax(np.abs(pivot.to_numpy())))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)
    matrix = np.ma.masked_invalid(pivot.to_numpy())
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, norm=norm, interpolation="none")
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels(metrics, weight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=6.4)
    ax.tick_params(axis="both", length=0)
    for i in range(len(pivot)):
        for j, metric in enumerate(metrics):
            value = pivot.iloc[i, j]
            if pd.isna(value):
                ax.text(j, i, "-", ha="center", va="center", color="#C6CBD0", fontsize=6.0)
            else:
                color = "white" if abs(value) > max_abs * 0.58 else INK
                ax.text(j, i, f"{value:.3f}", ha="center", va="center", color=color, fontsize=5.8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("SBP-associated white-matter tracts among participants not using BP medication", loc="left", pad=24)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Standardized difference per 10-mm Hg higher SBP", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.5)
    fig.text(0.01, 0.006, f"Cells show all {len(d)} globally FDR-significant results. Positive values are adverse-direction effects for FA/MD/ICVF/ISOVF; OD retains its original signed beta. Dashes indicate no FDR-significant result.", fontsize=6.7, color=MUTED)
    save(fig, FIG_SUPP / "Figure_S6_all_significant_tract_IDPs")
    record("Figure S6", "all", file_name, "All globally FDR-significant unweighted SBP tract results in untreated participants; adverse-direction standardized effects displayed as a complete tract-by-metric heatmap without manual tract selection.")


def figure_s7_sex_specific() -> None:
    file_name = "Table_P12_sex_specific_slopes.csv"
    d = pd.read_csv(source(file_name))
    d = d[d["analysis"].eq("Total LE8 per 10 points") & (~is_true(d["selection_IPW"])) & d["modality"].isin(["combined", "t1", "diffusion"])].copy()
    order = ["combined", "t1", "diffusion"]
    y = np.arange(3)[::-1]
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    shade_rows(ax, 3)
    for sex, off, color, marker in [("Female", 0.13, PURPLE, "o"), ("Male", -0.13, TEAL, "s")]:
        sub = d[d["sex"].eq(sex)].set_index("modality")
        for yy, mod in zip(y, order):
            r = sub.loc[mod]
            draw_ci(ax, r["beta"], r["CI_low"], r["CI_high"], yy + off, color, marker=marker)
    ax.axvline(0, color=AXIS, lw=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(["Combined", "T1", "Diffusion"])
    ax.set_xlabel("BAG difference, years\nper 10-point higher LE8")
    ax.set_title("Sex-specific LE8–BAG associations", loc="left", pad=25)
    ax.legend(handles=[legend_marker(PURPLE, "Female", "o"), legend_marker(TEAL, "Male", "s")], loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    clean_axis(ax)
    save(fig, FIG_SUPP / "Figure_S7_sex_specific_LE8_BAG")
    record("Figure S7", "all", file_name, "Unweighted female- and male-specific total-LE8 slopes after sex-specific brain-age calibration.")


def write_provenance() -> None:
    with PROVENANCE_FILE.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=["figure", "panel", "source_file", "selection_or_transformation"])
        writer.writeheader()
        writer.writerows(PROVENANCE)


def write_main_figure_legends() -> None:
    legends = """Figure 1. Study flow for the full WebQ cohort and the brain MRI cohort. (A) Full WebQ cohort for the prospective analysis of Life's Essential 8 and incident chronic kidney disease. (B) MRI cohort for analyses of Life's Essential 8 and multimodal brain age.

Figure 2. Development and evaluation of the combined multimodal brain-age model. (A) Candidate-model selection. (B) LASSO tuning. (C) Predicted versus chronological age in the healthy internal-validation sample. (D) Predicted versus chronological age in the out-of-sample nonhealthy test sample. (E) Distribution of chronological age at MRI. (F) Distribution of raw predicted brain age. (G) Distribution of bias-corrected out-of-sample brain age gap.

Figure 3. Mutually adjusted associations of individual Life's Essential 8 components with brain and kidney outcomes. (A) Combined multimodal BAG. (B) Incident CKD.

Figure 4. Modality-specific, white matter microstructural, and sex-specific brain imaging associations. (A) Modality-specific component associations. (B) Systolic blood pressure and diffusion MRI metric composites. (C) Sex-specific LE8–BAG associations.
"""
    MAIN_LEGENDS_FILE.write_text(legends, encoding="utf-8")


def main() -> None:
    configure_paths()
    ensure_dirs()
    figure_1_flow()
    figure_2_components()
    figure_3_imaging()
    figure_s1_performance()
    figure_s2_absolute_risk()
    figure_s3_ckd_sensitivity()
    figure_s4_components_by_modality()
    figure_s5_continuous_exposures()
    figure_s6_tract_heatmap()
    figure_s7_sex_specific()
    write_provenance()
    write_main_figure_legends()
    print(f"Created 3 association main figures and 7 supplementary figures in: {OUT}")
    print("Run 05_make_brain_age_development_figure.py to create main Figure 2.")


if __name__ == "__main__":
    main()
