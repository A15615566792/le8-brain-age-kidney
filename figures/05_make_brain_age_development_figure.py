from __future__ import annotations

import csv
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


PKG = Path(__file__).resolve().parent.parent
OUT_DIR = PKG / "figures_v2" / "main"
DATA_DIR = PKG / "tables" / "supplement"

SOURCE_DIR = PKG / "figures" / "data"
MODEL_RANKING_SOURCE = SOURCE_DIR / "Figure_2A_model_selection.csv"
TUNING_SOURCE = SOURCE_DIR / "Figure_2B_lasso_tuning.csv"
PREDICTION_SOURCE = PKG / "restricted" / "tables" / "candidate_predictions_all.csv"

FIGURE_STEM = OUT_DIR / "Figure_2_brain_age_model_development"
PROVENANCE = OUT_DIR / "Figure_2_brain_age_model_development_provenance.csv"


def configure_paths() -> None:
    global PKG, OUT_DIR, DATA_DIR, SOURCE_DIR
    global MODEL_RANKING_SOURCE, TUNING_SOURCE, PREDICTION_SOURCE
    global FIGURE_STEM, PROVENANCE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=PKG)
    parser.add_argument("--source-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    PKG = args.package_root.resolve()
    DATA_DIR = PKG / "tables" / "supplement"
    SOURCE_DIR = args.source_root.resolve() if args.source_root else PKG / "restricted"
    OUT_DIR = (
        args.output_root.resolve() if args.output_root
        else PKG / "figures_v2" / "main"
    )
    MODEL_RANKING_SOURCE = PKG / "figures" / "data" / "Figure_2A_model_selection.csv"
    TUNING_SOURCE = PKG / "figures" / "data" / "Figure_2B_lasso_tuning.csv"
    PREDICTION_SOURCE = SOURCE_DIR / "tables" / "candidate_predictions_all.csv"
    FIGURE_STEM = OUT_DIR / "Figure_2_brain_age_model_development"
    PROVENANCE = OUT_DIR / "Figure_2_brain_age_model_development_provenance.csv"


INK = "#20262E"
MUTED = "#5D6670"
AXIS = "#89939D"
GRID = "#E3E7EA"
BLUE = "#3C78A8"
TEAL = "#2A9D8F"
ORANGE = "#E76F51"
PALE_BLUE = "#DCE9F3"


mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.0,
        "axes.titlesize": 8.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": AXIS,
        "xtick.color": INK,
        "ytick.color": INK,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def validate_sources() -> None:
    for path in (MODEL_RANKING_SOURCE, TUNING_SOURCE):
        if not path.exists():
            raise FileNotFoundError(f"Required public summary file not found: {path}")
    if not PREDICTION_SOURCE.exists():
        raise FileNotFoundError(
            "Figure 2 panels C-G require the restricted local UK Biobank-derived "
            f"prediction table: {PREDICTION_SOURCE}"
        )


def prepare_figure_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

    ranking = pd.read_csv(MODEL_RANKING_SOURCE)
    tuning = pd.read_csv(TUNING_SOURCE)

    raw = pd.read_csv(
        PREDICTION_SOURCE,
        usecols=[
            "age",
            "split",
            "pred_raw__lasso__featurewiz",
            "pred_corrected__lasso__featurewiz",
            "bag__lasso__featurewiz",
        ],
    )
    predictions = raw.rename(
        columns={
            "pred_raw__lasso__featurewiz": "predicted_age_raw",
            "pred_corrected__lasso__featurewiz": "predicted_age_bias_corrected",
            "bag__lasso__featurewiz": "bag_bias_corrected",
        }
    )
    predictions["out_of_sample"] = predictions["split"].ne("Healthy train")

    return ranking, tuning, predictions


def clean_axis(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(AXIS)
    ax.spines["bottom"].set_color(AXIS)
    ax.spines["left"].set_linewidth(0.75)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.tick_params(length=3, width=0.7)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, lw=0.6, zorder=0)


def panel_label(ax: plt.Axes, letter: str, title: str) -> None:
    ax.text(-0.16, 1.12, letter, transform=ax.transAxes, fontsize=10.5, weight="bold", va="top")
    ax.text(0.00, 1.12, title, transform=ax.transAxes, fontsize=8.5, weight="bold", va="top")


def performance_stats(data: pd.DataFrame) -> tuple[int, float, float]:
    d = data[["age", "predicted_age_raw"]].dropna()
    mae = np.mean(np.abs(d["predicted_age_raw"] - d["age"]))
    r = np.corrcoef(d["age"], d["predicted_age_raw"])[0, 1]
    return len(d), float(mae), float(r)


def prediction_panel(ax: plt.Axes, data: pd.DataFrame, letter: str, title: str) -> None:
    d = data[["age", "predicted_age_raw"]].dropna()
    n, mae, r = performance_stats(d)
    lower = np.floor(d["age"].min() - 2)
    upper = np.ceil(d["age"].max() + 5)
    outside = ((d["predicted_age_raw"] < lower) | (d["predicted_age_raw"] > upper)).sum()
    ax.scatter(
        d["age"],
        d["predicted_age_raw"],
        s=5,
        alpha=0.075 if n > 5000 else 0.17,
        color=BLUE,
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )
    slope, intercept = np.polyfit(d["age"], d["predicted_age_raw"], 1)
    xx = np.array([lower, upper])
    ax.plot(xx, xx, color=INK, lw=0.85, ls=(0, (4, 3)), zorder=2)
    ax.plot(xx, intercept + slope * xx, color=ORANGE, lw=1.5, zorder=3)
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Chronological age, years")
    ax.set_ylabel("Predicted brain age, years")
    annotation = f"n = {n:,}\nr = {r:.3f}\nMAE = {mae:.2f} years"
    if outside:
        annotation += f"\n{outside:,} predictions outside display"
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        va="top",
        fontsize=7.0,
        linespacing=1.25,
    )
    panel_label(ax, letter, title)
    clean_axis(ax)


def histogram_panel(
    ax: plt.Axes,
    values: pd.Series,
    letter: str,
    title: str,
    xlabel: str,
    zero_line: bool = False,
    display_range: tuple[float, float] | None = None,
) -> None:
    values = pd.to_numeric(values, errors="coerce").dropna()
    outside = 0
    shown = values
    if display_range is not None:
        outside = ((values < display_range[0]) | (values > display_range[1])).sum()
        shown = values[(values >= display_range[0]) & (values <= display_range[1])]
        bins = np.linspace(display_range[0], display_range[1], 25)
    else:
        bins = np.histogram_bin_edges(values, bins=24)
    ax.hist(shown, bins=bins, color=BLUE, alpha=0.88, edgecolor="white", linewidth=0.45)
    if display_range is not None:
        ax.set_xlim(*display_range)
    if zero_line:
        ax.axvline(0, color=INK, lw=0.85, ls=(0, (4, 3)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Participants")
    note = f"n = {len(values):,}"
    if outside:
        note += f"\n{outside:,} outside display"
    ax.text(0.97, 0.92, note, transform=ax.transAxes, ha="right", va="top", fontsize=6.6, color=MUTED)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.08)
    panel_label(ax, letter, title)
    clean_axis(ax, "y")


def make_figure(ranking: pd.DataFrame, tuning: pd.DataFrame, predictions: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.35, 7.05))
    outer = fig.add_gridspec(3, 1, height_ratios=[0.92, 1.38, 0.92], hspace=0.57)
    top = outer[0].subgridspec(1, 2, wspace=0.40)
    middle = outer[1].subgridspec(1, 2, wspace=0.34)
    bottom = outer[2].subgridspec(1, 3, wspace=0.43)

    ax_a = fig.add_subplot(top[0, 0])
    rank_plot = ranking.sort_values("MAE", ascending=False).copy()
    colors = [BLUE if x == "LASSO" else PALE_BLUE for x in rank_plot["algorithm"]]
    edge = [BLUE if x == "LASSO" else AXIS for x in rank_plot["algorithm"]]
    ypos = np.arange(len(rank_plot))
    ax_a.barh(ypos, rank_plot["MAE"], color=colors, edgecolor=edge, linewidth=0.8, height=0.55, zorder=2)
    for y, (_, row) in zip(ypos, rank_plot.iterrows()):
        ax_a.text(row["MAE"] + 0.012, y, f"{row['MAE']:.2f}", va="center", fontsize=7.2)
    ax_a.set_yticks(ypos)
    ax_a.set_yticklabels(rank_plot["algorithm"])
    ax_a.set_xlim(0, max(rank_plot["MAE"]) + 0.42)
    ax_a.set_xlabel("Validation MAE, years (lower is better)")
    panel_label(ax_a, "a", "Candidate-model selection")
    clean_axis(ax_a, "x")

    ax_b = fig.add_subplot(top[0, 1])
    x = np.log10(tuning["alpha"].to_numpy())
    y = tuning["cv_mae"].to_numpy()
    ax_b.plot(x, y, color=BLUE, lw=1.2, zorder=2)
    ax_b.scatter(x, y, s=22, color=BLUE, edgecolor="white", linewidth=0.45, zorder=3)
    best = tuning.loc[tuning["selected"]].iloc[0]
    best_x = np.log10(best["alpha"])
    ax_b.axvline(best_x, color=ORANGE, lw=0.95, ls=(0, (3, 3)))
    ax_b.scatter([best_x], [best["cv_mae"]], s=44, color=ORANGE, edgecolor="white", linewidth=0.55, zorder=4)
    ax_b.annotate(
        f"Selected alpha = {best['alpha']:.4f}\nCV MAE = {best['cv_mae']:.2f}",
        xy=(best_x, best["cv_mae"]),
        xytext=(best_x + 0.35, best["cv_mae"] + 0.65),
        fontsize=6.6,
        color=MUTED,
        arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.65),
    )
    ax_b.set_xlabel("log10(alpha)")
    ax_b.set_ylabel("Five-fold CV MAE, years")
    panel_label(ax_b, "b", "LASSO tuning")
    clean_axis(ax_b, "y")

    validation = predictions.loc[predictions["split"].eq("Healthy validation")]
    nonhealthy = predictions.loc[predictions["split"].eq("Non-healthy test")]
    prediction_panel(fig.add_subplot(middle[0, 0]), validation, "c", "Internal validation")
    prediction_panel(fig.add_subplot(middle[0, 1]), nonhealthy, "d", "Out-of-sample non-healthy test")

    histogram_panel(
        fig.add_subplot(bottom[0, 0]),
        predictions["age"],
        "e",
        "Chronological age",
        "Age at MRI, years",
    )
    histogram_panel(
        fig.add_subplot(bottom[0, 1]),
        predictions["predicted_age_raw"],
        "f",
        "Predicted brain age",
        "Predicted age, years",
        display_range=(40, 90),
    )
    histogram_panel(
        fig.add_subplot(bottom[0, 2]),
        predictions.loc[predictions["out_of_sample"], "bag_bias_corrected"],
        "g",
        "Out-of-sample BAG",
        "Bias-corrected BAG, years",
        zero_line=True,
        display_range=(-20, 20),
    )

    fig.subplots_adjust(left=0.095, right=0.985, top=0.965, bottom=0.075)
    png = FIGURE_STEM.with_suffix(".png")
    fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(FIGURE_STEM.with_suffix(".pdf"), dpi=600, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)
    with Image.open(png) as image:
        rgb = image.convert("RGB")
        rgb.save(str(FIGURE_STEM.with_suffix(".tif")), format="TIFF", dpi=(600, 600))
        rgb.save(str(FIGURE_STEM.with_suffix(".jpg")), quality=95, optimize=True, dpi=(600, 600))


def write_provenance() -> None:
    rows = [
        {
            "panel": "A",
            "source_file": str(MODEL_RANKING_SOURCE),
            "packaged_plot_data": str(MODEL_RANKING_SOURCE),
            "rule": "Raw performance in the prespecified healthy validation split; all three completed candidate algorithms.",
        },
        {
            "panel": "B",
            "source_file": str(TUNING_SOURCE),
            "packaged_plot_data": str(TUNING_SOURCE),
            "rule": "All 10 completed LASSO alpha trials; minimum five-fold CV MAE marked as selected.",
        },
        {
            "panel": "C-G",
            "source_file": str(PREDICTION_SOURCE),
            "packaged_plot_data": "not packaged; restricted participant-level input",
            "rule": "Frozen LASSO plus FeatureWiz predictions; raw predictions for performance panels and age distributions; bias-corrected BAG outside the healthy training split for panel G.",
        },
    ]
    with PROVENANCE.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=["panel", "source_file", "packaged_plot_data", "rule"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    configure_paths()
    validate_sources()
    ranking, tuning, predictions = prepare_figure_data()
    make_figure(ranking, tuning, predictions)
    write_provenance()
    print(f"Created: {FIGURE_STEM.with_suffix('.png')}")
    print(f"Created: {FIGURE_STEM.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
