from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DIET_POSITIVE = [
    "Oliveoilcat", "greenvegetablecat", "othervegetablecat", "Berriscat",
    "Otherfruitcat", "Fish_24hcat", "Beanscat", "Wholegrainscat",
    "Nuts_24hcat",
]
DIET_NEGATIVE = [
    "Meat_24hcat", "Chicken_24hcat", "Cheese_24hcat", "Butter_24hcat",
    "Sweets_24hcat",
]
DIET_OTHER = ["Fastfoodcat", "Alcohol_24hcat"]

SOURCE_COLUMNS = [
    "participant_id", "age", "sex", "TPA_min_weekly",
    "Quit_Smoke_Yearcat", "Passivesmoking_home", "sleep_duration", "BMI",
    "Non_HDL", "FBG_mgdl", "HbA1c_percent", "DMPMH_SR", "SBP", "DBP",
    "BPM", "Essential8",
] + DIET_POSITIVE + DIET_NEGATIVE + DIET_OTHER


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def normalize_id(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)


def score_by_bins(x: pd.Series, rules: list[tuple[pd.Series, float]]) -> pd.Series:
    score = pd.Series(np.nan, index=x.index, dtype=float)
    for condition, value in rules:
        score.loc[condition & x.notna()] = value
    return score


def safe_qcut(series: pd.Series, q: int = 3) -> pd.Series:
    result = pd.Series(np.nan, index=series.index, dtype=float)
    valid = series.notna()
    if valid.sum() >= q:
        result.loc[valid] = (
            pd.qcut(series.loc[valid].rank(method="first"), q=q, labels=False)
            .astype(float) + 1
        )
    return result


def read_source(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".dta":
        data = pd.read_stata(path, columns=SOURCE_COLUMNS, convert_categoricals=False)
    else:
        data = pd.read_csv(path, usecols=SOURCE_COLUMNS, low_memory=False)
    data["participant_id"] = normalize_id(data["participant_id"])
    if data["participant_id"].duplicated().any():
        raise ValueError("participant_id must be unique in the LE8 source table")
    return data


def score_le8(data: pd.DataFrame) -> pd.DataFrame:
    d = data.copy()

    diet_items = pd.DataFrame(index=d.index)
    for column in DIET_POSITIVE:
        diet_items[column] = numeric(d[column]).fillna(0).clip(0, 1)
    for column in DIET_NEGATIVE:
        diet_items[column] = numeric(d[column]).fillna(1).clip(0, 1)
    fast_food = numeric(d["Fastfoodcat"])
    diet_items["Fastfood_healthy"] = (
        (1 - fast_food).where(fast_food.notna(), 0).clip(0, 1)
    )
    diet_items["Alcohol_healthy"] = (
        numeric(d["Alcohol_24hcat"]).fillna(0).clip(0, 1)
    )
    d["diet_le8_webq16"] = diet_items.sum(axis=1)
    diet = d["diet_le8_webq16"]
    d["le8_diet_primary"] = score_by_bins(diet, [
        (diet.between(0, 3), 0), (diet.between(4, 7), 25),
        (diet.between(8, 11), 50), (diet.between(12, 14), 80),
        (diet.between(15, 16), 100),
    ])

    pa = numeric(d["TPA_min_weekly"])
    d["le8_pa"] = score_by_bins(pa, [
        (pa.eq(0), 0), (pa.gt(0) & pa.lt(30), 20),
        (pa.ge(30) & pa.lt(60), 40), (pa.ge(60) & pa.lt(90), 60),
        (pa.ge(90) & pa.lt(120), 80), (pa.ge(120) & pa.lt(150), 90),
        (pa.ge(150), 100),
    ])

    quit_category = numeric(d["Quit_Smoke_Yearcat"])
    d["le8_nicotine"] = quit_category.map(
        {0: 100, 1: 75, 2: 50, 3: 25, 4: 0}
    ).astype(float)
    passive = numeric(d["Passivesmoking_home"])
    exposed = passive.gt(0) & d["le8_nicotine"].gt(0)
    d.loc[exposed, "le8_nicotine"] -= 20
    d["le8_nicotine"] = d["le8_nicotine"].clip(0, 100)

    sleep = numeric(d["sleep_duration"])
    d["le8_sleep"] = score_by_bins(sleep, [
        (sleep.le(3), 0), (sleep.eq(4), 20), (sleep.eq(5), 40),
        (sleep.eq(6), 70), (sleep.between(7, 8), 100),
        (sleep.eq(9), 90), (sleep.ge(10), 40),
    ])

    bmi = numeric(d["BMI"])
    d["le8_bmi"] = score_by_bins(bmi, [
        (bmi.lt(25), 100), (bmi.ge(25) & bmi.lt(30), 70),
        (bmi.ge(30) & bmi.lt(35), 30), (bmi.ge(35) & bmi.lt(40), 15),
        (bmi.ge(40), 0),
    ])

    non_hdl = numeric(d["Non_HDL"])
    d["le8_nonhdl"] = score_by_bins(non_hdl, [
        (non_hdl.lt(130), 100), (non_hdl.ge(130) & non_hdl.lt(160), 60),
        (non_hdl.ge(160) & non_hdl.lt(190), 40),
        (non_hdl.ge(190) & non_hdl.lt(220), 20), (non_hdl.ge(220), 0),
    ])

    fasting_glucose = numeric(d["FBG_mgdl"])
    hba1c = numeric(d["HbA1c_percent"])
    diabetes = numeric(d["DMPMH_SR"]).eq(1)
    glycaemia = pd.Series(np.nan, index=d.index, dtype=float)
    no_diabetes = ~diabetes
    glycaemia.loc[no_diabetes & (fasting_glucose.lt(100) | hba1c.lt(5.7))] = 100
    glycaemia.loc[no_diabetes & (fasting_glucose.ge(100) | hba1c.ge(5.7))] = 60
    glycaemia.loc[no_diabetes & (fasting_glucose.ge(126) | hba1c.ge(6.5))] = 40
    glycaemia.loc[diabetes & hba1c.lt(7)] = 40
    glycaemia.loc[diabetes & hba1c.ge(7) & hba1c.lt(8)] = 30
    glycaemia.loc[diabetes & hba1c.ge(8) & hba1c.lt(9)] = 20
    glycaemia.loc[diabetes & hba1c.ge(9) & hba1c.lt(10)] = 10
    glycaemia.loc[diabetes & hba1c.ge(10)] = 0
    glycaemia.loc[diabetes & hba1c.isna() & fasting_glucose.notna()] = 40
    glycaemia.loc[fasting_glucose.isna() & hba1c.isna()] = np.nan
    d["le8_glucose"] = glycaemia

    sbp = numeric(d["SBP"])
    dbp = numeric(d["DBP"])
    blood_pressure = pd.Series(np.nan, index=d.index, dtype=float)
    valid = sbp.notna() & dbp.notna()
    blood_pressure.loc[valid & sbp.lt(120) & dbp.lt(80)] = 100
    blood_pressure.loc[valid & sbp.between(120, 129) & dbp.lt(80)] = 75
    blood_pressure.loc[valid & (sbp.between(130, 139) | dbp.between(80, 89))] = 50
    blood_pressure.loc[valid & (sbp.between(140, 159) | dbp.between(90, 99))] = 25
    blood_pressure.loc[valid & (sbp.ge(160) | dbp.ge(100))] = 0
    treated = numeric(d["BPM"]).eq(1) & blood_pressure.notna()
    blood_pressure.loc[treated] = (blood_pressure.loc[treated] - 20).clip(lower=0)
    d["le8_bp"] = blood_pressure

    components = [
        "le8_diet_primary", "le8_pa", "le8_nicotine", "le8_sleep",
        "le8_bmi", "le8_nonhdl", "le8_glucose", "le8_bp",
    ]
    d["le8_n_components"] = d[components].notna().sum(axis=1)
    d["LE8"] = d[components].mean(axis=1).where(
        d["le8_n_components"].eq(8)
    )
    behavior = components[:4]
    biological = components[4:]
    d["LE8_behavior"] = d[behavior].mean(axis=1).where(d[behavior].notna().all(axis=1))
    d["LE8_biological"] = d[biological].mean(axis=1).where(d[biological].notna().all(axis=1))
    d["LE8_category"] = pd.cut(
        d["LE8"], [-np.inf, 49.999, 79.999, np.inf],
        labels=["low", "moderate", "high"],
    ).astype("string")
    d["LE8_tertile"] = safe_qcut(d["LE8"])

    d = d.rename(columns={
        "Essential8": "LE8_legacy_do",
    })
    keep = [
        "participant_id", "age", "sex", "LE8_legacy_do", "LE8",
        "LE8_category", "LE8_tertile", "LE8_behavior", "LE8_biological",
        "le8_n_components", *components, "diet_le8_webq16",
        "TPA_min_weekly", "sleep_duration", "BMI",
        "Non_HDL", "FBG_mgdl", "HbA1c_percent", "DMPMH_SR", "SBP", "DBP",
        "BPM",
    ]
    return d[keep].copy()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Prepared .dta or .csv source table")
    parser.add_argument("output", type=Path, help="Output .csv file")
    args = parser.parse_args()
    result = score_le8(read_source(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"Wrote {len(result):,} participants to {args.output}")


if __name__ == "__main__":
    main()
