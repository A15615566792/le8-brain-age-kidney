from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def egfr_2021_creatinine_cystatin_c(
    creatinine: pd.Series,
    cystatin_c_mg_l: pd.Series,
    age_years: pd.Series,
    sex: pd.Series,
    creatinine_unit: str,
) -> pd.Series:
    scr = pd.to_numeric(creatinine, errors="coerce").astype(float)
    if creatinine_unit == "umol/L":
        scr = scr / 88.4
    scys = pd.to_numeric(cystatin_c_mg_l, errors="coerce").astype(float)
    age = pd.to_numeric(age_years, errors="coerce").astype(float)
    sex_text = sex.astype("string").str.strip().str.lower().str.replace(r"\.0$", "", regex=True)
    female = sex_text.isin(["female", "f", "0"])
    male = sex_text.isin(["male", "m", "1"])

    valid = scr.gt(0) & scys.gt(0) & age.ge(18) & (female | male)
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.219, -0.144)
    scr_ratio = scr / kappa
    cys_ratio = scys / 0.8
    result = (
        135.0
        * np.minimum(scr_ratio, 1.0) ** alpha
        * np.maximum(scr_ratio, 1.0) ** -0.544
        * np.minimum(cys_ratio, 1.0) ** -0.323
        * np.maximum(cys_ratio, 1.0) ** -0.778
        * 0.9961**age
        * np.where(female, 0.963, 1.0)
    )
    return pd.Series(np.where(valid, result, np.nan), index=scr.index, name="egfr_2021_ckd_epi_cr_cys")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--creatinine-column", default="Creatinine_0")
    parser.add_argument("--cystatin-column", default="CYSC_0")
    parser.add_argument("--age-column", default="Age_at_recruitment")
    parser.add_argument("--sex-column", default="Sex")
    parser.add_argument("--creatinine-unit", choices=["umol/L", "mg/dL"], default="umol/L")
    parser.add_argument("--urine-albumin-column", default="malb_0")
    parser.add_argument("--urine-creatinine-column", default="ucr_0")
    args = parser.parse_args()

    data = pd.read_csv(args.input_csv)
    data["egfr_2021_ckd_epi_cr_cys"] = egfr_2021_creatinine_cystatin_c(
        data[args.creatinine_column],
        data[args.cystatin_column],
        data[args.age_column],
        data[args.sex_column],
        args.creatinine_unit,
    )
    albumin = pd.to_numeric(data[args.urine_albumin_column], errors="coerce")
    urine_creatinine = pd.to_numeric(data[args.urine_creatinine_column], errors="coerce")
    uacr = 8840.0 * albumin / urine_creatinine
    uacr = uacr.where(np.isfinite(uacr) & uacr.ge(0))
    data["uacr_mg_g"] = uacr
    data["egfr_cys_do"] = data["egfr_2021_ckd_epi_cr_cys"]
    data["uacr"] = data["uacr_mg_g"]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
