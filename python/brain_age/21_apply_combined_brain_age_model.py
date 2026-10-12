from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


def read_table(path: Path, columns: list[str]) -> pd.DataFrame:
    if path.suffix.lower() == ".dta":
        return pd.read_stata(path, columns=columns, convert_categoricals=False)
    return pd.read_csv(path, usecols=columns, low_memory=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--model-dir", type=Path, default=Path(__file__).resolve().parents[2] / "models")
    args = parser.parse_args()

    prep = pd.read_csv(args.model_dir / "combined_preprocessing.csv")
    metadata = json.loads((args.model_dir / "metadata_combined.json").read_text(encoding="utf-8"))
    features = prep["feature"].astype(str).tolist()
    data = read_table(args.input, ["participant_id", "age", *features])
    x = data[features].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    eligible = np.isfinite(pd.to_numeric(data["age"], errors="coerce")) & (
        np.mean(~np.isfinite(x), axis=1) <= metadata["maximum_mri_missing_proportion"]
    )
    medians = prep["median"].to_numpy(dtype=np.float32)
    centers = prep["center"].to_numpy(dtype=np.float32)
    scales = prep["scale"].to_numpy(dtype=np.float32)
    x = np.where(np.isfinite(x), x, medians)
    x = (x - centers) / scales
    selected = prep.loc[prep["selected_order"].notna()].sort_values("selected_order")
    positions = prep.index.get_indexer(selected.index)
    model = joblib.load(args.model_dir / "lasso_featurewiz_combined.joblib")
    raw = np.full(len(data), np.nan)
    raw[eligible] = model.predict(x[eligible][:, positions])
    corrected = (raw - metadata["age_bias_intercept"]) / metadata["age_bias_slope"]
    age = pd.to_numeric(data["age"], errors="coerce").to_numpy()
    result = pd.DataFrame({
        "participant_id": data["participant_id"],
        "age": age,
        "predicted_age_raw_combined": raw,
        "predicted_age_corrected_combined": corrected,
        "brain_age_gap_oos_candidate": corrected - age,
        "mri_model_eligible": eligible,
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
