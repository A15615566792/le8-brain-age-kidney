from __future__ import annotations

import gc
import importlib.util
import json
import math
import os
import sys
import types
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_score


if len(sys.argv) < 4:
    raise SystemExit(
        "Usage: python 22_train_modality_brain_age.py combined_mri.dta "
        "candidate_predictions_all.csv output_dir"
    )

input_dta = Path(sys.argv[1])
combined_predictions_path = Path(sys.argv[2])
out_dir = Path(sys.argv[3])
out_dir.mkdir(parents=True, exist_ok=True)
(out_dir / "models").mkdir(exist_ok=True)
(out_dir / "tables").mkdir(exist_ok=True)
(out_dir / "logs").mkdir(exist_ok=True)

SEED = 20260714
MAX_MISSING_PROP = 0.20
CV_FOLDS = 5
ALPHAS = np.geomspace(0.003, 0.20, 10)


def install_featurewiz_tensorflow_stub() -> None:
    if "tensorflow" in sys.modules or importlib.util.find_spec("tensorflow") is not None:
        return
    tf = types.ModuleType("tensorflow")
    keras = types.ModuleType("tensorflow.keras")
    layers = types.ModuleType("tensorflow.keras.layers")
    layers.Layer = type("Layer", (object,), {})
    keras.layers = layers
    tf.keras = keras
    sys.modules["tensorflow"] = tf
    sys.modules["tensorflow.keras"] = keras
    sys.modules["tensorflow.keras.layers"] = layers


install_featurewiz_tensorflow_stub()


def metrics(observed: np.ndarray, predicted: np.ndarray) -> dict:
    ok = np.isfinite(observed) & np.isfinite(predicted)
    y = observed[ok]
    p = predicted[ok]
    if y.size < 3:
        return {"N": int(y.size), "MAE": np.nan, "RMSE": np.nan,
                "R2": np.nan, "Pearson_r": np.nan}
    return {
        "N": int(y.size),
        "MAE": float(mean_absolute_error(y, p)),
        "RMSE": float(math.sqrt(mean_squared_error(y, p))),
        "R2": float(r2_score(y, p)),
        "Pearson_r": float(np.corrcoef(y, p)[0, 1]),
    }


def fit_age_bias(age: np.ndarray, predicted: np.ndarray) -> tuple[float, float]:
    slope, intercept = np.polyfit(age, predicted, 1)
    if not np.isfinite(slope) or abs(slope) < 0.05:
        raise RuntimeError(f"Invalid age-bias slope: {slope}")
    return float(intercept), float(slope)


def correct_age_bias(predicted: np.ndarray, intercept: float, slope: float) -> np.ndarray:
    return (predicted - intercept) / slope


def evaluate_splits(meta: pd.DataFrame, raw: np.ndarray, corrected: np.ndarray,
                    modality: str, n_features: int, alpha: float,
                    cv_mae: float) -> pd.DataFrame:
    rows = []
    for split in meta["split"].dropna().unique():
        idx = meta["split"].eq(split).to_numpy()
        for correction, pred in (("Raw", raw), ("Bias-corrected", corrected)):
            row = metrics(meta.loc[idx, "age"].to_numpy(float), pred[idx])
            row.update({"modality": modality, "split": split,
                        "correction": correction, "n_selected_features": n_features,
                        "alpha": alpha, "best_cv_mae": cv_mae})
            rows.append(row)
    return pd.DataFrame(rows)


def fit_modality(modality: str, feature_names: list[str], split_meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    checkpoint = out_dir / "tables" / f"predictions_{modality.lower()}.csv"
    performance_file = out_dir / "tables" / f"performance_{modality.lower()}.csv"
    if checkpoint.exists() and performance_file.exists():
        print(f"Loading completed {modality} checkpoint", flush=True)
        saved_prediction = pd.read_csv(checkpoint)
        saved_prediction["participant_id"] = (
            saved_prediction["participant_id"].astype("Int64").astype(str)
        )
        return saved_prediction, pd.read_csv(performance_file)

    print(f"Reading {modality}: {len(feature_names)} source features", flush=True)
    dat = pd.read_stata(
        input_dta,
        columns=["participant_id", "age", "healthy"] + feature_names,
        convert_categoricals=False,
        preserve_dtypes=False,
    )
    dat["participant_id"] = dat["participant_id"].astype("Int64").astype(str)
    dat = split_meta.merge(
        dat.drop(columns=["age", "healthy"]),
        on="participant_id", how="left", validate="one_to_one"
    )
    x_raw = dat[feature_names].to_numpy(dtype=np.float32, copy=True)
    missing_prop = np.mean(~np.isfinite(x_raw), axis=1)
    if np.any(missing_prop > MAX_MISSING_PROP):
        raise RuntimeError(
            f"{modality}: modality-specific missingness >20% occurred inside the fixed "
            "combined eligible sample; same-participant comparison would be violated"
        )

    train_idx = dat["split"].eq("Healthy train").to_numpy()
    validation_idx = dat["split"].eq("Healthy validation").to_numpy()
    if train_idx.sum() < 500 or validation_idx.sum() < 100:
        raise RuntimeError(
            f"Too few fixed-split participants for {modality}: "
            f"train={train_idx.sum()}, validation={validation_idx.sum()}"
        )

    x_train = x_raw[train_idx].astype(np.float64, copy=True)
    medians = np.nanmedian(np.where(np.isfinite(x_train), x_train, np.nan), axis=0)
    medians[~np.isfinite(medians)] = 0.0
    bad = ~np.isfinite(x_train)
    if bad.any():
        x_train[bad] = np.take(medians, np.where(bad)[1])
    centers = x_train.mean(axis=0)
    scales = x_train.std(axis=0, ddof=1)
    usable = np.isfinite(scales) & (scales > 1e-8)
    if not usable.all():
        feature_names = [x for x, keep in zip(feature_names, usable) if keep]
        medians = medians[usable]
        centers = centers[usable]
        scales = scales[usable]
        x_raw = x_raw[:, usable]

    for j in range(x_raw.shape[1]):
        col = x_raw[:, j]
        miss = ~np.isfinite(col)
        if miss.any():
            col[miss] = medians[j]
        x_raw[:, j] = (col - centers[j]) / scales[j]
    x_scaled = x_raw.astype(np.float32, copy=False)
    del x_train, x_raw, bad
    gc.collect()

    print(f"Running training-only FeatureWiz for {modality}", flush=True)
    from featurewiz import FeatureWiz
    fw = FeatureWiz(
        corr_limit=0.90,
        feature_engg="",
        nrows=None,
        transform_target=False,
        scalers=None,
        category_encoders="",
        add_missing=False,
        verbose=0,
        imbalanced=False,
    )
    x_train_df = pd.DataFrame(x_scaled[train_idx], columns=feature_names)
    y_train_series = pd.Series(dat.loc[train_idx, "age"].to_numpy(float),
                               index=x_train_df.index, name="age")
    fw.fit_transform(x_train_df, y_train_series)
    selected = [str(x) for x in fw.features if str(x) in feature_names]
    if not selected:
        raise RuntimeError(f"FeatureWiz selected no features for {modality}")
    selected_index = np.array([feature_names.index(x) for x in selected], dtype=int)
    pd.DataFrame({"feature": selected}).to_csv(
        out_dir / "tables" / f"selected_features_{modality.lower()}.csv", index=False
    )
    del fw, x_train_df, y_train_series
    gc.collect()

    x_train_selected = x_scaled[train_idx][:, selected_index]
    y_train = dat.loc[train_idx, "age"].to_numpy(np.float32)
    cv = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=SEED)
    tuning_rows = []
    print(f"Tuning fixed LASSO grid for {modality}: {len(ALPHAS)} x {CV_FOLDS}", flush=True)
    for alpha in ALPHAS:
        model = Lasso(alpha=float(alpha), max_iter=20000, tol=1e-4,
                      random_state=SEED, selection="cyclic")
        scores = cross_val_score(
            model, x_train_selected, y_train, cv=cv,
            scoring="neg_mean_absolute_error", n_jobs=1
        )
        row = {"modality": modality, "alpha": float(alpha),
               "cv_mae": float(-scores.mean()), "cv_mae_sd": float(scores.std(ddof=1))}
        tuning_rows.append(row)
        print(f"{modality} alpha={alpha:.6f}, CV MAE={row['cv_mae']:.4f}", flush=True)
    tuning = pd.DataFrame(tuning_rows).sort_values(["cv_mae", "alpha"])
    tuning.to_csv(out_dir / "tables" / f"lasso_tuning_{modality.lower()}.csv", index=False)
    best_alpha = float(tuning.iloc[0]["alpha"])
    best_cv_mae = float(tuning.iloc[0]["cv_mae"])

    model = Lasso(alpha=best_alpha, max_iter=30000, tol=1e-5,
                  random_state=SEED, selection="cyclic")
    model.fit(x_train_selected, y_train)
    raw_pred = np.empty(len(dat), dtype=np.float64)
    batch = 5000
    for start in range(0, len(dat), batch):
        stop = min(start + batch, len(dat))
        raw_pred[start:stop] = model.predict(x_scaled[start:stop][:, selected_index])
    intercept, slope = fit_age_bias(dat.loc[train_idx, "age"].to_numpy(float),
                                    raw_pred[train_idx])
    corrected = correct_age_bias(raw_pred, intercept, slope)
    bag = corrected - dat["age"].to_numpy(float)

    prediction = dat[["participant_id", "age", "healthy", "split"]].copy()
    prediction[f"pred_raw_{modality.lower()}"] = raw_pred
    prediction[f"pred_corrected_{modality.lower()}"] = corrected
    prediction[f"bag_{modality.lower()}"] = bag
    prediction[f"bag_oos_eligible_{modality.lower()}"] = prediction["split"].ne("Healthy train")
    prediction.to_csv(checkpoint, index=False)

    performance = evaluate_splits(dat, raw_pred, corrected, modality,
                                  len(selected), best_alpha, best_cv_mae)
    performance.to_csv(performance_file, index=False)
    nonzero = int(np.sum(np.abs(model.coef_) > 0))
    joblib.dump(model, out_dir / "models" / f"lasso_featurewiz_{modality.lower()}.joblib")
    joblib.dump(
        {"modality": modality, "source_features": feature_names,
         "selected_features": selected, "medians": medians,
         "centers": centers, "scales": scales, "alpha": best_alpha,
         "calibration_intercept": intercept, "calibration_slope": slope,
            "n_nonzero_coefficients": nonzero, "seed": SEED},
        out_dir / "models" / f"metadata_{modality.lower()}.joblib",
        compress=3,
    )
    print(
        f"Completed {modality}: selected={len(selected)}, nonzero={nonzero}, "
        f"alpha={best_alpha:.6f}, validation raw MAE="
        f"{performance.loc[(performance['split']=='Healthy validation') & (performance['correction']=='Raw'),'MAE'].iloc[0]:.4f}",
        flush=True,
    )
    del dat, x_scaled, x_train_selected, model, prediction
    gc.collect()
    saved_prediction = pd.read_csv(checkpoint)
    saved_prediction["participant_id"] = (
        saved_prediction["participant_id"].astype("Int64").astype(str)
    )
    return saved_prediction, performance


if not input_dta.exists():
    raise FileNotFoundError(input_dta)
if not combined_predictions_path.exists():
    raise FileNotFoundError(combined_predictions_path)

reader = pd.io.stata.StataReader(input_dta)
all_columns = list(reader.variable_labels().keys())
meta_columns = ["participant_id", "age", "healthy"]
feature_columns = [x for x in all_columns if x not in meta_columns]
first_diffusion = "wt_fa_tar_l_2"
if first_diffusion not in feature_columns:
    raise RuntimeError(f"Validated first diffusion field not found: {first_diffusion}")
boundary = feature_columns.index(first_diffusion)
t1_features = feature_columns[:boundary]
diffusion_features = feature_columns[boundary:]
if len(feature_columns) != 1816 or len(t1_features) != 1142 or len(diffusion_features) != 674:
    raise RuntimeError(
        f"Feature-map audit failed: total={len(feature_columns)}, "
        f"T1={len(t1_features)}, Diffusion={len(diffusion_features)}"
    )

feature_map = pd.DataFrame({
    "feature": feature_columns,
    "modality": ["T1"] * len(t1_features) + ["Diffusion"] * len(diffusion_features),
    "source_order": np.arange(1, len(feature_columns) + 1),
})
feature_map.to_csv(out_dir / "tables" / "modality_feature_map.csv", index=False)

combined = pd.read_csv(combined_predictions_path)
combined["participant_id"] = combined["participant_id"].astype("Int64").astype(str)
split_meta = combined[["participant_id", "age", "healthy", "split"]].copy()
if split_meta.empty:
    raise RuntimeError("No participants were available in the combined prediction table")
required_splits = {"Healthy train", "Healthy validation", "Non-healthy test"}
observed_splits = set(split_meta["split"].dropna().astype(str))
if not required_splits.issubset(observed_splits):
    raise RuntimeError(f"Required split labels are missing: {sorted(required_splits - observed_splits)}")


def audit_modality_missingness(modality: str, feature_names: list[str]) -> pd.DataFrame:
    print(f"Auditing {modality} missingness for the common-sample definition", flush=True)
    dd = pd.read_stata(
        input_dta,
        columns=["participant_id"] + feature_names,
        convert_categoricals=False,
        preserve_dtypes=False,
    )
    dd["participant_id"] = dd["participant_id"].astype("Int64").astype(str)
    x = dd[feature_names].to_numpy(dtype=np.float32, copy=False)
    ans = dd[["participant_id"]].copy()
    ans[f"missing_prop_{modality.lower()}"] = np.mean(~np.isfinite(x), axis=1)
    del dd, x
    gc.collect()
    return ans


t1_missing = audit_modality_missingness("T1", t1_features)
dwi_missing = audit_modality_missingness("Diffusion", diffusion_features)
common_qc = split_meta.merge(t1_missing, on="participant_id", how="left", validate="one_to_one")
common_qc = common_qc.merge(dwi_missing, on="participant_id", how="left", validate="one_to_one")
common_qc["common_modality_eligible"] = (
    common_qc["missing_prop_t1"].le(MAX_MISSING_PROP)
    & common_qc["missing_prop_diffusion"].le(MAX_MISSING_PROP)
)
common_qc.to_csv(out_dir / "tables" / "modality_common_sample_QC.csv", index=False)
split_meta_original = split_meta.copy()
split_meta = common_qc.loc[
    common_qc["common_modality_eligible"],
    ["participant_id", "age", "healthy", "split"],
].copy()
print(
    f"Common modality sample: {len(split_meta)} / {len(split_meta_original)}; "
    f"healthy train={split_meta['split'].eq('Healthy train').sum()}, "
    f"validation={split_meta['split'].eq('Healthy validation').sum()}",
    flush=True,
)

t1_pred, t1_perf = fit_modality("T1", t1_features, split_meta)
dwi_pred, dwi_perf = fit_modality("Diffusion", diffusion_features, split_meta)

combined_small = combined[combined["participant_id"].isin(split_meta["participant_id"])][["participant_id", "age", "healthy", "split",
                           "pred_raw__lasso__featurewiz",
                           "pred_corrected__lasso__featurewiz",
                           "bag__lasso__featurewiz"]].copy()
combined_small = combined_small.rename(columns={
    "pred_raw__lasso__featurewiz": "pred_raw_combined",
    "pred_corrected__lasso__featurewiz": "pred_corrected_combined",
    "bag__lasso__featurewiz": "bag_combined",
})
combined_small["bag_oos_eligible_combined"] = combined_small["split"].ne("Healthy train")

pred_all = combined_small.merge(
    t1_pred.drop(columns=["age", "healthy", "split"]), on="participant_id",
    how="left", validate="one_to_one"
).merge(
    dwi_pred.drop(columns=["age", "healthy", "split"]), on="participant_id",
    how="left", validate="one_to_one"
)
pred_all.to_csv(out_dir / "tables" / "modality_brain_age_predictions.csv", index=False)

combined_perf_rows = []
for split in combined_small["split"].dropna().unique():
    idx = combined_small["split"].eq(split).to_numpy()
    for correction, col in (("Raw", "pred_raw_combined"),
                            ("Bias-corrected", "pred_corrected_combined")):
        row = metrics(combined_small.loc[idx, "age"].to_numpy(float),
                      combined_small.loc[idx, col].to_numpy(float))
        row.update({"modality": "Combined", "split": split,
                    "correction": correction, "n_selected_features": 823,
                    "alpha": 0.0122258374448315,
                    "best_cv_mae": 3.01748563990987})
        combined_perf_rows.append(row)
performance_all = pd.concat(
    [t1_perf, dwi_perf, pd.DataFrame(combined_perf_rows)], ignore_index=True
)
performance_all.to_csv(out_dir / "tables" / "modality_model_performance.csv", index=False)

correlations = []
for split in ["All eligible", "Healthy validation", "Non-healthy test",
              "Outside healthy train"]:
    if split == "All eligible":
        dd = pred_all
    elif split == "Outside healthy train":
        dd = pred_all[pred_all["split"].ne("Healthy train")]
    else:
        dd = pred_all[pred_all["split"].eq(split)]
    for a, b in (("bag_t1", "bag_diffusion"),
                 ("bag_t1", "bag_combined"),
                 ("bag_diffusion", "bag_combined")):
        ok = dd[a].notna() & dd[b].notna()
        correlations.append({"split": split, "bag_x": a, "bag_y": b,
                             "N": int(ok.sum()),
                             "Pearson_r": float(dd.loc[ok, [a, b]].corr().iloc[0, 1])})
pd.DataFrame(correlations).to_csv(
    out_dir / "tables" / "modality_BAG_correlations.csv", index=False
)

audit = {
    "prepared_file_total_columns": len(all_columns),
    "metadata_columns": len(meta_columns),
    "actual_MRI_predictors": len(feature_columns),
    "T1_predictors": len(t1_features),
    "Diffusion_predictors": len(diffusion_features),
    "original_combined_eligible_participants": len(split_meta_original),
    "fixed_common_modality_participants": len(split_meta),
    "healthy_train": int(split_meta["split"].eq("Healthy train").sum()),
    "healthy_validation": int(split_meta["split"].eq("Healthy validation").sum()),
    "nonhealthy_test": int(split_meta["split"].eq("Non-healthy test").sum()),
    "health_unknown": int(split_meta["split"].eq("Health status unknown").sum()),
    "algorithm": "LASSO",
    "selector": "FeatureWiz fitted independently within each modality on healthy training only",
    "cv": f"{len(ALPHAS)} alpha values x {CV_FOLDS} folds, MAE selection",
    "combined_model": "Existing LASSO+FeatureWiz candidate retained without refitting",
}
with open(out_dir / "modality_brain_age_audit.json", "w", encoding="utf-8") as f:
    json.dump(audit, f, ensure_ascii=False, indent=2)

print("Modality-specific brain-age training completed", flush=True)
