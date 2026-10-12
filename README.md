# LE8, Incident CKD, and Multimodal Brain Aging

Code for the UK Biobank analyses reported in the accompanying manuscript.

## Contents

```text
R/analysis/       statistical analyses
R/brain_age/      fixed brain-age split assignment
python/cohort/    LE8 scoring
python/kidney/    eGFR and kidney-function analyses
python/brain_age/ T1 and diffusion brain-age models
figures/          main and supplementary figures
models/           fitted combined, T1, and diffusion model objects
config/           path template
data_dictionary/  required input variables
```

UK Biobank participant-level data, identifiers, predictions, clinical records, and MRI feature matrices are not included.

## Environment

- R 4.x with `survival` and `splines`
- Python 3.12
- Python packages in `requirements.txt`

## Run

Copy `config/paths.example.json` to `config/paths.local.json`, update the local paths, and run:

```powershell
.\run_pipeline.ps1 -Mode Analysis -Config config\paths.local.json -OutputRoot results
.\run_pipeline.ps1 -Mode ScoreLE8 -Config config\paths.local.json -OutputRoot results
.\run_pipeline.ps1 -Mode AssignSplits -Config config\paths.local.json -OutputRoot results
.\run_pipeline.ps1 -Mode BrainAge -Config config\paths.local.json -OutputRoot results
.\run_pipeline.ps1 -Mode BuildKidney -Config config\paths.local.json -OutputRoot results
.\run_pipeline.ps1 -Mode Figures -Config config\paths.local.json -OutputRoot results
```

The analysis scripts use prepared cohort tables. Raw UK Biobank data must be obtained through an approved application.

The diet component was operationalized from the prespecified 16-item Oxford WebQ-based dietary score used in this study and mapped to the 0–100 LE8 component scale.

`python/kidney/egfr_ckd_epi_2021.py` implements the 2021 race-free CKD-EPI creatinine–cystatin C equation.

`python/kidney/23_build_composite_ckd_endpoint.py` combines code-qualified primary-care, hospital, and death-record events with repeated eGFR and UACR abnormalities into one incident CKD endpoint. It outputs only the unified endpoint, event date, follow-up time, and death-competing indicator.

The fitted combined, T1, and diffusion model objects are provided in `models/`. `python/brain_age/21_apply_combined_brain_age_model.py` reproduces the combined candidate-model predictions used to derive the primary out-of-sample BAG.

The BrainAge pipeline applies the released combined model before fitting the modality-specific models. The modality training script uses a minimal in-process TensorFlow compatibility module when TensorFlow is unavailable because FeatureWiz imports TensorFlow interfaces; TensorFlow is not used for model fitting. Figure 2 panels A and B use the summary tables in `figures/data/`. Panels C through G require the training-stage `candidate_predictions_all.csv` in the approved UK Biobank environment. The output of `21_apply_combined_brain_age_model.py` is used for downstream BAG analyses and is not a replacement for that training-stage file because it does not contain the original split labels. The required Figure 2 columns are documented in `data_dictionary/input_schema.md`. Participant-level predictions are not included in this release. XGBoost and LightGBM were comparator algorithms, whereas the released fitted models are the selected LASSO models.

## Randomness

- Brain-age split assignment and LASSO model fitting: seed 20260714.
- Absolute-risk coefficient simulation: seed 20260721.
- Standardized modality comparison: seed 20260722 with 500 bootstrap repetitions.
- Sex-specific calibration: seed 20260714, supplied by `run_pipeline.ps1`.
