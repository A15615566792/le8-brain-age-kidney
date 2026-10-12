# Brain-age model card

## Included artifacts

- `lasso_featurewiz_combined.joblib`
- `metadata_combined.json`
- `combined_preprocessing.csv`
- `lasso_featurewiz_t1.joblib`
- `metadata_t1.joblib`
- `lasso_featurewiz_diffusion.joblib`
- `metadata_diffusion.joblib`
- `SHA256SUMS.txt`

The combined candidate model is the fitted LASSO model used to generate the primary out-of-sample combined BAG. It uses 823 FeatureWiz-selected features from the 1,816-feature T1 and diffusion input set. The preprocessing table contains the training-sample medians, centers, scales, and selected-feature order. The metadata file contains the model, calibration, and software-version parameters.

The T1 and diffusion artifacts contain the fitted modality-specific LASSO models and their preprocessing and calibration metadata. No participant identifiers, predictions, clinical records, or UK Biobank feature matrices are included.

## Intended use

Use `python/brain_age/21_apply_combined_brain_age_model.py` with the ordered UK Biobank imaging-derived phenotypes to reproduce the combined candidate-model predictions and age-bias-corrected BAG. The T1 and diffusion artifacts support the modality-specific analyses reported in the manuscript.

## Compatibility

The combined model was trained with Python 3.12, NumPy 1.26.4, scikit-learn 1.5.2, and joblib 1.5.3. Model artifacts should only be loaded from this trusted release. The models were developed in UK Biobank and are not intended for clinical use.
