param(
    [ValidateSet('Analysis','ScoreLE8','AssignSplits','BrainAge','BuildKidney','Figures')]
    [string]$Mode = 'Analysis',
    [string]$Config = (Join-Path $PSScriptRoot 'config\paths.local.json'),
    [string]$OutputRoot = (Join-Path $PSScriptRoot 'results')
)

$ErrorActionPreference = 'Stop'
$release = $PSScriptRoot
$logs = Join-Path $OutputRoot 'logs'
New-Item -ItemType Directory -Force -Path $OutputRoot,$logs | Out-Null

function Resolve-Executable([string]$preferred, [string]$fallback) {
    if ($preferred -and (Test-Path -LiteralPath $preferred -PathType Leaf)) { return $preferred }
    $cmd = Get-Command $fallback -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    throw "Executable not found: $fallback"
}

function Assert-File([string]$path, [string]$label) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Missing ${label}: $path" }
}

function Invoke-Step([string]$name, [string]$exe, [string[]]$arguments) {
    $log = Join-Path $logs ($name + '.log')
    "[$(Get-Date -Format o)] START $name" | Tee-Object -FilePath $log
    & $exe @arguments *>&1 | Tee-Object -FilePath $log -Append
    if ($LASTEXITCODE -ne 0) { throw "$name failed; see $log" }
    "[$(Get-Date -Format o)] END $name" | Tee-Object -FilePath $log -Append
}

if (-not (Test-Path -LiteralPath $Config -PathType Leaf)) {
    throw "Create config\paths.local.json from config\paths.example.json first."
}
$cfg = Get-Content -LiteralPath $Config -Raw -Encoding UTF8 | ConvertFrom-Json
$Rscript = Resolve-Executable $cfg.rscript 'Rscript'
$Python = Resolve-Executable $cfg.python 'python'
$env:PYTHONUTF8 = '1'

if ($Mode -eq 'ScoreLE8') {
    Assert-File $cfg.le8_source_table 'prepared LE8 source table'
    Invoke-Step 'score_le8' $Python @(
        (Join-Path $release 'python\cohort\01_score_le8.py'),
        $cfg.le8_source_table,(Join-Path $OutputRoot 'cohort\le8_scores.csv')
    )
    exit 0
}

if ($Mode -eq 'AssignSplits') {
    Assert-File $cfg.brain_age_health_source 'brain-age health-status source table'
    Invoke-Step 'assign_brain_age_splits' $Rscript @(
        (Join-Path $release 'R\brain_age\20_assign_brain_age_splits.R'),
        $cfg.brain_age_health_source,(Join-Path $OutputRoot 'brain_age\split_labels.csv'),'20260714'
    )
    exit 0
}

if ($Mode -eq 'BrainAge') {
    Assert-File $cfg.combined_mri_features 'combined MRI feature table'
    Assert-File $cfg.combined_candidate_predictions 'combined-model split/prediction metadata'
    Invoke-Step 'brain_age_combined_predictions' $Python @(
        (Join-Path $release 'python\brain_age\21_apply_combined_brain_age_model.py'),
        $cfg.combined_mri_features,(Join-Path $OutputRoot 'brain_age\combined_predictions.csv')
    )
    Invoke-Step 'brain_age_modality_models' $Python @(
        (Join-Path $release 'python\brain_age\22_train_modality_brain_age.py'),
        $cfg.combined_mri_features,$cfg.combined_candidate_predictions,
        (Join-Path $OutputRoot 'brain_age')
    )
    exit 0
}

if ($Mode -eq 'BuildKidney') {
    Assert-File $cfg.ckd_participants 'CKD participant index/censoring table'
    Assert-File $cfg.ckd_coded_events 'code-qualified CKD event table'
    Assert-File $cfg.ckd_laboratory_measurements 'longitudinal kidney laboratory table'
    Invoke-Step 'build_composite_ckd_endpoint' $Python @(
        (Join-Path $release 'python\kidney\23_build_composite_ckd_endpoint.py'),
        $cfg.ckd_participants,$cfg.ckd_coded_events,$cfg.ckd_laboratory_measurements,
        (Join-Path $OutputRoot 'kidney\composite_ckd_endpoint.csv')
    )
    exit 0
}

if ($Mode -eq 'Figures') {
    $figureInputs = @(
        'cohort_flow_detailed.csv',
        'Table_3_LE8_component_profiles.csv',
        'Table_P9_standardized_modality_specific_slopes.csv',
        'Table_P10_continuous_BP_glycaemia_modality_slopes.csv',
        'Table_P11_diffusion_metric_composites.csv',
        'Table_P11_significant_tracts_FDR134.csv',
        'Table_P12_sex_specific_slopes.csv',
        'modality_model_performance.csv',
        'Table_P4_CKD_adjusted_FG_CIF.csv',
        'Table_P2_total_LE8_CKD_sensitivities.csv'
    )
    foreach ($file in $figureInputs) {
        Assert-File (Join-Path $cfg.figure_table_root "tables\supplement\$file") "figure input $file"
    }
    Assert-File (Join-Path $cfg.brain_age_figure_source 'tables\candidate_predictions_all.csv') 'restricted participant-level brain-age prediction table'
    Invoke-Step 'association_figures' $Python @(
        (Join-Path $release 'figures\02_make_figures_v2.py'),
        '--package-root',$cfg.figure_table_root,
        '--output-root',(Join-Path $OutputRoot 'figures')
    )
    Invoke-Step 'brain_age_figure' $Python @(
        (Join-Path $release 'figures\05_make_brain_age_development_figure.py'),
        '--package-root',$release,
        '--source-root',$cfg.brain_age_figure_source,
        '--output-root',(Join-Path $OutputRoot 'figures\main')
    )
    exit 0
}

foreach ($item in @(
    @($cfg.frozen_full_cohort,'full cohort'),
    @($cfg.frozen_mri_cohort,'MRI cohort'),
    @($cfg.cohort_flow,'cohort flow'),
    @($cfg.modality_predictions,'modality predictions'),
    @($cfg.egfr_progression_subset,'eGFR progression subset'),
    @($cfg.diffusion_tract_features,'diffusion tract features'),
    @($cfg.diffusion_tract_map,'diffusion tract map')
)) { Assert-File $item[0] $item[1] }

$prepared=Join-Path $OutputRoot 'prepared'; New-Item -ItemType Directory -Force -Path $prepared | Out-Null
$full=Join-Path $prepared 'full_cohort_kidney_covariates.csv'
$mri=Join-Path $prepared 'mri_cohort_kidney_covariates.csv'
Invoke-Step 'prepare_full_kidney_covariates' $Python @((Join-Path $release 'python\kidney\egfr_ckd_epi_2021.py'),$cfg.frozen_full_cohort,$full)
Invoke-Step 'prepare_mri_kidney_covariates' $Python @((Join-Path $release 'python\kidney\egfr_ckd_epi_2021.py'),$cfg.frozen_mri_cohort,$mri)
$p1=Join-Path $OutputRoot 'P1_core'; $p2=Join-Path $OutputRoot 'P2_sensitivity'
$p3=Join-Path $OutputRoot 'P3_targeted'; $p4=Join-Path $OutputRoot 'P4_tables_absolute_risk'
Invoke-Step 'P1_core' $Rscript @((Join-Path $release 'R\analysis\18_run_p1_core_profiles.R'),$full,$mri,$p1)
Invoke-Step 'P2_sensitivity' $Rscript @((Join-Path $release 'R\analysis\19_run_p2_sensitivity.R'),$full,$mri,$p2)
Invoke-Step 'P3_competing_risk' $Rscript @((Join-Path $release 'R\analysis\20_run_competing_risk.R'),$full,$p3)
Invoke-Step 'P4_tables_absolute_risk' $Rscript @((Join-Path $release 'R\analysis\21_run_ph_absolute_tables.R'),$full,$mri,$p1,$p2,$p3,$cfg.cohort_flow,$p4)
Invoke-Step 'P5_modality_associations' $Rscript @((Join-Path $release 'R\analysis\23_run_modality_associations.R'),$full,$mri,$cfg.modality_predictions,(Join-Path $OutputRoot 'P5_modality_associations'))
Invoke-Step 'P8_kidney_progression' $Rscript @((Join-Path $release 'R\analysis\25_run_kidney_progression.R'),$full,$cfg.egfr_progression_subset,(Join-Path $OutputRoot 'P8_kidney_progression'))
Invoke-Step 'P9_P10_standardized_modality' $Rscript @((Join-Path $release 'R\analysis\26_run_standardized_modality_continuous.R'),$full,$mri,$cfg.modality_predictions,(Join-Path $OutputRoot 'P9_P10_standardized_modality'))
Invoke-Step 'P11_white_matter' $Rscript @((Join-Path $release 'R\analysis\27_run_white_matter_localization.R'),$full,$mri,$cfg.diffusion_tract_features,$cfg.diffusion_tract_map,(Join-Path $OutputRoot 'P11_white_matter'))
Invoke-Step 'P12_sex_calibration' $Rscript @((Join-Path $release 'R\analysis\28_run_sex_calibration.R'),$full,$mri,$cfg.modality_predictions,(Join-Path $OutputRoot 'P12_sex_calibration'),'20260714')

Write-Host "Analysis completed: $OutputRoot"
