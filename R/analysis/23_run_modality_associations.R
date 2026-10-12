args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: Rscript 23_run_modality_associations.R full.csv mri.csv modality_predictions.csv out_dir")
}
full_path <- args[1]
mri_path <- args[2]
pred_path <- args[3]
out_dir <- args[4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
suppressPackageStartupMessages(library(splines))

save_csv <- function(x, n) write.csv(x, file.path(out_dir, n), row.names = FALSE,
                                      na = "", fileEncoding = "UTF-8")
z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
eth4 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  w <- x %in% c("british", "irish", "white", "any other white background")
  a <- x %in% c("indian", "pakistani", "bangladeshi", "chinese",
                "asian or asian british", "any other asian background")
  b <- x %in% c("caribbean", "african", "black or black british",
                "any other black background")
  factor(ifelse(w, "White", ifelse(a, "Asian", ifelse(b, "Black", "Other/unknown"))),
         levels = c("White", "Asian", "Black", "Other/unknown"))
}

vcov_hcw <- function(f) {
  X <- model.matrix(f)
  k <- !is.na(coef(f))
  X <- X[, k, drop = FALSE]
  u <- residuals(f)
  h <- pmin(hatvalues(f), 0.999999)
  w <- weights(f)
  if (is.null(w)) w <- rep(1, nrow(X))
  B <- tryCatch(solve(crossprod(X, X * w)), error = function(e) NULL)
  if (is.null(B)) return(vcov(f)[k, k, drop = FALSE])
  S <- X * (w * u / (1 - h))
  V <- B %*% crossprod(S) %*% B
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

tidy_lm <- function(f, term, modality, model, component, weighted = FALSE) {
  V <- vcov_hcw(f)
  b <- coef(f)[term]
  se <- sqrt(V[term, term])
  data.frame(modality = modality, model = model, component = component,
             term = term, beta_years = unname(b), se = unname(se),
             CI_low = unname(b - 1.96 * se), CI_high = unname(b + 1.96 * se),
             p = unname(2 * pnorm(-abs(b / se))), n = nobs(f),
             selection_IPW = weighted, stringsAsFactors = FALSE)
}

cat("Reading frozen cohorts and modality predictions...\n")
full <- read.csv(full_path, check.names = FALSE, na.strings = c("", "NA"))
mri <- read.csv(mri_path, check.names = FALSE, na.strings = c("", "NA"))
pred <- read.csv(pred_path, check.names = FALSE, na.strings = c("", "NA"))
full$participant_id <- as.character(full$participant_id)
mri$participant_id <- as.character(mri$participant_id)
pred$participant_id <- as.character(pred$participant_id)
stopifnot(!anyDuplicated(full$participant_id), !anyDuplicated(mri$participant_id),
          !anyDuplicated(pred$participant_id))

full$baseline_date_r <- as.Date(substr(full$baseline_date, 1, 10))
full$time0_full_r <- as.Date(substr(full$time0_full, 1, 10))
full$age_webq <- as.numeric(full$Age_at_recruitment) +
  as.numeric(full$time0_full_r - full$baseline_date_r) / 365.25
full$sex_factor <- relevel(factor(full$Sex), ref = "Female")
full$ethnicity4 <- eth4(full$Ethnic_0)
full$tdi <- as.numeric(full$TDI_at_recruitment)
full$edu_z <- as.numeric(full$area_education_deprivation_z)
full$edu_missing <- factor(ifelse(is.na(full$edu_z), "Missing", "Observed"),
                           levels = c("Observed", "Missing"))
full$edu_imp <- ifelse(is.na(full$edu_z), 0, full$edu_z)
full$baseline_center_factor <- factor(full$baseline_center)
full$energy_z <- z(full$energy_kcal_primary)
full$rounds <- as.numeric(full$webq_n_dates)
full$LE8_per10 <- as.numeric(full$LE8) / 10
full$egfr <- as.numeric(full$egfr_cys_do)

mri <- merge(mri, pred, by = "participant_id", all = FALSE, sort = FALSE,
             suffixes = c("", "_pred"))
mri$sex_factor <- relevel(factor(mri$Sex), ref = "Female")
mri$ethnicity4 <- eth4(mri$Ethnic_0)
mri$tdi <- as.numeric(mri$TDI_at_recruitment)
mri$edu_z <- as.numeric(mri$area_education_deprivation_z)
mri$edu_missing <- factor(ifelse(is.na(mri$edu_z), "Missing", "Observed"),
                          levels = c("Observed", "Missing"))
mri$edu_imp <- ifelse(is.na(mri$edu_z), 0, mri$edu_z)
mri$center <- factor(mri$v4)
mri$energy_z <- z(mri$energy_kcal_pre_mri_mean)
mri$gap_years <- as.numeric(mri$webq_mri_gap_days) / 365.25
mri$v_scaling_z <- z(mri$v_scaling_t1_to_standard_2)
mri$scanner_x_z <- z(mri$scanner_lateral_x_brain_2)
mri$scanner_y_z <- z(mri$scanner_transverse_y_brain_2)
mri$scanner_z_z <- z(mri$scanner_longitudinal_z_brain_2)
mri$scanner_table_z <- z(mri$scanner_table_position_2)
mri$t2flair_qc <- factor(mri$t2flair_run_freesurfer_2)
mri$LE8_per10 <- as.numeric(mri$LE8) / 10

if (!all(c("bag_t1", "bag_diffusion", "bag_combined") %in% names(mri))) {
  stop("Missing modality BAG variables after merge")
}
outside_train <- mri$split != "Healthy train"
if (any(is.na(outside_train))) stop("Missing split label in formal MRI merge")
mri_oos <- mri[outside_train, , drop = FALSE]

delta <- mri_oos$bag_combined - mri_oos$brain_age_gap_oos_candidate
bag_validation <- data.frame(
  formal_MRI_merged = nrow(mri),
  outside_healthy_train = nrow(mri_oos),
  all_three_BAG_available = sum(complete.cases(mri_oos[, c("bag_t1", "bag_diffusion", "bag_combined")])),
  combined_vs_frozen_correlation = cor(mri_oos$bag_combined,
                                       mri_oos$brain_age_gap_oos_candidate,
                                       use = "complete.obs"),
  combined_vs_frozen_max_abs_difference = max(abs(delta), na.rm = TRUE),
  stringsAsFactors = FALSE)
save_csv(bag_validation, "Table_P5_modality_merge_validation.csv")
if (bag_validation$combined_vs_frozen_correlation < 0.999999) {
  stop("Combined candidate BAG does not reproduce the frozen OOS candidate BAG")
}

full$in_mri <- as.integer(full$participant_id %in% mri$participant_id)
selection_fit <- glm(
  in_mri ~ ns(age_webq, df = 4) + sex_factor + ethnicity4 + tdi + edu_imp +
    edu_missing + baseline_center_factor + LE8_per10 + egfr,
  data = full, family = binomial())
full$psel <- predict(selection_fit, newdata = full, type = "response")
pm <- mean(full$in_mri)
wmap <- full[, c("participant_id", "psel")]
mri_oos <- merge(mri_oos, wmap, by = "participant_id", all.x = TRUE, sort = FALSE)
mri_oos$sw <- pm / mri_oos$psel
q <- quantile(mri_oos$sw, c(0.01, 0.99), na.rm = TRUE)
mri_oos$sw_trim <- pmin(pmax(mri_oos$sw, q[1]), q[2])

bag_cov <- paste(c("age", "sex_factor", "ethnicity4", "tdi", "edu_imp",
                   "edu_missing", "center", "energy_z", "gap_years",
                   "v_scaling_z", "scanner_x_z", "scanner_y_z", "scanner_z_z",
                   "scanner_table_z", "t2flair_qc"), collapse = "+")
modalities <- c(T1 = "bag_t1", Diffusion = "bag_diffusion", Combined = "bag_combined")
components <- c("le8_diet_primary", "le8_pa", "le8_nicotine", "le8_sleep",
                "le8_bmi", "le8_nonhdl", "le8_glucose", "le8_bp")
labels <- setNames(c("Diet", "Physical activity", "Nicotine exposure", "Sleep",
                     "BMI", "Non-HDL cholesterol", "Blood glucose",
                     "Blood pressure"), components)
for (v in components) mri_oos[[paste0(v, "_per10")]] <- as.numeric(mri_oos[[v]]) / 10
component_terms <- paste0(components, "_per10")

cat("Running total LE8 and component models for each BAG modality...\n")
total_res <- list()
component_res <- list()
difference_res <- list()
kt <- kc <- kd <- 1L
for (modality in names(modalities)) {
  outcome <- modalities[[modality]]
  f_unw <- lm(as.formula(paste(outcome, "~LE8_per10+", bag_cov)), data = mri_oos)
  f_ipw <- lm(as.formula(paste(outcome, "~LE8_per10+", bag_cov)),
              data = mri_oos, weights = sw_trim)
  total_res[[kt]] <- tidy_lm(f_unw, "LE8_per10", modality,
                             "Total LE8", "Total LE8", FALSE); kt <- kt + 1L
  total_res[[kt]] <- tidy_lm(f_ipw, "LE8_per10", modality,
                             "Total LE8", "Total LE8", TRUE); kt <- kt + 1L

  for (i in seq_along(components)) {
    term <- component_terms[i]
    fs <- lm(as.formula(paste(outcome, "~", term, "+", bag_cov)), data = mri_oos)
    fw <- lm(as.formula(paste(outcome, "~", term, "+", bag_cov)),
             data = mri_oos, weights = sw_trim)
    component_res[[kc]] <- tidy_lm(fs, term, modality, "Separate component",
                                   labels[components[i]], FALSE); kc <- kc + 1L
    component_res[[kc]] <- tidy_lm(fw, term, modality, "Separate component",
                                   labels[components[i]], TRUE); kc <- kc + 1L
  }
  fm <- lm(as.formula(paste(outcome, "~", paste(component_terms, collapse = "+"),
                            "+", bag_cov)), data = mri_oos)
  for (i in seq_along(components)) {
    component_res[[kc]] <- tidy_lm(fm, component_terms[i], modality,
                                   "Mutually adjusted 8 components",
                                   labels[components[i]], FALSE); kc <- kc + 1L
  }
}
total_table <- do.call(rbind, total_res)
total_table$p_FDR <- ave(total_table$p, total_table$selection_IPW,
                         FUN = function(x) p.adjust(x, "BH"))
component_table <- do.call(rbind, component_res)
component_table$p_FDR <- ave(component_table$p,
  interaction(component_table$modality, component_table$model,
              component_table$selection_IPW, drop = TRUE),
  FUN = function(x) p.adjust(x, "BH"))
save_csv(total_table, "Table_P5_total_LE8_by_BAG_modality.csv")
save_csv(component_table, "Table_P5_components_by_BAG_modality.csv")

mri_oos$bag_diffusion_minus_t1 <- mri_oos$bag_diffusion - mri_oos$bag_t1
for (i in seq_along(components)) {
  term <- component_terms[i]
  fs <- lm(as.formula(paste("bag_diffusion_minus_t1~", term, "+", bag_cov)),
           data = mri_oos)
  difference_res[[kd]] <- tidy_lm(fs, term, "Diffusion minus T1",
                                  "Separate component difference",
                                  labels[components[i]], FALSE)
  kd <- kd + 1L
}
fd <- lm(as.formula(paste("bag_diffusion_minus_t1~",
                          paste(component_terms, collapse = "+"), "+", bag_cov)),
         data = mri_oos)
for (i in seq_along(components)) {
  difference_res[[kd]] <- tidy_lm(fd, component_terms[i], "Diffusion minus T1",
                                  "Mutually adjusted difference",
                                  labels[components[i]], FALSE)
  kd <- kd + 1L
}
difference_table <- do.call(rbind, difference_res)
difference_table$p_FDR <- ave(difference_table$p, difference_table$model,
                              FUN = function(x) p.adjust(x, "BH"))
save_csv(difference_table, "Table_P5_T1_vs_Diffusion_component_differences.csv")

cor_table <- data.frame(
  bag_x = c("T1", "T1", "Diffusion"),
  bag_y = c("Diffusion", "Combined", "Combined"),
  n = c(sum(complete.cases(mri_oos[, c("bag_t1", "bag_diffusion")])),
        sum(complete.cases(mri_oos[, c("bag_t1", "bag_combined")])),
        sum(complete.cases(mri_oos[, c("bag_diffusion", "bag_combined")]))),
  Pearson_r = c(cor(mri_oos$bag_t1, mri_oos$bag_diffusion, use = "complete.obs"),
                cor(mri_oos$bag_t1, mri_oos$bag_combined, use = "complete.obs"),
                cor(mri_oos$bag_diffusion, mri_oos$bag_combined, use = "complete.obs")))
save_csv(cor_table, "Table_P5_BAG_correlations_formal_MRI.csv")
cat("Modality association stage completed.\n")
