args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: Rscript 26_run_standardized_modality_continuous.R full.csv mri.csv modality_predictions.csv out_dir")
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
zs <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm = TRUE)
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

cluster_vcov <- function(f, cluster) {
  X <- model.matrix(f)
  keep <- !is.na(coef(f))
  X <- X[, keep, drop = FALSE]
  u <- residuals(f)
  w <- weights(f)
  if (is.null(w)) w <- rep(1, nrow(X))
  B <- tryCatch(solve(crossprod(X, X * w)), error = function(e) NULL)
  if (is.null(B)) return(vcov(f)[keep, keep, drop = FALSE])
  score <- X * as.numeric(w * u)
  grouped <- rowsum(score, group = as.character(cluster), reorder = FALSE)
  G <- nrow(grouped); n <- nrow(X); p <- ncol(X)
  V <- B %*% crossprod(grouped) %*% B
  if (G > 1 && n > p) V <- V * (G / (G - 1)) * ((n - 1) / (n - p))
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

interaction_name <- function(f, exposure) {
  nms <- names(coef(f))
  hit <- nms[(grepl("modalityDiffusion", nms, fixed = TRUE) &
              grepl(exposure, nms, fixed = TRUE))]
  hit <- hit[vapply(strsplit(hit, ":", fixed = TRUE), length, integer(1)) == 2L]
  if (length(hit) != 1L) stop("Could not uniquely identify modality interaction for ", exposure,
                              ": ", paste(hit, collapse = ", "))
  hit
}

tidy_interaction <- function(f, cluster, exposure, model, component, weighted) {
  V <- cluster_vcov(f, cluster)
  term <- interaction_name(f, exposure)
  b <- coef(f)[term]; se <- sqrt(V[term, term])
  data.frame(model = model, component = component, exposure = exposure,
             term = term, beta_Diffusion_minus_T1_SD = unname(b),
             cluster_robust_SE = unname(se), CI_low = unname(b - 1.96 * se),
             CI_high = unname(b + 1.96 * se),
             p = unname(2 * pnorm(-abs(b / se))),
             n_participants = length(unique(cluster)), selection_IPW = weighted,
             stringsAsFactors = FALSE)
}

tidy_slopes <- function(f, cluster, exposure, model, component, weighted) {
  V <- cluster_vcov(f, cluster)
  int <- interaction_name(f, exposure)
  b <- coef(f)
  bt <- b[exposure]; set <- sqrt(V[exposure, exposure])
  bd <- b[exposure] + b[int]
  sed <- sqrt(V[exposure, exposure] + V[int, int] + 2 * V[exposure, int])
  data.frame(model = model, component = component, exposure = exposure,
             modality = c("T1", "Diffusion"), beta_BAG_SD = c(bt, bd),
             cluster_robust_SE = c(set, sed),
             CI_low = c(bt - 1.96 * set, bd - 1.96 * sed),
             CI_high = c(bt + 1.96 * set, bd + 1.96 * sed),
             p = c(2 * pnorm(-abs(bt / set)), 2 * pnorm(-abs(bd / sed))),
             n_participants = length(unique(cluster)), selection_IPW = weighted,
             stringsAsFactors = FALSE)
}

make_long <- function(d) {
  a <- d; a$modality <- factor("T1", levels = c("T1", "Diffusion")); a$zbag <- a$z_bag_t1
  b <- d; b$modality <- factor("Diffusion", levels = c("T1", "Diffusion")); b$zbag <- b$z_bag_diffusion
  rbind(a, b)
}

cat("Reading frozen cohorts and modality predictions...\n")
full <- read.csv(full_path, check.names = FALSE, na.strings = c("", "NA"))
mri <- read.csv(mri_path, check.names = FALSE, na.strings = c("", "NA"))
pred <- read.csv(pred_path, check.names = FALSE, na.strings = c("", "NA"))
for (nm in c("full", "mri", "pred")) {
  d <- get(nm)
  if (!"participant_id" %in% names(d) && grepl("participant_id$", names(d)[1])) names(d)[1] <- "participant_id"
  assign(nm, d)
}
for (dname in c("full", "mri", "pred")) {
  d <- get(dname); d$participant_id <- as.character(d$participant_id); assign(dname, d)
}
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
full$LE8_per10 <- as.numeric(full$LE8) / 10
full$egfr <- as.numeric(full$egfr_cys_do)

mri <- merge(mri, pred, by = "participant_id", all = FALSE, sort = FALSE,
             suffixes = c("", "_pred"))
formal_merged_n <- nrow(mri)
mri$sex_factor <- relevel(factor(mri$Sex), ref = "Female")
mri$ethnicity4 <- eth4(mri$Ethnic_0)
mri$tdi <- as.numeric(mri$TDI_at_recruitment)
mri$edu_z <- as.numeric(mri$area_education_deprivation_z)
mri$edu_missing <- factor(ifelse(is.na(mri$edu_z), "Missing", "Observed"),
                          levels = c("Observed", "Missing"))
mri$edu_imp <- ifelse(is.na(mri$edu_z), 0, mri$edu_z)
mri$center <- factor(mri$v4)
mri$energy_z <- zs(mri$energy_kcal_pre_mri_mean)
mri$gap_years <- as.numeric(mri$webq_mri_gap_days) / 365.25
mri$v_scaling_z <- zs(mri$v_scaling_t1_to_standard_2)
mri$scanner_x_z <- zs(mri$scanner_lateral_x_brain_2)
mri$scanner_y_z <- zs(mri$scanner_transverse_y_brain_2)
mri$scanner_z_z <- zs(mri$scanner_longitudinal_z_brain_2)
mri$scanner_table_z <- zs(mri$scanner_table_position_2)
mri$t2flair_qc <- factor(mri$t2flair_run_freesurfer_2)
mri$LE8_per10 <- as.numeric(mri$LE8) / 10
mri <- mri[mri$split != "Healthy train" & !is.na(mri$split), , drop = FALSE]
outside_train_n <- nrow(mri)
mri$z_bag_t1 <- zs(mri$bag_t1)
mri$z_bag_diffusion <- zs(mri$bag_diffusion)

full$in_mri <- as.integer(full$participant_id %in% mri$participant_id)
sel_fit <- glm(in_mri ~ ns(age_webq, df = 4) + sex_factor + ethnicity4 + tdi +
                 edu_imp + edu_missing + baseline_center_factor + LE8_per10 + egfr,
               data = full, family = binomial())
full$psel <- predict(sel_fit, newdata = full, type = "response")
pm <- mean(full$in_mri)
mri <- merge(mri, full[, c("participant_id", "psel")], by = "participant_id",
             all.x = TRUE, sort = FALSE)
mri$sw <- pm / mri$psel
qlim <- quantile(mri$sw, c(0.01, 0.99), na.rm = TRUE)
mri$sw_trim <- pmin(pmax(mri$sw, qlim[1]), qlim[2])

components <- c("le8_diet_primary", "le8_pa", "le8_nicotine", "le8_sleep",
                "le8_bmi", "le8_nonhdl", "le8_glucose", "le8_bp")
labels <- setNames(c("Diet", "Physical activity", "Nicotine exposure", "Sleep",
                     "BMI", "Non-HDL cholesterol", "Blood glucose", "Blood pressure"),
                   components)
for (v in components) mri[[paste0(v, "_per10")]] <- as.numeric(mri[[v]]) / 10
component_terms <- paste0(components, "_per10")
covars <- c("age", "sex_factor", "ethnicity4", "tdi", "edu_imp", "edu_missing",
            "center", "energy_z", "gap_years", "v_scaling_z", "scanner_x_z",
            "scanner_y_z", "scanner_z_z", "scanner_table_z", "t2flair_qc")

cat("Running standardized modality interaction models...\n")
standardized_results <- list(); standardized_slopes <- list(); kr <- ks <- 1L
for (weighted in c(FALSE, TRUE)) {
  for (i in seq_along(component_terms)) {
    exposure <- component_terms[i]
    vars <- unique(c("participant_id", "z_bag_t1", "z_bag_diffusion", "sw_trim",
                     exposure, covars))
    d <- mri[complete.cases(mri[, vars]), vars, drop = FALSE]
    long <- make_long(d)
    form <- as.formula(paste("zbag ~ modality * (", exposure, "+",
                             paste(covars, collapse = "+"), ")"))
    f <- if (weighted) lm(form, data = long, weights = long$sw_trim) else lm(form, data = long)
    standardized_results[[kr]] <- tidy_interaction(
      f, long$participant_id, exposure, "Separate component", labels[components[i]], weighted)
    standardized_slopes[[ks]] <- tidy_slopes(
      f, long$participant_id, exposure, "Separate component", labels[components[i]], weighted)
    kr <- kr + 1L; ks <- ks + 1L
  }
  vars <- unique(c("participant_id", "z_bag_t1", "z_bag_diffusion", "sw_trim",
                   component_terms, covars))
  d <- mri[complete.cases(mri[, vars]), vars, drop = FALSE]
  long <- make_long(d)
  form <- as.formula(paste("zbag ~ modality * (",
                           paste(c(component_terms, covars), collapse = "+"), ")"))
  f <- if (weighted) lm(form, data = long, weights = long$sw_trim) else lm(form, data = long)
  for (i in seq_along(component_terms)) {
    standardized_results[[kr]] <- tidy_interaction(
      f, long$participant_id, component_terms[i], "Mutually adjusted 8 components",
      labels[components[i]], weighted)
    standardized_slopes[[ks]] <- tidy_slopes(
      f, long$participant_id, component_terms[i], "Mutually adjusted 8 components",
      labels[components[i]], weighted)
    kr <- kr + 1L; ks <- ks + 1L
  }
}
standardized_table <- do.call(rbind, standardized_results)
standardized_table$p_FDR_8 <- ave(
  standardized_table$p,
  interaction(standardized_table$model, standardized_table$selection_IPW, drop = TRUE),
  FUN = function(x) p.adjust(x, "BH"))
standardized_table$p_FDR_prespecified_glucose_BP <- NA_real_
for (model in unique(standardized_table$model)) {
  for (weighted in c(FALSE, TRUE)) {
    ii <- standardized_table$model == model &
      standardized_table$selection_IPW == weighted &
      standardized_table$component %in% c("Blood glucose", "Blood pressure")
    standardized_table$p_FDR_prespecified_glucose_BP[ii] <- p.adjust(standardized_table$p[ii], "BH")
  }
}
save_csv(standardized_table, "Table_P9_standardized_T1_Diffusion_interactions.csv")
save_csv(do.call(rbind, standardized_slopes), "Table_P9_standardized_modality_specific_slopes.csv")

cat("Bootstrapping glucose and blood-pressure modality differences...\n")
vars <- unique(c("participant_id", "z_bag_t1", "z_bag_diffusion", component_terms, covars))
bd <- mri[complete.cases(mri[, vars]), vars, drop = FALSE]
bd$zdiff <- bd$z_bag_diffusion - bd$z_bag_t1
bf <- lm(as.formula(paste("zdiff ~", paste(c(component_terms, covars), collapse = "+"))), data = bd)
X <- model.matrix(bf); y <- bd$zdiff
targets <- c(le8_glucose_per10 = "Blood glucose", le8_bp_per10 = "Blood pressure")
B <- 500L; set.seed(20260722)
boot <- matrix(NA_real_, nrow = B, ncol = length(targets), dimnames = list(NULL, names(targets)))
for (b in seq_len(B)) {
  ii <- sample.int(nrow(X), nrow(X), replace = TRUE)
  fitb <- lm.fit(X[ii, , drop = FALSE], y[ii])
  boot[b, ] <- fitb$coefficients[names(targets)]
}
boot_rows <- lapply(names(targets), function(term) {
  x <- boot[, term]; x <- x[is.finite(x)]
  data.frame(component = targets[[term]], term = term,
             original_beta_Diffusion_minus_T1_SD = coef(bf)[term],
             bootstrap_repetitions_valid = length(x),
             bootstrap_SE = sd(x), CI_low_percentile = quantile(x, 0.025),
             CI_high_percentile = quantile(x, 0.975),
             bootstrap_two_sided_p = 2 * min(mean(x <= 0), mean(x >= 0)),
             proportion_negative = mean(x < 0), stringsAsFactors = FALSE)
})
save_csv(do.call(rbind, boot_rows), "Table_P9_glucose_BP_modality_bootstrap.csv")

cat("Running continuous blood-pressure and glycaemia validation...\n")
mri$sbp_per10 <- as.numeric(mri$sbp) / 10
mri$dbp_per10 <- as.numeric(mri$dbp) / 10
mri$hba1c_per5 <- as.numeric(mri$HbA1c_0) / 5
mri$glucose_per1 <- as.numeric(mri$Glucose_0)
mri$bp_med_factor <- factor(mri$bp_med)
mri$diabetes_factor <- factor(mri$diabetes_sr)

specs <- list(
  list(term = "sbp_per10", label = "Systolic BP per 10 mmHg", omit = "le8_bp_per10",
       adjust = "bp_med_factor", restrict = "bp_med", restrict_value = 0,
       restriction = "No BP medication"),
  list(term = "dbp_per10", label = "Diastolic BP per 10 mmHg", omit = "le8_bp_per10",
       adjust = "bp_med_factor", restrict = "bp_med", restrict_value = 0,
       restriction = "No BP medication"),
  list(term = "hba1c_per5", label = "HbA1c per 5 mmol/mol", omit = "le8_glucose_per10",
       adjust = "diabetes_factor", restrict = "diabetes_sr", restrict_value = 0,
       restriction = "No self-reported diabetes/insulin"),
  list(term = "glucose_per1", label = "Glucose per 1 mmol/L", omit = "le8_glucose_per10",
       adjust = "diabetes_factor", restrict = "diabetes_sr", restrict_value = 0,
       restriction = "No self-reported diabetes/insulin")
)
continuous_interactions <- list(); continuous_slopes <- list(); ci <- cs <- 1L
for (sp in specs) {
  other_components <- setdiff(component_terms, sp$omit)
  for (restriction in c("Full sample adjusted for treatment/disease", sp$restriction)) {
    dd <- mri
    adjust_extra <- sp$adjust
    if (restriction != "Full sample adjusted for treatment/disease") {
      dd <- dd[as.numeric(dd[[sp$restrict]]) == sp$restrict_value, , drop = FALSE]
      adjust_extra <- character(0)
    }
    vars <- unique(c("participant_id", "z_bag_t1", "z_bag_diffusion", sp$term,
                     other_components, covars, adjust_extra))
    dd <- dd[complete.cases(dd[, vars]), vars, drop = FALSE]
    long <- make_long(dd)
    rhs <- c(sp$term, other_components, covars, adjust_extra)
    f <- lm(as.formula(paste("zbag ~ modality * (", paste(rhs, collapse = "+"), ")")),
            data = long)
    continuous_interactions[[ci]] <- tidy_interaction(
      f, long$participant_id, sp$term, restriction, sp$label, FALSE)
    continuous_slopes[[cs]] <- tidy_slopes(
      f, long$participant_id, sp$term, restriction, sp$label, FALSE)
    ci <- ci + 1L; cs <- cs + 1L
  }
}
continuous_table <- do.call(rbind, continuous_interactions)
continuous_table$p_FDR_within_sample <- ave(continuous_table$p, continuous_table$model,
                                             FUN = function(x) p.adjust(x, "BH"))
save_csv(continuous_table, "Table_P10_continuous_BP_glycaemia_modality_interactions.csv")
save_csv(do.call(rbind, continuous_slopes), "Table_P10_continuous_BP_glycaemia_modality_slopes.csv")

audit <- data.frame(
  item = c("Formal MRI merged", "Outside Healthy train", "All modality BAG complete",
           "Standardized mutual model complete", "Bootstrap repetitions"),
  n = c(formal_merged_n, outside_train_n, sum(complete.cases(mri[, c("z_bag_t1", "z_bag_diffusion")])),
        nrow(bd), B)
)
save_csv(audit, "P9_P10_analysis_sample_audit.csv")

notes <- c(
  "T1 and Diffusion BAG were standardized separately in the common outside-Healthy-train sample.",
  "Formal modality interactions use two records per participant and participant-cluster robust covariance.",
  "All modality-by-covariate interactions are included, allowing nuisance associations to differ by modality.",
  "Glucose and BP bootstrap uses 500 participant-level resamples of the mutually adjusted standardized difference model.",
  "Continuous BP and glycaemia models adjust the other seven LE8 components and treatment/disease status; restricted untreated/non-diabetic models are also reported.",
  "UK Biobank glucose is not guaranteed fasting and is labelled simply as glucose."
)
writeLines(notes, file.path(out_dir, "P9_P10_analysis_notes.txt"), useBytes = TRUE)
cat("P9 standardized modality and P10 continuous validation completed.\n")
