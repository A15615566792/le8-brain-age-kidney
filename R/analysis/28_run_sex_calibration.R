args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) stop("Usage: Rscript 28_run_sex_calibration.R full.csv mri.csv predictions.csv out_dir seed")
full_path <- args[1]; mri_path <- args[2]; pred_path <- args[3]; out_dir <- args[4]; seed <- as.integer(args[5])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
suppressPackageStartupMessages(library(splines))

read_bom_csv <- function(p) {
  d <- read.csv(p, check.names = FALSE, na.strings = c("", "NA"))
  if (!"participant_id" %in% names(d) && grepl("participant_id$", names(d)[1])) names(d)[1] <- "participant_id"
  d
}
save_csv <- function(x, n) write.csv(x, file.path(out_dir, n), row.names = FALSE, na = "", fileEncoding = "UTF-8")
zs <- function(x) { x <- as.numeric(x); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(NA_real_, length(x)) else (x - mean(x, na.rm = TRUE)) / s }
eth4 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  w <- x %in% c("british", "irish", "white", "any other white background")
  a <- x %in% c("indian", "pakistani", "bangladeshi", "chinese", "asian or asian british", "any other asian background")
  b <- x %in% c("caribbean", "african", "black or black british", "any other black background")
  factor(ifelse(w, "White", ifelse(a, "Asian", ifelse(b, "Black", "Other/unknown"))), levels = c("White", "Asian", "Black", "Other/unknown"))
}
robust_vcov <- function(f, cluster = NULL, hc3 = TRUE) {
  X <- model.matrix(f); keep <- !is.na(coef(f)); X <- X[, keep, drop = FALSE]
  u <- residuals(f); w <- weights(f); if (is.null(w)) w <- rep(1, nrow(X))
  B <- solve(crossprod(X, X * w)); h <- rowSums((X %*% B) * X) * w
  if (is.null(cluster)) {
    adj <- if (hc3) pmax(1 - h, 1e-6) else rep(1, length(h))
    score <- X * as.numeric(w * u / adj)
    V <- B %*% crossprod(score) %*% B
  } else {
    score <- X * as.numeric(w * u)
    grouped <- rowsum(score, group = as.character(cluster), reorder = FALSE)
    G <- nrow(grouped); n <- nrow(X); p <- ncol(X)
    V <- B %*% crossprod(grouped) %*% B
    if (G > 1 && n > p) V <- V * (G / (G - 1)) * ((n - 1) / (n - p))
  }
  dimnames(V) <- list(colnames(X), colnames(X)); V
}
find_term <- function(f, pieces) {
  nms <- names(coef(f)); hit <- nms[vapply(nms, function(z) all(vapply(pieces, grepl, logical(1), x = z, fixed = TRUE)), logical(1))]
  hit <- hit[vapply(strsplit(hit, ":", fixed = TRUE), length, integer(1)) == length(pieces)]
  if (length(hit) != 1L) stop("Could not identify term: ", paste(pieces, collapse = " x "), "; hits=", paste(hit, collapse = ","))
  hit
}
tidy_sex <- function(f, exposure, modality, analysis, weighted, n) {
  V <- robust_vcov(f); int <- find_term(f, c(exposure, "sex_factorMale")); b <- coef(f)
  bf <- b[exposure]; bm <- b[exposure] + b[int]
  sef <- sqrt(V[exposure, exposure]); sem <- sqrt(V[exposure, exposure] + V[int, int] + 2 * V[exposure, int])
  sei <- sqrt(V[int, int]); pi <- 2 * pnorm(-abs(b[int] / sei))
  list(
    interaction = data.frame(analysis = analysis, modality = modality, exposure = exposure, beta_interaction_male_minus_female = unname(b[int]), HC3_SE = sei, CI_low = unname(b[int] - 1.96 * sei), CI_high = unname(b[int] + 1.96 * sei), p = pi, n = n, selection_IPW = weighted),
    slopes = data.frame(analysis = analysis, modality = modality, exposure = exposure, sex = c("Female", "Male"), beta = c(bf, bm), HC3_SE = c(sef, sem), CI_low = c(bf - 1.96 * sef, bm - 1.96 * sem), CI_high = c(bf + 1.96 * sef, bm + 1.96 * sem), p = c(2 * pnorm(-abs(bf / sef)), 2 * pnorm(-abs(bm / sem))), n = n, selection_IPW = weighted)
  )
}
make_long <- function(d) {
  a <- d; a$modality <- factor("T1", levels = c("T1", "Diffusion")); a$zbag <- a$zbag_t1
  b <- d; b$modality <- factor("Diffusion", levels = c("T1", "Diffusion")); b$zbag <- b$zbag_diffusion
  rbind(a, b)
}

cat("Reading cohorts and raw brain-age predictions...\n")
full <- read_bom_csv(full_path); mri <- read_bom_csv(mri_path); pred <- read_bom_csv(pred_path)
for (nm in c("full", "mri", "pred")) { d <- get(nm); d$participant_id <- as.character(d$participant_id); assign(nm, d) }
if (!"Sex" %in% names(pred)) pred <- merge(pred, full[, c("participant_id", "Sex")], by = "participant_id", all.x = TRUE, sort = FALSE)
modalities <- c("combined", "t1", "diffusion")
cal_rows <- list(); perf_rows <- list(); kc <- kp <- 1L
for (mod in modalities) {
  raw_col <- paste0("pred_raw_", mod); out_col <- paste0("bag_sexcal_", mod)
  pred[[out_col]] <- NA_real_
  for (sx in c("Female", "Male")) {
    tr <- pred$split == "Healthy train" & !is.na(pred$Sex) & !is.na(pred$Sex) & pred$Sex == sx & complete.cases(pred[, c("age", raw_col)])
    if (sum(tr) < 300) stop("Too few calibration participants for ", mod, " ", sx)
    cf <- coef(lm(pred[[raw_col]][tr] ~ pred$age[tr])); intercept <- unname(cf[1]); slope <- unname(cf[2])
    use <- !is.na(pred$Sex) & pred$Sex == sx & complete.cases(pred[, c("age", raw_col)])
    corrected <- (pred[[raw_col]][use] - intercept) / slope
    pred[[out_col]][use] <- corrected - pred$age[use]
    cal_rows[[kc]] <- data.frame(modality = mod, sex = sx, n_train = sum(tr), intercept = intercept, slope = slope); kc <- kc + 1L
    va <- pred$split == "Healthy validation" & !is.na(pred$Sex) & pred$Sex == sx & complete.cases(pred[, c("age", raw_col)])
    raw_bag <- pred[[raw_col]][va] - pred$age[va]
    corr_bag <- (pred[[raw_col]][va] - intercept) / slope - pred$age[va]
    for (correction in c("Raw", "Sex-specific calibrated")) {
      bag <- if (correction == "Raw") raw_bag else corr_bag
      perf_rows[[kp]] <- data.frame(modality = mod, sex = sx, correction = correction, n_validation = sum(va), mean_BAG = mean(bag), SD_BAG = sd(bag), MAE = mean(abs(bag)), BAG_age_correlation = cor(bag, pred$age[va])); kp <- kp + 1L
    }
  }
}
save_csv(do.call(rbind, cal_rows), "Table_P12_sex_specific_calibration_coefficients.csv")
save_csv(do.call(rbind, perf_rows), "Table_P12_sex_calibration_validation_performance.csv")

keep_pred <- c("participant_id", paste0("bag_sexcal_", modalities))
mri <- merge(mri, pred[, keep_pred], by = "participant_id", all = FALSE, sort = FALSE)
mri <- mri[mri$split != "Healthy train" & !is.na(mri$split), , drop = FALSE]
outside_train_n <- nrow(mri)
mri$sex_factor <- relevel(factor(mri$Sex), ref = "Female")
mri$ethnicity4 <- eth4(mri$Ethnic_0); mri$tdi <- as.numeric(mri$TDI_at_recruitment)
mri$edu_z <- as.numeric(mri$area_education_deprivation_z); mri$edu_missing <- factor(ifelse(is.na(mri$edu_z), "Missing", "Observed"), levels = c("Observed", "Missing")); mri$edu_imp <- ifelse(is.na(mri$edu_z), 0, mri$edu_z)
mri$center <- factor(mri$v4); mri$energy_z <- zs(mri$energy_kcal_pre_mri_mean); mri$gap_years <- as.numeric(mri$webq_mri_gap_days) / 365.25
mri$v_scaling_z <- zs(mri$v_scaling_t1_to_standard_2); mri$scanner_x_z <- zs(mri$scanner_lateral_x_brain_2); mri$scanner_y_z <- zs(mri$scanner_transverse_y_brain_2); mri$scanner_z_z <- zs(mri$scanner_longitudinal_z_brain_2); mri$scanner_table_z <- zs(mri$scanner_table_position_2); mri$t2flair_qc <- factor(mri$t2flair_run_freesurfer_2)
mri$LE8_per10 <- as.numeric(mri$LE8) / 10; mri$sbp_per10 <- as.numeric(mri$sbp) / 10; mri$hba1c_per5 <- as.numeric(mri$HbA1c_0) / 5
mri$bp_med_factor <- factor(mri$bp_med); mri$diabetes_factor <- factor(mri$diabetes_sr)
components <- c("le8_diet_primary", "le8_pa", "le8_nicotine", "le8_sleep", "le8_bmi", "le8_nonhdl", "le8_glucose", "le8_bp")
for (v in components) mri[[paste0(v, "_per10")]] <- as.numeric(mri[[v]]) / 10
component_terms <- paste0(components, "_per10")

full$baseline_date_r <- as.Date(substr(full$baseline_date, 1, 10)); full$time0_full_r <- as.Date(substr(full$time0_full, 1, 10)); full$age_webq <- as.numeric(full$Age_at_recruitment) + as.numeric(full$time0_full_r - full$baseline_date_r) / 365.25
full$sex_factor <- relevel(factor(full$Sex), ref = "Female"); full$ethnicity4 <- eth4(full$Ethnic_0); full$tdi <- as.numeric(full$TDI_at_recruitment); full$edu_z <- as.numeric(full$area_education_deprivation_z); full$edu_missing <- factor(ifelse(is.na(full$edu_z), "Missing", "Observed"), levels = c("Observed", "Missing")); full$edu_imp <- ifelse(is.na(full$edu_z), 0, full$edu_z); full$baseline_center_factor <- factor(full$baseline_center); full$LE8_per10 <- as.numeric(full$LE8) / 10; full$egfr <- as.numeric(full$egfr_cys_do)
full$in_mri <- as.integer(full$participant_id %in% mri$participant_id)
sel_fit <- glm(in_mri ~ ns(age_webq, df = 4) + sex_factor + ethnicity4 + tdi + edu_imp + edu_missing + baseline_center_factor + LE8_per10 + egfr, data = full, family = binomial())
full$psel <- predict(sel_fit, newdata = full, type = "response"); pm <- mean(full$in_mri)
mri <- merge(mri, full[, c("participant_id", "psel")], by = "participant_id", all.x = TRUE, sort = FALSE)
mri$sw <- pm / mri$psel; qlim <- quantile(mri$sw, c(0.01, 0.99), na.rm = TRUE); mri$sw_trim <- pmin(pmax(mri$sw, qlim[1]), qlim[2])
mri$zbag_t1 <- zs(mri$bag_sexcal_t1); mri$zbag_diffusion <- zs(mri$bag_sexcal_diffusion)

covars_no_sex <- c("ns(age, df = 4)", "ethnicity4", "tdi", "edu_imp", "edu_missing", "center", "energy_z", "gap_years", "v_scaling_z", "scanner_x_z", "scanner_y_z", "scanner_z_z", "scanner_table_z", "t2flair_qc")
specs <- list(
  list(name = "Total LE8 per 10 points", term = "LE8_per10", rhs = c("LE8_per10", covars_no_sex)),
  list(name = "LE8 glucose component per 10 points", term = "le8_glucose_per10", rhs = c(component_terms, covars_no_sex)),
  list(name = "LE8 BP component per 10 points", term = "le8_bp_per10", rhs = c(component_terms, covars_no_sex)),
  list(name = "Systolic BP per 10 mmHg", term = "sbp_per10", rhs = c("sbp_per10", setdiff(component_terms, "le8_bp_per10"), "bp_med_factor", covars_no_sex)),
  list(name = "HbA1c per 5 mmol/mol", term = "hba1c_per5", rhs = c("hba1c_per5", setdiff(component_terms, "le8_glucose_per10"), "diabetes_factor", covars_no_sex))
)

cat("Running sex-interaction sensitivity models after sex-specific calibration...\n")
int_rows <- list(); slope_rows <- list(); three_rows <- list(); ki <- ks <- k3 <- 1L
for (sp in specs) {
  for (mod in modalities) {
    outcome <- paste0("bag_sexcal_", mod)
    vars <- unique(c(outcome, sp$term, "sex_factor", "sw_trim", gsub("ns\\(age, df = 4\\)", "age", sp$rhs)))
    d <- mri[complete.cases(mri[, vars]), , drop = FALSE]
    form <- as.formula(paste(outcome, "~", sp$term, "* sex_factor +", paste(setdiff(sp$rhs, sp$term), collapse = "+")))
    for (weighted in c(FALSE, TRUE)) {
      f <- if (weighted) lm(form, data = d, weights = d$sw_trim) else lm(form, data = d)
      z <- tidy_sex(f, sp$term, mod, sp$name, weighted, nrow(d)); int_rows[[ki]] <- z$interaction; slope_rows[[ks]] <- z$slopes; ki <- ki + 1L; ks <- ks + 1L
    }
  }
  vars <- unique(c("participant_id", "zbag_t1", "zbag_diffusion", sp$term, "sex_factor", "sw_trim", gsub("ns\\(age, df = 4\\)", "age", sp$rhs)))
  d <- mri[complete.cases(mri[, vars]), , drop = FALSE]; long <- make_long(d)
  other_rhs <- setdiff(sp$rhs, sp$term)
  form3 <- as.formula(paste("zbag ~ modality * (sex_factor *", sp$term, "+", paste(other_rhs, collapse = "+"), ")"))
  for (weighted in c(FALSE, TRUE)) {
    f3 <- if (weighted) lm(form3, data = long, weights = long$sw_trim) else lm(form3, data = long)
    V <- robust_vcov(f3, cluster = long$participant_id, hc3 = FALSE)
    term3 <- find_term(f3, c("modalityDiffusion", "sex_factorMale", sp$term)); b <- coef(f3)[term3]; se <- sqrt(V[term3, term3])
    three_rows[[k3]] <- data.frame(analysis = sp$name, exposure = sp$term, term = term3, beta_threeway_Diffusion_vs_T1_male_difference = unname(b), cluster_robust_SE = se, CI_low = unname(b - 1.96 * se), CI_high = unname(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), n = nrow(d), selection_IPW = weighted); k3 <- k3 + 1L
  }
}
interactions <- do.call(rbind, int_rows); interactions$p_FDR_5 <- ave(interactions$p, interaction(interactions$modality, interactions$selection_IPW, drop = TRUE), FUN = function(x) p.adjust(x, "BH"))
three <- do.call(rbind, three_rows); three$p_FDR_5 <- ave(three$p, three$selection_IPW, FUN = function(x) p.adjust(x, "BH"))
save_csv(interactions, "Table_P12_sex_interactions_by_modality.csv")
save_csv(do.call(rbind, slope_rows), "Table_P12_sex_specific_slopes.csv")
save_csv(three, "Table_P12_modality_sex_threeway_interactions.csv")
save_csv(data.frame(item = c("Formal MRI outside Healthy train", "Female", "Male"), n = c(nrow(mri), sum(mri$Sex == "Female"), sum(mri$Sex == "Male"))), "P12_analysis_sample_audit.csv")
notes <- c(
  "Sex-specific age-bias calibration was estimated separately in female and male Healthy-train participants from raw predicted age and then applied outside Healthy train.",
  "Healthy validation performance is reported separately by sex to audit calibration transport.",
  "Sex interactions use sex-specifically calibrated BAG in years with HC3 covariance and flexible MRI age adjustment.",
  "The T1-versus-Diffusion three-way interaction uses separately standardized sex-calibrated BAG and participant-cluster robust covariance.",
  "Selection-IPW estimates are sensitivity analyses; FDR is applied across five prespecified exposures."
)
writeLines(notes, file.path(out_dir, "P12_analysis_notes.txt"), useBytes = TRUE)
cat("P12 sex calibration completed.\n")
