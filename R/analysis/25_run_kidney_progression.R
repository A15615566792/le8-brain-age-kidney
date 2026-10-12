args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript 25_run_kidney_progression.R full.csv egfr_subset.csv out_dir")
}

full_path <- args[1]
egfr_path <- args[2]
out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
suppressPackageStartupMessages(library(survival))

save_csv <- function(x, name) {
  write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

eth4 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  white <- x %in% c("british", "irish", "white", "any other white background")
  asian <- x %in% c("indian", "pakistani", "bangladeshi", "chinese",
                    "asian or asian british", "any other asian background")
  black <- x %in% c("caribbean", "african", "black or black british",
                    "any other black background")
  factor(ifelse(white, "White", ifelse(asian, "Asian",
         ifelse(black, "Black", "Other/unknown"))),
         levels = c("White", "Asian", "Black", "Other/unknown"))
}

robust_lm_term <- function(fit, term, analysis, outcome, weighted = FALSE) {
  X <- model.matrix(fit)
  e <- residuals(fit)
  w <- weights(fit)
  if (is.null(w)) w <- rep(1, length(e))
  bread <- solve(crossprod(X, X * w))
  score <- X * as.numeric(w * e)
  meat <- crossprod(score)
  n <- nrow(X); p <- ncol(X)
  V <- bread %*% meat %*% bread * n / (n - p)
  b <- coef(fit)[term]
  se <- sqrt(V[term, term])
  data.frame(
    analysis = analysis, outcome = outcome, exposure = term,
    estimate = unname(b), robust_SE = unname(se),
    CI_low = unname(b - 1.96 * se), CI_high = unname(b + 1.96 * se),
    p = 2 * pnorm(-abs(b / se)), n = nobs(fit), weighted = weighted,
    stringsAsFactors = FALSE
  )
}

cat("Reading frozen cohort and audited repeat-eGFR subset...\n")
full <- read.csv(full_path, check.names = FALSE, na.strings = c("", "NA"))
egfr <- read.csv(egfr_path, check.names = FALSE, na.strings = c("", "NA"))
full$participant_id <- as.character(full$participant_id)
egfr$participant_id <- as.character(egfr$participant_id)
stopifnot(!anyDuplicated(full$participant_id), !anyDuplicated(egfr$participant_id))

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
full$edu_z_imp <- ifelse(is.na(full$edu_z), 0, full$edu_z)
full$baseline_center_factor <- factor(full$baseline_center)
full$energy_z <- zscore(full$energy_kcal_primary)
full$webq_rounds <- as.numeric(full$webq_n_dates)
full$LE8_per10 <- as.numeric(full$LE8) / 10
full$egfr0 <- as.numeric(full$egfr_cys_do)

covars <- c("age_webq", "sex_factor", "ethnicity4", "tdi", "edu_z_imp",
            "edu_missing", "baseline_center_factor", "energy_z", "webq_rounds")

cat("Running repeat-eGFR progression models...\n")
keep_full <- c("participant_id", "age_webq", "sex_factor", "ethnicity4", "tdi",
               "edu_z_imp", "edu_missing", "baseline_center_factor", "energy_z",
               "webq_rounds", "LE8_per10", "egfr0")
prog <- merge(egfr, full[, keep_full], by = "participant_id", all.x = TRUE, sort = FALSE)

progression_vars <- c("egfr1_rebuilt", "egfr0_rebuilt", "egfr_annual_change",
                      "baseline_to_repeat_years", "time0_to_repeat_years", "LE8_per10",
                      covars)
prog_cc <- prog[complete.cases(prog[, progression_vars]), , drop = FALSE]
if (nrow(prog_cc) < 1000L) stop("Unexpectedly small progression complete-case sample")

rhs_prog <- paste(c("LE8_per10", "egfr0_rebuilt", "baseline_to_repeat_years",
                    "time0_to_repeat_years", covars), collapse = " + ")
fit_ancova <- lm(as.formula(paste("egfr1_rebuilt ~", rhs_prog)), data = prog_cc)
fit_slope <- lm(as.formula(paste("egfr_annual_change ~", rhs_prog)), data = prog_cc)

full$repeat_egfr_eligible <- as.integer(full$participant_id %in% egfr$participant_id)
sel_vars <- c("repeat_egfr_eligible", "LE8_per10", "egfr0", covars)
sel <- full[complete.cases(full[, sel_vars]), , drop = FALSE]
sel_fit <- glm(as.formula(paste("repeat_egfr_eligible ~",
                                paste(c("LE8_per10", "egfr0", covars), collapse = " + "))),
               family = binomial(), data = sel)
sel$p_repeat <- pmin(pmax(predict(sel_fit, type = "response"), 0.001), 0.999)
prev <- mean(sel$repeat_egfr_eligible)
sel$sw <- ifelse(sel$repeat_egfr_eligible == 1, prev / sel$p_repeat, NA_real_)
eligible_w <- sel$sw[sel$repeat_egfr_eligible == 1]
lims <- quantile(eligible_w, c(0.01, 0.99), na.rm = TRUE)
sel$sw_trim <- pmin(pmax(sel$sw, lims[1]), lims[2])
prog_cc <- merge(prog_cc, sel[, c("participant_id", "sw_trim")],
                 by = "participant_id", all.x = TRUE, sort = FALSE)
fit_ancova_ipw <- lm(as.formula(paste("egfr1_rebuilt ~", rhs_prog)),
                     data = prog_cc, weights = sw_trim)
fit_slope_ipw <- lm(as.formula(paste("egfr_annual_change ~", rhs_prog)),
                    data = prog_cc, weights = sw_trim)

progression_results <- rbind(
  robust_lm_term(fit_ancova, "LE8_per10", "Primary unweighted ANCOVA",
                 "Follow-up eGFR", FALSE),
  robust_lm_term(fit_ancova_ipw, "LE8_per10", "IPW ANCOVA sensitivity",
                 "Follow-up eGFR", TRUE),
  robust_lm_term(fit_slope, "LE8_per10", "Two-point slope sensitivity",
                 "Annual eGFR change", FALSE),
  robust_lm_term(fit_slope_ipw, "LE8_per10", "IPW two-point slope sensitivity",
                 "Annual eGFR change", TRUE)
)
progression_results$q_FDR <- p.adjust(progression_results$p, method = "BH")
save_csv(progression_results, "Table_P8_eGFR_progression_models.csv")

weight_diag <- data.frame(
  metric = c("full_n", "repeat_eligible_n", "progression_complete_case_n",
             "selection_model_n", "selection_model_events", "IPW_mean", "IPW_SD",
             "IPW_min_after_trim", "IPW_max_after_trim", "IPW_effective_sample_size"),
  value = c(nrow(full), nrow(egfr), nrow(prog_cc), nrow(sel), sum(sel$repeat_egfr_eligible),
            mean(prog_cc$sw_trim), sd(prog_cc$sw_trim), min(prog_cc$sw_trim),
            max(prog_cc$sw_trim), sum(prog_cc$sw_trim)^2 / sum(prog_cc$sw_trim^2))
)
save_csv(weight_diag, "Table_P8_eGFR_progression_IPW_diagnostics.csv")

cat("Progression results:\n")
print(progression_results)
cat("Completed P8 output:", normalizePath(out_dir), "\n")
