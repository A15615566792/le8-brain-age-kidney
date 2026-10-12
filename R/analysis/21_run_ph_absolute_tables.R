args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7L) {
  stop("Usage: Rscript 21_run_ph_absolute_tables.R full.csv mri.csv p1_dir p2_dir p3_dir flow.csv out_dir")
}

full_path <- args[1]
mri_path <- args[2]
p1_dir <- args[3]
p2_dir <- args[4]
p3_dir <- args[5]
flow_path <- args[6]
out_dir <- args[7]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
suppressPackageStartupMessages(library(survival))

save_csv <- function(x, name) {
  write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

read_result <- function(dir, name) {
  p <- file.path(dir, name)
  if (!file.exists(p)) stop("Missing prior result: ", p)
  read.csv(p, check.names = FALSE, na.strings = c("", "NA"))
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

le8_category <- function(x) {
  factor(ifelse(x < 50, "Low (<50)", ifelse(x < 80, "Moderate (50-79)", "High (>=80)")),
         levels = c("Low (<50)", "Moderate (50-79)", "High (>=80)"))
}

fmt_num <- function(x, digits = 1L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits, big.mark = ","), "")
}

fmt_p <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3)))
}

cat("Reading frozen cohorts...\n")
full <- read.csv(full_path, check.names = FALSE, na.strings = c("", "NA"))
mri <- read.csv(mri_path, check.names = FALSE, na.strings = c("", "NA"))
stopifnot(!anyDuplicated(full$participant_id), !anyDuplicated(mri$participant_id))

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
full$LE8_category <- le8_category(as.numeric(full$LE8))
full$egfr <- as.numeric(full$egfr_cys_do)
full$loguacr <- log1p(as.numeric(full$uacr))

mri$LE8_category <- le8_category(as.numeric(mri$LE8))
mri$ethnicity4 <- eth4(mri$Ethnic_0)
mri$gap_years <- as.numeric(mri$webq_mri_gap_days) / 365.25

components <- c("le8_diet_primary", "le8_pa", "le8_nicotine", "le8_sleep",
                "le8_bmi", "le8_nonhdl", "le8_glucose", "le8_bp")
component_labels <- setNames(c("Diet", "Physical activity", "Nicotine exposure",
                               "Sleep", "BMI", "Non-HDL cholesterol",
                               "Blood glucose", "Blood pressure"), components)
for (v in components) full[[paste0(v, "_per10")]] <- as.numeric(full[[v]]) / 10

full_cov <- paste(c("age_webq", "sex_factor", "ethnicity4", "tdi", "edu_z_imp",
                    "edu_missing", "baseline_center_factor", "energy_z",
                    "webq_rounds"), collapse = "+")


make_split5 <- function(d) {
  early <- d[d$followup_years > 0 & is.finite(d$followup_years), , drop = FALSE]
  early$tstart <- 0
  early$tstop <- pmin(early$followup_years, 5)
  early$event_split <- as.integer(early$kidney_event == 1 & early$followup_years <= 5)
  early$after5 <- 0
  late <- d[d$followup_years > 5 & is.finite(d$followup_years), , drop = FALSE]
  late$tstart <- 5
  late$tstop <- late$followup_years
  late$event_split <- as.integer(late$kidney_event == 1)
  late$after5 <- 1
  out <- rbind(early, late)
  out$participant_cluster <- as.character(out$participant_id)
  out
}

tidy_piecewise <- function(fit, main_term, interaction_term, model, component) {
  b <- coef(fit)
  V <- vcov(fit)
  b0 <- unname(b[main_term])
  se0 <- sqrt(V[main_term, main_term])
  b1 <- unname(b[main_term] + b[interaction_term])
  se1 <- sqrt(V[main_term, main_term] + V[interaction_term, interaction_term] +
              2 * V[main_term, interaction_term])
  pint <- 2 * pnorm(-abs(b[interaction_term] / sqrt(V[interaction_term, interaction_term])))
  omitted_rows <- if (is.null(fit$na.action)) integer(0) else as.integer(fit$na.action)
  used_rows <- setdiff(seq_len(nrow(split5)), omitted_rows)
  n_participants <- length(unique(split5$participant_cluster[used_rows]))
  data.frame(
    model = model,
    component = component,
    period = c("0-5 years", ">5 years"),
    logHR = c(b0, b1),
    HR = exp(c(b0, b1)),
    CI_low = exp(c(b0 - 1.96 * se0, b1 - 1.96 * se1)),
    CI_high = exp(c(b0 + 1.96 * se0, b1 + 1.96 * se1)),
    p_period_effect = c(2 * pnorm(-abs(b0 / se0)), 2 * pnorm(-abs(b1 / se1))),
    p_time_interaction = pint,
    n_participants = n_participants,
    events = sum(fit$y[, "status"]),
    stringsAsFactors = FALSE
  )
}

cat("Running piecewise proportional-hazards models...\n")
split5 <- make_split5(full)
split5$pa_after5 <- split5$le8_pa_per10 * split5$after5
split5$nonhdl_after5 <- split5$le8_nonhdl_per10 * split5$after5

ph_results <- list()
counter <- 1L
for (spec in list(
  list(main = "le8_pa_per10", int = "pa_after5", label = "Physical activity"),
  list(main = "le8_nonhdl_per10", int = "nonhdl_after5", label = "Non-HDL cholesterol")
)) {
  form_sep <- as.formula(paste0(
    "Surv(tstart,tstop,event_split)~", spec$main, "+", spec$int, "+", full_cov,
    "+cluster(participant_cluster)"))
  fit_sep <- coxph(form_sep, data = split5, model = TRUE, x = TRUE)
  ph_results[[counter]] <- tidy_piecewise(fit_sep, spec$main, spec$int,
                                           "Separate component", spec$label)
  counter <- counter + 1L
}

mutual_terms <- paste0(components, "_per10")
form_mut <- as.formula(paste0(
  "Surv(tstart,tstop,event_split)~", paste(mutual_terms, collapse = "+"),
  "+pa_after5+nonhdl_after5+", full_cov, "+cluster(participant_cluster)"))
fit_mut <- coxph(form_mut, data = split5, model = TRUE, x = TRUE)
ph_results[[counter]] <- tidy_piecewise(fit_mut, "le8_pa_per10", "pa_after5",
                                         "Mutually adjusted 8 components",
                                         "Physical activity")
counter <- counter + 1L
ph_results[[counter]] <- tidy_piecewise(fit_mut, "le8_nonhdl_per10", "nonhdl_after5",
                                         "Mutually adjusted 8 components",
                                         "Non-HDL cholesterol")
ph_table <- do.call(rbind, ph_results)
ph_table$p_FDR_time_interaction <- ave(ph_table$p_time_interaction,
                                        ph_table$model,
                                        FUN = function(x) p.adjust(x, "BH"))
save_csv(ph_table, "Table_P4_PH_piecewise_0_5_vs_after5.csv")


horizons <- c(5, 10)
cat("Fitting category-based Fine-Gray model for standardized absolute risks...\n")
full$status_fg <- factor(ifelse(full$kidney_event == 1, "CKD",
                         ifelse(full$death_competing == 1, "Death", "Censor")),
                         levels = c("Censor", "CKD", "Death"))
fg_vars <- c("followup_years", "status_fg", "LE8_category", "age_webq",
             "sex_factor", "ethnicity4", "tdi", "edu_z_imp", "edu_missing",
             "baseline_center_factor", "energy_z", "webq_rounds")
fg_source <- full[complete.cases(full[, fg_vars]), fg_vars, drop = FALSE]
fg_source$id_num <- seq_len(nrow(fg_source))
fg_expanded <- finegray(
  Surv(followup_years, status_fg) ~ LE8_category + age_webq + sex_factor +
    ethnicity4 + tdi + edu_z_imp + edu_missing + baseline_center_factor +
    energy_z + webq_rounds + id_num,
  data = fg_source, etype = "CKD", id = id_num)
fg_fit <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~ LE8_category + age_webq + sex_factor +
    ethnicity4 + tdi + edu_z_imp + edu_missing + baseline_center_factor +
    energy_z + webq_rounds + cluster(id_num),
  data = fg_expanded, weights = fgwt, x = TRUE, model = TRUE)

bh <- basehaz(fg_fit, centered = FALSE)
get_h0 <- function(t) {
  ii <- which(bh$time <= t)
  if (!length(ii)) 0 else bh$hazard[max(ii)]
}

pred_formula <- ~ LE8_category + age_webq + sex_factor + ethnicity4 + tdi +
  edu_z_imp + edu_missing + baseline_center_factor + energy_z + webq_rounds

beta <- coef(fg_fit)
V <- vcov(fg_fit)
coef_names <- names(beta)

make_x <- function(nd) {
  X <- model.matrix(pred_formula, nd)
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  missing_cols <- setdiff(coef_names, colnames(X))
  if (length(missing_cols)) {
    add <- matrix(0, nrow(X), length(missing_cols),
                  dimnames = list(NULL, missing_cols))
    X <- cbind(X, add)
  }
  X[, coef_names, drop = FALSE]
}

set.seed(20260721)
sim_n <- min(10000L, nrow(fg_source))
sim_idx <- sort(sample.int(nrow(fg_source), sim_n))
eig <- eigen((V + t(V)) / 2, symmetric = TRUE)
eig$values[eig$values < 0] <- 0
rootV <- eig$vectors %*% diag(sqrt(eig$values), nrow = length(eig$values))
n_sim <- 500L
beta_draws <- matrix(rnorm(length(beta) * n_sim), nrow = n_sim) %*% t(rootV)
beta_draws <- sweep(beta_draws, 2, beta, "+")

adjusted_list <- list()
sim_store <- list()
for (g in levels(fg_source$LE8_category)) {
  nd <- fg_source
  nd$LE8_category <- factor(g, levels = levels(fg_source$LE8_category))
  Xfull <- make_x(nd)
  lp_full <- drop(Xfull %*% beta)
  nd_sim <- nd[sim_idx, , drop = FALSE]
  Xsim <- make_x(nd_sim)
  sim_store[[g]] <- list()
  for (tt in horizons) {
    H0 <- get_h0(tt)
    risk <- mean(1 - exp(-H0 * exp(lp_full)))
    lp_draws <- Xsim %*% t(beta_draws)
    risk_draws <- colMeans(1 - exp(-H0 * exp(lp_draws)))
    sim_store[[g]][[as.character(tt)]] <- risk_draws
    adjusted_list[[length(adjusted_list) + 1L]] <- data.frame(
      LE8_category = g, horizon_years = tt,
      adjusted_cumulative_incidence = risk,
      CI_low_approx = as.numeric(quantile(risk_draws, 0.025)),
      CI_high_approx = as.numeric(quantile(risk_draws, 0.975)),
      adjusted_risk_per_1000 = 1000 * risk,
      source_n = nrow(fg_source), standardization_n = nrow(fg_source),
      CI_method = "500 coefficient simulations; baseline hazard held fixed",
      stringsAsFactors = FALSE)
  }
  rm(nd, Xfull, lp_full)
  gc()
}
adjusted_cif <- do.call(rbind, adjusted_list)
save_csv(adjusted_cif, "Table_P4_CKD_adjusted_FG_CIF.csv")

contrast_list <- list()
for (tt in horizons) {
  for (cmp in list(c("Moderate (50-79)", "Low (<50)"),
                   c("High (>=80)", "Low (<50)"),
                   c("High (>=80)", "Moderate (50-79)"))) {
    r1 <- adjusted_cif$adjusted_cumulative_incidence[
      adjusted_cif$LE8_category == cmp[1] & adjusted_cif$horizon_years == tt]
    r0 <- adjusted_cif$adjusted_cumulative_incidence[
      adjusted_cif$LE8_category == cmp[2] & adjusted_cif$horizon_years == tt]
    draws <- sim_store[[cmp[1]]][[as.character(tt)]] -
      sim_store[[cmp[2]]][[as.character(tt)]]
    contrast_list[[length(contrast_list) + 1L]] <- data.frame(
      comparison = paste(cmp[1], "vs", cmp[2]), horizon_years = tt,
      adjusted_risk_difference = r1 - r0,
      CI_low_approx = as.numeric(quantile(draws, 0.025)),
      CI_high_approx = as.numeric(quantile(draws, 0.975)),
      events_difference_per_1000 = 1000 * (r1 - r0),
      CI_low_per_1000 = 1000 * as.numeric(quantile(draws, 0.025)),
      CI_high_per_1000 = 1000 * as.numeric(quantile(draws, 0.975)),
      stringsAsFactors = FALSE)
  }
}
risk_contrasts <- do.call(rbind, contrast_list)
save_csv(risk_contrasts, "Table_P4_CKD_adjusted_absolute_risk_differences.csv")


continuous_rows <- function(d, vars, labels, skewed, group_var, cohort) {
  groups <- c("Overall", levels(d[[group_var]]))
  out <- list()
  k <- 1L
  for (i in seq_along(vars)) {
    v <- vars[i]
    for (g in groups) {
      x <- if (g == "Overall") as.numeric(d[[v]]) else
        as.numeric(d[d[[group_var]] == g, v])
      nmiss <- sum(!is.finite(x))
      x <- x[is.finite(x)]
      if (v %in% skewed) {
        stat <- if (length(x)) paste0(fmt_num(median(x), 1), " [",
                                      fmt_num(quantile(x, 0.25), 1), ", ",
                                      fmt_num(quantile(x, 0.75), 1), "]") else ""
        summary_type <- "Median [IQR]"
      } else {
        stat <- if (length(x)) paste0(fmt_num(mean(x), 1), " (", fmt_num(sd(x), 1), ")") else ""
        summary_type <- "Mean (SD)"
      }
      out[[k]] <- data.frame(cohort = cohort, variable = labels[i], level = "",
                             group = g, statistic = stat,
                             summary_type = summary_type, missing_n = nmiss)
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

categorical_rows <- function(d, vars, labels, group_var, cohort) {
  groups <- c("Overall", levels(d[[group_var]]))
  out <- list()
  k <- 1L
  for (i in seq_along(vars)) {
    v <- vars[i]
    lv <- levels(factor(d[[v]]))
    for (lev in lv) {
      for (g in groups) {
        x <- if (g == "Overall") d[[v]] else d[d[[group_var]] == g, v]
        denom <- sum(!is.na(x))
        n <- sum(as.character(x) == lev, na.rm = TRUE)
        stat <- if (denom) paste0(format(n, big.mark = ","), " (", fmt_num(100 * n / denom, 1), ")") else ""
        out[[k]] <- data.frame(cohort = cohort, variable = labels[i], level = lev,
                               group = g, statistic = stat,
                               summary_type = "n (%)", missing_n = sum(is.na(x)))
        k <- k + 1L
      }
    }
  }
  do.call(rbind, out)
}

wide_table <- function(long) {
  id <- unique(long[, c("cohort", "variable", "level", "summary_type")])
  groups <- c("Overall", "Low (<50)", "Moderate (50-79)", "High (>=80)")
  for (g in groups) {
    z <- long[long$group == g, c("cohort", "variable", "level", "summary_type", "statistic")]
    names(z)[5] <- g
    id <- merge(id, z, by = c("cohort", "variable", "level", "summary_type"), all.x = TRUE, sort = FALSE)
  }
  id
}

cat("Building Table 1 for full and MRI cohorts...\n")
full$sex_table <- factor(full$Sex)
full$smoking_table <- factor(ifelse(is.na(full$smoking_status_0), NA,
                                    as.character(full$smoking_status_0)))
full$diabetes_table <- factor(ifelse(as.numeric(full$diabetes_sr) == 1, "Yes", "No"),
                              levels = c("No", "Yes"))
full$bp_med_table <- factor(ifelse(as.numeric(full$bp_med) == 1, "Yes", "No"),
                            levels = c("No", "Yes"))
full$lipid_med_table <- factor(ifelse(as.numeric(full$lipid_lowering_med) == 1, "Yes", "No"),
                               levels = c("No", "Yes"))

full_cont_vars <- c("age_webq", "LE8", "LE8_behavior", "LE8_biological",
                    "BMI_0", "sbp", "dbp", "Glucose_0", "HbA1c_0",
                    "non_hdl_mgdl", "sleep_duration_0", "pa_minutes_week",
                    "egfr_cys_do", "uacr", "TDI_at_recruitment",
                    "area_education_deprivation_z", "energy_kcal_primary",
                    "webq_n_dates", "followup_years")
full_cont_labels <- c("Age at first valid WebQ, years", "LE8 total score",
                      "LE8 behavior score", "LE8 biological score",
                      "BMI, kg/m2", "Systolic blood pressure, mmHg",
                      "Diastolic blood pressure, mmHg", "Glucose, mmol/L",
                      "HbA1c, mmol/mol", "Non-HDL cholesterol, mg/dL",
                      "Sleep duration, hours/day", "Physical activity, min/week",
                      "eGFR, mL/min/1.73m2", "UACR", "Townsend deprivation index",
                      "Area education deprivation z-score", "Energy intake, kcal/day",
                      "Valid WebQ rounds", "Follow-up, years")
full_cat_vars <- c("sex_table", "ethnicity4", "smoking_table", "diabetes_table",
                   "bp_med_table", "lipid_med_table")
full_cat_labels <- c("Sex", "Ethnicity", "Smoking status", "Diabetes",
                     "Blood-pressure medication", "Lipid-lowering medication")
full_long <- rbind(
  continuous_rows(full, full_cont_vars, full_cont_labels,
                  skewed = c("pa_minutes_week", "uacr", "webq_n_dates"),
                  group_var = "LE8_category", cohort = "Full WebQ cohort"),
  categorical_rows(full, full_cat_vars, full_cat_labels,
                   group_var = "LE8_category", cohort = "Full WebQ cohort"))
save_csv(full_long, "Table_1A_full_cohort_long.csv")
save_csv(wide_table(full_long), "Table_1A_full_cohort_formatted.csv")

mri$sex_table <- factor(mri$Sex)
mri$smoking_table <- factor(ifelse(is.na(mri$smoking_status_0), NA,
                                   as.character(mri$smoking_status_0)))
mri$diabetes_table <- factor(ifelse(as.numeric(mri$diabetes_sr) == 1, "Yes", "No"),
                             levels = c("No", "Yes"))
mri$lipid_med_table <- factor(ifelse(as.numeric(mri$lipid_lowering_med) == 1, "Yes", "No"),
                              levels = c("No", "Yes"))
mri$mri_center_table <- factor(mri$v4)

mri_cont_vars <- c("age", "LE8", "LE8_behavior", "LE8_biological",
                   "brain_age_gap_oos_candidate", "BMI_0", "sbp", "dbp",
                   "Glucose_0", "HbA1c_0", "non_hdl_mgdl", "sleep_duration_0",
                   "pa_minutes_week", "egfr_cys_do", "uacr", "TDI_at_recruitment",
                   "area_education_deprivation_z", "energy_kcal_pre_mri_mean",
                   "gap_years")
mri_cont_labels <- c("Age at MRI, years", "LE8 total score", "LE8 behavior score",
                     "LE8 biological score", "Leakage-controlled BAG, years",
                     "BMI, kg/m2", "Systolic blood pressure, mmHg",
                     "Diastolic blood pressure, mmHg", "Glucose, mmol/L",
                     "HbA1c, mmol/mol", "Non-HDL cholesterol, mg/dL",
                     "Sleep duration, hours/day", "Physical activity, min/week",
                     "eGFR, mL/min/1.73m2", "UACR", "Townsend deprivation index",
                     "Area education deprivation z-score", "Pre-MRI energy intake, kcal/day",
                     "WebQ-to-MRI interval, years")
mri_cat_vars <- c("sex_table", "ethnicity4", "smoking_table", "diabetes_table",
                  "lipid_med_table", "mri_center_table")
mri_cat_labels <- c("Sex", "Ethnicity", "Smoking status", "Diabetes",
                    "Lipid-lowering medication", "MRI assessment centre")
mri_long <- rbind(
  continuous_rows(mri, mri_cont_vars, mri_cont_labels,
                  skewed = c("pa_minutes_week", "uacr", "gap_years"),
                  group_var = "LE8_category", cohort = "MRI cohort"),
  categorical_rows(mri, mri_cat_vars, mri_cat_labels,
                   group_var = "LE8_category", cohort = "MRI cohort"))
save_csv(mri_long, "Table_1B_MRI_cohort_long.csv")
save_csv(wide_table(mri_long), "Table_1B_MRI_cohort_formatted.csv")

category_counts <- do.call(rbind, list(
  transform(as.data.frame(table(full$LE8_category)), cohort = "Full WebQ cohort"),
  transform(as.data.frame(table(mri$LE8_category)), cohort = "MRI cohort")))
names(category_counts)[1:2] <- c("LE8_category", "n")
category_counts$events <- c(
  sapply(levels(full$LE8_category), function(g) sum(full$kidney_event[full$LE8_category == g], na.rm = TRUE)),
  sapply(levels(mri$LE8_category), function(g) sum(mri$kidney_event_mri[mri$LE8_category == g], na.rm = TRUE)))
save_csv(category_counts, "Table_P4_LE8_category_counts_events.csv")


cat("Consolidating prior model outputs into manuscript tables...\n")
p1_total <- read_result(p1_dir, "Table_P1_total_LE8_models.csv")
p1_components <- read_result(p1_dir, "Table_P1_component_association_profiles.csv")
p1_domains <- read_result(p1_dir, "Table_P1_behavior_biological_domains.csv")
p2_total <- read_result(p2_dir, "Table_P2_total_LE8_CKD_sensitivities.csv")
p2_egfr <- read_result(p2_dir, "Table_P2_components_CKD_plus_eGFR.csv")
p2_ipw_total <- read_result(p2_dir, "Table_P2_MRI_selection_IPW_total_LE8.csv")
p3_fg <- read_result(p3_dir, "Table_P3_FineGray_total_LE8.csv")

bind_fill <- function(xs) {
  alln <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(x) {
    miss <- setdiff(alln, names(x))
    for (v in miss) x[[v]] <- NA
    x[, alln, drop = FALSE]
  })
  do.call(rbind, xs)
}

table2 <- bind_fill(list(
  cbind(source = "P1_total", p1_total),
  cbind(source = "P2_CKD_sensitivity", p2_total),
  cbind(source = "P2_IPW_total", p2_ipw_total),
  cbind(source = "P3_FineGray", p3_fg)))
save_csv(table2, "Table_2_total_LE8_BAG_CKD_models.csv")

table3 <- bind_fill(list(
  cbind(source = "P1_primary_components", p1_components[p1_components$scale == "per_10_points", ]),
  cbind(source = "P2_CKD_plus_eGFR", p2_egfr)))
save_csv(table3, "Table_3_LE8_component_profiles.csv")

table4 <- bind_fill(list(
  cbind(source = "Behavior_biological_domains", p1_domains),
  cbind(source = "PH_piecewise", ph_table),
  cbind(source = "Adjusted_absolute_risk", adjusted_cif),
  cbind(source = "Adjusted_risk_difference", risk_contrasts)))
save_csv(table4, "Table_4_sensitivity_and_clinical_interpretation.csv")

flow <- read.csv(flow_path, check.names = FALSE, na.strings = c("", "NA"))
save_csv(flow, "Figure_1_cohort_flow_counts.csv")

notes <- c(
  "PA and non-HDL proportional-hazards departures were handled with prespecified 0-5-year and >5-year coefficients.",
  "Crude CKD cumulative incidence uses an Aalen-Johansen estimator with death as a competing event.",
  "Adjusted absolute risks use direct standardization from a Fine-Gray subdistribution model.",
  "LE8 clinical categories are Low <50, Moderate 50-79, and High >=80."
)
writeLines(notes, file.path(out_dir, "P4_analysis_notes.txt"), useBytes = TRUE)

manifest <- data.frame(file = list.files(out_dir),
                       bytes = file.info(file.path(out_dir, list.files(out_dir)))$size,
                       stringsAsFactors = FALSE)
save_csv(manifest, "P4_output_manifest.csv")
cat("P4 PH, absolute-risk, and formal-table stage completed.\n")
