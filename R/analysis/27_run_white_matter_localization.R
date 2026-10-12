args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) stop("Usage: Rscript 27_run_white_matter_localization.R full.csv mri.csv feature.csv map.csv out_dir")
full_path <- args[1]; mri_path <- args[2]; feature_path <- args[3]; map_path <- args[4]; out_dir <- args[5]
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

cat("Reading frozen cohorts and core diffusion features...\n")
full <- read_bom_csv(full_path); mri <- read_bom_csv(mri_path); feat <- read_bom_csv(feature_path)
map <- read.csv(map_path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
for (nm in c("full", "mri", "feat")) { d <- get(nm); d$participant_id <- as.character(d$participant_id); assign(nm, d) }
stopifnot(!anyDuplicated(full$participant_id), !anyDuplicated(mri$participant_id), !anyDuplicated(feat$participant_id))
feature_names <- as.character(map$feature)
stopifnot(length(feature_names) == 134L, all(feature_names %in% names(feat)))
mri <- merge(mri, feat[, c("participant_id", feature_names)], by = "participant_id", all = FALSE, sort = FALSE)
formal_feature_merged_n <- nrow(mri)

mri$sex_factor <- relevel(factor(mri$Sex), ref = "Female")
mri$ethnicity4 <- eth4(mri$Ethnic_0)
mri$tdi <- as.numeric(mri$TDI_at_recruitment)
mri$edu_z <- as.numeric(mri$area_education_deprivation_z)
mri$edu_missing <- factor(ifelse(is.na(mri$edu_z), "Missing", "Observed"), levels = c("Observed", "Missing"))
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
mri$sbp_per10 <- as.numeric(mri$sbp) / 10
mri$hba1c_per5 <- as.numeric(mri$HbA1c_0) / 5
mri$bp_med_factor <- factor(mri$bp_med)
mri$diabetes_factor <- factor(mri$diabetes_sr)

components <- c("le8_diet_primary", "le8_pa", "le8_nicotine", "le8_sleep", "le8_bmi", "le8_nonhdl", "le8_glucose", "le8_bp")
for (v in components) mri[[paste0(v, "_per10")]] <- as.numeric(mri[[v]]) / 10
component_terms <- paste0(components, "_per10")

full$baseline_date_r <- as.Date(substr(full$baseline_date, 1, 10)); full$time0_full_r <- as.Date(substr(full$time0_full, 1, 10))
full$age_webq <- as.numeric(full$Age_at_recruitment) + as.numeric(full$time0_full_r - full$baseline_date_r) / 365.25
full$sex_factor <- relevel(factor(full$Sex), ref = "Female"); full$ethnicity4 <- eth4(full$Ethnic_0)
full$tdi <- as.numeric(full$TDI_at_recruitment); full$edu_z <- as.numeric(full$area_education_deprivation_z)
full$edu_missing <- factor(ifelse(is.na(full$edu_z), "Missing", "Observed"), levels = c("Observed", "Missing"))
full$edu_imp <- ifelse(is.na(full$edu_z), 0, full$edu_z); full$baseline_center_factor <- factor(full$baseline_center)
full$LE8_per10 <- as.numeric(full$LE8) / 10; full$egfr <- as.numeric(full$egfr_cys_do)
full$in_mri <- as.integer(full$participant_id %in% mri$participant_id)
sel_fit <- glm(in_mri ~ ns(age_webq, df = 4) + sex_factor + ethnicity4 + tdi + edu_imp + edu_missing + baseline_center_factor + LE8_per10 + egfr, data = full, family = binomial())
full$psel <- predict(sel_fit, newdata = full, type = "response"); pm <- mean(full$in_mri)
mri <- merge(mri, full[, c("participant_id", "psel")], by = "participant_id", all.x = TRUE, sort = FALSE)
mri$sw <- pm / mri$psel; qlim <- quantile(mri$sw, c(0.01, 0.99), na.rm = TRUE)
mri$sw_trim <- pmin(pmax(mri$sw, qlim[1]), qlim[2])

covars <- c("ns(age, df = 4)", "sex_factor", "ethnicity4", "tdi", "edu_imp", "edu_missing", "center", "energy_z", "gap_years", "v_scaling_z", "scanner_x_z", "scanner_y_z", "scanner_z_z", "scanner_table_z", "t2flair_qc")
specs <- list(
  list(name = "SBP full sample", term = "sbp_per10", label = "Systolic BP per 10 mmHg", omit = "le8_bp_per10", adjust = "bp_med_factor", subset = rep(TRUE, nrow(mri))),
  list(name = "SBP no BP medication", term = "sbp_per10", label = "Systolic BP per 10 mmHg", omit = "le8_bp_per10", adjust = character(0), subset = as.numeric(mri$bp_med) == 0),
  list(name = "HbA1c full sample", term = "hba1c_per5", label = "HbA1c per 5 mmol/mol", omit = "le8_glucose_per10", adjust = "diabetes_factor", subset = rep(TRUE, nrow(mri))),
  list(name = "HbA1c no diabetes/insulin", term = "hba1c_per5", label = "HbA1c per 5 mmol/mol", omit = "le8_glucose_per10", adjust = character(0), subset = as.numeric(mri$diabetes_sr) == 0)
)

centers <- vapply(mri[, feature_names], mean, numeric(1), na.rm = TRUE)
sds <- vapply(mri[, feature_names], sd, numeric(1), na.rm = TRUE)
stopifnot(all(is.finite(sds) & sds > 0))
Zall <- sweep(sweep(as.matrix(mri[, feature_names]), 2, centers, "-"), 2, sds, "/")
colnames(Zall) <- feature_names
metric <- setNames(as.character(map$metric), as.character(map$feature))
description <- setNames(as.character(map$description), as.character(map$feature))
adverse_sign <- ifelse(metric %in% c("FA", "ICVF"), -1, ifelse(metric %in% c("MD", "ISOV"), 1, NA_real_)); names(adverse_sign) <- names(metric)

fit_many <- function(dd, Y, rhs, exposure, w) {
  X <- model.matrix(as.formula(paste("~", paste(rhs, collapse = "+"))), data = dd)
  j <- match(exposure, colnames(X)); if (is.na(j)) stop("Exposure absent from model matrix: ", exposure)
  sw <- as.numeric(w); Xw <- X * sqrt(sw); Yw <- Y * sqrt(sw)
  fit <- lm.fit(Xw, Yw); cf <- fit$coefficients
  bread <- solve(crossprod(X, X * sw))
  resid <- Y - X %*% cf
  h <- rowSums((X %*% bread) * X) * sw
  h <- pmin(h, 0.999999)
  influence <- as.numeric(X %*% bread[, j]) * sw
  se <- sqrt(colSums((resid * (influence / (1 - h)))^2))
  b <- as.numeric(cf[j, ]); p <- 2 * pnorm(-abs(b / se))
  data.frame(outcome = colnames(Y), beta = b, HC3_SE = se, CI_low = b - 1.96 * se, CI_high = b + 1.96 * se, p = p, stringsAsFactors = FALSE)
}

cat("Running diffusion metric composites and tract localization...\n")
tract_rows <- list(); composite_rows <- list(); audit_rows <- list(); kt <- kc <- ka <- 1L
for (sp in specs) {
  other_components <- setdiff(component_terms, sp$omit)
  rhs <- c(sp$term, other_components, covars, sp$adjust)
  needed <- unique(c(sp$term, other_components, "age", "sex_factor", "ethnicity4", "tdi", "edu_imp", "edu_missing", "center", "energy_z", "gap_years", "v_scaling_z", "scanner_x_z", "scanner_y_z", "scanner_z_z", "scanner_table_z", "t2flair_qc", sp$adjust, "sw_trim"))
  keep <- sp$subset & complete.cases(mri[, needed]) & complete.cases(Zall)
  dd <- mri[keep, , drop = FALSE]; Y <- Zall[keep, , drop = FALSE]
  comp <- sapply(c("FA", "MD", "ICVF", "ISOV"), function(mm) {
    jj <- which(metric[colnames(Y)] == mm); rowMeans(sweep(Y[, jj, drop = FALSE], 2, adverse_sign[colnames(Y)[jj]], "*"))
  })
  od_jj <- which(metric[colnames(Y)] == "OD")
  comp <- cbind(comp, `OD original direction` = rowMeans(Y[, od_jj, drop = FALSE]), `Global adverse 4-metric` = rowMeans(comp))
  colnames(comp)[1:4] <- c("FA lower-direction", "MD higher-direction", "ICVF lower-direction", "ISOVF higher-direction")
  for (weighted in c(FALSE, TRUE)) {
    ww <- if (weighted) dd$sw_trim else rep(1, nrow(dd))
    tr <- fit_many(dd, Y, rhs, sp$term, ww)
    tr$analysis <- sp$name; tr$exposure_label <- sp$label; tr$selection_IPW <- weighted; tr$n <- nrow(dd)
    tr$description <- unname(description[tr$outcome]); tr$metric <- unname(metric[tr$outcome]); tr$adverse_sign <- unname(adverse_sign[tr$outcome])
    tr$adverse_beta <- tr$beta * tr$adverse_sign
    tr$adverse_CI_low <- ifelse(tr$adverse_sign == 1, tr$CI_low, -tr$CI_high)
    tr$adverse_CI_high <- ifelse(tr$adverse_sign == 1, tr$CI_high, -tr$CI_low)
    tract_rows[[kt]] <- tr; kt <- kt + 1L
    cr <- fit_many(dd, comp, rhs, sp$term, ww)
    cr$analysis <- sp$name; cr$exposure_label <- sp$label; cr$selection_IPW <- weighted; cr$n <- nrow(dd)
    composite_rows[[kc]] <- cr; kc <- kc + 1L
  }
  audit_rows[[ka]] <- data.frame(analysis = sp$name, n = nrow(dd), exposure = sp$term, stringsAsFactors = FALSE); ka <- ka + 1L
}

tract <- do.call(rbind, tract_rows)
grp <- interaction(tract$analysis, tract$selection_IPW, drop = TRUE)
tract$p_FDR_134 <- ave(tract$p, grp, FUN = function(x) p.adjust(x, "BH"))
tract$p_FDR_within_metric <- ave(tract$p, interaction(grp, tract$metric, drop = TRUE), FUN = function(x) p.adjust(x, "BH"))
tract <- tract[, c("analysis", "exposure_label", "selection_IPW", "n", "outcome", "description", "metric", "beta", "HC3_SE", "CI_low", "CI_high", "p", "p_FDR_134", "p_FDR_within_metric", "adverse_beta", "adverse_CI_low", "adverse_CI_high")]
save_csv(tract, "Table_P11_diffusion_tract_localization.csv")
save_csv(tract[tract$p_FDR_134 < 0.05, ], "Table_P11_significant_tracts_FDR134.csv")

comp <- do.call(rbind, composite_rows)
comp$p_FDR_6 <- ave(comp$p, interaction(comp$analysis, comp$selection_IPW, drop = TRUE), FUN = function(x) p.adjust(x, "BH"))
comp <- comp[, c("analysis", "exposure_label", "selection_IPW", "n", "outcome", "beta", "HC3_SE", "CI_low", "CI_high", "p", "p_FDR_6")]
save_csv(comp, "Table_P11_diffusion_metric_composites.csv")

sig_summary <- aggregate(as.integer(tract$p_FDR_134 < 0.05), by = list(analysis = tract$analysis, selection_IPW = tract$selection_IPW), FUN = sum)
names(sig_summary)[3] <- "n_FDR134_significant"
sig_summary$n_tested <- 134L
save_csv(sig_summary, "Table_P11_localization_significance_summary.csv")
save_csv(do.call(rbind, audit_rows), "P11_analysis_sample_audit.csv")

notes <- c(
  "White-matter localization is a mechanistic follow-up using raw weighted-mean diffusion tract IDPs, not feature-importance attribution from the brain-age model.",
  "Five prespecified core metrics were used: FA, MD, ICVF, ISOVF, and OD.",
  "The local UKB code-list mapping was read directly. wt_fa_tfm_2 was excluded because only 1,843 of 16,535 formal OOS MRI participants had data and its local code-list Instance label was misaligned.",
  "wt_md_tml_r_2 was restored using the observed symmetric left/right naming across the other metrics; the analysis includes 134 complete core tract metrics.",
  "FA and ICVF were sign-reversed; MD and ISOVF retained their higher-is-adverse direction. OD was reported in its original direction but excluded from the global adverse composite because it has no uniform adverse direction across tracts.",
  "Models adjust the other seven LE8 components, flexible MRI age, sex, ethnicity, deprivation, education, MRI centre, WebQ energy and timing, and scanner/image covariates.",
  "Primary tract multiplicity control is Benjamini-Hochberg FDR across all 134 outcomes per analysis and weighting specification; within-metric FDR is secondary.",
  "Selection-IPW results are sensitivity analyses."
)
writeLines(notes, file.path(out_dir, "P11_analysis_notes.txt"), useBytes = TRUE)
cat("P11 white-matter localization completed.\n")
