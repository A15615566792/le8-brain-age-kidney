args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript 20_assign_brain_age_splits.R input.csv output.csv [seed]")
input_path <- args[1]
output_path <- args[2]
seed <- if (length(args) >= 3L) as.integer(args[3]) else 20260714L

dat <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("participant_id", "age")
if (!all(required %in% names(dat))) stop("Missing required columns: ", paste(setdiff(required, names(dat)), collapse = ", "))

if ("healthy" %in% names(dat)) {
  healthy <- suppressWarnings(as.numeric(dat$healthy))
} else {
  health_columns <- c("mri_date", "earliest_inpatient_icd10_date", "long_standing_illness", "self_rated_health")
  if (!all(health_columns %in% names(dat))) stop("Missing health-definition columns: ", paste(setdiff(health_columns, names(dat)), collapse = ", "))
  mri_date <- as.Date(dat$mri_date)
  diagnosis_date <- as.Date(dat$earliest_inpatient_icd10_date)
  illness_text <- tolower(trimws(as.character(dat$long_standing_illness)))
  illness_num <- suppressWarnings(as.numeric(dat$long_standing_illness))
  illness <- rep(NA, nrow(dat))
  illness[illness_num == 0 | illness_text %in% c("no", "false")] <- 0
  illness[illness_num == 1 | illness_text %in% c("yes", "true")] <- 1
  health_text <- tolower(trimws(as.character(dat$self_rated_health)))
  health_num <- suppressWarnings(as.numeric(dat$self_rated_health))
  favorable <- rep(NA, nrow(dat))
  favorable[health_num %in% c(1, 2) | health_text %in% c("excellent", "good")] <- 1
  favorable[health_num %in% c(3, 4) | health_text %in% c("fair", "poor")] <- 0
  complete <- !is.na(dat$age) & !is.na(mri_date) & !is.na(illness) & !is.na(favorable)
  no_pre_mri_diagnosis <- is.na(diagnosis_date) | diagnosis_date > mri_date
  healthy <- rep(NA_real_, nrow(dat))
  healthy[complete] <- as.numeric(no_pre_mri_diagnosis[complete] & illness[complete] == 0 & favorable[complete] == 1)
}

make_age_strata <- function(age, groups = 10L) {
  breaks <- unique(quantile(age, probs = seq(0, 1, length.out = groups + 1L), na.rm = TRUE, type = 2))
  if (length(breaks) < 3L) return(rep(1L, length(age)))
  as.integer(cut(age, breaks = breaks, include.lowest = TRUE, labels = FALSE))
}

stratified_validation_split <- function(age, validation_prop, seed) {
  set.seed(seed)
  strata <- make_age_strata(age)
  validation <- unlist(lapply(split(seq_along(age), strata), function(ii) {
    if (length(ii) < 2L) return(integer())
    n_take <- max(1L, round(length(ii) * validation_prop))
    sample(ii, min(n_take, length(ii) - 1L))
  }), use.names = FALSE)
  list(train = setdiff(seq_along(age), validation), validation = sort(validation))
}

healthy_idx <- which(healthy == 1)
split_local <- stratified_validation_split(as.numeric(dat$age[healthy_idx]), 0.20, seed)
train_idx <- healthy_idx[split_local$train]
validation_idx <- healthy_idx[split_local$validation]
nonhealthy_idx <- which(healthy == 0)
unknown_idx <- which(is.na(healthy))

split_label <- rep("Other", nrow(dat))
split_label[train_idx] <- "Healthy train"
split_label[validation_idx] <- "Healthy validation"
split_label[nonhealthy_idx] <- "Non-healthy test"
split_label[unknown_idx] <- "Health status unknown"

result <- data.frame(
  participant_id = as.character(dat$participant_id),
  age = as.numeric(dat$age),
  healthy = healthy,
  split = split_label,
  stringsAsFactors = FALSE
)
observed <- table(result$split)
expected <- c("Health status unknown" = 966L, "Healthy train" = 5507L, "Healthy validation" = 1380L, "Non-healthy test" = 32216L)
if (!identical(as.integer(observed[names(expected)]), as.integer(expected))) {
  stop("Split-count audit failed: ", paste(names(observed), observed, collapse = "; "))
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(result, output_path, row.names = FALSE, na = "")
print(observed)
