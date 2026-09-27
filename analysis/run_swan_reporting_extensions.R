# Aggregate reporting extensions for Supplementary Tables S32-S37.
# Sourced by run_swan_longitudinal_analysis.R after `analysis` and `long_all`
# have been constructed. Participant-level data are never written.

report_dir <- file.path(root, "outputs", "swan_reporting_extensions")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

clean_number <- function(x) {
  value <- as.numeric(zap_missing(x))
  value[value < 0] <- NA_real_
  value
}

baseline_raw <- read_dta(file.path(input_dir, "28762-0001-Data.dta"))
baseline_extra <- tibble(
  swanid = as.character(zap_missing(baseline_raw$SWANID)),
  income = clean_number(baseline_raw$INCOME0),
  recreational_activity = clean_number(baseline_raw$RECACTI0),
  alcohol_servings_day = clean_number(baseline_raw$ALCHSRV0)
)

extended_analysis <- analysis %>% left_join(baseline_extra, by = "swanid") %>% mutate(
  income_factor = factor(income),
  recreational_activity_factor = factor(recreational_activity),
  log1p_alcohol = log1p(alcohol_servings_day),
  diabetes_med_factor = factor(diabetes_med),
  lipid_med_factor = factor(lipid_med)
)

fit_extended <- function(outcome, exposure, stage) {
  rhs <- c(paste0(exposure, c("_between", "_within")), "visit_factor", "age_centered",
           "race_factor", "status_factor", "smoking_baseline")
  if (stage == "bmi") rhs <- c(rhs, "bmi_between", "bmi_within")
  if (stage == "waist") rhs <- c(rhs, "waist_between", "waist_within")
  rhs <- c(rhs, "income_factor", "recreational_activity_factor", "log1p_alcohol",
           "diabetes_med_factor", "lipid_med_factor")
  d <- extended_analysis[complete.cases(extended_analysis[, c(outcome, rhs, "swanid")]), , drop = FALSE]
  fit <- lme(reformulate(rhs, outcome), random = ~1 | swanid, data = d, method = "REML",
             control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
  term <- paste0(exposure, "_within")
  tt <- summary(fit)$tTable
  tibble(outcome = outcome, exposure = exposure, adiposity_stage = stage,
         adjustment = "Primary covariates plus baseline income, recreational activity, alcohol intake, and time-varying diabetes/lipid-lowering medication use",
         estimate = tt[term, "Value"], std_error = tt[term, "Std.Error"],
         conf_low = tt[term, "Value"] - 1.96 * tt[term, "Std.Error"],
         conf_high = tt[term, "Value"] + 1.96 * tt[term, "Std.Error"],
         p_value = tt[term, "p-value"], observations = nrow(d), participants = n_distinct(d$swanid))
}

extended_models <- expand_grid(
  outcome = c("metabolic_score", "z_log_homa"),
  exposure = c("testosterone", "dheas", "shbg", "fai"),
  stage = c("bmi", "waist")
) %>% pmap_dfr(fit_extended)
write_csv(extended_models, file.path(report_dir, "S32_extended_adjustment_models.csv"))

baseline <- long_all %>% filter(visit == 0L) %>% left_join(baseline_extra, by = "swanid") %>% mutate(
  included = swanid %in% unique(analysis$swanid),
  homa_ir = glucose * insulin / 405,
  map = (systolic_bp + 2 * diastolic_bp) / 3,
  fai = 100 * (testosterone * 0.0347) / shbg,
  current_smoker = smoking_baseline,
  natural_post = menopause_status == 2
)

format_mean <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
format_median <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x, na.rm = TRUE),
  quantile(x, .25, na.rm = TRUE, type = 2), quantile(x, .75, na.rm = TRUE, type = 2))
format_category <- function(x, condition) sprintf("%d (%.1f%%)", sum(condition & !is.na(x)),
  100 * mean(condition[!is.na(x)]))

make_baseline_table <- function(selected, label) {
  d <- baseline[selected, , drop = FALSE]
  continuous <- list(
    c("Age, years", "age", "mean"), c("Body mass index, kg/m2", "bmi", "mean"),
    c("Waist circumference, cm", "waist", "mean"), c("HOMA-IR", "homa_ir", "median"),
    c("Triglycerides, mg/dL", "triglycerides", "median"), c("HDL cholesterol, mg/dL", "hdl", "mean"),
    c("Mean arterial pressure, mmHg", "map", "mean"), c("Testosterone, ng/dL", "testosterone", "median"),
    c("DHEAS, µg/dL", "dheas", "median"), c("SHBG, nmol/L", "shbg", "median"),
    c("Free androgen index", "fai", "median"), c("Alcohol servings/day", "alcohol_servings_day", "median")
  )
  continuous_rows <- bind_rows(lapply(continuous, function(spec) {
    x <- d[[spec[2]]]
    tibble(group = label, characteristic = spec[1],
           value = if (spec[3] == "mean") format_mean(x) else format_median(x),
           nonmissing_n = sum(!is.na(x)))
  }))
  race_rows <- bind_rows(lapply(list(c("Black/African American", 1), c("Chinese/Chinese American", 2),
    c("Japanese/Japanese American", 3), c("White, non-Hispanic", 4), c("Hispanic", 5)), function(spec) {
      x <- d$race_baseline
      tibble(group = label, characteristic = spec[1], value = format_category(x, x == as.numeric(spec[2])),
             nonmissing_n = sum(!is.na(x)))
  }))
  categorical <- list(
    list("Current smoker", "current_smoker", function(x) x),
    list("Natural postmenopause", "natural_post", function(x) x),
    list("Diabetes medication use", "diabetes_med", function(x) x),
    list("Lipid-lowering medication use", "lipid_med", function(x) x),
    list("Family income <20,000 USD", "income", function(x) x == 1),
    list("Recreational activity at least same as peers", "recreational_activity", function(x) x >= 3)
  )
  category_rows <- bind_rows(lapply(categorical, function(spec) {
    x <- d[[spec[[2]]]]
    tibble(group = label, characteristic = spec[[1]], value = format_category(x, spec[[3]](x)),
           nonmissing_n = sum(!is.na(x)))
  }))
  bind_rows(continuous_rows %>% slice(1), race_rows, continuous_rows %>% slice(-1), category_rows)
}

included <- make_baseline_table(baseline$included, "Included (N=2746)")
not_included <- make_baseline_table(!baseline$included, "Not included (N=556)")
comparison <- bind_rows(included, not_included) %>%
  pivot_wider(names_from = group, values_from = c(value, nonmissing_n))
write_csv(comparison, file.path(report_dir, "S33_included_vs_not_included.csv"))

candidate <- long_all %>% mutate(
  structural_eligible = visit %in% metabolic_visits & fasting & menopause_status %in% 2:5 &
    !steroid_use & is.finite(age) & is.finite(race_baseline),
  map = (systolic_bp + 2 * diastolic_bp) / 3
) %>% filter(structural_eligible)

variables <- c("testosterone", "dheas", "shbg", "glucose", "insulin", "triglycerides", "hdl", "map", "bmi", "waist")
missingness <- bind_rows(lapply(variables, function(variable) {
  x <- candidate[[variable]]
  tibble(variable = variable, candidate_observations = nrow(candidate),
         missing_or_nonpositive = sum(!is.finite(x) | x <= 0),
         percent_missing_or_nonpositive = 100 * mean(!is.finite(x) | x <= 0))
}))
write_csv(missingness, file.path(report_dir, "S34_variable_missingness.csv"))

flow <- tibble(
  step = c("SWAN baseline enrollment", "At least one observation at a selected metabolic visit",
    "At least one structurally eligible observation",
    "At least one complete eligible hormone-metabolic observation",
    "At least two complete eligible observations (longitudinal cohort)"),
  participants = c(n_distinct(long_all$swanid),
    n_distinct(long_all$swanid[long_all$visit %in% metabolic_visits]),
    n_distinct(candidate$swanid), n_distinct(long_all$swanid[long_all$primary_eligible]),
    n_distinct(analysis$swanid))
)
write_csv(flow, file.path(report_dir, "S35_participant_flow.csv"))

visit_distribution <- analysis %>% distinct(swanid, n_primary_visits) %>%
  count(n_primary_visits, name = "participants") %>% mutate(percent = 100 * participants / sum(participants))
write_csv(visit_distribution, file.path(report_dir, "S36_visit_contribution_distribution.csv"))

# Refit SHBG on its natural-log scale so the coefficient can be translated to a two-fold contrast.
log_analysis <- analysis %>% mutate(raw_log_shbg = log(shbg)) %>%
  group_by(swanid) %>% mutate(raw_log_shbg_between = mean(raw_log_shbg, na.rm = TRUE),
    raw_log_shbg_within = raw_log_shbg - raw_log_shbg_between) %>% ungroup()
fit_doubling <- function(stage) {
  rhs <- c("raw_log_shbg_between", "raw_log_shbg_within", "visit_factor", "age_centered",
           "race_factor", "status_factor", "smoking_baseline")
  if (stage == "bmi") rhs <- c(rhs, "bmi_between", "bmi_within")
  if (stage == "waist") rhs <- c(rhs, "waist_between", "waist_within")
  d <- log_analysis[complete.cases(log_analysis[, c("metabolic_score", rhs, "swanid")]), , drop = FALSE]
  fit <- lme(reformulate(rhs, "metabolic_score"), random = ~1 | swanid, data = d, method = "REML",
             control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
  tt <- summary(fit)$tTable; term <- "raw_log_shbg_within"
  log_contrast_value <- base::log(2)
  beta_value <- as.double(unlist(tt[term, "Value"], use.names = FALSE))
  se_value <- as.double(unlist(tt[term, "Std.Error"], use.names = FALSE))
  tibble(stage = stage, contrast = "Two-fold higher SHBG within person", log_contrast = log_contrast_value,
    estimate_for_doubling = beta_value * log_contrast_value,
    conf_low_for_doubling = (beta_value - 1.96 * se_value) * log_contrast_value,
    conf_high_for_doubling = (beta_value + 1.96 * se_value) * log_contrast_value,
    p_value = as.numeric(tt[term, "p-value"]), observations = nrow(d), participants = n_distinct(d$swanid))
}
shbg_doubling <- bind_rows(lapply(c("base", "bmi", "waist"), fit_doubling))
write_csv(shbg_doubling, file.path(report_dir, "S37_SHBG_doubling_interpretation.csv"))

message("Reporting extensions complete: aggregate Tables S32-S37 written to ", report_dir)
