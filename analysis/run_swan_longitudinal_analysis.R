suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(nlme)
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) stop("Run this script from the repository root.")
input_dir <- Sys.getenv("SWAN_DIR", unset = file.path(root, "data", "restricted", "swan"))
output_dir <- file.path(root, "outputs", "swan_longitudinal_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

study_visit <- c(
  `28762` = 0L, `29221` = 1L, `29401` = 2L, `29701` = 3L,
  `30142` = 4L, `30501` = 5L, `31181` = 6L, `31901` = 7L,
  `32122` = 8L, `32721` = 9L, `32961` = 10L
)
metabolic_visits <- c(0L, 1L, 3L, 4L, 5L, 6L, 7L)

get_num <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) return(rep(NA_real_, nrow(data)))
  as.numeric(zap_missing(data[[hit[[1]]]]))
}

is_yes_any <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) return(rep(FALSE, nrow(data)))
  Reduce(`|`, lapply(hit, function(v) get_num(data, v) == 2))
}

row_mean_available <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) return(rep(NA_real_, nrow(data)))
  values <- as.data.frame(lapply(data[hit], function(x) as.numeric(zap_missing(x))))
  out <- rowMeans(values, na.rm = TRUE)
  out[rowSums(is.finite(as.matrix(values))) == 0] <- NA_real_
  out
}

extract_visit <- function(file) {
  study <- sub("-.*", "", basename(file))
  visit <- unname(study_visit[[study]])
  s <- as.character(visit)
  data <- read_dta(file)

  steroid_use <- if (visit == 0L) {
    get_num(data, "STEROID0") == 2
  } else {
    is_yes_any(data, c(paste0("STERTW1", s), paste0("STERTW2", s)))
  }
  diabetes_med <- if (visit == 0L) {
    get_num(data, "INSUYS0") == 2 | get_num(data, "INSULIN0") == 2
  } else {
    is_yes_any(data, c(paste0("INSUTW1", s), paste0("INSUTW2", s)))
  }
  lipid_med <- if (visit == 0L) {
    get_num(data, "CHOLYS0") == 2 | get_num(data, "CHOLEST0") == 2
  } else {
    is_yes_any(data, c(paste0("CHOLTW1", s), paste0("CHOLTW2", s)))
  }

  tibble(
    swanid = as.character(zap_missing(data$SWANID)),
    study = study,
    visit = visit,
    age = get_num(data, paste0("AGE", s)),
    race = get_num(data, "RACE"),
    menopause_status = get_num(data, paste0("STATUS", s)),
    fasting = get_num(data, paste0("EATDRIN", s)) == 1,
    steroid_use = replace_na(steroid_use, FALSE),
    diabetes_med = replace_na(diabetes_med, FALSE),
    lipid_med = replace_na(lipid_med, FALSE),
    testosterone = get_num(data, paste0("T", s)),
    dheas = get_num(data, paste0("DHAS", s)),
    shbg = get_num(data, paste0("SHBG", s)),
    glucose = get_num(data, c(paste0("GLUCRES", s), paste0("GLUCRE", s))),
    insulin = get_num(data, paste0("INSURES", s)),
    triglycerides = get_num(data, paste0("TRIGRES", s)),
    hdl = get_num(data, paste0("HDLRESU", s)),
    bmi = get_num(data, paste0("BMI", s)),
    waist = get_num(data, paste0("WAIST", s)),
    systolic_bp = row_mean_available(data, c(paste0("SYSBP1", s), paste0("SYSBP2", s), paste0("SYSBP3", s))),
    diastolic_bp = row_mean_available(data, c(paste0("DIABP1", s), paste0("DIABP2", s), paste0("DIABP3", s)))
  )
}

files <- Sys.glob(file.path(input_dir, "*-Data.dta"))
if (!length(files)) stop("No SWAN visit files found. Set SWAN_DIR to the local directory containing *-Data.dta files.")
long_all <- map_dfr(files, extract_visit) %>% arrange(swanid, visit)

baseline_covariates <- long_all %>%
  filter(visit == 0L) %>%
  transmute(
    swanid,
    race_baseline = race,
    smoking_baseline = {
      raw <- read_dta(file.path(input_dir, "28762-0001-Data.dta"), col_select = c(SWANID, SMOKENO0))
      as.numeric(zap_missing(raw$SMOKENO0)) == 2
    }
  )

long_all <- long_all %>%
  select(-race) %>%
  left_join(baseline_covariates, by = "swanid") %>%
  mutate(
    natural_stage = menopause_status %in% 2:5,
    positive_hormones = testosterone > 0 & dheas > 0 & shbg > 0,
    positive_metabolic = glucose > 0 & insulin > 0 & triglycerides > 0 & hdl > 0,
    primary_eligible = visit %in% metabolic_visits & fasting & natural_stage & !steroid_use &
      positive_hormones & positive_metabolic & is.finite(age) & is.finite(race_baseline),
    homa_ir = glucose * insulin / 405,
    map = (systolic_bp + 2 * diastolic_bp) / 3,
    fai = 100 * (testosterone * 0.0347) / shbg
  )

z_within_visit <- function(x, visit) {
  ave(x, visit, FUN = function(v) {
    ok <- is.finite(v)
    out <- rep(NA_real_, length(v))
    if (sum(ok) > 1 && sd(v[ok]) > 0) out[ok] <- as.numeric(scale(v[ok]))
    out
  })
}

analysis <- long_all %>%
  filter(primary_eligible) %>%
  mutate(
    z_log_testosterone = z_within_visit(log(testosterone), visit),
    z_log_dheas = z_within_visit(log(dheas), visit),
    z_log_shbg = z_within_visit(log(shbg), visit),
    z_log_fai = z_within_visit(log(fai), visit),
    z_log_homa = z_within_visit(log(homa_ir), visit),
    z_log_triglycerides = z_within_visit(log(triglycerides), visit),
    z_inverse_hdl = -z_within_visit(hdl, visit),
    z_map = z_within_visit(map, visit),
    z_bmi = z_within_visit(bmi, visit),
    z_waist = z_within_visit(waist, visit),
    metabolic_score = rowMeans(across(c(z_log_homa, z_log_triglycerides, z_inverse_hdl, z_map)), na.rm = FALSE),
    age_centered = age - 46,
    race_factor = factor(race_baseline),
    status_factor = factor(menopause_status),
    visit_factor = factor(visit)
  )

decompose <- function(data, variable, prefix) {
  mean_name <- paste0(prefix, "_between")
  within_name <- paste0(prefix, "_within")
  data %>% group_by(swanid) %>% mutate(
    "{mean_name}" := mean(.data[[variable]], na.rm = TRUE),
    "{within_name}" := .data[[variable]] - .data[[mean_name]]
  ) %>% ungroup()
}

analysis <- analysis %>%
  decompose("z_log_testosterone", "testosterone") %>%
  decompose("z_log_dheas", "dheas") %>%
  decompose("z_log_shbg", "shbg") %>%
  decompose("z_log_fai", "fai") %>%
  decompose("z_bmi", "bmi") %>%
  decompose("z_waist", "waist") %>%
  group_by(swanid) %>% mutate(n_primary_visits = n()) %>% ungroup() %>%
  filter(n_primary_visits >= 2)

fit_lme <- function(data, outcome, exposure, stage = "base", sensitivity = "primary") {
  between <- paste0(exposure, "_between")
  within <- paste0(exposure, "_within")
  adiposity <- switch(stage,
    base = character(),
    bmi = c("bmi_between", "bmi_within"),
    waist = c("waist_between", "waist_within")
  )
  rhs <- c(between, within, "visit_factor", "age_centered", "race_factor", "status_factor", "smoking_baseline", adiposity)
  needed <- c(outcome, rhs, "swanid")
  d <- data[complete.cases(data[, needed]), , drop = FALSE]
  form <- reformulate(rhs, outcome)
  fit <- try(lme(form, random = ~1 | swanid, data = d, method = "REML", na.action = na.omit,
                 control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(tibble())
  tt <- summary(fit)$tTable
  terms <- c(between, within)
  bind_rows(lapply(terms, function(term) tibble(
    outcome = outcome, exposure = exposure, stage = stage, sensitivity = sensitivity,
    component = ifelse(term == within, "within-person", "between-person"),
    estimate = tt[term, "Value"], std_error = tt[term, "Std.Error"],
    df = tt[term, "DF"], p_value = tt[term, "p-value"],
    conf_low = tt[term, "Value"] - 1.96 * tt[term, "Std.Error"],
    conf_high = tt[term, "Value"] + 1.96 * tt[term, "Std.Error"],
    n_observations = nrow(d), n_participants = n_distinct(d$swanid), AIC = AIC(fit)
  )))
}

exposures <- c("testosterone", "dheas", "shbg", "fai")
outcomes <- c("metabolic_score", "z_log_homa")
stages <- c("base", "bmi", "waist")

primary_models <- expand_grid(outcome = outcomes, exposure = exposures, stage = stages) %>%
  pmap_dfr(~fit_lme(analysis, ..1, ..2, ..3, "primary")) %>%
  group_by(outcome, component, stage) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>% ungroup()

fit_joint <- function(data, outcome, stage = "base") {
  adiposity <- switch(stage, base = character(), bmi = c("bmi_between", "bmi_within"), waist = c("waist_between", "waist_within"))
  hormones <- unlist(lapply(c("testosterone", "dheas", "shbg"), function(x) paste0(x, c("_between", "_within"))))
  rhs <- c(hormones, "visit_factor", "age_centered", "race_factor", "status_factor", "smoking_baseline", adiposity)
  needed <- c(outcome, rhs, "swanid")
  d <- data[complete.cases(data[, needed]), , drop = FALSE]
  fit <- lme(reformulate(rhs, outcome), random = ~1 | swanid, data = d, method = "REML", na.action = na.omit,
             control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
  tt <- summary(fit)$tTable
  bind_rows(lapply(hormones, function(term) tibble(
    outcome = outcome, stage = stage, exposure = sub("_(between|within)$", "", term),
    component = ifelse(grepl("within$", term), "within-person", "between-person"),
    estimate = tt[term, "Value"], std_error = tt[term, "Std.Error"], df = tt[term, "DF"],
    p_value = tt[term, "p-value"], conf_low = tt[term, "Value"] - 1.96 * tt[term, "Std.Error"],
    conf_high = tt[term, "Value"] + 1.96 * tt[term, "Std.Error"],
    n_observations = nrow(d), n_participants = n_distinct(d$swanid), AIC = AIC(fit)
  )))
}

joint_models <- expand_grid(outcome = outcomes, stage = stages) %>%
  pmap_dfr(~fit_joint(analysis, ..1, ..2)) %>%
  group_by(outcome, component, stage) %>% mutate(p_fdr = p.adjust(p_value, method = "BH")) %>% ungroup()

sensitivity_sets <- list(
  no_metabolic_medications = analysis %>% filter(!diabetes_med, !lipid_med),
  menopause_transition_only = analysis %>% filter(menopause_status %in% 3:5)
)
sensitivity_models <- imap_dfr(sensitivity_sets, function(d, label) {
  expand_grid(outcome = outcomes, exposure = c("testosterone", "dheas", "shbg"), stage = c("base", "bmi")) %>%
    pmap_dfr(~fit_lme(d, ..1, ..2, ..3, label))
}) %>% group_by(sensitivity, outcome, component, stage) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>% ungroup()

fit_car1 <- function(data, outcome, exposure, stage = "base") {
  between <- paste0(exposure, "_between")
  within <- paste0(exposure, "_within")
  adiposity <- if (stage == "bmi") c("bmi_between", "bmi_within") else character()
  rhs <- c(between, within, "visit_factor", "age_centered", "race_factor", "status_factor", "smoking_baseline", adiposity)
  needed <- c(outcome, rhs, "swanid", "visit")
  d <- data[complete.cases(data[, needed]), , drop = FALSE] %>% arrange(swanid, visit)
  fit <- try(lme(reformulate(rhs, outcome), random = ~1 | swanid,
                 correlation = corCAR1(form = ~visit | swanid), data = d,
                 method = "REML", na.action = na.omit,
                 control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(tibble())
  tt <- summary(fit)$tTable
  tibble(
    outcome = outcome, exposure = exposure, stage = stage, component = "within-person",
    correlation = "continuous AR(1)", estimate = tt[within,"Value"], std_error = tt[within,"Std.Error"],
    df = tt[within,"DF"], p_value = tt[within,"p-value"],
    conf_low = tt[within,"Value"] - 1.96*tt[within,"Std.Error"],
    conf_high = tt[within,"Value"] + 1.96*tt[within,"Std.Error"],
    n_observations = nrow(d), n_participants = n_distinct(d$swanid),
    residual_correlation = as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  )
}

car1_models <- expand_grid(outcome = outcomes, exposure = c("testosterone", "dheas", "shbg"), stage = c("base", "bmi")) %>%
  pmap_dfr(~fit_car1(analysis, ..1, ..2, ..3)) %>%
  group_by(outcome, stage) %>% mutate(p_fdr = p.adjust(p_value, method = "BH")) %>% ungroup()

lag_data <- analysis %>%
  arrange(swanid, visit) %>% group_by(swanid) %>%
  mutate(next_visit = lead(visit), next_metabolic_score = lead(metabolic_score), next_z_log_homa = lead(z_log_homa),
         time_gap = next_visit - visit) %>% ungroup() %>%
  filter(is.finite(next_metabolic_score), time_gap > 0)

fit_lag <- function(data, outcome_next, outcome_current, exposure, adjust_current = FALSE) {
  rhs <- c(paste0(exposure, "_between"), paste0(exposure, "_within"),
           "visit_factor", "time_gap", "age_centered", "race_factor", "status_factor", "smoking_baseline")
  if (adjust_current) rhs <- c(rhs, outcome_current)
  needed <- c(outcome_next, rhs, "swanid")
  d <- data[complete.cases(data[, needed]), , drop = FALSE]
  fit <- lme(reformulate(rhs, outcome_next), random = ~1 | swanid, data = d, method = "REML", na.action = na.omit,
             control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
  tt <- summary(fit)$tTable
  term <- paste0(exposure, "_within")
  tibble(
    outcome = outcome_next, exposure = exposure,
    model_type = ifelse(adjust_current, "conditional change (current outcome adjusted)", "next outcome (current outcome unadjusted)"),
    component = "within-person lagged", estimate = tt[term,"Value"],
    std_error = tt[term,"Std.Error"], df = tt[term,"DF"], p_value = tt[term,"p-value"],
    conf_low = tt[term,"Value"] - 1.96*tt[term,"Std.Error"], conf_high = tt[term,"Value"] + 1.96*tt[term,"Std.Error"],
    n_observations = nrow(d), n_participants = n_distinct(d$swanid)
  )
}

lag_models <- bind_rows(lapply(c("testosterone", "dheas", "shbg"), function(e) {
  bind_rows(
    fit_lag(lag_data, "next_metabolic_score", "metabolic_score", e, FALSE),
    fit_lag(lag_data, "next_z_log_homa", "z_log_homa", e, FALSE),
    fit_lag(lag_data, "next_metabolic_score", "metabolic_score", e, TRUE),
    fit_lag(lag_data, "next_z_log_homa", "z_log_homa", e, TRUE)
  )
})) %>% group_by(outcome, model_type) %>% mutate(p_fdr = p.adjust(p_value, method = "BH")) %>% ungroup()

visit_counts <- analysis %>% group_by(visit) %>% summarise(
  observations = n(), participants = n_distinct(swanid),
  mean_age = mean(age), natural_post_percent = mean(menopause_status == 2) * 100,
  mean_homa = mean(homa_ir), median_homa = median(homa_ir),
  .groups = "drop"
)

write_csv(primary_models, file.path(output_dir, "primary_mixed_models.csv"))
write_csv(joint_models, file.path(output_dir, "joint_hormone_mixed_models.csv"))
write_csv(sensitivity_models, file.path(output_dir, "sensitivity_mixed_models.csv"))
write_csv(car1_models, file.path(output_dir, "car1_sensitivity_models.csv"))
write_csv(lag_models, file.path(output_dir, "lagged_models.csv"))
write_csv(visit_counts, file.path(output_dir, "analytic_visit_counts.csv"))

palette <- c(testosterone = "#4C78A8", dheas = "#8A6FB6", shbg = "#D55E5E", fai = "#2A9D8F")
plot_data <- primary_models %>%
  filter(outcome == "metabolic_score", component == "within-person", stage %in% c("base", "bmi"), exposure != "fai") %>%
  mutate(
    exposure_label = recode(exposure, testosterone = "Testosterone", dheas = "DHEAS", shbg = "SHBG"),
    stage_label = recode(stage, base = "Core adjusted", bmi = "+ BMI"),
    exposure_label = factor(exposure_label, levels = c("DHEAS", "Testosterone", "SHBG"))
  )

theme_pub <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 9),
        plot.tag = element_text(face = "bold", size = 10))

p_a <- ggplot(plot_data, aes(estimate, exposure_label, colour = exposure, shape = stage_label)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.14, orientation = "y",
                position = position_dodge(width = 0.45), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.45), size = 2.1) +
  scale_colour_manual(values = palette, labels = c(dheas="DHEAS", shbg="SHBG", testosterone="Testosterone")) +
  labs(title = "Separate hormone models", x = "Within-person association with metabolic score\n(SD per 1-SD hormone deviation)", y = NULL, shape = NULL) + theme_pub

joint_plot_data <- joint_models %>%
  filter(outcome == "metabolic_score", component == "within-person", stage %in% c("base", "bmi")) %>%
  mutate(
    exposure_label = recode(exposure, testosterone = "Testosterone", dheas = "DHEAS", shbg = "SHBG"),
    stage_label = recode(stage, base = "Core adjusted", bmi = "+ BMI"),
    exposure_label = factor(exposure_label, levels = c("DHEAS", "Testosterone", "SHBG"))
  )

p_b <- ggplot(joint_plot_data, aes(estimate, exposure_label, colour = exposure, shape = stage_label)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.14, orientation = "y",
                position = position_dodge(width = 0.45), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.45), size = 2.1) +
  scale_colour_manual(values = palette, labels = c(dheas="DHEAS", shbg="SHBG", testosterone="Testosterone")) +
  labs(title = "Mutually adjusted hormone model", x = "Within-person association with metabolic score\n(SD per 1-SD hormone deviation)", y = NULL, shape = NULL) + theme_pub

fig <- (p_a | p_b) + plot_layout(widths = c(1.25, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
ggsave(file.path(output_dir, "Figure_SWAN_longitudinal.png"), fig, width = 183/25.4, height = 82/25.4, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Figure_SWAN_longitudinal.pdf"), fig, width = 183/25.4, height = 82/25.4, device = cairo_pdf, family = "Arial", bg = "white")
ggsave(file.path(output_dir, "Figure_SWAN_longitudinal.tiff"), fig, width = 183/25.4, height = 82/25.4, dpi = 600, compression = "lzw", bg = "white")

cat("Primary analysis observations:", nrow(analysis), "\n")
cat("Primary analysis participants:", n_distinct(analysis$swanid), "\n")
print(primary_models %>% filter(component == "within-person", outcome == "metabolic_score"), n = Inf, width = Inf)
print(lag_models, n = Inf, width = Inf)
