prepare_common_sample <- function(data) {
  required <- c(
    "z_testosterone", "z_androstenedione", "z_dheas", "z_amh", "z_shbg",
    "z_fai", "z_homa_ir", "z_insulin", "z_glucose", "z_triglycerides",
    "z_hba1c", "z_inverse_hdl", "z_map", "nonfasting_score", "fasting_score",
    "waist_height_ratio", "RIDAGEYR", "RIDRETH3", "INDFMPIR", "current_smoking"
  )
  data$common_complete <- data$strict_eligible & complete.cases(data[, required]) &
    positive_weight(data$fasting_weight)
  common <- data[data$common_complete, , drop = FALSE]
  weights <- common$fasting_weight

  rest <- c(
    testosterone = "z_testosterone", androstenedione = "z_androstenedione",
    dheas = "z_dheas", amh = "z_amh", shbg = "z_shbg", fai = "z_fai",
    homa_ir = "z_homa_ir", insulin = "z_insulin", glucose = "z_glucose",
    triglycerides = "z_triglycerides", hba1c = "z_hba1c",
    inverse_hdl = "z_inverse_hdl", map = "z_map",
    nonfasting_score = "nonfasting_score", fasting_score = "fasting_score",
    adiposity = "waist_height_ratio"
  )
  for (target in names(rest)) {
    common[[paste0("c_", target)]] <- weighted_center_scale(common[[rest[[target]]]], weights)
  }
  common$c_inverse_shbg <- -common$c_shbg
  common$c_three_androgens <- rowMeans(common[, c("c_testosterone", "c_androstenedione", "c_dheas")])
  common$c_three_plus_amh <- rowMeans(common[, c("c_testosterone", "c_androstenedione", "c_dheas", "c_amh")])
  common$c_three_plus_inverse_shbg <- rowMeans(common[, c("c_testosterone", "c_androstenedione", "c_dheas", "c_inverse_shbg")])
  common$c_all_five <- rowMeans(common[, c("c_testosterone", "c_androstenedione", "c_dheas", "c_amh", "c_inverse_shbg")])

  common
}

make_common_design <- function(common) {
  svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~fasting_weight,
    nest = TRUE,
    data = common
  )
}

run_operationalization_matrix <- function(common, root = project_root()) {
  design <- make_common_design(common)
  exposures <- c(
    testosterone = "c_testosterone", androstenedione = "c_androstenedione",
    dheas = "c_dheas", amh = "c_amh", shbg = "c_shbg", fai = "c_fai",
    three_androgens = "c_three_androgens", three_plus_amh = "c_three_plus_amh",
    three_plus_inverse_shbg = "c_three_plus_inverse_shbg", all_five = "c_all_five"
  )
  outcomes <- c(
    nonfasting_score = "c_nonfasting_score", fasting_score = "c_fasting_score",
    homa_ir = "c_homa_ir", insulin = "c_insulin", glucose = "c_glucose",
    triglycerides = "c_triglycerides", hba1c = "c_hba1c",
    inverse_hdl = "c_inverse_hdl", map = "c_map"
  )
  adjustment <- list(
    base = c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking"),
    adiposity = c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking", "c_adiposity")
  )
  rows <- list()
  index <- 0L
  for (stage in names(adjustment)) {
    for (exposure_name in names(exposures)) {
      for (outcome_name in names(outcomes)) {
        exposure <- exposures[[exposure_name]]
        outcome <- outcomes[[outcome_name]]
        fit <- svyglm(model_formula(outcome, exposure, adjustment[[stage]]), design = design)
        index <- index + 1L
        rows[[index]] <- tidy_svyglm(
          fit, exposure, "common_sample_matrix", stage, exposure_name, outcome_name
        )
      }
    }
  }
  result <- bind_rows(rows) |>
    group_by(stage) |>
    mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
    ungroup()
  write_csv_safe(result, file.path(root, "outputs", "tables", "operationalization_construct_matrix.csv"))
  result
}

run_incremental_fit <- function(common, root = project_root()) {
  design <- make_common_design(common)
  outcomes <- c("c_nonfasting_score", "c_fasting_score", "c_homa_ir")
  base_terms <- c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking", "c_adiposity")
  orders <- list(
    androgens_then_amh = c("c_three_androgens", "c_amh"),
    amh_then_androgens = c("c_amh", "c_three_androgens"),
    androgens_then_inverse_shbg = c("c_three_androgens", "c_inverse_shbg"),
    inverse_shbg_then_androgens = c("c_inverse_shbg", "c_three_androgens")
  )
  rows <- list()
  index <- 0L
  for (outcome in outcomes) {
    base <- svyglm(model_formula(outcome, base_terms[1], base_terms[-1]), design = design, y = TRUE)
    base_r2 <- weighted_r_squared(base)
    for (order_name in names(orders)) {
      cumulative <- base_terms
      prior_r2 <- base_r2
      for (step in seq_along(orders[[order_name]])) {
        term <- orders[[order_name]][step]
        cumulative <- c(cumulative, term)
        fit <- svyglm(reformulate(cumulative, response = outcome), design = design, y = TRUE)
        current_r2 <- weighted_r_squared(fit)
        index <- index + 1L
        rows[[index]] <- tibble(
          outcome = sub("^c_", "", outcome), order = order_name, step = step,
          added_term = sub("^c_", "", term), n = nobs(fit),
          r_squared = current_r2, delta_r_squared = current_r2 - prior_r2
        )
        prior_r2 <- current_r2
      }
    }
  }
  result <- bind_rows(rows)
  write_csv_safe(result, file.path(root, "outputs", "tables", "incremental_construct_fit.csv"))
  result
}

