primary_complete_indicator <- function(data) {
  data$strict_eligible & complete.cases(data[, c(
    "androgen_score", "z_amh", "nonfasting_score", "RIDAGEYR", "RIDRETH3",
    "INDFMPIR", "current_smoking", "waist_height_ratio"
  )])
}

run_primary_models <- function(data, root = project_root()) {
  design <- make_nhanes_design(data, "hormone")
  data$primary_complete <- primary_complete_indicator(data)
  design <- make_nhanes_design(data, "hormone") |>
    subset(primary_complete)

  stages <- list(
    exposure_only = character(),
    demographic = c("RIDAGEYR", "factor(RIDRETH3)"),
    social_behavioral = c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking"),
    adiposity = c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking", "waist_height_ratio")
  )
  exposures <- c(androgen_score = "androgen_score", amh = "z_amh")

  estimates <- list()
  index <- 0L
  for (exposure_name in names(exposures)) {
    exposure <- exposures[[exposure_name]]
    for (stage in names(stages)) {
      fit <- svyglm(model_formula("nonfasting_score", exposure, stages[[stage]]), design = design)
      index <- index + 1L
      estimates[[index]] <- tidy_svyglm(
        fit, exposure, "discovery_primary", stage, exposure_name, "nonfasting_score"
      )
    }
    interaction_formula <- as.formula(paste(
      "nonfasting_score ~", exposure, "* waist_height_ratio + RIDAGEYR +",
      "factor(RIDRETH3) + INDFMPIR + current_smoking"
    ))
    fit_interaction <- svyglm(interaction_formula, design = design)
    term <- paste0(exposure, ":waist_height_ratio")
    index <- index + 1L
    estimates[[index]] <- tidy_svyglm(
      fit_interaction, term, "discovery_primary", "interaction", exposure_name, "nonfasting_score"
    )
  }
  result <- bind_rows(estimates) |>
    mutate(p_holm = p.adjust(p_value, method = "holm"))

  core_complete <- data$strict_eligible & complete.cases(data[, c(
    "androgen_score", "z_amh", "nonfasting_score", "waist_height_ratio"
  )])
  flow <- tibble(
    milestone = c("source_records", "women_age_20_44", "strict_eligible",
                  "nonfasting_core_complete", "fully_adjusted_primary"),
    n = c(
      nrow(data),
      sum(data$female & data$age_eligible, na.rm = TRUE),
      sum(data$strict_eligible, na.rm = TRUE),
      sum(core_complete, na.rm = TRUE),
      sum(data$primary_complete & positive_weight(data$hormone_weight), na.rm = TRUE)
    )
  )
  write_csv_safe(result, file.path(root, "outputs", "tables", "main_models.csv"))
  write_csv_safe(flow, file.path(root, "outputs", "audit", "cohort_flow.csv"))

  correlation_data <- data[data$primary_complete & positive_weight(data$hormone_weight), ]
  correlations <- weighted_cor(
    correlation_data,
    c("androgen_score", "z_amh", "nonfasting_score", "waist_height_ratio"),
    correlation_data$hormone_weight
  )
  write_csv_safe(
    as.data.frame(as.table(correlations)) |>
      rename(variable_1 = Var1, variable_2 = Var2, correlation = Freq),
    file.path(root, "outputs", "tables", "weighted_correlations.csv")
  )
  list(models = result, flow = flow, correlations = correlations, design = design)
}

run_fasting_primary <- function(data, root = project_root()) {
  complete <- data$strict_eligible & complete.cases(data[, c(
    "androgen_score", "z_amh", "fasting_score", "RIDAGEYR", "RIDRETH3",
    "INDFMPIR", "current_smoking", "waist_height_ratio"
  )])
  data$fasting_primary_complete <- complete
  design <- make_nhanes_design(data, "fasting") |> subset(fasting_primary_complete)
  exposures <- c(androgen_score = "androgen_score", amh = "z_amh")
  output <- lapply(names(exposures), function(label) {
    exposure <- exposures[[label]]
    model <- svyglm(
      model_formula(
        "fasting_score", exposure,
        c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking")
      ),
      design = design
    )
    tidy_svyglm(model, exposure, "discovery_primary", "social_behavioral", label, "fasting_score")
  }) |> bind_rows() |> mutate(p_holm = p.adjust(p_value, method = "holm"))
  write_csv_safe(output, file.path(root, "outputs", "tables", "fasting_primary_models.csv"))
  output
}

run_primary_sensitivities <- function(data, root = project_root()) {
  definitions <- list(
    all = rep(TRUE, nrow(data)),
    no_prescriptions = !is.na(data$prescription_any) & !data$prescription_any,
    age_20_39 = data$RIDAGEYR <= 39,
    no_diabetes = is.na(data$diabetes) | !data$diabetes,
    central_99_percent = data$LBXTST <= quantile(data$LBXTST, 0.995, na.rm = TRUE) &
      data$LBXAND <= quantile(data$LBXAND, 0.995, na.rm = TRUE) &
      data$LBXDHE <= quantile(data$LBXDHE, 0.995, na.rm = TRUE)
  )
  substitutions <- list(
    waist_height_ratio = "waist_height_ratio",
    bmi = "BMXBMI",
    waist_circumference = "BMXWAIST"
  )
  rows <- list()
  index <- 0L
  for (subset_name in names(definitions)) {
    for (adiposity_name in names(substitutions)) {
      adiposity <- substitutions[[adiposity_name]]
      data$.sensitivity_complete <- primary_complete_indicator(data) & definitions[[subset_name]] &
        is.finite(data[[adiposity]])
      design <- make_nhanes_design(data, "hormone") |> subset(.sensitivity_complete)
      for (exposure in c("androgen_score", "z_amh", "z_fai", "z_inverse_shbg")) {
        if (all(is.na(data[[exposure]]))) next
        formula <- model_formula(
          "nonfasting_score", exposure,
          c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking", adiposity)
        )
        model <- svyglm(formula, design = design)
        index <- index + 1L
        rows[[index]] <- tidy_svyglm(
          model, exposure, "sensitivity", paste(subset_name, adiposity_name, sep = ":"),
          exposure, "nonfasting_score"
        )
      }
    }
  }
  result <- bind_rows(rows)
  write_csv_safe(result, file.path(root, "outputs", "tables", "core_sensitivities.csv"))
  result
}

leave_one_psu_out <- function(data, root = project_root()) {
  data$primary_complete <- primary_complete_indicator(data)
  full_design <- make_nhanes_design(data, "hormone") |> subset(primary_complete)
  psus <- unique(data.frame(
    stratum = full_design$variables$SDMVSTRA,
    psu = full_design$variables$SDMVPSU
  ))
  rows <- lapply(seq_len(nrow(psus)), function(index) {
    keep <- !(full_design$variables$SDMVSTRA == psus$stratum[index] &
              full_design$variables$SDMVPSU == psus$psu[index])
    design <- full_design[keep, ]
    bind_rows(lapply(c("androgen_score", "z_amh"), function(exposure) {
      fit <- svyglm(
        model_formula(
          "nonfasting_score", exposure,
          c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking")
        ),
        design = design
      )
      tidy_svyglm(fit, exposure, "leave_one_psu_out",
                  paste(psus$stratum[index], psus$psu[index], sep = ":"),
                  exposure, "nonfasting_score")
    }))
  }) |> bind_rows()
  write_csv_safe(rows, file.path(root, "outputs", "tables", "leave_one_psu_out.csv"))
  rows
}

