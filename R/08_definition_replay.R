definition_replay_spec <- function() {
  tibble::tribble(
    ~definition, ~label, ~expression, ~rationale,
    "high_total_testosterone", "Total testosterone", "tt_high",
    "First-line concentration marker recommended for biochemical androgen assessment.",
    "high_fai", "Free androgen index", "fai_high",
    "Calculated free-testosterone proxy that incorporates total testosterone and SHBG.",
    "testosterone_or_fai", "Testosterone or FAI", "tt_high | fai_high",
    "Combined total- and free-testosterone approach used in guideline-based assessment.",
    "high_androstenedione", "Androstenedione", "a4_high",
    "Additional ovarian androgen marker used when testosterone measures are not elevated.",
    "high_dheas", "DHEAS", "dheas_high",
    "Adrenal androgen marker used as an additional biochemical measure.",
    "testosterone_or_androstenedione", "Testosterone or androstenedione", "tt_high | a4_high",
    "Concentration-based two-marker definition evaluated in LC-MS/MS studies.",
    "testosterone_fai_or_androstenedione", "Testosterone, FAI, or androstenedione",
    "tt_high | fai_high | a4_high",
    "Expanded definition combining first-line testosterone measures with androstenedione.",
    "any_four_marker", "Any of four markers", "tt_high | fai_high | a4_high | dheas_high",
    "Expanded biochemical definition that also includes an adrenal androgen marker."
  )
}

prepare_definition_replay <- function(common) {
  weights <- common$fasting_weight
  thresholds <- c(
    testosterone = weighted_quantile(common$c_testosterone, weights, 0.95),
    fai = weighted_quantile(common$c_fai, weights, 0.95),
    androstenedione = weighted_quantile(common$c_androstenedione, weights, 0.95),
    dheas = weighted_quantile(common$c_dheas, weights, 0.95)
  )
  common <- common |>
    dplyr::mutate(
      tt_high = c_testosterone >= thresholds[["testosterone"]],
      fai_high = c_fai >= thresholds[["fai"]],
      a4_high = c_androstenedione >= thresholds[["androstenedione"]],
      dheas_high = c_dheas >= thresholds[["dheas"]]
    )
  specification <- definition_replay_spec()
  for (index in seq_len(nrow(specification))) {
    common[[specification$definition[index]]] <- as.numeric(with(
      common,
      eval(parse(text = specification$expression[index]))
    ))
  }
  list(data = common, thresholds = thresholds, specification = specification)
}

weighted_binary_agreement <- function(x, y, weights) {
  ok <- !is.na(x) & !is.na(y) & positive_weight(weights)
  x <- as.logical(x[ok])
  y <- as.logical(y[ok])
  weights <- weights[ok]
  weights <- weights / sum(weights)
  observed <- sum(weights * (x == y))
  px <- sum(weights * x)
  py <- sum(weights * y)
  expected <- px * py + (1 - px) * (1 - py)
  kappa <- if (expected < 1) (observed - expected) / (1 - expected) else NA_real_
  intersection <- sum(weights * (x & y))
  union <- sum(weights * (x | y))
  tibble::tibble(
    agreement = observed,
    kappa = kappa,
    jaccard = if (union > 0) intersection / union else NA_real_
  )
}

run_definition_replay <- function(common, root = project_root()) {
  prepared <- prepare_definition_replay(common)
  data <- prepared$data
  specification <- prepared$specification
  design <- make_common_design(data)

  prevalence_rows <- lapply(seq_len(nrow(specification)), function(index) {
    definition <- specification$definition[index]
    estimate <- survey::svyciprop(
      stats::as.formula(paste0("~I(", definition, ")")),
      design,
      method = "logit",
      level = 0.95
    )
    interval <- stats::confint(estimate)
    tibble::tibble(
      definition = definition,
      label = specification$label[index],
      unweighted_n = sum(data[[definition]], na.rm = TRUE),
      total_n = sum(!is.na(data[[definition]])),
      weighted_proportion = as.numeric(stats::coef(estimate)),
      conf_low = as.numeric(interval[1]),
      conf_high = as.numeric(interval[2])
    )
  })
  prevalence <- dplyr::bind_rows(prevalence_rows)

  outcomes <- c(homa_ir = "c_homa_ir", fasting_score = "c_fasting_score")
  adjustment <- list(
    base = c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking"),
    adiposity = c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking", "c_adiposity")
  )
  association_rows <- list()
  row_index <- 0L
  for (stage in names(adjustment)) {
    for (definition in specification$definition) {
      for (outcome_name in names(outcomes)) {
        outcome <- outcomes[[outcome_name]]
        fit <- survey::svyglm(
          reformulate(c(definition, adjustment[[stage]]), response = outcome),
          design = design
        )
        row_index <- row_index + 1L
        association_rows[[row_index]] <- tidy_svyglm(
          fit, definition, "definition_replay", stage, definition, outcome_name
        )
      }
    }
  }
  associations <- dplyr::bind_rows(association_rows) |>
    dplyr::group_by(stage, outcome) |>
    dplyr::mutate(p_fdr = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::left_join(specification |> dplyr::select(definition, label), by = c("exposure" = "definition"))

  agreement_rows <- list()
  row_index <- 0L
  for (i in seq_len(nrow(specification))) {
    for (j in i:nrow(specification)) {
      row_index <- row_index + 1L
      agreement_rows[[row_index]] <- weighted_binary_agreement(
        data[[specification$definition[i]]],
        data[[specification$definition[j]]],
        data$fasting_weight
      ) |>
        dplyr::mutate(
          definition_1 = specification$definition[i],
          definition_2 = specification$definition[j],
          label_1 = specification$label[i],
          label_2 = specification$label[j],
          .before = 1
        )
    }
  }
  agreement <- dplyr::bind_rows(agreement_rows)

  threshold_table <- tibble::tibble(
    marker = names(prepared$thresholds),
    common_sample_weighted_95th_percentile = as.numeric(prepared$thresholds),
    scale = "Age-residualized and survey-weighted standardized value",
    note = paste(
      "Internal assay-specific reference limit used for structural replay;",
      "not a diagnostic cutoff and not a PCOS/PMOS prevalence estimate."
    )
  )
  definition_table <- specification |>
    dplyr::mutate(
      threshold_rule = "Marker-specific survey-weighted 95th percentile in the common NHANES sample",
      interpretation = "Biochemical androgen classification only; not a syndrome diagnosis"
    )

  write_csv_safe(threshold_table, file.path(root, "outputs", "tables", "definition_replay_thresholds.csv"))
  write_csv_safe(definition_table, file.path(root, "outputs", "tables", "definition_replay_specifications.csv"))
  write_csv_safe(prevalence, file.path(root, "outputs", "tables", "definition_replay_proportions.csv"))
  write_csv_safe(associations, file.path(root, "outputs", "tables", "definition_replay_associations.csv"))
  write_csv_safe(agreement, file.path(root, "outputs", "tables", "definition_replay_agreement.csv"))

  list(
    thresholds = threshold_table,
    specifications = definition_table,
    proportions = prevalence,
    associations = associations,
    agreement = agreement,
    data = data
  )
}
