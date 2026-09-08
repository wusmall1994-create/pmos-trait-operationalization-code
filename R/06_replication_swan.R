run_temporal_replication <- function(data, root = project_root()) {
  exposures <- c(
    androstenedione = "z_androstenedione",
    amh = "z_amh",
    shbg = "z_shbg"
  )
  data$replication_complete <- data$strict_eligible &
    (is.na(data$hormone_or_fertility_drug) | !data$hormone_or_fertility_drug) &
    complete.cases(data[, c(
      unname(exposures), "nonfasting_score", "RIDAGEYR", "RIDRETH3", "INDFMPIR",
      "current_smoking", "waist_height_ratio"
    )])
  design <- make_nhanes_design(data, "hormone") |> subset(replication_complete)
  result <- bind_rows(lapply(names(exposures), function(label) {
    exposure <- exposures[[label]]
    fit <- svyglm(
      model_formula(
        "nonfasting_score", exposure,
        c("RIDAGEYR", "factor(RIDRETH3)", "INDFMPIR", "current_smoking")
      ),
      design = design
    )
    tidy_svyglm(fit, exposure, "temporal_replication", "social_behavioral", label, "nonfasting_score")
  })) |>
    mutate(p_holm = p.adjust(p_value, method = "holm"))
  write_csv_safe(result, file.path(root, "outputs", "tables", "temporal_replication_models.csv"))
  result
}

read_tabular_file <- function(path) {
  extension <- tolower(tools::file_ext(path))
  switch(
    extension,
    dta = haven::read_dta(path),
    sav = haven::read_sav(path),
    sas7bdat = haven::read_sas(path),
    xpt = haven::read_xpt(path),
    rds = readRDS(path),
    csv = readr::read_csv(path, show_col_types = FALSE),
    stop("Unsupported SWAN file extension: ", extension)
  )
}

read_swan_mapping <- function(root = project_root()) {
  mapping <- readr::read_csv(
    file.path(root, "config", "swan_variable_map.csv"), show_col_types = FALSE
  )
  required_rows <- mapping$required == "yes"
  if (any(is.na(mapping$source[required_rows]) | mapping$source[required_rows] == "")) {
    stop("Complete the required 'source' cells in config/swan_variable_map.csv.")
  }
  mapping
}

prepare_swan <- function(path, mapping) {
  raw <- read_tabular_file(path)
  used <- mapping |> filter(!is.na(source), source != "")
  assert_columns(raw, used$source, "SWAN baseline file")
  output <- raw[, used$source, drop = FALSE]
  names(output) <- used$canonical
  output <- as.data.frame(output)

  eligible <- output$menopause_stage %in% c("premenopausal", "early perimenopausal", 1, 2) &
    output$reproductive_hormone_use %in% c(FALSE, 0, 2) &
    output$systemic_steroid_use %in% c(FALSE, 0, 2)
  output <- output[!is.na(eligible) & eligible, , drop = FALSE]
  output$waist_height_ratio <- output$waist_cm / output$height_cm
  output$homa_ir <- output$glucose * output$insulin / 405
  output$map <- (output$systolic_bp + 2 * output$diastolic_bp) / 3
  output$fai <- 100 * (output$testosterone * 0.0347) / output$shbg

  z <- function(x, log_transform = FALSE) {
    if (log_transform) x <- log(x)
    as.numeric(scale(x))
  }
  output$z_testosterone <- z(output$testosterone, TRUE)
  output$z_dheas <- z(output$dheas, TRUE)
  output$z_shbg <- z(output$shbg, TRUE)
  output$z_fai <- z(output$fai, TRUE)
  output$z_homa_ir <- z(output$homa_ir, TRUE)
  output$z_triglycerides <- z(output$triglycerides, TRUE)
  output$z_inverse_hdl <- -z(output$hdl)
  output$z_map <- z(output$map)
  metabolic_columns <- c("z_homa_ir", "z_triglycerides", "z_inverse_hdl", "z_map")
  if ("hba1c" %in% names(output)) {
    output$z_hba1c <- z(output$hba1c)
    metabolic_columns <- c(metabolic_columns, "z_hba1c")
  }
  output$fasting_score <- as.numeric(scale(row_mean_min(output, metabolic_columns, length(metabolic_columns))))
  output$androgen_composite <- rowMeans(output[, c("z_testosterone", "z_dheas")])
  output
}

run_swan_models <- function(root = project_root()) {
  path <- Sys.getenv("SWAN_FILE", unset = "")
  if (!nzchar(path)) {
    message("SWAN_FILE is not set; restricted-data analysis was skipped.")
    return(NULL)
  }
  if (!file.exists(path)) stop("SWAN_FILE does not exist: ", path)
  mapping <- read_swan_mapping(root)
  data <- prepare_swan(path, mapping)
  exposures <- c(
    testosterone = "z_testosterone", dheas = "z_dheas",
    androgen_composite = "androgen_composite", shbg = "z_shbg", fai = "z_fai"
  )
  stages <- list(
    base = c("age", "factor(race_ethnicity)", "factor(menopause_stage)", "current_smoking"),
    adiposity = c("age", "factor(race_ethnicity)", "factor(menopause_stage)",
                  "current_smoking", "waist_height_ratio")
  )
  rows <- list()
  index <- 0L
  for (stage in names(stages)) {
    for (label in names(exposures)) {
      exposure <- exposures[[label]]
      fit <- lm(model_formula("fasting_score", exposure, stages[[stage]]), data = data)
      index <- index + 1L
      rows[[index]] <- tidy_hc3(
        fit, exposure, "swan", stage, label, "fasting_score"
      )
    }
  }
  result <- bind_rows(rows) |> mutate(p_fdr = p.adjust(p_value, method = "BH"))
  write_csv_safe(result, file.path(root, "outputs", "tables", "swan_baseline_models.csv"))
  result
}

