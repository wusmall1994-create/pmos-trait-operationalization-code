nhanes_file_spec <- function() {
  tribble(
    ~period, ~stem, ~role,
    "2021-2023", "DEMO_L", "demographics",
    "2021-2023", "RHQ_L", "reproductive_health",
    "2021-2023", "BMX_L", "anthropometry",
    "2021-2023", "BPXO_L", "blood_pressure",
    "2021-2023", "TST_L", "sex_steroids",
    "2021-2023", "GHB_L", "glycohemoglobin",
    "2021-2023", "HDL_L", "hdl",
    "2021-2023", "TRIGLY_L", "triglycerides",
    "2021-2023", "GLU_L", "glucose",
    "2021-2023", "INS_L", "insulin",
    "2021-2023", "SMQ_L", "smoking",
    "2021-2023", "DIQ_L", "diabetes",
    "2021-2023", "RXQ_RX_L", "prescriptions",
    "2017-2020", "P_DEMO", "demographics",
    "2017-2020", "P_RHQ", "reproductive_health",
    "2017-2020", "P_BMX", "anthropometry",
    "2017-2020", "P_BPXO", "blood_pressure",
    "2017-2020", "P_TST", "sex_steroids",
    "2017-2020", "P_GHB", "glycohemoglobin",
    "2017-2020", "P_HDL", "hdl",
    "2017-2020", "P_TRIGLY", "triglycerides",
    "2017-2020", "P_GLU", "glucose",
    "2017-2020", "P_INS", "insulin",
    "2017-2020", "P_SMQ", "smoking",
    "2017-2020", "P_DIQ", "diabetes",
    "2017-2020", "P_RXQ_RX", "prescriptions"
  ) |>
    mutate(
      path_year = if_else(period == "2021-2023", "2021", "2017"),
      url = paste0(
        "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", path_year,
        "/DataFiles/", stem, ".xpt"
      )
    )
}

download_nhanes <- function(root = project_root(), overwrite = FALSE) {
  destination <- dir_create(root, "data", "public", "nhanes")
  specification <- nhanes_file_spec()
  records <- vector("list", nrow(specification))

  for (index in seq_len(nrow(specification))) {
    item <- specification[index, ]
    output <- file.path(destination, paste0(item$stem, ".xpt"))
    if (!file.exists(output) || overwrite) {
      message("Downloading ", item$stem)
      download.file(item$url, output, mode = "wb", quiet = TRUE)
    }
    if (file.info(output)$size <= 0) stop("Empty download: ", output)
    records[[index]] <- tibble(
      period = item$period,
      role = item$role,
      file = basename(output),
      source_url = item$url,
      bytes = file.info(output)$size,
      sha256 = digest::digest(output, algo = "sha256", file = TRUE)
    )
  }

  manifest <- bind_rows(records)
  write_csv_safe(manifest, file.path(root, "outputs", "audit", "nhanes_download_manifest.csv"))
  invisible(manifest)
}

read_nhanes_component <- function(root, stem) {
  path <- file.path(root, "data", "public", "nhanes", paste0(stem, ".xpt"))
  if (!file.exists(path)) stop("Missing NHANES component: ", path)
  haven::read_xpt(path) |> mutate(SEQN = as.numeric(SEQN))
}

select_existing <- function(data, columns) {
  data |> select(any_of(unique(c("SEQN", columns))))
}

aggregate_prescriptions <- function(data, period) {
  if (period == "2021-2023") {
    return(data |>
      transmute(
        SEQN,
        prescription_any = RXQ033 == 1,
        prescription_count = if ("RXQ050" %in% names(data)) RXQ050 else NA_real_,
        hormone_or_fertility_drug = NA
      ))
  }

  drug_name <- if ("RXDDRUG" %in% names(data)) tolower(as.character(data$RXDDRUG)) else rep("", nrow(data))
  hormone_pattern <- paste(
    c("estradiol", "estrogen", "progesterone", "progestin", "testosterone",
      "etonogestrel", "levonorgestrel", "norgestrel", "norgestimate", "norethindrone",
      "drospirenone", "desogestrel", "gestodene", "ethynodiol", "norelgestromin",
      "spironolactone",
      "clomiphene", "letrozole", "gonadotropin", "human chorionic", "follitropin",
      "menotropin", "leuprolide", "ganirelix", "cetrorelix", "fertility"),
    collapse = "|"
  )
  data |>
    mutate(.hormone_or_fertility = grepl(hormone_pattern, drug_name)) |>
    group_by(SEQN) |>
    summarise(
      prescription_any = any(RXDUSE == 1, na.rm = TRUE),
      prescription_count = suppressWarnings(max(RXDCOUNT, na.rm = TRUE)),
      hormone_or_fertility_drug = any(.hormone_or_fertility, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(prescription_count = if_else(is.infinite(prescription_count), NA_real_, prescription_count))
}

merge_nhanes_period <- function(root = project_root(), period = c("2021-2023", "2017-2020")) {
  period <- match.arg(period)
  stems <- nhanes_file_spec() |> filter(.data$period == .env$period)
  components <- setNames(lapply(stems$stem, function(x) read_nhanes_component(root, x)), stems$role)

  medication <- aggregate_prescriptions(components$prescriptions, period)
  components$prescriptions <- medication

  keep <- list(
    demographics = c("RIAGENDR", "RIDAGEYR", "RIDRETH3", "RIDEXPRG", "RIDSTATR", "INDFMPIR",
                     "SDMVSTRA", "SDMVPSU", "WTMEC2YR", "WTMECPRP"),
    reproductive_health = c("RHD143", "RHQ200", "RHD280", "RHQ305", "RHQ332"),
    anthropometry = c("BMXHT", "BMXBMI", "BMXWAIST"),
    blood_pressure = c("BPXOSY1", "BPXOSY2", "BPXOSY3", "BPXODI1", "BPXODI2", "BPXODI3"),
    sex_steroids = c("WTPH2YR", "WTTSTPP", "LBXTST", "LBDTSTSI", "LBXAND", "LBXAMH", "LBXDHE", "LBXSHBG"),
    glycohemoglobin = c("WTPH2YR", "LBXGH"),
    hdl = c("WTPH2YR", "LBDHDD"),
    triglycerides = c("WTSAF2YR", "WTSAFPRP", "LBXTLG", "LBXTR"),
    glucose = c("WTSAF2YR", "WTSAFPRP", "LBXGLU"),
    insulin = c("WTSAF2YR", "WTSAFPRP", "LBXIN", "LBDINSI"),
    smoking = c("SMQ020", "SMQ040"),
    diabetes = c("DIQ010"),
    prescriptions = c("prescription_any", "prescription_count", "hormone_or_fertility_drug")
  )

  selected <- Map(select_existing, components, keep[names(components)])
  merged <- Reduce(function(x, y) full_join(x, y, by = "SEQN"), selected)

  duplicated_names <- names(merged)[duplicated(names(merged))]
  if (length(duplicated_names)) stop("Duplicated columns after merge: ", paste(duplicated_names, collapse = ", "))
  merged
}

coalesce_named <- function(data, candidates) {
  available <- candidates[candidates %in% names(data)]
  if (!length(available)) return(rep(NA_real_, nrow(data)))
  Reduce(dplyr::coalesce, lapply(available, function(x) data[[x]]))
}

derive_nhanes <- function(data, period) {
  weights_hormone <- if (period == "2021-2023") {
    coalesce_named(data, c("WTPH2YR", "WTPH2YR.x", "WTPH2YR.y"))
  } else data$WTTSTPP
  weights_fasting <- if (period == "2021-2023") {
    coalesce_named(data, c("WTSAF2YR", "WTSAF2YR.x", "WTSAF2YR.y"))
  } else coalesce_named(data, c("WTSAFPRP", "WTSAFPRP.x", "WTSAFPRP.y"))
  insulin <- if (period == "2021-2023") data$LBXIN else data$LBDINSI / 6.0
  triglycerides <- if (period == "2021-2023") data$LBXTLG else data$LBXTR
  testosterone_si <- coalesce_named(data, "LBDTSTSI")

  systolic_columns <- intersect(c("BPXOSY1", "BPXOSY2", "BPXOSY3"), names(data))
  diastolic_columns <- intersect(c("BPXODI1", "BPXODI2", "BPXODI3"), names(data))
  data$systolic_bp <- row_mean_min(data, systolic_columns, 1)
  data$diastolic_bp <- row_mean_min(data, diastolic_columns, 1)
  period_medication_eligible <- if (period == "2017-2020") {
    is.na(data$hormone_or_fertility_drug) | !data$hormone_or_fertility_drug
  } else rep(TRUE, nrow(data))

  data |>
    mutate(
      period = period,
      hormone_weight = weights_hormone,
      fasting_weight = weights_fasting,
      female = RIAGENDR == 2,
      age_eligible = between(RIDAGEYR, 20, 44),
      mec_examined = RIDSTATR == 2,
      not_pregnant = RIDEXPRG == 2 & (is.na(RHD143) | RHD143 == 2),
      not_breastfeeding = is.na(RHQ200) | RHQ200 != 1,
      no_hysterectomy = RHD280 == 2,
      no_bilateral_oophorectomy = RHQ305 == 2,
      strict_eligible = female & age_eligible & mec_examined & not_pregnant & not_breastfeeding &
        no_hysterectomy & no_bilateral_oophorectomy & period_medication_eligible,
      current_smoking = case_when(
        SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ 1,
        SMQ020 %in% c(1, 2) ~ 0,
        TRUE ~ NA_real_
      ),
      mean_arterial_pressure = (systolic_bp + 2 * diastolic_bp) / 3,
      waist_height_ratio = BMXWAIST / BMXHT,
      insulin_uuml = insulin,
      triglycerides_mgdl = triglycerides,
      homa_ir = LBXGLU * insulin_uuml / 405,
      free_androgen_index = if_else(
        is.finite(testosterone_si) & is.finite(LBXSHBG) & LBXSHBG > 0,
        100 * testosterone_si / LBXSHBG,
        NA_real_
      ),
      diabetes = DIQ010 == 1
    )
}

prepare_nhanes_period <- function(root = project_root(), period = c("2021-2023", "2017-2020")) {
  period <- match.arg(period)
  data <- merge_nhanes_period(root, period) |> derive_nhanes(period)
  domain <- !is.na(data$strict_eligible) & data$strict_eligible & positive_weight(data$hormone_weight)
  domain_weight <- ifelse(domain, data$hormone_weight, 0)

  hormone_map <- c(
    testosterone = "LBXTST",
    androstenedione = "LBXAND",
    dheas = "LBXDHE",
    amh = "LBXAMH",
    shbg = "LBXSHBG",
    fai = "free_androgen_index"
  )
  for (canonical in names(hormone_map)) {
    source <- hormone_map[[canonical]]
    if (source %in% names(data) && any(is.finite(data[[source]]) & data[[source]] > 0 & domain, na.rm = TRUE)) {
      data[[paste0("z_", canonical)]] <- weighted_age_residual_z(
        data[[source]], data$RIDAGEYR, domain_weight
      )
    } else data[[paste0("z_", canonical)]] <- NA_real_
  }

  data$androgen_score <- row_mean_min(
    data, c("z_testosterone", "z_androstenedione", "z_dheas"), minimum = 2
  )
  data$z_inverse_shbg <- -data$z_shbg
  metabolic_sources <- c(
    hba1c = "LBXGH",
    hdl = "LBDHDD",
    map = "mean_arterial_pressure"
  )
  for (canonical in names(metabolic_sources)) {
    source <- metabolic_sources[[canonical]]
    data[[paste0("z_", canonical)]] <- weighted_center_scale(data[[source]], domain_weight)
  }
  data$z_inverse_hdl <- -data$z_hdl
  data$nonfasting_score <- row_mean_min(data, c("z_hba1c", "z_inverse_hdl", "z_map"), 3)
  data$nonfasting_score <- weighted_center_scale(data$nonfasting_score, domain_weight)

  fasting_domain <- !is.na(data$strict_eligible) & data$strict_eligible & positive_weight(data$fasting_weight)
  fasting_weight <- ifelse(fasting_domain, data$fasting_weight, 0)
  data$z_homa_ir <- weighted_center_scale(log(data$homa_ir), fasting_weight)
  data$z_insulin <- weighted_center_scale(log(data$insulin_uuml), fasting_weight)
  data$z_glucose <- weighted_center_scale(data$LBXGLU, fasting_weight)
  data$z_triglycerides <- weighted_center_scale(log(data$triglycerides_mgdl), fasting_weight)
  data$fasting_score <- row_mean_min(
    data, c("z_hba1c", "z_inverse_hdl", "z_map", "z_homa_ir", "z_triglycerides"), 5
  )
  data$fasting_score <- weighted_center_scale(data$fasting_score, fasting_weight)
  data
}

make_nhanes_design <- function(data, weight = c("hormone", "fasting")) {
  weight <- match.arg(weight)
  column <- if (weight == "hormone") "hormone_weight" else "fasting_weight"
  sampled <- positive_weight(data[[column]]) & is.finite(data$SDMVSTRA) & is.finite(data$SDMVPSU)
  survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = as.formula(paste0("~", column)),
    nest = TRUE,
    data = data[sampled, , drop = FALSE]
  )
}
