suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) stop("Run this script from the repository root.")
input_dir <- Sys.getenv("SWAN_DIR", unset = file.path(root, "data", "restricted", "swan"))
output_dir <- file.path(root, "outputs", "swan_longitudinal_assessment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

study_visit <- c(
  `28762` = 0L, `29221` = 1L, `29401` = 2L, `29701` = 3L,
  `30142` = 4L, `30501` = 5L, `31181` = 6L, `31901` = 7L,
  `32122` = 8L, `32721` = 9L, `32961` = 10L
)

positive_value <- function(x) {
  x <- zap_missing(x)
  is.finite(x) & x > 0
}

get_var <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) return(rep(NA_real_, nrow(data)))
  as.numeric(zap_missing(data[[hit[[1]]]]))
}

extract_visit <- function(file) {
  study <- sub("-.*", "", basename(file))
  visit <- unname(study_visit[[study]])
  suffix <- as.character(visit)
  data <- read_dta(file)

  id <- as.character(zap_missing(data$SWANID))
  testosterone <- get_var(data, paste0("T", suffix))
  dheas <- get_var(data, paste0("DHAS", suffix))
  shbg <- get_var(data, paste0("SHBG", suffix))
  glucose <- get_var(data, c(paste0("GLUCRES", suffix), paste0("GLUCRE", suffix)))
  insulin <- get_var(data, paste0("INSURES", suffix))
  triglycerides <- get_var(data, paste0("TRIGRES", suffix))
  hdl <- get_var(data, paste0("HDLRESU", suffix))
  bmi <- get_var(data, paste0("BMI", suffix))
  waist <- get_var(data, paste0("WAIST", suffix))

  tibble(
    study = study,
    visit = visit,
    swanid = id,
    testosterone = testosterone,
    dheas = dheas,
    shbg = shbg,
    glucose = glucose,
    insulin = insulin,
    triglycerides = triglycerides,
    hdl = hdl,
    bmi = bmi,
    waist = waist
  ) %>%
    mutate(
      hormone3_complete = positive_value(testosterone) & positive_value(dheas) & positive_value(shbg),
      homa_complete = positive_value(glucose) & positive_value(insulin),
      metabolic_lab_complete = homa_complete & positive_value(triglycerides) & positive_value(hdl),
      hormone_homa_complete = hormone3_complete & homa_complete,
      hormone_metabolic_complete = hormone3_complete & metabolic_lab_complete,
      hormone_homa_bmi_complete = hormone_homa_complete & positive_value(bmi),
      hormone_homa_waist_complete = hormone_homa_complete & positive_value(waist)
    )
}

files <- Sys.glob(file.path(input_dir, "*-Data.dta"))
if (!length(files)) stop("No SWAN visit files found. Set SWAN_DIR to the local directory containing *-Data.dta files.")
long <- map_dfr(files, extract_visit) %>% arrange(visit, swanid)

visit_summary <- long %>%
  group_by(study, visit) %>%
  summarise(
    records = n(),
    unique_ids = n_distinct(swanid),
    duplicate_ids = records - unique_ids,
    testosterone_valid = sum(positive_value(testosterone)),
    dheas_valid = sum(positive_value(dheas)),
    shbg_valid = sum(positive_value(shbg)),
    hormone3_complete = sum(hormone3_complete),
    glucose_valid = sum(positive_value(glucose)),
    insulin_valid = sum(positive_value(insulin)),
    triglycerides_valid = sum(positive_value(triglycerides)),
    hdl_valid = sum(positive_value(hdl)),
    homa_complete = sum(homa_complete),
    metabolic_lab_complete = sum(metabolic_lab_complete),
    hormone_homa_complete = sum(hormone_homa_complete),
    hormone_metabolic_complete = sum(hormone_metabolic_complete),
    hormone_homa_bmi_complete = sum(hormone_homa_bmi_complete),
    hormone_homa_waist_complete = sum(hormone_homa_waist_complete),
    .groups = "drop"
  ) %>% arrange(visit)

participant_summary <- long %>%
  group_by(swanid) %>%
  summarise(
    visits_observed = n(),
    hormone3_visits = sum(hormone3_complete),
    hormone_homa_visits = sum(hormone_homa_complete),
    hormone_metabolic_visits = sum(hormone_metabolic_complete),
    hormone_homa_bmi_visits = sum(hormone_homa_bmi_complete),
    hormone_homa_waist_visits = sum(hormone_homa_waist_complete),
    .groups = "drop"
  )

repeat_summary <- tibble(
  criterion = c("hormone3", "hormone_homa", "hormone_metabolic", "hormone_homa_bmi", "hormone_homa_waist"),
  n_at_least_1 = c(
    sum(participant_summary$hormone3_visits >= 1),
    sum(participant_summary$hormone_homa_visits >= 1),
    sum(participant_summary$hormone_metabolic_visits >= 1),
    sum(participant_summary$hormone_homa_bmi_visits >= 1),
    sum(participant_summary$hormone_homa_waist_visits >= 1)
  ),
  n_at_least_2 = c(
    sum(participant_summary$hormone3_visits >= 2),
    sum(participant_summary$hormone_homa_visits >= 2),
    sum(participant_summary$hormone_metabolic_visits >= 2),
    sum(participant_summary$hormone_homa_bmi_visits >= 2),
    sum(participant_summary$hormone_homa_waist_visits >= 2)
  ),
  n_at_least_3 = c(
    sum(participant_summary$hormone3_visits >= 3),
    sum(participant_summary$hormone_homa_visits >= 3),
    sum(participant_summary$hormone_metabolic_visits >= 3),
    sum(participant_summary$hormone_homa_bmi_visits >= 3),
    sum(participant_summary$hormone_homa_waist_visits >= 3)
  ),
  n_at_least_4 = c(
    sum(participant_summary$hormone3_visits >= 4),
    sum(participant_summary$hormone_homa_visits >= 4),
    sum(participant_summary$hormone_metabolic_visits >= 4),
    sum(participant_summary$hormone_homa_bmi_visits >= 4),
    sum(participant_summary$hormone_homa_waist_visits >= 4)
  ),
  max_visits = c(
    max(participant_summary$hormone3_visits),
    max(participant_summary$hormone_homa_visits),
    max(participant_summary$hormone_metabolic_visits),
    max(participant_summary$hormone_homa_bmi_visits),
    max(participant_summary$hormone_homa_waist_visits)
  )
)

write_csv(visit_summary, file.path(output_dir, "visit_variable_completeness.csv"))
write_csv(participant_summary, file.path(output_dir, "participant_repeat_counts.csv"))
write_csv(repeat_summary, file.path(output_dir, "repeat_measurement_eligibility.csv"))

print(visit_summary, n = Inf, width = Inf)
print(repeat_summary, n = Inf, width = Inf)

