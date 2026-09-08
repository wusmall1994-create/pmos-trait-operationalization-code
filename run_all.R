root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the repository root.")
}

source(file.path(root, "R", "00_dependencies.R"))
assert_dependencies()
for (script in list.files(file.path(root, "R"), pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE)) {
  if (!grepl("00_dependencies\\.R$", script)) source(script)
}
quiet_library()

quick <- identical(tolower(Sys.getenv("PMOS_QUICK", unset = "false")), "true")
download_nhanes(root)

message("Preparing NHANES August 2021-August 2023...")
discovery <- prepare_nhanes_period(root, "2021-2023")
primary <- run_primary_models(discovery, root)
fasting <- run_fasting_primary(discovery, root)
sensitivities <- run_primary_sensitivities(discovery, root)
influence <- leave_one_psu_out(discovery, root)

message("Running common-sample construct analyses...")
common <- prepare_common_sample(discovery)
construct_matrix <- run_operationalization_matrix(common, root)
incremental_fit <- run_incremental_fit(common, root)
pca <- run_weighted_pca(common, root, quick = quick)
discordance <- run_discordance_benchmarks(discovery, root, quick = quick)

message("Preparing NHANES 2017-March 2020...")
replication_data <- prepare_nhanes_period(root, "2017-2020")
replication <- run_temporal_replication(replication_data, root)
swan <- run_swan_models(root)
make_code_release_figures(root)

session <- capture.output(sessionInfo())
dir.create(file.path(root, "outputs", "audit"), recursive = TRUE, showWarnings = FALSE)
writeLines(session, file.path(root, "outputs", "audit", "sessionInfo.txt"), useBytes = TRUE)
message("Analysis complete. Generated files are under outputs/ (ignored by Git).")
