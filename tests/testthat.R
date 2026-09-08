source(file.path("R", "00_dependencies.R"))
assert_dependencies(include_tests = TRUE)
for (script in list.files("R", pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE)) {
  if (!grepl("00_dependencies\\.R$", script)) source(script)
}
quiet_library()
testthat::test_dir("tests/testthat", reporter = "summary")

