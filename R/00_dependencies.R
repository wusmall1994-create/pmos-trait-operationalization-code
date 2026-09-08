required_packages <- c(
  "digest", "dplyr", "ggplot2", "haven", "jsonlite", "lmtest", "readr",
  "sandwich", "survey", "tibble", "tidyr", "yaml"
)

install_dependencies <- function(repos = "https://cloud.r-project.org") {
  missing <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    install.packages(missing, repos = repos)
  }
  invisible(required_packages)
}

assert_dependencies <- function(include_tests = FALSE) {
  packages <- c(required_packages, if (include_tests) "testthat")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Run source('R/00_dependencies.R'); install_dependencies().")
  }
  invisible(TRUE)
}
