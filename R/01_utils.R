quiet_library <- function() {
  suppressPackageStartupMessages({
    library(dplyr)
    library(haven)
    library(readr)
    library(survey)
    library(tibble)
    library(tidyr)
  })
  options(survey.lonely.psu = "adjust")
}

project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        file.exists(file.path(current, "run_all.R"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Run from the repository or a subdirectory.")
    current <- parent
  }
}

dir_create <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_config <- function(root = project_root()) {
  yaml::read_yaml(file.path(root, "config", "analysis.yml"))
}

write_csv_safe <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, path, na = "")
  invisible(path)
}

assert_columns <- function(data, columns, label = deparse(substitute(data))) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

positive_weight <- function(x) is.finite(x) & x > .Machine$double.xmin

row_mean_min <- function(data, columns, minimum = length(columns)) {
  values <- as.matrix(data[, columns, drop = FALSE])
  observed <- rowSums(is.finite(values))
  output <- rowMeans(values, na.rm = TRUE)
  output[observed < minimum] <- NA_real_
  output
}

weighted_center_scale <- function(x, weights) {
  ok <- is.finite(x) & positive_weight(weights)
  if (sum(ok) < 2L) return(rep(NA_real_, length(x)))
  center <- weighted.mean(x[ok], weights[ok])
  spread <- sqrt(weighted.mean((x[ok] - center)^2, weights[ok]))
  if (!is.finite(spread) || spread == 0) stop("Cannot standardize a constant variable.")
  (x - center) / spread
}

weighted_age_residual_z <- function(x, age, weights, transform = log) {
  transformed <- rep(NA_real_, length(x))
  valid_x <- is.finite(x) & x > 0
  transformed[valid_x] <- transform(x[valid_x])
  ok <- is.finite(transformed) & is.finite(age) & positive_weight(weights)
  fit <- lm(transformed ~ age, weights = weights, subset = ok, na.action = na.exclude)
  residual <- rep(NA_real_, length(x))
  residual[which(ok)] <- residuals(fit)
  weighted_center_scale(residual, weights)
}

weighted_cor <- function(data, variables, weights) {
  values <- as.matrix(data[, variables, drop = FALSE])
  ok <- complete.cases(values) & positive_weight(weights)
  values <- values[ok, , drop = FALSE]
  w <- weights[ok] / sum(weights[ok])
  centers <- colSums(values * w)
  centered <- sweep(values, 2, centers, "-")
  covariance <- crossprod(centered * sqrt(w), centered * sqrt(w))
  scales <- sqrt(diag(covariance))
  correlation <- covariance / outer(scales, scales)
  dimnames(correlation) <- list(variables, variables)
  correlation
}

weighted_quantile <- function(x, weights, probs) {
  ok <- is.finite(x) & positive_weight(weights)
  x <- x[ok]
  weights <- weights[ok]
  order_index <- order(x)
  x <- x[order_index]
  cumulative <- cumsum(weights[order_index]) / sum(weights)
  vapply(probs, function(p) x[which(cumulative >= p)[1]], numeric(1))
}

weighted_r_squared <- function(model) {
  response <- model$y
  if (is.null(response)) response <- model$model[[1]]
  fitted <- fitted(model)
  weights <- weights(model$survey.design, "sampling")
  used <- as.numeric(rownames(model$model))
  if (length(weights) != length(response)) weights <- weights[used]
  mean_y <- weighted.mean(response, weights)
  1 - sum(weights * (response - fitted)^2) / sum(weights * (response - mean_y)^2)
}

tidy_svyglm <- function(model, term, analysis, stage, exposure, outcome) {
  coefficient <- coef(model)[term]
  interval <- confint(model)[term, ]
  p_value <- summary(model)$coefficients[term, "Pr(>|t|)"]
  tibble(
    analysis = analysis,
    stage = stage,
    exposure = exposure,
    outcome = outcome,
    term = term,
    n = stats::nobs(model),
    design_df = survey::degf(model$survey.design),
    estimate = unname(coefficient),
    conf_low = unname(interval[1]),
    conf_high = unname(interval[2]),
    p_value = unname(p_value)
  )
}

tidy_hc3 <- function(model, term, analysis, stage, exposure, outcome) {
  test <- lmtest::coeftest(model, vcov. = sandwich::vcovHC(model, type = "HC3"))
  estimate <- test[term, 1]
  standard_error <- test[term, 2]
  critical <- qnorm(0.975)
  tibble(
    analysis = analysis,
    stage = stage,
    exposure = exposure,
    outcome = outcome,
    term = term,
    n = stats::nobs(model),
    estimate = unname(estimate),
    conf_low = unname(estimate - critical * standard_error),
    conf_high = unname(estimate + critical * standard_error),
    p_value = unname(test[term, 4])
  )
}

model_formula <- function(outcome, exposure, covariates = character(), interaction = NULL) {
  terms <- c(exposure, covariates)
  if (!is.null(interaction)) terms <- c(terms, paste0(exposure, "*", interaction))
  reformulate(unique(terms), response = outcome)
}

with_seed <- function(seed, expression) {
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

