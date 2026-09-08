pca_from_weighted_cor <- function(data, variables, weights) {
  correlation <- weighted_cor(data, variables, weights)
  decomposition <- eigen(correlation, symmetric = TRUE)
  list(
    correlation = correlation,
    values = decomposition$values,
    vectors = decomposition$vectors,
    explained = decomposition$values / sum(decomposition$values)
  )
}

permute_within_psu <- function(values, psu) {
  output <- values
  for (group in unique(psu)) {
    index <- which(psu == group)
    output[index] <- sample(values[index], length(index), replace = FALSE)
  }
  output
}

run_parallel_analysis <- function(common, variables, replicates, seed) {
  observed <- pca_from_weighted_cor(common, variables, common$fasting_weight)
  null_values <- with_seed(seed, {
    replicate(replicates, {
      permuted <- common
      for (variable in variables) {
        permuted[[variable]] <- permute_within_psu(permuted[[variable]], permuted$SDMVPSU)
      }
      pca_from_weighted_cor(permuted, variables, permuted$fasting_weight)$values
    })
  })
  thresholds <- apply(null_values, 1, quantile, probs = 0.95, na.rm = TRUE)
  list(observed = observed, thresholds = thresholds, null_values = null_values)
}

bootstrap_psu_weights <- function(strata, psu, base_weights) {
  multiplier <- numeric(length(base_weights))
  for (stratum in unique(strata)) {
    stratum_index <- which(strata == stratum)
    psus <- unique(psu[stratum_index])
    sampled <- sample(psus, length(psus), replace = TRUE)
    counts <- table(sampled)
    for (cluster in names(counts)) {
      multiplier[stratum_index[psu[stratum_index] == type.convert(cluster, as.is = TRUE)]] <- counts[[cluster]]
    }
  }
  base_weights * multiplier
}

align_two_axes <- function(reference, candidate) {
  direct <- sum(abs(diag(crossprod(reference[, 1:2, drop = FALSE], candidate[, 1:2, drop = FALSE]))))
  swapped_candidate <- candidate[, c(2, 1), drop = FALSE]
  swapped <- sum(abs(diag(crossprod(reference[, 1:2, drop = FALSE], swapped_candidate))))
  aligned <- if (swapped > direct) swapped_candidate else candidate[, 1:2, drop = FALSE]
  for (axis in 1:2) {
    if (sum(reference[, axis] * aligned[, axis]) < 0) aligned[, axis] <- -aligned[, axis]
  }
  aligned
}

minimum_canonical_correlation <- function(reference, candidate) {
  reference_q <- qr.Q(qr(reference[, 1:2, drop = FALSE]))
  candidate_q <- qr.Q(qr(candidate[, 1:2, drop = FALSE]))
  min(svd(crossprod(reference_q, candidate_q))$d)
}

run_weighted_pca <- function(common, root = project_root(), config = read_config(root), quick = FALSE) {
  variables <- c(
    "c_testosterone", "c_androstenedione", "c_dheas", "c_amh", "c_inverse_shbg",
    "c_homa_ir", "c_triglycerides", "c_hba1c", "c_inverse_hdl", "c_map", "c_adiposity"
  )
  parallel_n <- if (quick) 20L else as.integer(config$parallel_permutations)
  bootstrap_n <- if (quick) 30L else as.integer(config$pca_bootstrap_replicates)
  parallel <- run_parallel_analysis(common, variables, parallel_n, config$seed_parallel)
  observed <- parallel$observed

  variance <- tibble(
    component = seq_along(observed$values), eigenvalue = observed$values,
    explained_variance = observed$explained,
    cumulative_variance = cumsum(observed$explained),
    parallel_threshold_95 = parallel$thresholds,
    retained = observed$values > parallel$thresholds
  )
  loadings <- as.data.frame(observed$vectors) |>
    mutate(variable = variables, .before = 1) |>
    pivot_longer(-variable, names_to = "component", values_to = "loading") |>
    mutate(component = as.integer(sub("V", "", component)))

  bootstrap <- with_seed(config$seed_bootstrap, lapply(seq_len(bootstrap_n), function(replicate_id) {
    replicate_weights <- bootstrap_psu_weights(
      common$SDMVSTRA, common$SDMVPSU, common$fasting_weight
    )
    if (sum(replicate_weights > 0) < length(variables) + 1L) return(NULL)
    pca <- pca_from_weighted_cor(common, variables, replicate_weights)
    axes <- align_two_axes(observed$vectors, pca$vectors)
    tibble(
      replicate = replicate_id,
      pc1_explained = pca$explained[1],
      pc2_explained = pca$explained[2],
      cumulative_pc12 = sum(pca$explained[1:2]),
      pc1_retained = pca$values[1] > parallel$thresholds[1],
      pc2_retained = pca$values[2] > parallel$thresholds[2],
      pc3_retained = pca$values[3] > parallel$thresholds[3],
      minimum_canonical_correlation = minimum_canonical_correlation(observed$vectors, pca$vectors),
      aligned_loadings = list(axes)
    )
  })) |> bind_rows()

  loading_bootstrap <- bind_rows(lapply(seq_len(nrow(bootstrap)), function(index) {
    as.data.frame(bootstrap$aligned_loadings[[index]]) |>
      mutate(variable = variables, replicate = bootstrap$replicate[index], .before = 1) |>
      rename(pc1 = V1, pc2 = V2) |>
      pivot_longer(c(pc1, pc2), names_to = "component", values_to = "loading")
  }))
  reference_assignment <- ifelse(abs(observed$vectors[, 1]) >= abs(observed$vectors[, 2]), "pc1", "pc2")
  assignment <- loading_bootstrap |>
    pivot_wider(names_from = component, values_from = loading) |>
    mutate(assignment = if_else(abs(pc1) >= abs(pc2), "pc1", "pc2")) |>
    group_by(variable) |>
    summarise(
      reference_axis = reference_assignment[match(first(variable), variables)],
      retention_frequency = mean(assignment == reference_axis),
      .groups = "drop"
    )
  loading_intervals <- loading_bootstrap |>
    group_by(variable, component) |>
    summarise(
      median = median(loading),
      conf_low = quantile(loading, 0.025),
      conf_high = quantile(loading, 0.975),
      .groups = "drop"
    )
  bootstrap_summary <- bootstrap |>
    select(-aligned_loadings) |>
    summarise(across(
      c(pc1_explained, pc2_explained, cumulative_pc12, minimum_canonical_correlation),
      list(median = median, conf_low = ~quantile(.x, 0.025), conf_high = ~quantile(.x, 0.975))
    ))

  write_csv_safe(variance, file.path(root, "outputs", "tables", "weighted_pca_variance.csv"))
  write_csv_safe(loadings, file.path(root, "outputs", "tables", "weighted_pca_loadings.csv"))
  write_csv_safe(loading_intervals, file.path(root, "outputs", "tables", "pca_bootstrap_loadings.csv"))
  write_csv_safe(assignment, file.path(root, "outputs", "tables", "pca_bootstrap_assignment_stability.csv"))
  write_csv_safe(bootstrap |> select(-aligned_loadings), file.path(root, "outputs", "tables", "pca_bootstrap_replicates.csv"))
  write_csv_safe(bootstrap_summary, file.path(root, "outputs", "tables", "pca_bootstrap_global_stability.csv"))
  list(variance = variance, loadings = loadings, bootstrap = bootstrap,
       loading_intervals = loading_intervals, assignment = assignment)
}

directional_statistics <- function(exposure, outcome, weights) {
  quartiles <- weighted_quantile(exposure, weights, c(0.25, 0.75))
  median_outcome <- weighted_quantile(outcome, weights, 0.5)
  total <- sum(weights)
  c(
    low_high = sum(weights[exposure <= quartiles[1] & outcome >= median_outcome]) / total,
    high_low = sum(weights[exposure >= quartiles[2] & outcome < median_outcome]) / total
  )
}

run_discordance_benchmarks <- function(data, root = project_root(), config = read_config(root), quick = FALSE) {
  replicates <- if (quick) 50L else as.integer(config$discordance_permutations)
  analysis_rows <- primary_complete_indicator(data) & positive_weight(data$hormone_weight)
  analysis_data <- data[analysis_rows, , drop = FALSE]
  contrasts <- list(androgen = "androgen_score", amh = "z_amh")
  rows <- list()
  index <- 0L
  for (label in names(contrasts)) {
    exposure <- analysis_data[[contrasts[[label]]]]
    observed <- directional_statistics(
      exposure, analysis_data$nonfasting_score, analysis_data$hormone_weight
    )
    null <- with_seed(config$seed_discordance, replicate(replicates, {
      permuted <- permute_within_psu(exposure, analysis_data$SDMVPSU)
      directional_statistics(
        permuted, analysis_data$nonfasting_score, analysis_data$hormone_weight
      )
    }))
    for (direction in names(observed)) {
      null_values <- null[direction, ]
      index <- index + 1L
      rows[[index]] <- tibble(
        exposure = label,
        direction = direction,
        n = nrow(analysis_data),
        observed_proportion = observed[[direction]],
        analytical_independence = 0.125,
        permutation_mean = mean(null_values),
        observed_minus_null = observed[[direction]] - mean(null_values),
        permutation_p = (1 + sum(abs(null_values - mean(null_values)) >=
          abs(observed[[direction]] - mean(null_values)))) / (replicates + 1)
      )
    }
  }
  result <- bind_rows(rows) |> mutate(p_fdr = p.adjust(permutation_p, method = "BH"))
  write_csv_safe(result, file.path(root, "outputs", "tables", "discordance_independence_benchmark.csv"))
  result
}
