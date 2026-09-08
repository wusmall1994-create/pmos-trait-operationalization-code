make_code_release_figures <- function(root = project_root()) {
  dir.create(file.path(root, "outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
  matrix_path <- file.path(root, "outputs", "tables", "operationalization_construct_matrix.csv")
  if (file.exists(matrix_path)) {
    matrix <- readr::read_csv(matrix_path, show_col_types = FALSE)
    plot <- ggplot2::ggplot(
      matrix,
      ggplot2::aes(x = outcome, y = exposure, fill = estimate)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~stage) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
      ggplot2::labs(x = NULL, y = NULL, fill = "Standardized beta") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(
      file.path(root, "outputs", "figures", "operationalization_matrix.png"),
      plot, width = 12, height = 7, dpi = 300
    )
  }
  invisible(TRUE)
}
