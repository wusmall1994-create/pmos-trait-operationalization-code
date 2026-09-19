suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) stop("Run this script from the repository root.")
out_dir <- Sys.getenv("PMOS_FIGURE_DIR", unset = file.path(root, "outputs", "figures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

box <- function(xmin, xmax, ymin, ymax, fill) {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           fill = fill, colour = "#667085", linewidth = 0.75)
}

p <- ggplot() +
  box(0.4, 4.6, 1.2, 8.8, "#EAF4F8") +
  box(5.9, 10.6, 5.5, 9.0, "#F4F5F7") +
  box(5.9, 10.6, 0.9, 4.9, "#FFF8EC") +
  annotate("segment", x = 4.7, xend = 5.65, y = 5.9, yend = 7.0,
           colour = "#667085", linewidth = 0.8,
           arrow = arrow(length = grid::unit(0.16, "in"), type = "open")) +
  annotate("segment", x = 4.7, xend = 5.65, y = 4.3, yend = 3.4,
           colour = "#667085", linewidth = 0.8,
           arrow = arrow(length = grid::unit(0.16, "in"), type = "open")) +
  annotate("text", x = 2.5, y = 7.8, label = "NHANES 2021-2023",
           family = "Arial", fontface = "bold", size = 5.2, colour = "#101828") +
  annotate("text", x = 2.5, y = 6.7, label = "Cross-sectional discovery",
           family = "Arial", fontface = "bold", size = 4.1, colour = "#155E75") +
  annotate("text", x = 2.5, y = 4.8,
           label = paste("Strict eligibility: n = 824",
                         "Fully adjusted primary model: n = 643",
                         "Common fasting sample: n = 316",
                         "Testosterone, androstenedione, DHEAS,",
                         "AMH, and SHBG", sep = "\n"),
           family = "Arial", size = 3.15, lineheight = 1.35, colour = "#344054") +
  annotate("text", x = 8.25, y = 8.25, label = "NHANES 2017-March 2020",
           family = "Arial", fontface = "bold", size = 4.55, colour = "#101828") +
  annotate("text", x = 8.25, y = 7.55, label = "Temporal comparison",
           family = "Arial", fontface = "bold", size = 3.8, colour = "#155E75") +
  annotate("text", x = 8.25, y = 6.55,
           label = paste("Strict eligibility: n = 1,180",
                         "AMH, androstenedione,", "and SHBG available", sep = "\n"),
           family = "Arial", size = 3.05, lineheight = 1.3, colour = "#344054") +
  annotate("text", x = 8.25, y = 4.1, label = "SWAN repeated visits",
           family = "Arial", fontface = "bold", size = 4.8, colour = "#101828") +
  annotate("text", x = 8.25, y = 3.35, label = "Longitudinal within-person analysis",
           family = "Arial", fontface = "bold", size = 3.55, colour = "#155E75") +
  annotate("text", x = 8.25, y = 2.25,
           label = paste("7 visits; 2,746 women; 13,051 observations",
                         "Testosterone, DHEAS, and SHBG",
                         "Fasting metabolic measures and adiposity", sep = "\n"),
           family = "Arial", size = 2.95, lineheight = 1.3, colour = "#344054") +
  coord_cartesian(xlim = c(0, 11), ylim = c(0.5, 9.5), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 8, 8, 8))

ggsave(file.path(out_dir, "Figure_1_study_cohorts_and_analytical_roles.pdf"), p,
       width = 183 / 25.4, height = 105 / 25.4, device = cairo_pdf)
ragg::agg_tiff(file.path(out_dir, "Figure_1_study_cohorts_and_analytical_roles.tiff"),
               width = 183 / 25.4, height = 105 / 25.4, units = "in", res = 600,
               compression = "lzw")
print(p)
dev.off()
ragg::agg_png(file.path(out_dir, "Figure_1_study_cohorts_and_analytical_roles.png"),
              width = 183 / 25.4, height = 105 / 25.4, units = "in", res = 300)
print(p)
dev.off()

