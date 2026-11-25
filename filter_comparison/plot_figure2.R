script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    normalizePath(getwd())
  } else {
    file_path <- sub("^--file=", "", file_arg[[1L]])
    file_path <- gsub("~+~", " ", file_path, fixed = TRUE)
    normalizePath(dirname(file_path))
  }
})

packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "ragg")
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Please install ", paste(missing_packages, collapse = ", "), ".")
}

read_results <- function(results_directory, scenario_name) {
  result_files <- list.files(
    results_directory,
    pattern = "-results[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(result_files) == 0L) {
    stop("No result files were found under ", results_directory)
  }
  dplyr::bind_rows(lapply(result_files, function(file_name) {
    value <- utils::read.csv(file_name, check.names = FALSE)
    if (names(value)[1L] %in% c("", "X", "X1")) {
      value <- value[, -1L, drop = FALSE]
    }
    value$state <- as.integer(
      sub("^result-state", "", basename(dirname(file_name)))
    )
    value$p <- as.integer(
      sub("-results[.]csv$", "", basename(file_name))
    )
    value$scenario <- scenario_name
    value
  }))
}

admixture_directory <- file.path(script_dir, "admixture", "results_bic")
family_directory <- file.path(script_dir, "family", "results_bic")

results <- dplyr::bind_rows(
  read_results(admixture_directory, "Low-rank admixture"),
  read_results(family_directory, "Family structure")
) |>
  dplyr::mutate(
    method = dplyr::case_when(
      Method == "BH" ~ "LMM-BH",
      Method == "LMM" ~ "LMM-Lasso",
      Method == "PCA-O" ~ "PCA-Lasso (O)",
      Method == "PCA-E" ~ "PCA-Lasso (E)",
      grepl("DDMM-Trim", Method, fixed = TRUE) ~ "DDMM-Trim",
      grepl("DDMM-LMM", Method, fixed = TRUE) ~ "DDMM-LMM",
      grepl("DDMM-PCA", Method, fixed = TRUE) ~ "DDMM-PCA",
      grepl("DDMM-O", Method, fixed = TRUE) ~ "DDMM-PCA (O)",
      grepl("DDMM-E", Method, fixed = TRUE) ~ "DDMM-PCA (E)",
      TRUE ~ Method
    )
  )

admixture_methods <- c(
  "DDMM-PCA (O)",
  "DDMM-PCA (E)",
  "PCA-Lasso (O)",
  "PCA-Lasso (E)",
  "ggmix",
  "plmmr"
)
family_methods <- c(
  "DDMM-Trim",
  "DDMM-PCA",
  "LMM-Lasso",
  "LMM-BH",
  "ggmix",
  "plmmr"
)

plot_data <- results |>
  dplyr::filter(
    (scenario == "Low-rank admixture" & method %in% admixture_methods) |
      (scenario == "Family structure" & method %in% family_methods)
  ) |>
  tidyr::pivot_longer(
    cols = c("FDR", "Power"),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::group_by(scenario, method, p, metric) |>
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    sd = stats::sd(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    lower = pmax(mean - sd, 0),
    upper = pmin(mean + sd, 1)
  )

estimation_data <- results |>
  dplyr::filter(
    (scenario == "Low-rank admixture" & method %in% admixture_methods) |
      (scenario == "Family structure" & method %in% family_methods)
  ) |>
  tidyr::pivot_longer(
    cols = c("L1_Error", "MSPE"),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::group_by(scenario, method, p, metric) |>
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    sd = stats::sd(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    lower = pmax(mean - sd, 0),
    upper = mean + sd
  )

method_colors <- c(
  "LMM-BH" = "#777777",
  "ggmix" = "#CC79A7",
  "plmmr" = "#7A3E9D",
  "LMM-Lasso" = "#56B4E9",
  "PCA-Lasso (O)" = "#E69F00",
  "PCA-Lasso (E)" = "#E69F00",
  "DDMM-PCA (O)" = "#009E73",
  "DDMM-PCA (E)" = "#009E73",
  "DDMM-Trim" = "#0072B2",
  "DDMM-PCA" = "#009E73"
)
method_shapes <- c(
  "LMM-BH" = 4,
  "ggmix" = 18,
  "plmmr" = 8,
  "LMM-Lasso" = 17,
  "PCA-Lasso (O)" = 16,
  "PCA-Lasso (E)" = 17,
  "DDMM-PCA (O)" = 16,
  "DDMM-PCA (E)" = 17,
  "DDMM-Trim" = 15,
  "DDMM-PCA" = 16
)
method_lines <- c(
  "LMM-BH" = "dotted",
  "ggmix" = "dashed",
  "plmmr" = "dashed",
  "LMM-Lasso" = "dashed",
  "PCA-Lasso (O)" = "dashed",
  "PCA-Lasso (E)" = "dotdash",
  "DDMM-PCA (O)" = "solid",
  "DDMM-PCA (E)" = "longdash",
  "DDMM-Trim" = "solid",
  "DDMM-PCA" = "solid"
)

make_panel <- function(scenario_name, metric_name, panel_label, method_order) {
  panel_data <- plot_data |>
    dplyr::filter(
      scenario == scenario_name,
      metric == metric_name
    ) |>
    dplyr::mutate(
      method = factor(method, levels = method_order)
    )

  panel <- ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = p,
      y = mean,
      color = method,
      shape = method,
      linetype = method
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 140,
      linewidth = 0.65,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::geom_point(size = 3.6, stroke = 0.9) +
    ggplot2::scale_color_manual(
      values = method_colors,
      breaks = method_order,
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = method_shapes,
      breaks = method_order,
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = method_lines,
      breaks = method_order,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1000, 2000, 4000, 6000, 8000),
      labels = function(x) {
        format(x, big.mark = ",", scientific = FALSE)
      }
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(mult = c(0.01, 0.035))
    ) +
    ggplot2::labs(
      title = paste0(
        panel_label,
        "  ",
        if (metric_name == "FDR") "Empirical FDR" else "Power"
      ),
      x = "Number of tested variants, p",
      y = if (metric_name == "FDR") "Empirical FDR" else "Power",
      color = NULL,
      shape = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_classic(base_family = "sans", base_size = 17) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 20,
        hjust = 0
      ),
      axis.title = ggplot2::element_text(size = 17),
      axis.text = ggplot2::element_text(
        color = "black",
        size = 15
      ),
      axis.line = ggplot2::element_line(linewidth = 0.8),
      axis.ticks = ggplot2::element_line(linewidth = 0.7),
      panel.grid.major.y = ggplot2::element_line(
        color = "#E3E3E3",
        linewidth = 0.55
      ),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 14.5),
      legend.key.width = grid::unit(1.55, "lines"),
      legend.spacing.x = grid::unit(0.35, "cm"),
      plot.margin = ggplot2::margin(8, 12, 3, 8)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    )

  if (metric_name == "FDR") {
    panel <- panel + ggplot2::geom_hline(
      yintercept = 0.1,
      color = "#333333",
      linewidth = 0.75,
      linetype = "dotted",
      inherit.aes = FALSE
    )
  }
  panel
}

make_estimation_panel <- function(
    scenario_name,
    metric_name,
    panel_label,
    method_order) {
  panel_data <- estimation_data |>
    dplyr::filter(
      scenario == scenario_name,
      metric == metric_name
    ) |>
    dplyr::mutate(
      method = factor(method, levels = method_order)
    )

  ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = p,
      y = mean,
      color = method,
      shape = method,
      linetype = method
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 140,
      linewidth = 0.65,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::geom_point(size = 3.6, stroke = 0.9) +
    ggplot2::scale_color_manual(
      values = method_colors,
      breaks = method_order,
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = method_shapes,
      breaks = method_order,
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = method_lines,
      breaks = method_order,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1000, 2000, 4000, 6000, 8000),
      labels = function(x) {
        format(x, big.mark = ",", scientific = FALSE)
      }
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.04, 0.10))
    ) +
    ggplot2::labs(
      title = paste0(
        panel_label,
        "  ",
        if (metric_name == "L1_Error") {
          "L1 estimation error"
        } else {
          "Mean squared prediction error"
        }
      ),
      x = "Number of tested variants, p",
      y = if (metric_name == "L1_Error") {
        "L1 estimation error"
      } else {
        "MSPE"
      },
      color = NULL,
      shape = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_classic(base_family = "sans", base_size = 17) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 20,
        hjust = 0
      ),
      axis.title = ggplot2::element_text(size = 17),
      axis.text = ggplot2::element_text(
        color = "black",
        size = 15
      ),
      axis.line = ggplot2::element_line(linewidth = 0.8),
      axis.ticks = ggplot2::element_line(linewidth = 0.7),
      panel.grid.major.y = ggplot2::element_line(
        color = "#E3E3E3",
        linewidth = 0.55
      ),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 14.5),
      legend.key.width = grid::unit(1.55, "lines"),
      legend.spacing.x = grid::unit(0.35, "cm"),
      plot.margin = ggplot2::margin(8, 12, 3, 8)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    )
}

admixture_row <- (
  make_panel(
    "Low-rank admixture",
    "FDR",
    "A",
    admixture_methods
  ) +
    make_panel(
      "Low-rank admixture",
      "Power",
      "B",
      admixture_methods
    ) +
    patchwork::plot_layout(guides = "collect")
) & ggplot2::theme(legend.position = "bottom")

family_row <- (
  make_panel(
    "Family structure",
    "FDR",
    "A",
    family_methods
  ) +
    make_panel(
      "Family structure",
      "Power",
      "B",
      family_methods
    ) +
    patchwork::plot_layout(guides = "collect")
) & ggplot2::theme(legend.position = "bottom")

admixture_figure <- admixture_row +
  patchwork::plot_annotation(
    title = "Low-rank admixture with observed and external PCA filters",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 21,
        hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      )
    )
  )

family_figure <- family_row +
  patchwork::plot_annotation(
    title = "Family structure with PCA and Trim filters",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 21,
        hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      )
    )
  )

admixture_estimation_figure <- (
  make_estimation_panel(
    "Low-rank admixture",
    "L1_Error",
    "A",
    admixture_methods
  ) +
    make_estimation_panel(
      "Low-rank admixture",
      "MSPE",
      "B",
      admixture_methods
    ) +
    patchwork::plot_layout(guides = "collect")
) & ggplot2::theme(legend.position = "bottom")

admixture_estimation_figure <- admixture_estimation_figure +
  patchwork::plot_annotation(
    title = "Estimation and prediction error under low-rank admixture",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 21,
        hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      )
    )
  )

family_estimation_figure <- (
  make_estimation_panel(
    "Family structure",
    "L1_Error",
    "A",
    family_methods
  ) +
    make_estimation_panel(
      "Family structure",
      "MSPE",
      "B",
      family_methods
    ) +
    patchwork::plot_layout(guides = "collect")
) & ggplot2::theme(legend.position = "bottom")

family_estimation_figure <- family_estimation_figure +
  patchwork::plot_annotation(
    title = "Estimation and prediction error under family structure",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 21,
        hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      )
    )
  )

output_directory <- file.path(script_dir, "figure")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

save_publication_figure <- function(plot, file_stem) {
  width <- 12
  height <- 5.8

  ggplot2::ggsave(
    paste0(file_stem, ".pdf"),
    plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    paste0(file_stem, ".png"),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    background = "white"
  )
  ggplot2::ggsave(
    paste0(file_stem, ".tiff"),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = ragg::agg_tiff,
    compression = "lzw",
    background = "white"
  )
}

save_publication_figure(
  admixture_figure,
  file.path(output_directory, "admixture_external_filter_comparison")
)
save_publication_figure(
  family_figure,
  file.path(output_directory, "family_filter_comparison")
)
save_publication_figure(
  admixture_estimation_figure,
  file.path(output_directory, "supplementary_admixture_l1_mspe")
)
save_publication_figure(
  family_estimation_figure,
  file.path(output_directory, "supplementary_family_l1_mspe")
)

message("Saved the main and supplementary simulation figures under ", output_directory)
