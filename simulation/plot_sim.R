library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(patchwork)
library(cowplot)
library(scales)
library(ragg)
library(latex2exp)

results_dir <- "results"
output_dir <- "figures"

model_pal <- c(
  "GSFA" = "#984EA3",
  "PerturbVI" = "#FF7F00"
)
fill_pal <- c(
  model_pal,
  "Mis-specified GSFA" = "#CBA6D6",
  "Mis-specified PerturbVI" = "#FDBE85"
)
legend_levels <- names(fill_pal)

read_sweep <- function(prefix, suffixes, values, parameter, method) {
  map2_dfr(suffixes, values, function(suffix, value) {
    path <- file.path(results_dir, paste0(prefix, suffix, ".csv"))
    dat <- read_csv(path, show_col_types = FALSE, progress = FALSE)
    dat %>% mutate(
      !!parameter := as.numeric(value),
      method = method
    )
  })
}

# ----------------------------- data assembly -----------------------------

b_values <- c(0.05, 0.10, 0.15, 0.20)
g_values <- c(50, 100, 150, 200)
l_values <- c(100, 125, 150, 175, 200)
z_values <- c(3, 4, 5, 6)
missing_values <- c(10, 20, 30, 40, 50)

pvi_b <- read_sweep(
  "perturbvi_", paste0("b_", b_values), b_values,
  "b_sparsity", "PerturbVI"
)
gsfa_b <- read_sweep(
  "gsfa_", paste0("b_", b_values), b_values,
  "b_sparsity", "GSFA"
)
b_data <- bind_rows(pvi_b, gsfa_b)

pvi_g <- read_sweep(
  "perturbvi_", paste0("g_", g_values), g_values,
  "g_dim", "PerturbVI"
)
gsfa_g <- read_sweep(
  "gsfa_", paste0("g_", g_values), g_values,
  "g_dim", "GSFA"
)
g_data <- bind_rows(pvi_g, gsfa_g)

# GSFA's L sweep is data sparsity (not a misspecified-L grid), so the L panels
# show PerturbVI only.
pvi_l <- read_sweep(
  "perturbvi_", paste0("l_", l_values), l_values,
  "l_dim", "PerturbVI"
)

pvi_z <- read_sweep(
  "perturbvi_", paste0("z_", z_values), z_values,
  "z_dim", "PerturbVI"
)
gsfa_z <- read_sweep(
  "gsfa_", paste0("z_", z_values), z_values,
  "z_dim", "GSFA"
)
z_data <- bind_rows(pvi_z, gsfa_z)

pvi_missing <- read_sweep(
  "perturbvi_", paste0("misg_", missing_values),
  missing_values, "remove_prop", "PerturbVI"
)
gsfa_missing <- read_sweep(
  "gsfa_", paste0("misg_", missing_values),
  missing_values, "remove_prop", "GSFA"
)
missing_data <- bind_rows(pvi_missing, gsfa_missing)

# reference medians are the correctly specified b=.2 and K=4 runs.
reference_stats <- function(dat, parameter, value) {
  dat %>%
    filter(.data[[parameter]] == value) %>%
    group_by(method) %>%
    summarise(
      reference_sensitivity = median(sensitivity, na.rm = TRUE),
      reference_beta = median(beta_err, na.rm = TRUE),
      reference_w = median(w_err, na.rm = TRUE),
      .groups = "drop"
    )
}

refs_b <- reference_stats(b_data, "b_sparsity", 0.2)
refs_z <- reference_stats(z_data, "z_dim", 4)

# Keep the sensitivity panels in the main figure on one scale, based on all
# plotted values with a small cushion so observed extrema remain visible.
main_sensitivity_values <- c(
  g_data$sensitivity,
  pvi_l %>% filter(l_dim != 150) %>% pull(sensitivity),
  z_data %>% filter(z_dim != 4) %>% pull(sensitivity),
  missing_data$sensitivity,
  missing_data$sensitivity_correct
)
main_sensitivity_range <- range(main_sensitivity_values, finite = TRUE)
main_sensitivity_padding <- diff(main_sensitivity_range) * 0.04
main_sensitivity_limits <- main_sensitivity_range +
  c(-main_sensitivity_padding, main_sensitivity_padding)

# ------------------------------- plot theme -------------------------------

STYLE <- list(
  family = "sans",
  ink = "#342D38",
  muted = "#342D38",
  base_size = 8.5,
  tag_size = 10,
  axis_title_size = 8.5,
  axis_tick_size = 7,
  legend_size = 6,
  plot_title_size = 7.5,
  axis_title_gap_pt = 3,
  axis_tick_gap_pt = 1.5,
  axis_line_width = 0.30,
  axis_tick_width = 0.22,
  axis_tick_length_mm = 1.2,
  plot_margin = c(2.5, 3, 2.5, 3)
)

figure_theme <- theme_classic(base_size = STYLE$base_size, base_family = STYLE$family) +
  theme(
    text = element_text(colour = STYLE$ink),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(
      colour = STYLE$ink, fill = NA,
      linewidth = STYLE$axis_line_width
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "#77717C", linewidth = STYLE$axis_tick_width),
    axis.ticks.length = grid::unit(STYLE$axis_tick_length_mm, "mm"),
    axis.text = element_text(colour = STYLE$ink, size = STYLE$axis_tick_size),
    axis.title = element_text(
      face = "plain", colour = STYLE$ink,
      size = STYLE$axis_title_size
    ),
    axis.title.x.bottom = element_text(
      margin = margin(t = STYLE$axis_title_gap_pt)
    ),
    axis.title.y.left = element_text(
      angle = 90, margin = margin(r = STYLE$axis_title_gap_pt)
    ),
    axis.text.x.bottom = element_text(
      margin = margin(t = STYLE$axis_tick_gap_pt)
    ),
    axis.text.y.left = element_text(
      margin = margin(r = STYLE$axis_tick_gap_pt)
    ),
    plot.title = element_text(
      face = "plain", hjust = 0.5,
      size = STYLE$plot_title_size,
      colour = STYLE$ink
    ),
    plot.tag = element_text(
      face = "bold", size = STYLE$tag_size,
      colour = STYLE$ink
    ),
    legend.title = element_text(
      face = "plain", colour = STYLE$ink,
      size = STYLE$legend_size
    ),
    legend.text = element_text(colour = STYLE$ink, size = STYLE$legend_size),
    legend.position = "bottom",
    legend.key = element_rect(fill = "white", colour = NA),
    legend.key.width = grid::unit(0.18 * 0.76, "cm"),
    legend.key.height = grid::unit(0.65, "lines"),
    legend.key.spacing = grid::unit(1.2 * 0.76, "mm"),
    legend.key.spacing.x = grid::unit(3 * 0.76, "mm"),
    legend.spacing.x = grid::unit(0.08, "lines"),
    legend.box.spacing = grid::unit(0.08, "lines"),
    legend.margin = margin(t = 2 * 0.76, b = 0),
    plot.margin = do.call(margin, as.list(STYLE$plot_margin))
  )

box_panel <- function(data, x, y, title, x_label, y_label, x_levels,
                      fill = "method", references = NULL, y_limits = NULL, log_y = FALSE) {
  dat <- data %>% mutate(
    .x_value = factor(as.character(.data[[x]]), levels = as.character(x_levels)),
    .fill_value = as.character(.data[[fill]])
  )

  p <- ggplot(dat, aes(x = .x_value, y = .data[[y]], fill = .fill_value)) +
    geom_boxplot(
      width = 0.62,
      linewidth = 0.42,
      colour = STYLE$ink,
      outlier.shape = 21,
      outlier.size = 1.25,
      outlier.fill = "white",
      outlier.colour = STYLE$ink,
      outlier.stroke = 0.35
    ) +
    scale_fill_manual(
      values = fill_pal, limits = legend_levels,
      drop = FALSE, name = NULL
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(x = x_label, y = y_label) +
    figure_theme

  if (log_y) {
    p <- p + scale_y_log10(expand = expansion(mult = c(0.04, 0.08)))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.04, 0.08)))
  }

  if (!is.null(references)) {
    p <- p +
      geom_hline(
        data = references,
        aes(yintercept = reference, colour = method),
        linewidth = 0.45,
        show.legend = FALSE
      ) +
      scale_colour_manual(values = model_pal)
  }
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  p
}

beta_title <- expression("Procrustes error in " * bolditalic(B))
w_title <- expression("Procrustes error in " * bolditalic(W))
k_beta_title <- expression("Misspecified " * italic(K) *
  ": error in " * bolditalic(B))
k_w_title <- expression("Misspecified " * italic(K) *
  ": error in " * bolditalic(W))
l_beta_title <- expression("Misspecified " * italic(L) *
  ": error in " * bolditalic(B))
l_w_title <- expression("Misspecified " * italic(L) *
  ": error in " * bolditalic(W))
k_main_title <- expression("Misspecified " * italic(K) *
  " (true " * italic(K) * " = 4)")
l_main_title <- expression("Misspecified " * italic(L) *
  " (true " * italic(L) * " = 150)")
k_sensitivity_title <- expression("Misspecified " * italic(K) *
  ": sensitivity")
l_sensitivity_title <- expression("Misspecified " * italic(L) *
  ": sensitivity")
k_specificity_title <- expression("Specificity: " * italic(K) * " sweep")
l_specificity_title <- expression("Specificity: " * italic(L) * " sweep")
beta_label <- TeX("$||B-Q\\hat{B}||^2_P$", bold = TRUE)
w_label <- TeX("$||W-Q\\hat{W}||^2_P$", bold = TRUE)
sensitivity_label <- "Sensitivity"

make_shared_legend <- function(levels) {
  legend_data <- data.frame(
    x = seq_along(levels),
    y = 1,
    method = factor(levels, levels = levels)
  )
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, fill = method)) +
    geom_point(shape = 22, size = 2.6, colour = STYLE$ink, stroke = 0.3) +
    scale_fill_manual(
      values = fill_pal[levels], limits = levels,
      drop = FALSE, name = NULL
    ) +
    guides(fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(shape = 22, size = 2.6)
    )) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(
        face = "plain", family = STYLE$family,
        colour = STYLE$ink, size = STYLE$legend_size
      ),
      legend.text = element_text(
        family = STYLE$family, colour = STYLE$ink,
        size = STYLE$legend_size
      ),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.width = grid::unit(0.18 * 0.76, "cm"),
      legend.key.height = grid::unit(0.65, "lines"),
      legend.key.spacing = grid::unit(1.2 * 0.76, "mm"),
      legend.key.spacing.x = grid::unit(3 * 0.76, "mm"),
      legend.spacing.x = grid::unit(0.08, "lines"),
      legend.box.spacing = grid::unit(0.08, "lines"),
      legend.margin = margin(t = 2 * 0.76, b = 0)
    )
  wrap_elements(full = cowplot::get_legend(legend_plot))
}

add_shared_legend <- function(panels, levels = legend_levels, title = NULL,
                              caption = NULL, keep_panel_layout = FALSE) {
  panel_plot <- panels & theme(legend.position = "none")
  if (keep_panel_layout) panel_plot <- wrap_elements(full = panel_plot)
  combined <- panel_plot /
    make_shared_legend(levels) +
    plot_layout(heights = c(1, 0.045))
  combined + plot_annotation(
    title = title,
    caption = caption,
    theme = theme(
      plot.title = element_text(
        face = "plain", family = STYLE$family,
        size = STYLE$plot_title_size, hjust = 0.5, colour = STYLE$ink
      ),
      plot.caption = element_text(
        family = STYLE$family, size = 10,
        hjust = 0.5, colour = STYLE$muted
      )
    )
  )
}

# ------------------------------- main figure ------------------------------

g_beta <- box_panel(
  g_data, "g_dim", "beta_err", beta_title,
  "Number of perturbations", beta_label, g_values,
  log_y = TRUE
)
g_w <- box_panel(
  g_data, "g_dim", "w_err", w_title,
  "Number of perturbations", w_label, g_values,
  log_y = TRUE
)
g_sensitivity <- box_panel(
  g_data, "g_dim", "sensitivity", "Sensitivity in overall effects",
  "Number of perturbations", sensitivity_label, g_values,
  y_limits = main_sensitivity_limits
)
g_beta <- g_beta + labs(tag = "A")
g_w <- g_w + labs(tag = "B")
g_sensitivity <- g_sensitivity + labs(tag = "C")

l_main <- pvi_l %>%
  filter(l_dim != 150) %>%
  mutate(plot_method = "Mis-specified PerturbVI")
l_refs_sensitivity <- refs_b %>% transmute(method, reference = reference_sensitivity)
l_sensitivity_main <- box_panel(
  l_main, "l_dim", "sensitivity", l_main_title,
  "Number of single effects", sensitivity_label, c(100, 125, 175, 200),
  fill = "plot_method", references = l_refs_sensitivity,
  y_limits = main_sensitivity_limits
)
l_sensitivity_main <- l_sensitivity_main + labs(tag = "D")
l_sensitivity_main <- l_sensitivity_main +
  theme(axis.title.y.left = element_text(colour = NA))

z_mis <- z_data %>%
  filter(z_dim != 4) %>%
  mutate(plot_method = paste("Mis-specified", method))
z_refs_sensitivity <- refs_z %>% transmute(method, reference = reference_sensitivity)
z_sensitivity_main <- box_panel(
  z_mis, "z_dim", "sensitivity", k_main_title,
  "Number of components", sensitivity_label, c(3, 5, 6),
  fill = "plot_method", references = z_refs_sensitivity,
  y_limits = main_sensitivity_limits
)
z_sensitivity_main <- z_sensitivity_main + labs(tag = "E")
z_sensitivity_main <- z_sensitivity_main +
  theme(axis.title.y.left = element_text(colour = NA))

missing_long <- missing_data %>%
  pivot_longer(
    cols = c(sensitivity, sensitivity_correct),
    names_to = "sensitivity_type",
    values_to = "sensitivity_value"
  ) %>%
  mutate(
    plot_method = case_when(
      method == "GSFA" & sensitivity_type == "sensitivity" ~ "Mis-specified GSFA",
      method == "GSFA" & sensitivity_type == "sensitivity_correct" ~ "GSFA",
      method == "PerturbVI" & sensitivity_type == "sensitivity" ~ "Mis-specified PerturbVI",
      TRUE ~ "PerturbVI"
    )
  )
missing_main <- box_panel(
  missing_long, "remove_prop", "sensitivity_value", "Missing targets (true = 0%)",
  "Missed targets %", sensitivity_label, missing_values,
  fill = "plot_method", y_limits = main_sensitivity_limits
)
missing_main <- missing_main + labs(tag = "F")

main_panels <- patchwork::wrap_plots(
  g_beta, g_w, g_sensitivity,
  plot_spacer(), plot_spacer(), plot_spacer(),
  l_sensitivity_main, z_sensitivity_main, missing_main,
  ncol = 3,
  widths = c(1, 1, 1),
  heights = c(1, 0.03, 1)
)
main_figure <- add_shared_legend(
  main_panels
)
main_figure <- cowplot::ggdraw(main_figure) +
  cowplot::draw_label(
    sensitivity_label,
    x = 0.068, y = 0.33, angle = 90,
    fontfamily = STYLE$family, size = STYLE$axis_title_size,
    colour = STYLE$ink, hjust = 0.5, vjust = 0.5
  ) +
  cowplot::draw_label(
    sensitivity_label,
    x = 0.405, y = 0.33, angle = 90,
    fontfamily = STYLE$family, size = STYLE$axis_title_size,
    colour = STYLE$ink, hjust = 0.5, vjust = 0.5
  )

# -------------------------- supplementary figure 1 ------------------------

b_beta <- box_panel(
  b_data, "b_sparsity", "beta_err", beta_title,
  "Perturbation density", beta_label, b_values,
  log_y = TRUE
)
b_w <- box_panel(
  b_data, "b_sparsity", "w_err", w_title,
  "Perturbation density", w_label, b_values,
  log_y = TRUE
)
b_sensitivity <- box_panel(
  b_data, "b_sparsity", "sensitivity", "Sensitivity in overall effects",
  "Perturbation density", sensitivity_label, b_values
)
b_beta <- b_beta + labs(tag = "A")
b_w <- b_w + labs(tag = "B")
b_sensitivity <- b_sensitivity + labs(tag = "C")

supplementary_density_panels <- (b_beta + b_w + b_sensitivity) +
  plot_layout(ncol = 3)
supplementary_density <- add_shared_legend(
  supplementary_density_panels,
  c("GSFA", "PerturbVI")
)

# -------------------------- supplementary figure 2 ------------------------

z_beta <- box_panel(
  z_mis, "z_dim", "beta_err", k_beta_title,
  "Number of components", beta_label, c(3, 5, 6),
  fill = "plot_method", references = refs_z %>% transmute(method, reference = reference_beta),
  log_y = TRUE
)
z_w <- box_panel(
  z_mis, "z_dim", "w_err", k_w_title,
  "Number of components", w_label, c(3, 5, 6),
  fill = "plot_method", references = refs_z %>% transmute(method, reference = reference_w),
  log_y = TRUE
)
z_sensitivity <- box_panel(
  z_mis, "z_dim", "sensitivity", k_sensitivity_title,
  "Number of components", sensitivity_label, c(3, 5, 6),
  fill = "plot_method", references = z_refs_sensitivity
)
z_beta <- z_beta + labs(tag = "A")
z_w <- z_w + labs(tag = "B")
z_sensitivity <- z_sensitivity + labs(tag = "C")

l_robust <- pvi_l %>% mutate(plot_method = "Mis-specified PerturbVI")
l_beta <- box_panel(
  l_robust, "l_dim", "beta_err", l_beta_title,
  "Number of single effects", beta_label, l_values,
  fill = "plot_method", references = refs_b %>% transmute(method, reference = reference_beta),
  log_y = TRUE
)
l_w <- box_panel(
  l_robust, "l_dim", "w_err", l_w_title,
  "Number of single effects", w_label, l_values,
  fill = "plot_method", references = refs_b %>% transmute(method, reference = reference_w),
  log_y = TRUE
)
l_sensitivity <- box_panel(
  l_robust, "l_dim", "sensitivity", l_sensitivity_title,
  "Number of single effects", sensitivity_label, l_values,
  fill = "plot_method", references = l_refs_sensitivity
)
l_beta <- l_beta + labs(tag = "D")
l_w <- l_w + labs(tag = "E")
l_sensitivity <- l_sensitivity + labs(tag = "F")

supplementary_robustness_panels <- patchwork::wrap_plots(
  z_beta, z_w, z_sensitivity,
  plot_spacer(), plot_spacer(), plot_spacer(),
  l_beta, l_w, l_sensitivity,
  ncol = 3, heights = c(1, 0.12, 1)
)
supplementary_robustness <- add_shared_legend(
  supplementary_robustness_panels,
  c("Mis-specified GSFA", "Mis-specified PerturbVI")
)

# -------------------------- supplementary figure 3 ------------------------

with_specificity <- function(dat) {
  dat %>% mutate(specificity_value = 1 - specificity)
}

b_spec_data <- with_specificity(b_data)
g_spec_data <- with_specificity(g_data)
z_spec_data <- with_specificity(z_data)
l_spec_data <- with_specificity(pvi_l %>% mutate(method = "PerturbVI"))

missing_spec_long <- missing_data %>%
  pivot_longer(
    cols = c(specificity, specificity_correct),
    names_to = "specificity_type",
    values_to = "specificity_raw"
  ) %>%
  mutate(
    specificity_value = 1 - specificity_raw,
    plot_method = case_when(
      method == "GSFA" & specificity_type == "specificity" ~ "Mis-specified GSFA",
      method == "GSFA" & specificity_type == "specificity_correct" ~ "GSFA",
      method == "PerturbVI" & specificity_type == "specificity" ~ "Mis-specified PerturbVI",
      TRUE ~ "PerturbVI"
    )
  )

# Specificity is close to 1, so zoom the display, but use one data-derived
# range for every panel and include all observed values.
specificity_values <- c(
  b_spec_data$specificity_value,
  g_spec_data$specificity_value,
  z_spec_data$specificity_value,
  l_spec_data$specificity_value,
  missing_spec_long$specificity_value
)
specificity_range <- range(specificity_values, finite = TRUE)
specificity_padding <- max(diff(specificity_range) * 0.04, 0.0002)
specificity_limits <- specificity_range +
  c(-specificity_padding, specificity_padding)

b_specificity <- box_panel(
  b_spec_data, "b_sparsity", "specificity_value", "Specificity: density sweep",
  "Perturbation density", "Specificity", b_values,
  y_limits = specificity_limits
)
g_specificity <- box_panel(
  g_spec_data, "g_dim", "specificity_value", "Specificity: perturbation sweep",
  "Number of perturbations", "Specificity", g_values,
  y_limits = specificity_limits
)
z_specificity <- box_panel(
  z_spec_data, "z_dim", "specificity_value", k_specificity_title,
  "Number of components", "Specificity", z_values,
  y_limits = specificity_limits
)
l_specificity <- box_panel(
  l_spec_data, "l_dim", "specificity_value", l_specificity_title,
  "Number of single effects", "Specificity", l_values,
  y_limits = specificity_limits
)
missing_specificity <- box_panel(
  missing_spec_long, "remove_prop", "specificity_value", "Specificity: missing targets",
  "Missed targets %", "Specificity", missing_values,
  fill = "plot_method", y_limits = specificity_limits
)
b_specificity <- b_specificity + labs(tag = "A")
g_specificity <- g_specificity + labs(tag = "B")
z_specificity <- z_specificity + labs(tag = "C")
l_specificity <- l_specificity + labs(tag = "D")
missing_specificity <- missing_specificity + labs(tag = "E")

specificity_top <- patchwork::wrap_plots(
  b_specificity, g_specificity, z_specificity,
  ncol = 3, widths = c(1, 1, 1)
)
specificity_bottom <- patchwork::wrap_plots(
  plot_spacer(), l_specificity, missing_specificity, plot_spacer(),
  ncol = 4, widths = c(0.5, 1, 1, 0.5)
)
supplementary_specificity_panels <- specificity_top / specificity_bottom +
  plot_layout(heights = c(1, 1))
supplementary_specificity <- add_shared_legend(
  supplementary_specificity_panels,
  keep_panel_layout = TRUE
)

# ------------------------------- file output -------------------------------

save_figure <- function(plot, stem, width, height) {
  output_width <- 17.4 / 2.54
  ggsave(file.path(output_dir, paste0(stem, ".png")),
    plot = plot,
    device = ragg::agg_png, width = output_width,
    height = height * output_width / width, units = "in", dpi = 600,
    bg = "white", limitsize = FALSE
  )
}

save_figure(main_figure, "figure_1", 13.7, 10.0)
save_figure(supplementary_density, "suppl_figure_1_density", 10.3, 5.0)
save_figure(supplementary_robustness, "suppl_figure_2_robustness", 13.7, 10.0)
save_figure(supplementary_specificity, "suppl_figure_3_specificity", 13.7, 7.6)

message("Wrote regenerated figures to: ", output_dir)
