#!/usr/bin/env Rscript

# Generates Figure 2 and supplementary files.

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(scales)
library(RColorBrewer)
library(ragg)

# Keep layout measurement plots off the default Rplots.pdf device.
pdf(NULL)

STYLE <- list(
  family = "sans", ink = "#342D38", muted = "#766C7A",
  base_size = 8, tag_size = 12, tag_band = .17,
  axis_title_size = 9, axis_tick_size = 7.5,
  deg_tick_size = 6,
  marker_x_tick_size = 6.2,
  marker_tick_size = 5.2,
  axis_title_gap_pt = 6, axis_tick_gap_pt = 4,
  count_axis_lower_expansion = .03, count_label_gap_pt = .8,
  methods = setNames(
    RColorBrewer::brewer.pal(9, "Set1")[4:5],
    c("GSFA", "PerturbVI")
  ),
  enrichment = RColorBrewer::brewer.pal(9, "Reds")[3:8],
  factor_fdr_limits = c(1, 3), factor_fdr_ticks = c(1, 2, 3),
  effect = grDevices::colorRampPalette(
    c("#2166AC", "#67A9CF", "#FFFFFF", "#EF8A62", "#B2182B"),
    space = "rgb"
  )(257),
  heatmap_frame_color = "#555555", heatmap_frame_linewidth = .25,
  colorbar_length = .78, colorbar_thickness = .10,
  colorbar_title_size = 7, colorbar_tick_size = 6.5,
  effect_key_height = 1.40,
  width_cm = 17.4, row_height_cm = 5.8, dpi = 300
)

fit_dir <- "fit"
results_dir <- "results"
figures_dir <- "figures"

axis_typography <- function() {
  theme(
    axis.title = element_text(
      family = STYLE$family,
      size = STYLE$axis_title_size,
      face = "plain", color = STYLE$ink
    ),
    axis.title.x.bottom = element_text(
      margin = margin(t = STYLE$axis_title_gap_pt)
    ),
    axis.title.y.left = element_text(
      angle = 90, margin = margin(r = STYLE$axis_title_gap_pt)
    ),
    axis.text = element_text(
      family = STYLE$family,
      size = STYLE$axis_tick_size,
      color = STYLE$ink
    ),
    axis.text.x.bottom = element_text(
      margin = margin(t = STYLE$axis_tick_gap_pt)
    ),
    axis.text.y.left = element_text(
      margin = margin(r = STYLE$axis_tick_gap_pt)
    )
  )
}

theme_panel <- function() {
  theme_classic(base_size = STYLE$base_size, base_family = STYLE$family) +
    axis_typography() +
    theme(
      text = element_text(color = STYLE$ink),
      axis.line = element_line(color = "#55515A", linewidth = .25),
      axis.ticks = element_line(color = "#77717C", linewidth = .22),
      axis.ticks.length = grid::unit(1.2, "mm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      legend.background = element_blank(), legend.key = element_blank(),
      legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
      plot.margin = margin(4, 5, 4, 5)
    )
}

matrix_theme <- function() {
  theme_void(base_size = STYLE$base_size, base_family = STYLE$family) +
    axis_typography() +
    theme(
      axis.text.x.top = element_text(
        angle = 60, hjust = 0, vjust = 0,
        margin = margin(b = STYLE$axis_tick_gap_pt)
      ),
      axis.text.y.right = element_text(
        face = "italic", hjust = 0,
        margin = margin(l = STYLE$axis_tick_gap_pt)
      ),
      plot.margin = margin(4, 4, 4, 4)
    )
}

heatmap_frame <- function(n_columns, n_rows, x_offset = 0) {
  annotate("rect",
    xmin = .5 + x_offset, xmax = n_columns + .5 + x_offset,
    ymin = .5, ymax = n_rows + .5, fill = NA,
    color = STYLE$heatmap_frame_color,
    linewidth = STYLE$heatmap_frame_linewidth
  )
}

effect_scale <- function(limit) {
  scale_fill_gradientn(
    colors = STYLE$effect, limits = c(-limit, limit),
    guide = "none"
  )
}

# Panel A, copied from merged_luhmes.R.
make_factor_counts <- function(data, width, height) {
  ggplot(data, aes(Method, Count, color = Method, fill = Method)) +
    geom_boxplot(width = .44, alpha = .19, linewidth = .5, outlier.shape = NA) +
    geom_point(
      position = position_jitter(width = .08, height = 0, seed = 42),
      shape = 21, size = 1.7, stroke = .25, color = "white"
    ) +
    scale_color_manual(values = STYLE$methods) +
    scale_fill_manual(values = STYLE$methods) +
    scale_y_continuous(
      limits = c(0, 1080), breaks = seq(0, 1000, 200),
      labels = scales::label_comma(), expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Genes with PIP > 0.95") +
    theme_panel() +
    theme(legend.position = "none", axis.ticks.x = element_blank())
}

panel_frame <- function(body, tag, width, height) {
  ggdraw() +
    draw_plot(body,
      x = 0, y = 0, width = 1,
      height = (height - STYLE$tag_band) / height
    ) +
    draw_label(tag,
      x = .025 / width, y = (height - .015) / height,
      hjust = 0, vjust = 1, size = STYLE$tag_size,
      fontface = "bold", color = STYLE$ink
    )
}

place <- function(canvas, plot, left, bottom, width, height,
                  canvas_width, canvas_height) {
  canvas + draw_plot(plot,
    x = left / canvas_width, y = bottom / canvas_height,
    width = width / canvas_width, height = height / canvas_height
  )
}

load_factor_counts <- function() {
  read_mat <- function(path) {
    as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  }

  gsfa <- readRDS(file.path(fit_dir, "gsfa_fit.light.rds"))
  gsfa_pip <- gsfa$posterior_means$F_pm
  perturbvi_pip <- read_mat(file.path(fit_dir, "PIP_W.csv"))
  n_factors <- nrow(perturbvi_pip)

  counts <- rbind(
    data.frame(
      Factor = seq_len(n_factors), Method = "GSFA",
      Count = unname(colSums(gsfa_pip > 0.95))
    ),
    data.frame(
      Factor = seq_len(n_factors), Method = "PerturbVI",
      Count = unname(rowSums(perturbvi_pip > 0.95))
    )
  )
  counts$Method <- factor(counts$Method, levels = c("GSFA", "PerturbVI"))
  counts
}

load_deg_counts <- function() {
  read_mat <- function(path) {
    as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  }
  gsfa_fit <- readRDS(file.path(fit_dir, "gsfa_fit.light.rds"))
  gsfa_lfsr <- gsfa_fit$lfsr
  gsfa_counts <- colSums(gsfa_lfsr < 0.05, na.rm = TRUE)
  gsfa_counts <- gsfa_counts[!names(gsfa_counts) %in% c("Nontargeting", "offset")]

  perturbvi_lfsr <- read_mat(file.path(fit_dir, "LFSR_BW.csv"))
  perturbations <- rownames(perturbvi_lfsr)

  data.frame(
    perturbation = perturbations,
    GSFA = unname(gsfa_counts[perturbations]),
    PerturbVI = rowSums(perturbvi_lfsr < 0.05, na.rm = TRUE),
    check.names = FALSE
  )
}

load_enriched_go_counts <- function() {
  counts <- read.csv(
    file.path(
      results_dir,
      "go_bp_significant_term_counts.csv"
    ),
    check.names = FALSE
  )
  counts[
    rowSums(counts[c("GSFA", "PerturbVI")]) > 0,
    c("perturbation", "GSFA", "PerturbVI")
  ]
}

load_marker_effect_data <- function() {
  read_mat <- function(path) {
    as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  }

  settings <- list(
    perturbations = setdiff(
      load_enriched_go_counts()$perturbation,
      c("CTNND2", "DYRK1A", "HDAC5")
    ),
    color_limits = c(-0.34, 0.34)
  )
  markers <- read.csv(file.path(results_dir, "suppl_fig_2E_annotation.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  marker_colors <- c(
    "Neuronal differentiation and maturation" = "#66C2A5",
    "Regulation of neuronal differentiation" = "#FC8D62",
    "Neuron projection development" = "#E78AC3",
    "Regulation of neuron projection development" = "#8DA0CB"
  )
  markers$Color <- unname(marker_colors[markers$Annotation])
  markers$Row <- nrow(markers) + 1L - seq_len(nrow(markers))

  effects <- read_mat(file.path(fit_dir, "BW.csv"))
  data <- expand.grid(
    Gene = markers$Gene,
    Perturbation = settings$perturbations,
    stringsAsFactors = FALSE
  )
  ensembl <- markers$gene_ID[match(data$Gene, markers$Gene)]
  data$Overall_effect <- effects[cbind(data$Perturbation, ensembl)]

  # Match Figure 4's ordering: keep annotation groups contiguous, cluster genes
  # within each group, and cluster perturbations across all marker genes. Use
  # Ward.D2 with Euclidean distance and display no dendrograms.
  effect_profiles <- matrix(
    data$Overall_effect,
    nrow = nrow(markers),
    ncol = length(settings$perturbations)
  )
  annotation_order <- c(
    "Neuronal differentiation and maturation",
    "Regulation of neuronal differentiation",
    "Neuron projection development",
    "Regulation of neuron projection development"
  )
  annotation_order <- annotation_order[
    annotation_order %in% markers$Annotation
  ]
  gene_order <- unlist(lapply(annotation_order, function(annotation) {
    idx <- which(markers$Annotation == annotation)
    if (length(idx) < 2L) {
      return(idx)
    }
    idx[stats::hclust(
      stats::dist(effect_profiles[idx, , drop = FALSE]),
      method = "ward.D2"
    )$order]
  }), use.names = FALSE)
  markers <- markers[gene_order, , drop = FALSE]
  perturbation_order <- stats::hclust(
    stats::dist(t(effect_profiles)),
    method = "ward.D2"
  )$order
  settings$perturbations <- settings$perturbations[perturbation_order]
  markers$Row <- nrow(markers) + 1L - seq_len(nrow(markers))

  data$Column <- match(data$Perturbation, settings$perturbations)
  data$Row <- nrow(markers) + 1L - match(data$Gene, markers$Gene)
  list(data = data, markers = markers, settings = settings)
}

load_new_factor_enrichment <- function() {
  enrichment <- read.csv(
    file.path(results_dir, "suppl_fig_2F_perturbvi.csv"),
    check.names = FALSE, stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  neuronal_pattern <- paste0(
    "neuro|neural|nerve|nervous|axon|dendrit|brain|synap|glia|",
    "myelin|cerebr|cerebell|hippocamp"
  )
  enrichment <- enrichment[
    grepl(neuronal_pattern, enrichment$description, ignore.case = TRUE) &
      enrichment$fdr < .05, ,
    drop = FALSE
  ]
  enrichment$Factor <- as.integer(sub("^factor_", "", enrichment$factor)) + 1L
  enrichment$Neg_log10_FDR <- -log10(pmax(enrichment$fdr, .Machine$double.xmin))
  enrichment <- enrichment[order(
    enrichment$Factor, enrichment$fdr, enrichment$p_value, enrichment$term_id
  ), , drop = FALSE]
  enrichment$Display_rank <- ave(
    seq_len(nrow(enrichment)), enrichment$Factor,
    FUN = seq_along
  )
  enrichment
}

load_factor_effect_matrix <- function() {
  as.matrix(read.csv(file.path(fit_dir, "B.csv"),
    row.names = 1,
    check.names = FALSE
  ))
}

factor_effect_data <- function(effect_matrix, perturbations, factors) {
  factor_columns <- paste0("factor_", factors - 1L)

  data <- expand.grid(
    Perturbation = perturbations, Factor = factors,
    stringsAsFactors = FALSE
  )
  data$Effect <- effect_matrix[cbind(
    data$Perturbation,
    factor_columns[match(data$Factor, factors)]
  )]
  data$Column <- match(data$Factor, factors)
  data$Row <- length(perturbations) + 1L - match(data$Perturbation, perturbations)
  data
}

factor_dendrogram_segments <- function(clustering) {
  n <- length(clustering$order)
  segments <- data.frame(
    x = numeric(), xend = numeric(),
    y = numeric(), yend = numeric()
  )

  visit <- function(node) {
    if (node < 0) {
      leaf <- -node
      rank <- match(leaf, clustering$order)
      return(list(x = 0, y = n + 1L - rank))
    }

    left <- visit(clustering$merge[node, 1])
    right <- visit(clustering$merge[node, 2])
    parent_x <- clustering$height[node]
    segments <<- rbind(
      segments,
      data.frame(x = left$x, xend = parent_x, y = left$y, yend = left$y),
      data.frame(x = right$x, xend = parent_x, y = right$y, yend = right$y),
      data.frame(x = parent_x, xend = parent_x, y = left$y, yend = right$y)
    )
    list(x = parent_x, y = (left$y + right$y) / 2)
  }

  visit(n - 1L)
  segments
}

make_factor_effects <- function(data, settings, width, height,
                                row_dendrogram = NULL) {
  factors <- settings$selected_factor_numbers
  targets <- settings$selected_perturbations
  limit <- max(abs(settings$color_limits_from_full_figure))
  pretty_ticks <- pretty(c(-limit, limit), n = 5)
  tick_max <- max(pretty_ticks[pretty_ticks > 0 & pretty_ticks < limit])
  ticks <- c(-tick_max, 0, tick_max)
  dendrogram_start <- length(factors) + 2.8
  dendrogram_width <- 3.25
  x_upper <- if (is.null(row_dendrogram)) {
    NULL
  } else {
    dendrogram_start + dendrogram_width + 1.2
  }

  p <- ggplot(data, aes(Column, Row, fill = Effect)) +
    geom_tile(width = 1, height = 1) +
    heatmap_frame(length(factors), length(targets)) +
    effect_scale(limit) +
    scale_x_continuous(
      breaks = seq_along(factors),
      labels = paste("Factor", factors),
      position = "top", expand = c(0, 0),
      limits = if (is.null(x_upper)) NULL else c(.5, x_upper),
      sec.axis = dup_axis(name = "Factors", breaks = NULL)
    ) +
    scale_y_continuous(
      breaks = seq_along(targets), labels = rev(targets),
      position = "right", expand = c(0, 0),
      sec.axis = dup_axis(name = "Perturbations", breaks = NULL)
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    matrix_theme()

  if (!is.null(row_dendrogram)) {
    segments <- factor_dendrogram_segments(row_dendrogram)
    max_height <- max(row_dendrogram$height)
    segments$x <- dendrogram_start + segments$x / max_height * dendrogram_width
    segments$xend <- dendrogram_start + segments$xend / max_height * dendrogram_width
    row_labels <- data.frame(
      Perturbation = targets,
      Row = length(targets) + 1L - seq_along(targets)
    )
    p <- p +
      geom_text(
        data = row_labels,
        aes(x = length(factors) + .8, y = Row, label = Perturbation),
        inherit.aes = FALSE, hjust = 0, vjust = .5,
        size = STYLE$axis_tick_size / ggplot2::.pt,
        family = STYLE$family, fontface = "italic", color = STYLE$ink
      ) +
      geom_segment(
        data = segments,
        aes(x = x, xend = xend, y = y, yend = yend),
        inherit.aes = FALSE, color = STYLE$ink, linewidth = .35
      ) +
      theme(
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
      )
  }

  key_width <- .62
  plot_width <- width - key_width - .05
  canvas <- place(ggdraw(), p, 0, 0, plot_width, height, width, height)
  key_height <- min(STYLE$effect_key_height, height)
  key <- effect_key(limit, ticks, "Effect\nsize", key_width, key_height,
    title_gap = .25
  )
  place(
    canvas, key, width - key_width, (height - key_height) / 2,
    key_width, key_height, width, height
  )
}

marker_annotation_key <- function(markers, width, height = .78) {
  groups <- markers[!duplicated(markers$Annotation), ]
  labels <- groups$Annotation
  display_labels <- labels
  ncol <- 1L
  legend_size <- 5.8
  chip_size <- .065
  text_gap <- .025
  column_gap <- .10
  row_step <- .15
  title_space <- .14
  rows <- (seq_along(labels) - 1L) %/% ncol
  cols <- (seq_along(labels) - 1L) %% ncol

  text_width <- function(x) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(
        x,
        gp = grid::gpar(
          fontsize = legend_size, fontfamily = STYLE$family,
          lineheight = 1
        )
      )), "inches",
      valueOnly = TRUE
    )
  }
  col_widths <- vapply(0:(ncol - 1L), function(j) {
    x <- display_labels[cols == j]
    if (!length(x)) {
      return(0)
    }
    max(vapply(x, text_width, numeric(1)))
  }, numeric(1))
  column_lefts <- numeric(ncol)
  if (ncol > 1L) {
    for (j in 2:ncol) {
      column_lefts[j] <- column_lefts[j - 1L] + chip_size + text_gap +
        col_widths[j - 1L] + column_gap
    }
  }
  row_centers <- height - title_space - (rows * row_step + row_step / 2)
  p <- ggdraw() + draw_label("Marker annotation",
    x = 0, y = .99,
    hjust = 0, vjust = 1, size = 6.5, color = STYLE$ink
  )
  for (i in seq_along(labels)) {
    left <- column_lefts[cols[i] + 1L]
    center_y <- row_centers[i] / height
    p <- p +
      draw_grob(grid::rectGrob(
        x = grid::unit(left / width, "npc"),
        y = grid::unit(center_y, "npc"),
        width = grid::unit(chip_size / width, "npc"),
        height = grid::unit(chip_size / height, "npc"),
        just = c("left", "centre"),
        gp = grid::gpar(fill = groups$Color[i], col = NA)
      )) +
      draw_grob(grid::textGrob(
        display_labels[i],
        x = grid::unit((left + chip_size + text_gap) / width, "npc"),
        y = grid::unit(center_y, "npc") +
          grid::stringDescent(display_labels[i]) / 2,
        just = c("left", "centre"),
        gp = grid::gpar(
          fontsize = legend_size, fontfamily = STYLE$family,
          col = STYLE$ink, lineheight = 1
        )
      ))
  }
  p
}

make_marker_effects <- function(data, markers, settings, width, height) {
  n <- nrow(markers)
  targets <- settings$perturbations
  limit <- max(settings$color_limits)
  pretty_ticks <- pretty(c(-limit, limit), n = 5)
  tick_max <- max(pretty_ticks[pretty_ticks > 0 & pretty_ticks < limit])
  ticks <- c(-tick_max, 0, tick_max)
  p <- ggplot(data, aes(Column, Row)) +
    geom_tile(aes(fill = Overall_effect), width = 1, height = 1) +
    effect_scale(limit) +
    heatmap_frame(length(targets), n) +
    geom_rect(
      data = markers,
      aes(
        xmin = length(targets) + .63, xmax = length(targets) + .86,
        ymin = Row - .5, ymax = Row + .5
      ),
      inherit.aes = FALSE, fill = markers$Color
    ) +
    scale_x_continuous(
      limits = c(.5, length(targets) + .86),
      breaks = seq_along(targets),
      labels = targets, position = "top", expand = c(0, 0),
      sec.axis = dup_axis(name = "Perturbations", breaks = NULL)
    ) +
    scale_y_continuous(
      limits = c(.5, n + .5), breaks = 1:n,
      labels = rev(markers$Gene), position = "right", expand = c(0, 0),
      sec.axis = dup_axis(name = "Neuronal marker genes", breaks = NULL)
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    matrix_theme() +
    theme(
      axis.text.x.top = element_text(
        face = "italic",
        size = STYLE$marker_x_tick_size
      ),
      axis.text.y.right = element_text(
        size = STYLE$marker_tick_size,
        face = "italic"
      )
    )
  key_width <- .62
  annotation_height <- .78
  heat_bottom <- annotation_height + .12
  heat_height <- height - heat_bottom
  plot_width <- width - key_width - .05
  canvas <- place(
    ggdraw(), p, 0, heat_bottom, plot_width, heat_height,
    width, height
  )
  canvas <- place(
    canvas,
    effect_key(limit, ticks, "Overall\neffect",
      width = key_width,
      height = STYLE$effect_key_height, title_gap = .25
    ),
    width - key_width, heat_bottom +
      (heat_height - STYLE$effect_key_height) / 2,
    key_width, STYLE$effect_key_height, width, height
  )
  annotation_width <- width - .28
  place(
    canvas, marker_annotation_key(
      markers, annotation_width,
      annotation_height
    ),
    .23, .015, annotation_width, annotation_height, width, height
  )
}

make_factor_fdr_legend <- function(width, height) {
  bar_width <- .85
  bar_height <- STYLE$colorbar_thickness
  title_width <- .82
  group_width <- bar_width + .12 + title_width
  left <- min(
    width - group_width,
    max(0, .77 * width - group_width / 2)
  )
  bar_y <- .60 * height
  palette <- scales::gradient_n_pal(STYLE$enrichment)(
    seq(0, 1, length.out = 128)
  )
  legend <- ggdraw()
  for (i in seq_along(palette)) {
    legend <- legend + draw_grob(grid::rectGrob(
      x = grid::unit(left + (i - .5) * bar_width / length(palette), "in"),
      y = grid::unit(bar_y, "in"),
      width = grid::unit(bar_width / length(palette), "in"),
      height = grid::unit(bar_height, "in"),
      gp = grid::gpar(fill = palette[i], col = NA)
    ))
  }
  for (value in STYLE$factor_fdr_ticks) {
    x <- left + (value - min(STYLE$factor_fdr_limits)) /
      diff(STYLE$factor_fdr_limits) * bar_width
    tick_depth <- bar_height * .36
    legend <- legend + draw_grob(grid::segmentsGrob(
      x0 = grid::unit(rep(x, 2), "in"), x1 = grid::unit(rep(x, 2), "in"),
      y0 = grid::unit(c(
        bar_y - bar_height / 2,
        bar_y + bar_height / 2 - tick_depth
      ), "in"),
      y1 = grid::unit(c(
        bar_y - bar_height / 2 + tick_depth,
        bar_y + bar_height / 2
      ), "in"),
      gp = grid::gpar(col = "white", lwd = .65)
    )) + draw_grob(grid::textGrob(
      label = as.character(value),
      x = grid::unit(x, "in"),
      y = grid::unit(bar_y - bar_height / 2 - .045, "in"),
      just = "centre", gp = grid::gpar(
        fontsize = 5.5, col = STYLE$ink,
        fontfamily = STYLE$family
      )
    ))
  }
  legend + draw_grob(grid::textGrob(
    label = expression(-log[10](FDR)),
    x = grid::unit(left + bar_width + .12, "in"),
    y = grid::unit(bar_y, "in"), just = c("left", "centre"),
    gp = grid::gpar(
      fontsize = 6, col = STYLE$ink,
      fontfamily = STYLE$family
    )
  ))
}

make_factor_enrichment <- function(data, width, height) {
  factor_levels <- sort(unique(data$Factor))
  data$Factor <- factor(data$Factor,
    levels = factor_levels,
    labels = paste("Factor", factor_levels)
  )
  data <- data[order(
    data$Factor, -data$fold_enrichment, data$fdr,
    data$p_value, data$term_id
  ), , drop = FALSE]
  data$term_key <- paste(data$Factor, data$term_id, sep = "__")
  data$term_label <- vapply(data$description, function(label) {
    paste(strwrap(label, width = 28), collapse = "\n")
  }, character(1))
  data$term_key <- factor(data$term_key, levels = rev(unique(data$term_key)))
  label_map <- setNames(data$term_label, as.character(data$term_key))

  p <- ggplot(data, aes(fold_enrichment, term_key)) +
    geom_segment(
      aes(
        x = 0, xend = fold_enrichment, yend = term_key,
        color = Neg_log10_FDR
      ),
      linewidth = 3.8, lineend = "butt"
    ) +
    scale_x_continuous(
      limits = c(0, max(data$fold_enrichment, na.rm = TRUE) * 1.05),
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0, .02))
    ) +
    scale_y_discrete(labels = label_map) +
    scale_color_gradientn(
      colors = STYLE$enrichment,
      limits = STYLE$factor_fdr_limits,
      breaks = STYLE$factor_fdr_ticks,
      oob = scales::squish,
      guide = "none"
    ) +
    facet_grid(
      rows = vars(Factor), scales = "free_y", space = "free_y",
      switch = "y"
    ) +
    labs(x = "Fold enrichment", y = NULL) +
    theme_panel() +
    theme(
      axis.text.y = element_text(
        size = 5.5, lineheight = .9,
        color = STYLE$ink
      ),
      axis.text.x = element_text(size = 6, color = STYLE$ink),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        angle = 0, size = 6.2,
        color = STYLE$ink
      ),
      legend.position = "none",
      plot.margin = margin(3, 4, 3, 4)
    )
  legend_height <- .32
  legend <- make_factor_fdr_legend(width, legend_height)
  plot_grid(p, legend,
    ncol = 1,
    rel_heights = c(height - legend_height, legend_height)
  )
}

colorbar_key <- function(limits, ticks, title, colors, width, height,
                         orientation = "vertical", transform = identity,
                         title_size = STYLE$colorbar_title_size,
                         title_gap = .20) {
  horizontal <- identical(orientation, "horizontal")
  bar_length <- STYLE$colorbar_length
  thickness <- STYLE$colorbar_thickness
  bottom <- if (horizontal) .12 else .22
  palette <- scales::gradient_n_pal(colors)(seq(0, 1, length.out = 257))
  n <- length(palette)
  centers <- (seq_len(n) - .5) / n
  x <- thickness / (2 * width)
  y <- (bottom + centers * bar_length) / height
  dx <- thickness / width
  dy <- bar_length / (n * height)
  key <- ggdraw() + draw_grob(grid::rectGrob(
    x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
    width = grid::unit(dx, "npc"), height = grid::unit(dy, "npc"),
    gp = grid::gpar(fill = palette, col = NA)
  ))
  q <- (transform(ticks) - transform(limits[1])) / diff(transform(limits))
  tick_depth <- thickness * .36
  for (i in seq_along(ticks)) {
    tick_y <- (bottom + q[i] * bar_length) / height
    key <- key +
      draw_grob(grid::segmentsGrob(
        x0 = c(0, thickness) / width, x1 = c(tick_depth, thickness - tick_depth) / width,
        y0 = rep(tick_y, 2), y1 = rep(tick_y, 2),
        gp = grid::gpar(col = "white", lwd = .65)
      )) +
      draw_label(ticks[i],
        x = (thickness + .045) / width, y = tick_y,
        hjust = 0, vjust = .5, size = STYLE$colorbar_tick_size,
        color = STYLE$ink
      )
  }
  title_y <- if (horizontal) {
    .99
  } else {
    (bottom + bar_length + title_gap) / height
  }
  key + draw_label(title,
    x = 0, y = title_y, hjust = 0, vjust = 1,
    size = title_size, fontfamily = STYLE$family, color = STYLE$ink
  )
}

effect_key <- function(limit, ticks, title, width = .52,
                       height = STYLE$effect_key_height,
                       title_gap = .20) {
  colorbar_key(c(-limit, limit), ticks, title, STYLE$effect, width, height,
    orientation = "vertical", title_gap = title_gap
  )
}

# Count labels copied from the merged figure so crowded values are staggered
# while remaining positioned over their corresponding bars.
GeomStaggeredCountText <- ggproto("GeomStaggeredCountText", GeomText,
  draw_panel = function(data, panel_params, coord, na.rm = FALSE) {
    grid::gTree(
      data = coord$transform(data, panel_params), zero = data$y == 0,
      gap_pt = STYLE$count_label_gap_pt, cl = "staggered_count_labels"
    )
  }
)

makeContent.staggered_count_labels <- function(x) {
  d <- x$data
  panel_width <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  panel_height <- grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  gp <- grid::gpar(
    col = d$colour, fontsize = d$size * ggplot2::.pt,
    fontfamily = d$family, fontface = d$fontface, lineheight = d$lineheight
  )
  metrics <- vapply(seq_len(nrow(d)), function(i) {
    label <- grid::textGrob(d$label[i], gp = grid::gpar(
      fontsize = gp$fontsize[i], fontfamily = d$family[i], fontface = d$fontface[i]
    ))
    c(
      width = grid::convertWidth(grid::grobWidth(label), "inches", valueOnly = TRUE),
      height = grid::convertHeight(grid::grobHeight(label), "inches", valueOnly = TRUE)
    )
  }, numeric(2))
  label_width <- metrics["width", ] / panel_width
  label_height <- metrics["height", ] / panel_height
  gap_x <- x$gap_pt / (72.27 * panel_width)
  gap_y <- x$gap_pt / (72.27 * panel_height)
  bottom <- d$y + .45 * label_height
  placed <- integer()
  for (i in order(!x$zero, d$y, d$x)) {
    if (!x$zero[i]) {
      repeat {
        hits <- placed[
          abs(d$x[i] - d$x[placed]) < (label_width[i] + label_width[placed]) / 2 + gap_x &
            bottom[i] < bottom[placed] + label_height[placed] + gap_y &
            bottom[i] + label_height[i] + gap_y > bottom[placed]
        ]
        if (!length(hits)) break
        bottom[i] <- max(bottom[hits] + label_height[hits]) + gap_y
      }
    }
    placed <- c(placed, i)
  }
  grid::setChildren(x, grid::gList(grid::textGrob(
    d$label,
    x = d$x, y = bottom, hjust = .5, vjust = 0, gp = gp
  )))
}

# Panel C from merged_luhmes.R, using the two displayed methods.
make_count_bars <- function(data, methods, y_title, include_control, width, height) {
  if (!include_control) data <- data[data[[1]] != "NegCtrl", ]
  targets <- data[[1]]
  gene_col <- names(data)[1]
  data <- data[c(gene_col, methods)] |>
    pivot_longer(all_of(methods), names_to = "Method", values_to = "Count") |>
    mutate(
      Target = factor(.data[[gene_col]], levels = targets),
      Method = factor(Method, levels = methods)
    )
  multi <- grepl("DEGs", y_title, fixed = TRUE)
  dodge <- position_dodge(width = .70)
  max_count <- max(data$Count, na.rm = TRUE)
  top <- if (multi) 1600 else max(1, max_count) * 1.15
  y_breaks <- if (multi) {
    seq(0, 1500, 500)
  } else {
    scales::breaks_pretty(n = 4)(c(0, top))
  }
  count_labels <- if (multi) {
    layer(
      geom = GeomStaggeredCountText, stat = "identity", position = dodge,
      mapping = aes(label = scales::comma(Count)), show.legend = FALSE,
      params = list(
        size = 2.1, colour = STYLE$ink, family = STYLE$family,
        na.rm = FALSE
      )
    )
  } else {
    geom_text(aes(label = scales::comma(Count)),
      position = dodge,
      angle = 0, hjust = .5, vjust = -.45, size = 2.2, color = STYLE$ink
    )
  }
  ggplot(data, aes(Target, Count, fill = Method)) +
    geom_col(position = dodge, width = .70) +
    count_labels +
    scale_fill_manual(values = STYLE$methods, breaks = methods) +
    scale_x_discrete(expand = expansion(add = .55)) +
    scale_y_continuous(
      limits = c(0, top), breaks = y_breaks,
      labels = scales::label_comma(),
      expand = expansion(mult = c(STYLE$count_axis_lower_expansion, 0))
    ) +
    labs(x = "Perturbations", y = y_title) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_panel() +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1,
        face = "italic",
        size = if (multi) STYLE$deg_tick_size else STYLE$axis_tick_size
      ),
      axis.ticks.x = element_blank(), legend.position = "top",
      legend.title = element_blank(), legend.justification = "left",
      legend.key.width = grid::unit(3, "mm"),
      legend.key.height = grid::unit(3, "mm"),
      legend.margin = margin(0, 0, 2, 0),
      legend.spacing.x = grid::unit(2, "mm")
    )
}

main <- function() {
  canvas_width <- STYLE$width_cm / 2.54
  row_height <- STYLE$row_height_cm / 2.54
  marker_row_height_cm <- 12.5
  marker_row_height <- marker_row_height_cm / 2.54
  main_height <- 2 * row_height + marker_row_height
  slot_width <- canvas_width / 4

  # Main figure: A and B share the first row; C and D split the second row.
  panel_height <- row_height
  effect_matrix <- load_factor_effect_matrix()
  subset_settings <- list(
    selected_perturbations = setdiff(
      load_enriched_go_counts()$perturbation,
      c("HDAC5", "CTNND2", "DYRK1A", "ASH1L")
    ),
    selected_factor_numbers = c(2, 3, 4, 5, 7, 9, 10, 12, 13, 15, 17),
    color_limits_from_full_figure = c(-0.62, 0.62)
  )
  subset_matrix <- effect_matrix[
    subset_settings$selected_perturbations,
    paste0("factor_", subset_settings$selected_factor_numbers - 1L),
    drop = FALSE
  ]
  perturbation_order <- stats::hclust(
    stats::dist(subset_matrix),
    method = "ward.D2"
  )$order
  factor_order <- stats::hclust(
    stats::dist(t(subset_matrix)),
    method = "ward.D2"
  )$order
  subset_settings$selected_perturbations <-
    subset_settings$selected_perturbations[perturbation_order]
  subset_settings$selected_factor_numbers <-
    subset_settings$selected_factor_numbers[factor_order]
  subset_data <- factor_effect_data(
    effect_matrix, subset_settings$selected_perturbations,
    subset_settings$selected_factor_numbers
  )

  panel_a <- panel_frame(
    make_factor_counts(
      load_factor_counts(), slot_width,
      panel_height - STYLE$tag_band
    ),
    "A", slot_width, panel_height
  )
  panel_b_width <- canvas_width - slot_width
  panel_b <- panel_frame(
    make_count_bars(
      load_deg_counts(), c("GSFA", "PerturbVI"),
      "DEGs with LFSR < 0.05", FALSE,
      panel_b_width, panel_height - STYLE$tag_band
    ),
    "B", panel_b_width, panel_height
  )
  panel_c_width <- canvas_width / 2
  panel_c <- panel_frame(
    make_factor_effects(
      subset_data, subset_settings, panel_c_width,
      panel_height - STYLE$tag_band
    ),
    "C", panel_c_width, panel_height
  )
  panel_d <- panel_frame(
    make_count_bars(
      load_enriched_go_counts(), c("GSFA", "PerturbVI"),
      "Significant GO terms", FALSE,
      panel_c_width, panel_height - STYLE$tag_band
    ) +
      theme(legend.position = "none"),
    "D", panel_c_width, panel_height
  )
  marker_data <- load_marker_effect_data()
  panel_e_width <- canvas_width / 2
  panel_e <- panel_frame(
    make_marker_effects(
      marker_data$data, marker_data$markers,
      marker_data$settings, panel_e_width,
      marker_row_height - STYLE$tag_band
    ),
    "E", panel_e_width, marker_row_height
  )
  factor_enrichment_data <- load_new_factor_enrichment()
  panel_f_width <- canvas_width / 2
  panel_f <- panel_frame(
    make_factor_enrichment(
      factor_enrichment_data, panel_f_width,
      marker_row_height - STYLE$tag_band
    ),
    "F", panel_f_width, marker_row_height
  )

  figure <- place(
    ggdraw(), panel_a, 0, row_height + marker_row_height,
    slot_width, panel_height,
    canvas_width, main_height
  )
  figure <- place(
    figure, panel_b, slot_width, row_height + marker_row_height,
    panel_b_width, panel_height, canvas_width, main_height
  )
  figure <- place(
    figure, panel_c, 0, marker_row_height, panel_c_width, panel_height,
    canvas_width, main_height
  )
  figure <- place(
    figure, panel_d, panel_c_width, marker_row_height,
    panel_c_width, panel_height,
    canvas_width, main_height
  )
  figure <- place(
    figure, panel_e, 0, 0, panel_e_width, marker_row_height,
    canvas_width, main_height
  )
  figure <- place(
    figure, panel_f, panel_e_width, 0,
    panel_f_width, marker_row_height,
    canvas_width, main_height
  )
  ggsave(file.path(figures_dir, "figure_2.png"), figure,
    device = ragg::agg_png,
    width = STYLE$width_cm,
    height = 2 * STYLE$row_height_cm + marker_row_height_cm, units = "cm",
    dpi = STYLE$dpi, bg = "white"
  )

  message("Wrote ", file.path(figures_dir, "figure_2.png"))
}

main()
