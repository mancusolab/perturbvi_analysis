#!/usr/bin/env Rscript

# Merged LUHMES overview: seven panels from the patched no-control PerturbVI
# fit and the frozen WebGestalt annotation.
#
# Run:
#   Rscript merged_luhmes.R
#
# Outputs: figures/merged_luhmes.png

# 1. Configuration -------------------------------------------------------------
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(RColorBrewer)
library(ragg)

STYLE <- list(
  family = "sans", ink = "#342D38", muted = "#766C7A",
  base_size = 8, tag_size = 12, tag_band = .17,
  axis_title_size = 9, axis_tick_size = 7.5,
  axis_title_gap_pt = 6, axis_tick_gap_pt = 4,
  count_axis_lower_expansion = .03,
  count_label_gap_pt = .8,
  methods = setNames(RColorBrewer::brewer.pal(9, "Set1")[1:5],
                       c("scMAGeCK", "MAST", "DESeq2", "GSFA", "PerturbVI")),
  effect = grDevices::colorRampPalette(
    c("#2166AC", "#67A9CF", "#FFFFFF", "#EF8A62", "#B2182B"), space = "rgb")(257),
  enrichment = RColorBrewer::brewer.pal(9, "Reds")[3:8],
  go_comparison_colors = RColorBrewer::brewer.pal(9, "OrRd")[3:9],
  go_comparison_empty = "#F3F1F5",
  factor_bar_linewidth = 4.0, factor_title_size = 8,
  factor_color_limits = c(1, 3), factor_color_ticks = c(1, 2, 3),
  go_term_label_gap_pt = 4,
  heatmap_frame_color = "#555555", heatmap_frame_linewidth = .25,
  colorbar_length = 1.00, colorbar_thickness = .10,
  colorbar_title_size = 7, colorbar_tick_size = 6.5,
  effect_key_height = 1.40,
  width = 10, height = 14.14, dpi = 300, preview_dpi = 150
)

# Coordinates are inches from the top-left of the final canvas.
# Portrait, approximately A4 proportions: wide C; A/D, B, E on the left;
# compact factor-enrichment bars above the smaller GO-comparison tiles on the right.
# Titles are metadata for file names and the caption, never drawn on the figure.
LAYOUT <- tibble::tribble(
  ~tag, ~source, ~title, ~left, ~top, ~width, ~height,
  "A", "1",        "Genes per factor",                0.15, 3.31, 1.85, 2.10,
  "B", "2 subset", "Perturbation effects on factors", 0.15, 5.57, 4.64, 2.35,
  "C", "3 tools",  "DEGs across methods",             0.15, 0.15, 9.70, 3.00,
  "D", "4",        "Enriched GO terms",               2.16, 3.31, 2.63, 2.10,
  "E", "5 subset", "Effects on neuronal genes",       0.15, 8.08, 4.64, 5.91,
  "F", "6",        "Neuronal GO enrichment by factor",4.95, 3.31, 4.90, 4.80,
  "G", "7",        "GO enrichment across methods",    4.95, 8.27, 4.90, 5.40
)


load_data <- function() {
  read_mat <- function(path) as.matrix(read.csv(path, row.names = 1, check.names = FALSE))

  gsfa <- readRDS("input/gsfa_fit.light.rds")
  gsfa_pip <- gsfa$posterior_means$F_pm

  # A: genes per factor (PIP > 0.95)
  PIP_W <- read_mat("results/PIP_W.csv")
  n_factors <- nrow(PIP_W)
  a <- rbind(
    data.frame(Factor = seq_len(n_factors), Method = "GSFA",
               Count = unname(colSums(gsfa_pip > 0.95))),
    data.frame(Factor = seq_len(n_factors), Method = "PerturbVI",
               Count = unname(rowSums(PIP_W > 0.95)))
  )
  a$Method <- factor(a$Method, levels = c("GSFA", "PerturbVI"))

  # B: perturbation x factor effects (subset)
  B <- read_mat("results/B.csv")
  PIP_B <- read_mat("results/PIP_B.csv")
  b_settings <- list(
    selected_perturbations = c("ADNP", "ARID1B", "ASH1L", "CHD2", "PTEN", "SETD5"),
    selected_factor_numbers = c(2, 3, 4, 5, 7, 9, 10, 12, 13, 15, 17),
    color_limits_from_full_figure = c(-0.62, 0.62)
  )
  b_pert <- b_settings$selected_perturbations
  b_factors <- b_settings$selected_factor_numbers
  factor_cols <- paste0("factor_", b_factors - 1)
  b <- expand.grid(Perturbation = b_pert, Factor = b_factors, stringsAsFactors = FALSE)
  b$Effect <- B[cbind(b$Perturbation, factor_cols[match(b$Factor, b_factors)])]
  b$PIP <- PIP_B[cbind(b$Perturbation, factor_cols[match(b$Factor, b_factors)])]
  b$Column <- match(b$Factor, b_factors)
  b$Row <- length(b_pert) + 1L - match(b$Perturbation, b_pert)

  # C: DEGs across methods
  c <- read.csv("input/gsfa_author_reported_deg_counts.csv", check.names = FALSE)

  # D: enriched GO terms per perturbation
  d <- read.csv("input/deg_go_counts.csv", check.names = FALSE)
  d_omitted <- d[[1]][rowSums(d[c("GSFA", "PerturbVI")]) == 0]
  d <- d[rowSums(d[c("GSFA", "PerturbVI")]) > 0, ]

  # E: neuronal marker effects (subset)
  BW <- read_mat("results/BW.csv")
  LFSR_BW <- read_mat("results/LFSR_BW.csv")
  e_settings <- list(
    perturbations = c("ADNP", "ARID1B", "ASH1L", "CHD2", "PTEN", "SETD5"),
    color_limits = c(-0.34, 0.34)
  )
  markers <- read.csv("input/marker_annotations.csv", check.names = FALSE, stringsAsFactors = FALSE)
  marker_colors <- c(
    "neuron diff. / maturation" = "#66C2A5",
    "neg. reg. neuron differentiation" = "#FC8D62",
    "neg. reg. neuron projection dev." = "#8DA0CB",
    "neuron projection development" = "#E78AC3",
    "pos. reg. neuron differentiation" = "#A6D854",
    "pos. reg. neuron projection dev." = "#FFD92F",
    "reg. neuron differentiation" = "#E5C494"
  )
  markers$Color <- unname(marker_colors[markers$annotation])
  markers$Row <- nrow(markers) + 1L - seq_len(nrow(markers))
  e <- expand.grid(Gene = markers$gene_name, Perturbation = e_settings$perturbations,
                   stringsAsFactors = FALSE)
  ensembl <- markers$gene_ID[match(e$Gene, markers$gene_name)]
  e$Overall_effect <- BW[cbind(e$Perturbation, ensembl)]
  e$LFSR <- LFSR_BW[cbind(e$Perturbation, ensembl)]
  e$Column <- match(e$Perturbation, e_settings$perturbations)
  e$Row <- nrow(markers) + 1L - match(e$Gene, markers$gene_name)

  # F: neuronal factor enrichment
  fe <- read.csv("input/factor_enrichment.csv", stringsAsFactors = FALSE, check.names = FALSE)
  pattern <- paste0("neuro|neural|nerve|nervous|axon|dendrit|brain|synap|glia|",
                    "myelin|cerebr|cerebell|hippocamp")
  fe <- fe[grepl(pattern, fe$description, ignore.case = TRUE), ]
  fe$Factor <- as.integer(sub("^factor_", "", fe$group)) + 1L
  fe <- fe[order(fe$Factor, fe$fdr, fe$p_value, fe$term_id), ]
  fe$Display_rank <- ave(seq_len(nrow(fe)), fe$Factor, FUN = seq_along)
  f <- data.frame(
    Factor = fe$Factor,
    Description = fe$description,
    Fold_enrichment = fe$fold_enrichment,
    Neg_log10_FDR = -log10(pmax(fe$fdr, .Machine$double.xmin)),
    Display_rank = fe$Display_rank,
    stringsAsFactors = FALSE
  )
  f_settings <- list(color_limits = c(1, 3))

  # G: GO enrichment across methods
  ge <- read.csv("input/deg_go_enrichment.csv", stringsAsFactors = FALSE, check.names = FALSE)
  g_pert <- c("ADNP", "ARID1B", "ASH1L", "CHD2", "PTEN", "SETD5", "HDAC5")
  ge <- ge[ge$perturbation %in% g_pert, ]
  sig <- ge[ge$significant, ]
  stats <- do.call(rbind, lapply(split(sig, sig$term_id), function(z) {
    data.frame(term_id = z$term_id[1], description = z$description[1],
               recurrence = length(unique(z$perturbation)), min_fdr = min(z$fdr),
               stringsAsFactors = FALSE)
  }))
  rownames(stats) <- NULL
  stats$neuronal <- grepl(pattern, stats$description, ignore.case = TRUE)
  neuronal <- stats[stats$neuronal, ]
  broader <- stats[!stats$neuronal, ]
  broader <- broader[order(-broader$recurrence, broader$min_fdr, broader$term_id), ]
  broader <- head(broader, 28 - nrow(neuronal))
  terms <- rbind(neuronal, broader)
  terms <- terms[order(!terms$neuronal, -terms$recurrence, terms$min_fdr, terms$term_id), ]
  terms <- terms[c("term_id", "description")]
  names(terms) <- c("geneSet", "description")

  grid <- expand.grid(method = c("GSFA", "PerturbVI"), perturbation = g_pert,
                      term_id = terms$geneSet, stringsAsFactors = FALSE)
  grid <- merge(grid, ge[c("method", "perturbation", "term_id", "fold_enrichment", "fdr", "significant")],
                by = c("method", "perturbation", "term_id"), all.x = TRUE)
  section <- ifelse(grepl(pattern, terms$description[match(grid$term_id, terms$geneSet)],
                          ignore.case = TRUE), "Neuronal", "Broader")
  g <- data.frame(
    Method = factor(grid$method, levels = c("GSFA", "PerturbVI")),
    Perturbation = factor(grid$perturbation, levels = g_pert),
    geneSet = grid$term_id,
    enrichmentRatio = grid$fold_enrichment,
    Status = ifelse(!is.na(grid$significant) & grid$significant, "FDR < 0.05", "FDR >= 0.05"),
    Section = factor(section, levels = c("Neuronal", "Broader")),
    stringsAsFactors = FALSE
  )
  g$GO <- factor(g$geneSet, levels = rev(terms$geneSet))

  list(a = a, b = b, b_settings = b_settings, c = c, d = d, d_omitted = d_omitted,
       e = e, markers = markers, e_settings = e_settings,
       f = f, f_settings = f_settings, g = g, g_terms = terms)
}
# 3. Shared typography, legends, and panel placement ----------------------------
axis_typography <- function() {
  # Shared physical text sizes and gaps, including secondary heatmap axes.
  theme(axis.title = element_text(family = STYLE$family, size = STYLE$axis_title_size,
                                   face = "plain", color = STYLE$ink),
        axis.title.x.bottom = element_text(margin = margin(t = STYLE$axis_title_gap_pt)),
        axis.title.x.top = element_text(margin = margin(b = STYLE$axis_title_gap_pt)),
        axis.title.y.left = element_text(angle = 90, margin = margin(r = STYLE$axis_title_gap_pt)),
        axis.title.y.right = element_text(angle = 90, margin = margin(l = STYLE$axis_title_gap_pt)),
        axis.text = element_text(family = STYLE$family, size = STYLE$axis_tick_size,
                                  color = STYLE$ink),
        axis.text.x.bottom = element_text(margin = margin(t = STYLE$axis_tick_gap_pt)),
        axis.text.x.top = element_text(margin = margin(b = STYLE$axis_tick_gap_pt)),
        axis.text.y.left = element_text(margin = margin(r = STYLE$axis_tick_gap_pt)),
        axis.text.y.right = element_text(margin = margin(l = STYLE$axis_tick_gap_pt)))
}

theme_panel <- function() {
  theme_classic(base_size = STYLE$base_size, base_family = STYLE$family) +
    axis_typography() +
    theme(text = element_text(color = STYLE$ink),
          axis.line = element_line(color = "#55515A", linewidth = .25),
          axis.ticks = element_line(color = "#77717C", linewidth = .22),
          axis.ticks.length = grid::unit(1.2, "mm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.background = element_blank(), legend.key = element_blank(),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
          plot.margin = margin(4, 5, 4, 5))
}

place <- function(canvas, plot, left, bottom, width, height, canvas_width, canvas_height) {
  canvas + draw_plot(plot, x = left / canvas_width, y = bottom / canvas_height,
                     width = width / canvas_width, height = height / canvas_height)
}

panel_frame <- function(body, tag, width, height) {
  # Only the panel letter is drawn; descriptive titles belong in the caption.
  ggdraw() + draw_plot(body, x = 0, y = 0, width = 1,
                       height = (height - STYLE$tag_band) / height) +
    draw_label(tag, x = .025 / width, y = (height - .015) / height, hjust = 0, vjust = 1,
               size = STYLE$tag_size, fontface = "bold", color = STYLE$ink)
}

effect_scale <- function(limit) {
  scale_fill_gradientn(
    colors = STYLE$effect,
    limits = c(-limit, limit),
    guide = "none"
  )
}

colorbar_key <- function(limits, ticks, title, colors, width, height,
                         orientation = "horizontal", transform = identity,
                         title_size = STYLE$colorbar_title_size) {
  # One physical strip size and inward white tick style for B, E, F, and G.
  # The horizontal strip AND its title share the full panel midpoint.
  horizontal <- identical(orientation, "horizontal")
  bar_length <- STYLE$colorbar_length
  thickness <- STYLE$colorbar_thickness
  left <- if (horizontal) (width - bar_length) / 2 else 0
  bottom <- if (horizontal) .12 else .22
  palette <- scales::gradient_n_pal(colors)(seq(0, 1, length.out = 257))
  n <- length(palette)
  centers <- (seq_len(n) - .5) / n
  if (horizontal) {
    x <- (left + centers * bar_length) / width
    y <- (bottom + thickness / 2) / height
    dx <- bar_length / (n * width)
    dy <- thickness / height
  } else {
    x <- thickness / (2 * width)
    y <- (bottom + centers * bar_length) / height
    dx <- thickness / width
    dy <- bar_length / (n * height)
  }
  p <- ggdraw() + draw_grob(grid::rectGrob(
    x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
    width = grid::unit(dx, "npc"), height = grid::unit(dy, "npc"),
    gp = grid::gpar(fill = palette, col = NA)
  ))
  q <- (transform(ticks) - transform(limits[1])) / diff(transform(limits))
  stopifnot(all(q >= 0), all(q <= 1))
  tick_depth <- thickness * .36
  for (i in seq_along(ticks)) {
    if (horizontal) {
      tick_x <- (left + q[i] * bar_length) / width
      x0 <- x1 <- rep(tick_x, 2)
      y0 <- c(bottom, bottom + thickness) / height
      y1 <- c(bottom + tick_depth, bottom + thickness - tick_depth) / height
      label_x <- tick_x
      label_y <- .015 / height
    } else {
      tick_y <- (bottom + q[i] * bar_length) / height
      x0 <- c(0, thickness) / width
      x1 <- c(tick_depth, thickness - tick_depth) / width
      y0 <- y1 <- rep(tick_y, 2)
      label_x <- (thickness + .045) / width
      label_y <- tick_y
    }
    p <- p + draw_grob(grid::segmentsGrob(x0 = x0, x1 = x1, y0 = y0, y1 = y1,
                                           gp = grid::gpar(col = "white", lwd = .65))) +
      draw_label(ticks[i], x = label_x, y = label_y, hjust = if (horizontal) .5 else 0,
                   vjust = if (horizontal) 0 else .5,
                   size = STYLE$colorbar_tick_size, color = STYLE$ink)
  }
  p + draw_label(title, x = if (horizontal) .5 else 0, y = .99,
                   hjust = if (horizontal) .5 else 0, vjust = 1,
                   size = title_size, fontfamily = STYLE$family, color = STYLE$ink)
}

effect_key <- function(limit, ticks, title, width = .52,
                       height = STYLE$effect_key_height) {
  colorbar_key(
    c(-limit, limit), ticks, title, STYLE$effect,
    width, height,
    orientation = "vertical"
  )
}

marker_annotation_key <- function(markers, width, height = .78) {
  groups <- markers[!duplicated(markers$annotation), ]
  labels <- groups$annotation

  ncol <- 2L
  chip_size <- .08
  text_gap <- .03
  column_gap <- .30
  row_step <- .15
  title_space <- .17

  rows <- (seq_along(labels) - 1L) %/% ncol
  cols <- (seq_along(labels) - 1L) %% ncol

  text_width <- function(x) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(
        x,
        gp = grid::gpar(
          fontsize = 7,
          fontfamily = STYLE$family,
          lineheight = 1
        )
      )),
      "inches", valueOnly = TRUE
    )
  }

  # Dynamic column widths from actual rendered text
  col_widths <- vapply(0:(ncol - 1L), function(j) {
    x <- labels[cols == j]
    if (!length(x)) return(0)
    max(vapply(x, text_width, numeric(1)))
  }, numeric(1))

  # Column 2 begins after complete column 1 + fixed gap
  column_lefts <- numeric(ncol)

  for (j in 2:ncol) {
    column_lefts[j] <- column_lefts[j - 1L] +
      chip_size + text_gap + col_widths[j - 1L] + column_gap
  }

  row_centers <- height - title_space - (rows * row_step + row_step / 2)

  p <- ggdraw() +
    draw_label(
      "Marker annotation",
      x = 0, y = .99,
      hjust = 0, vjust = 1,
      size = 7.5,
      color = STYLE$ink
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
        labels[i],
        x = grid::unit((left + chip_size + text_gap) / width, "npc"),
        y = grid::unit(center_y, "npc") +
          grid::stringDescent(labels[i]) / 2,
        just = c("left", "centre"),
        gp = grid::gpar(
          fontsize = 7,
          fontfamily = STYLE$family,
          col = STYLE$ink,
          lineheight = 1
        )
      ))
  }

  p
}

heatmap_frame <- function(n_columns, n_rows) {
  annotate("rect", xmin = .5, xmax = n_columns + .5,
             ymin = .5, ymax = n_rows + .5, fill = NA,
             color = STYLE$heatmap_frame_color, linewidth = STYLE$heatmap_frame_linewidth)
}

factor_enrichment_key <- function(limits, width, height = .40) {
  ticks <- STYLE$factor_color_ticks

  # Physical dimensions
  label_width <- .72
  label_bar_gap <- .08
  bar_width <- STYLE$colorbar_length
  key_width <- label_width + label_bar_gap + bar_width

  # Colorbar without its own title
  bar <- colorbar_key(
    limits = limits,
    ticks = ticks,
    title = "",
    colors = STYLE$enrichment,
    width = bar_width,
    height = height
  )

  # FDR label immediately to the LEFT of the colorbar
  label <- ggdraw() +
    draw_label(
      expression(-log[10](FDR)),
      x = 1, y = .5,
      hjust = 1, vjust = .5,
      size = STYLE$colorbar_title_size,
      color = STYLE$ink
    )

  # One compact legend unit:
  # [-log10(FDR)] [gap] [colorbar]
  key <- plot_grid(
    label, NULL, bar,
    nrow = 1,
    rel_widths = c(label_width, label_bar_gap, bar_width)
  )

  ggdraw() +
    draw_plot(
      key,
      x = (width - key_width) / (2 * width),
      y = 0,
      width = key_width / width,
      height = 1
    )
}

matrix_theme <- function() {
  theme_void(base_size = STYLE$base_size, base_family = STYLE$family) +
    axis_typography() +
    theme(axis.text.x.top = element_text(angle = 60, hjust = 0, vjust = 0,
                                        margin = margin(b = STYLE$axis_tick_gap_pt)),
          axis.text.y.right = element_text(face = "italic", hjust = 0,
                                          margin = margin(l = STYLE$axis_tick_gap_pt)),
          plot.margin = margin(4, 4, 4, 4))
}

# 4. Native R panel builders ----------------------------------------------------

make_factor_counts <- function(data, width, height) {
  ggplot(data, aes(Method, Count, color = Method, fill = Method)) +
    geom_boxplot(width = .44, alpha = .19, linewidth = .5, outlier.shape = NA) +
    geom_point(position = position_jitter(width = .08, height = 0, seed = 42),
               shape = 21, size = 1.7, stroke = .25, color = "white") +
    scale_color_manual(values = STYLE$methods) + scale_fill_manual(values = STYLE$methods) +
    scale_y_continuous(limits = c(0, 1080), breaks = seq(0, 1000, 200),
                       labels = scales::label_comma(), expand = c(0, 0)) +
    labs(x = NULL, y = "Genes with PIP > 0.95") +
    theme_panel() + theme(legend.position = "none", axis.ticks.x = element_blank())
}

make_factor_effects <- function(data, settings, width, height) {
  factors <- settings$selected_factor_numbers
  targets <- settings$selected_perturbations
  limit <- max(settings$color_limits_from_full_figure)

  pretty_ticks <- pretty(c(-limit, limit), n = 5)
  tick_max <- max(pretty_ticks[pretty_ticks > 0 & pretty_ticks < limit])
  ticks <- c(-tick_max, 0, tick_max)

  p <- ggplot(data, aes(Column, Row, fill = Effect)) + geom_tile(width = 1, height = 1) +
    heatmap_frame(length(factors), length(targets)) +
    effect_scale(limit) +
    scale_x_continuous(breaks = seq_along(factors), labels = paste("Factor", factors),
                       position = "top", expand = c(0, 0),
                       sec.axis = dup_axis(name = "Factors", breaks = NULL)) +
    scale_y_continuous(breaks = seq_along(targets), labels = rev(targets),
                       position = "right", expand = c(0, 0),
                       sec.axis = dup_axis(name = "Perturbations", breaks = NULL)) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + matrix_theme()
  key_width <- .52
  plot_width <- width - key_width - .05
  canvas <- place(ggdraw(), p, 0, 0, plot_width, height, width, height)
  key_height <- min(STYLE$effect_key_height, height)
  canvas <- place(canvas, effect_key(limit, ticks, "Effect size",
           width = key_width, height = key_height),
                   width - key_width, (height - key_height) / 2, key_width, key_height, width, height)
  canvas
}

# measure and stagger labels when their panel is actually drawn, so collision
# spacing follows the font metrics and panel dimensions of each PNG/PDF device.
GeomStaggeredCountText <- ggproto("GeomStaggeredCountText", GeomText,
  draw_panel = function(data, panel_params, coord, na.rm = FALSE) {
    grid::gTree(data = coord$transform(data, panel_params), zero = data$y == 0,
                gap_pt = STYLE$count_label_gap_pt, cl = "staggered_count_labels")
  }
)

makeContent.staggered_count_labels <- function(x) {
  d <- x$data
  panel_width <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  panel_height <- grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  gp <- grid::gpar(col = d$colour, fontsize = d$size * ggplot2::.pt,
                   fontfamily = d$family, fontface = d$fontface, lineheight = d$lineheight)
  metrics <- vapply(seq_len(nrow(d)), function(i) {
    label <- grid::textGrob(d$label[i], gp = grid::gpar(
      fontsize = gp$fontsize[i], fontfamily = d$family[i], fontface = d$fontface[i]))
    c(width = grid::convertWidth(grid::grobWidth(label), "inches", valueOnly = TRUE),
      height = grid::convertHeight(grid::grobHeight(label), "inches", valueOnly = TRUE))
  }, numeric(2))
  label_width <- metrics["width", ] / panel_width
  label_height <- metrics["height", ] / panel_height
  gap_x <- x$gap_pt / (72.27 * panel_width)
  gap_y <- x$gap_pt / (72.27 * panel_height)
  bottom <- d$y + .45 * label_height
  placed <- integer()
  # Reserve all zero labels at their common level, then place other labels
  # from low to high. Only upward movement is allowed; x stays over its bar.
  for (i in order(!x$zero, d$y, d$x)) {
    if (!x$zero[i]) repeat {
      hits <- placed[
        abs(d$x[i] - d$x[placed]) < (label_width[i] + label_width[placed]) / 2 + gap_x &
          bottom[i] < bottom[placed] + label_height[placed] + gap_y &
          bottom[i] + label_height[i] + gap_y > bottom[placed]]
      if (!length(hits)) break
      bottom[i] <- max(bottom[hits] + label_height[hits]) + gap_y
    }
    placed <- c(placed, i)
  }
  grid::setChildren(x, grid::gList(grid::textGrob(
    d$label, x = d$x, y = bottom, hjust = .5, vjust = 0, gp = gp)))
}

make_count_bars <- function(data, methods, y_title, include_control, width, height) {
  if (!include_control) data <- data[data[[1]] != "NegCtrl", ]
  targets <- data[[1]]
  gene_col <- names(data)[1]
  data <- data[c(gene_col, methods)] |>
    pivot_longer(all_of(methods), names_to = "Method", values_to = "Count") |>
    mutate(Target = factor(.data[[gene_col]], levels = targets),
           Method = factor(Method, levels = methods))
  multi <- length(methods) > 2
  # Fill each gene's 0.88-wide block; leave a small gap between gene groups.
  dodge <- position_dodge(width = if (multi) .88 else .70)
  top <- if (multi) 1600 else 60
  count_labels <- if (multi) {
    layer(geom = GeomStaggeredCountText, stat = "identity", position = dodge,
            mapping = aes(label = scales::comma(Count)), show.legend = FALSE,
            params = list(size = 2.1, colour = STYLE$ink, family = STYLE$family, na.rm = FALSE))
  } else {
    geom_text(aes(label = scales::comma(Count)), position = dodge,
                 angle = 0, hjust = .5, vjust = -.45, size = 2.2, color = STYLE$ink)
  }
  p <- ggplot(data, aes(Target, Count, fill = Method)) +
    geom_col(position = dodge, width = if (multi) .88 else .70) +
    count_labels +
    scale_fill_manual(values = STYLE$methods, breaks = methods) +
    # Ordinary italic text matches the heatmap gene labels; plotmath expressions
    # introduce mathematical character spacing into alphanumeric gene symbols.
    scale_x_discrete(expand = expansion(add = .55)) +
    # Space below zero lifts the full data baseline and its tick off the x-axis.
    scale_y_continuous(limits = c(0, top), breaks = if (multi) seq(0, 1500, 500) else seq(0, 60, 20),
                       labels = scales::label_comma(),
                       expand = expansion(mult = c(STYLE$count_axis_lower_expansion, 0))) +
    labs(x = "Perturbations", y = y_title) +
    guides(fill = guide_legend(nrow = 1)) + theme_panel() +
    theme(axis.text.x = element_text(angle = if (multi) 0 else 45,
                                     hjust = if (multi) .5 else 1,
                                     face = "italic", size = STYLE$axis_tick_size),
          axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank(),
          legend.justification = "left", legend.key.width = grid::unit(3.5, "mm"),
          legend.key.height = grid::unit(2.5, "mm"), legend.margin = margin(0, 0, 2, 0),
          legend.spacing.x = grid::unit(2, "mm"))
  p
}

make_marker_effects <- function(data, markers, settings, width, height) {
  n <- nrow(markers)
  targets <- settings$perturbations
  limit <- max(settings$color_limits)

  pretty_ticks <- pretty(c(-limit, limit), n = 5)
  tick_max <- max(pretty_ticks[pretty_ticks > 0 & pretty_ticks < limit])
  ticks <- c(-tick_max, 0, tick_max)

  p <- ggplot(data, aes(Column, Row)) + geom_tile(aes(fill = Overall_effect), width = 1, height = 1) +
    effect_scale(limit) +
    heatmap_frame(length(targets), n) +
    # Fixed annotation colors do not share the signed-effect fill scale.
    geom_rect(data = markers, aes(xmin = 6.63, xmax = 6.86, ymin = Row - .5, ymax = Row + .5),
               inherit.aes = FALSE, fill = markers$Color) +
    scale_x_continuous(limits = c(.5, 6.86), breaks = 1:6, labels = targets,
                       position = "top", expand = c(0, 0),
                       sec.axis = dup_axis(name = "Perturbations", breaks = NULL)) +
    scale_y_continuous(limits = c(.5, n + .5), breaks = 1:n, labels = rev(markers$gene_name),
                       position = "right", expand = c(0, 0),
                       sec.axis = dup_axis(name = "Neuronal marker genes", breaks = NULL)) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + matrix_theme() +
    theme(axis.text.x.top = element_text(face = "italic"))
  key_width <- .62
  key_height <- STYLE$effect_key_height
  annotation_height <- .78
  heat_bottom <- annotation_height + .12
  heat_height <- height - heat_bottom
  plot_width <- width - key_width - .05
  canvas <- place(ggdraw(), p, 0, heat_bottom, plot_width, heat_height, width, height)
  canvas <- place(canvas, effect_key(limit, ticks, "Effect size",
           width = key_width, height = key_height),
                   width - key_width, heat_bottom + (heat_height - key_height) / 2,
                   key_width, key_height, width, height)
  annotation_width <- width - .28
  canvas <- place(canvas, marker_annotation_key(markers, annotation_width, annotation_height),
                   .23, .015, annotation_width, annotation_height, width, height)
  canvas
}

factor_grid_layout <- function(data, width) {
  n_columns <- if (width >= 8) 4L else 2L
  ids <- unique(data$Factor)
  weights <- vapply(ids, function(id) {
    .24 + .14 * sum(data$Factor == id)
  }, numeric(1))
  order_ids <- order(-weights, ids)
  packed <- data.frame(Factor = ids[order_ids], Weight = weights[order_ids], Column = 0L, Panel = 0L)
  column_heights <- numeric(n_columns)
  column_counts <- integer(n_columns)
  for (i in seq_len(nrow(packed))) {
    column <- which.min(round(column_heights, 10))
    packed$Column[i] <- column
    column_counts[column] <- column_counts[column] + 1L
    packed$Panel[i] <- column_counts[column]
    column_heights[column] <- column_heights[column] + packed$Weight[i]
  }
  packed[order(packed$Column, packed$Panel), ]
}

make_factor_go <- function(data, settings, width, height) {
  factors <- unique(data$Factor)
  wrap_width <- 26
  packing <- factor_grid_layout(data, width)

  font_size <- 7
  line_height <- 0.90
  em <- grid::convertHeight(grid::unit(font_size, "pt"), "inches", valueOnly = TRUE)

  plots <- setNames(lapply(factors, function(id) {
    d <- data[data$Factor == id, ]
    d <- d[order(d$Display_rank), ]

    labels <- stringr::str_wrap(d$Description, width = wrap_width)
    d$Label <- factor(labels, levels = rev(labels))
    x_max <- max(d$Fold_enrichment) * 1.10

    ggplot(d, aes(Fold_enrichment, Label)) +
      geom_segment(
        aes(x = 0, xend = Fold_enrichment, yend = Label, color = Neg_log10_FDR),
        linewidth = STYLE$factor_bar_linewidth, lineend = "butt"
      ) +
      scale_x_continuous(
        limits = c(0, x_max),
        breaks = scales::breaks_pretty(n = 3),
        expand = expansion(mult = 0)
      ) +
      scale_y_discrete(expand = expansion(add = .65)) +
      scale_color_gradientn(
        colors = STYLE$enrichment,
        limits = STYLE$factor_color_limits,
        breaks = STYLE$factor_color_ticks,
        name = expression(-log[10](FDR))
      ) +
      labs(title = paste("Factor", id), x = NULL, y = NULL) +
      theme_panel() +
      theme(
        plot.title = element_text(
          size = STYLE$factor_title_size, face = "plain",
          lineheight = 1.0, margin = margin(b = 2)
        ),
        plot.title.position = "panel",
        axis.text.y = element_text(
          size = font_size, lineheight = line_height,
          margin = margin(r = 0.6 * em, unit = "in")
        ),
        axis.text.x = element_text(
          size = STYLE$axis_tick_size, color = STYLE$ink
        ),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color = "black", linewidth = .25),
        legend.position = "none",
        plot.margin = margin(3, 3, 4, 4)
      )
  }), as.character(factors))

  # -------------------------------------------------------------------------
  # Preserve existing 2-column packing
  # -------------------------------------------------------------------------

  columns <- lapply(sort(unique(packing$Column)), function(column) {
    group <- packing[packing$Column == column, ]
    ids <- group$Factor

    plot_grid(
      plotlist = lapply(
        plots[as.character(ids)],
        function(p) p + theme(legend.position = "none")
      ),
      ncol = 1,
      rel_heights = group$Weight,
      align = "v",
      axis = "lr"
    )
  })

  body <- plot_grid(
    plotlist = columns,
    ncol = length(columns)
  )

  # -------------------------------------------------------------------------
  # Native ggplot legend
  # -------------------------------------------------------------------------

  legend_bar_width <- .2

  legend_source <- plots[[1]] +
    guides(
      color = guide_colorbar(title.position = "left")
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = grid::unit(legend_bar_width, "in"),
      legend.key.height = grid::unit(em, "in"),
      legend.title = element_text(
        size = font_size + 0.5,
        margin = margin(r = 4)
      ),
      legend.text = element_text(size = font_size),
      legend.margin = margin(0, 0, 0, 0)
    )

  legend_box <- get_legend(legend_source)

  guide_id <- which(legend_box$layout$name == "guides")[1]
  legend <- legend_box$grobs[[guide_id]]

  # Locate actual colorbar
  bar_cell <- legend$layout[
    legend$layout$name == "bar",
    , drop = FALSE
  ]

  if (nrow(bar_cell) != 1L) {
    stop("Could not locate the -log10(FDR) color strip in the Fig. F legend.")
  }

  # Locate native -log10(FDR) title
  title_idx <- which(grepl("title", legend$layout$name))

  if (length(title_idx) == 1L) {
    # Force title into exactly the same vertical row as the colorbar
    legend$layout$t[title_idx] <- bar_cell$t
    legend$layout$b[title_idx] <- bar_cell$b
  }

  # -------------------------------------------------------------------------
  # Add "Fold enrichment" above the ACTUAL colorbar
  # -------------------------------------------------------------------------

  legend <- gtable::gtable_add_rows(
    legend,
    grid::unit(2.6 * em, "in"),
    pos = 0
  )

  legend <- gtable::gtable_add_grob(
    legend,
    grid::textGrob(
      "Fold enrichment",
      x = .5,
      gp = grid::gpar(
        fontsize = STYLE$axis_title_size,
        fontfamily = STYLE$family,
        col = STYLE$ink
      )
    ),
    t = 1,
    l = bar_cell$l,
    r = bar_cell$r,
    clip = "off",
    name = "fold-enrichment-title"
  )

  # -------------------------------------------------------------------------
  # Compact body -> legend spacing
  # -------------------------------------------------------------------------

  legend_height <- grid::convertHeight(
    sum(legend$heights), "inches", valueOnly = TRUE
  )

  legend_gap <- .08

  plot_grid(
    body,
    NULL,
    ggdraw() + draw_grob(legend),
    ncol = 1,
    rel_heights = c(
      height - legend_height - legend_gap,
      legend_gap,
      legend_height
    )
  )
}

make_go_comparison <- function(data, terms, width, height) {
  significant <- data[data$Status == "FDR < 0.05", ]

  fill_limits <- c(1, 11)
  ticks <- c(2, 6, 10)

  # -------------------------------------------------------------------------
  # GO comparison grid
  # -------------------------------------------------------------------------

  p <- ggplot(data, aes(Perturbation, GO)) +
    geom_tile(fill = STYLE$go_comparison_empty, width = .94, height = .92) +
    geom_tile(
      data = significant, aes(fill = enrichmentRatio),
      width = .94, height = .92
    ) +
    facet_grid(
      rows = vars(Section), cols = vars(Method),
      scales = "free_y", space = "free_y"
    ) +
    scale_x_discrete(position = "top", expand = expansion(add = .5)) +
    scale_y_discrete(
      labels = setNames(terms$description, terms$geneSet),
      expand = expansion(add = .55)
    ) +
    scale_fill_gradientn(
      colors = STYLE$go_comparison_colors,
      limits = fill_limits,
      guide = "none"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8, base_family = STYLE$family) +
    axis_typography() +
    theme(
      panel.grid = element_blank(),
      panel.spacing.x = grid::unit(1.5, "mm"),
      panel.spacing.y = grid::unit(1.5, "mm"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y = element_blank(),
      strip.text.x = element_text(
        size = 9.5, face = "plain", color = STYLE$ink,
        margin = margin(b = 7)
      ),
      axis.text.x.top = element_text(
        face = "italic", angle = 55, hjust = 0, vjust = 0,
        margin = margin(b = STYLE$axis_tick_gap_pt)
      ),
      axis.text.y = element_text(
        size = 7.2, color = STYLE$ink,
        margin = margin(r = STYLE$go_term_label_gap_pt)
      ),
      axis.ticks = element_blank(),
      plot.margin = margin(6, 5, STYLE$axis_title_gap_pt, 15),
      legend.position = "none"
    )

  # -------------------------------------------------------------------------
  # Fold-enrichment legend
  # -------------------------------------------------------------------------

  legend_height <- .52
  bar_width <- .90
  bar_height <- STYLE$colorbar_thickness
  title_bar_gap <- .06

  # Physical region occupied by the actual grid columns
  grid_left <- 1.55
  grid_right <- .05

  grid_width <- width - grid_left - grid_right
  grid_center <- grid_left + grid_width / 2

  # Center colorbar on the grid region
  bar_left <- grid_center - bar_width / 2

  palette <- scales::gradient_n_pal(
    STYLE$go_comparison_colors
  )(seq(0, 1, length.out = 257))

  bar_bottom <- .12
  bar_top <- bar_bottom + bar_height
  title_y <- bar_top + title_bar_gap

  n <- length(palette)
  centers <- (seq_len(n) - .5) / n

  legend_plot <- ggdraw() +
    draw_label(
      "Fold enrichment",
      x = grid_center / width,
      y = title_y / legend_height,
      hjust = .5, vjust = 0,
      size = STYLE$axis_title_size,
      color = STYLE$ink
    ) +
    draw_grob(grid::rectGrob(
      x = grid::unit(
        (bar_left + centers * bar_width) / width, "npc"
      ),
      y = grid::unit(
        (bar_bottom + bar_height / 2) / legend_height, "npc"
      ),
      width = grid::unit(
        bar_width / (n * width), "npc"
      ),
      height = grid::unit(
        bar_height / legend_height, "npc"
      ),
      gp = grid::gpar(fill = palette, col = NA)
    ))

  # -------------------------------------------------------------------------
  # Inward white ticks: 1, 5, 10
  # -------------------------------------------------------------------------

  q <- (ticks - fill_limits[1]) / diff(fill_limits)
  tick_depth <- bar_height * .36

  for (i in seq_along(ticks)) {
    tick_x <- (bar_left + q[i] * bar_width) / width

    legend_plot <- legend_plot +
      draw_grob(grid::segmentsGrob(
        x0 = grid::unit(rep(tick_x, 2), "npc"),
        x1 = grid::unit(rep(tick_x, 2), "npc"),
        y0 = grid::unit(
          c(bar_bottom, bar_top) / legend_height, "npc"
        ),
        y1 = grid::unit(
          c(
            bar_bottom + tick_depth,
            bar_top - tick_depth
          ) / legend_height,
          "npc"
        ),
        gp = grid::gpar(col = "white", lwd = .65)
      )) +
      draw_label(
        ticks[i],
        x = tick_x,
        y = .015 / legend_height,
        hjust = .5, vjust = 0,
        size = STYLE$colorbar_tick_size,
        color = STYLE$ink
      )
  }

  # -------------------------------------------------------------------------
  # Grid + legend
  # -------------------------------------------------------------------------

  plot_grid(
    p,
    legend_plot,
    ncol = 1,
    rel_heights = c(height - legend_height, legend_height)
  )
}

# 5. Assemble, export, and record provenance ------------------------------------

build_panels <- function(data) {
  builders <- list(
    A = function(w, h) make_factor_counts(data$a, w, h),
    B = function(w, h) make_factor_effects(data$b, data$b_settings, w, h),
    C = function(w, h) make_count_bars(data$c, names(STYLE$methods), "Number of DEGs", FALSE, w, h),
    D = function(w, h) make_count_bars(data$d, c("GSFA", "PerturbVI"),
                                      "Enriched GO terms", FALSE, w, h),
    E = function(w, h) make_marker_effects(data$e, data$markers, data$e_settings, w, h),
    F = function(w, h) make_factor_go(data$f, data$f_settings, w, h),
    G = function(w, h) make_go_comparison(data$g, data$g_terms, w, h)
  )
  setNames(lapply(seq_len(nrow(LAYOUT)), function(i) {
    spec <- LAYOUT[i, ]
    body <- builders[[spec$tag]](spec$width, spec$height - STYLE$tag_band)
    panel_frame(body, spec$tag, spec$width, spec$height)
  }), LAYOUT$tag)
}

save_plots <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, device = ragg::agg_png, width = width, height = height,
           units = "in", dpi = STYLE$dpi, bg = "white")
  # ggsave(paste0(stem, ".pdf"), plot, device = grDevices::cairo_pdf, width = width, height = height,
  #          units = "in", bg = "white")
}

main <- function() {
  grDevices::pdf(NULL)  # dummy device
  data <- load_data()
  panels <- build_panels(data)
  plate <- ggdraw()
  for (i in seq_len(nrow(LAYOUT))) {
    spec <- LAYOUT[i, ]
    plate <- place(plate, panels[[spec$tag]], spec$left,
                   STYLE$height - spec$top - spec$height,
                   spec$width, spec$height, STYLE$width, STYLE$height)
  }
  save_plots(plate, "figures/merged_luhmes", STYLE$width, STYLE$height)
  invisible(grDevices::dev.off())
}

main()
