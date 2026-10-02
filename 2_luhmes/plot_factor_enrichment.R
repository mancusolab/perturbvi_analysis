#!/usr/bin/env Rscript
# LUHMES neuronal factor enrichment plot.
#
# Run:
#   Rscript plot_factor_enrichment.R
#
# Outputs: figures/6_neuronal_factor_enrichment.png

library(ggplot2)
library(cowplot)
library(RColorBrewer)

enrichment <- read.csv("input/factor_enrichment.csv", stringsAsFactors = FALSE)

plot_factor_go <- function(enrichment) {
  data <- subset(enrichment, fdr < 0.05)
  pattern <- paste0("neuro|neural|nerve|nervous|axon|dendrit|brain|synap|glia|",
                    "myelin|cerebr|cerebell|hippocamp")
  data <- data[grepl(pattern, data$description, ignore.case = TRUE), ]
  if (!nrow(data)) stop("No enriched terms in this selection.")

  data$factor_number <- as.integer(sub("^factor_", "", data$group))
  data$term <- data$description
  data$score <- -log10(pmax(data$fdr, .Machine$double.xmin))

  factor_numbers <- sort(unique(data$factor_number))
  groups <- lapply(factor_numbers, function(f) {
    d <- data[data$factor_number == f, , drop = FALSE]
    d <- d[order(-d$fold_enrichment, d$fdr, d$p_value, d$term_id), , drop = FALSE]
    d <- d[!duplicated(d$term_id), , drop = FALSE]
    d$panel_title <- paste0("Factor ", f + 1L)
    d
  })
  names(groups) <- paste0("factor_", factor_numbers)

  n_panels <- length(groups)
  n_columns <- 2L

  column_width <- 99
  figure_width_mm <- column_width * n_columns

  font_size <- 7
  line_height <- 1.2
  ink <- "#342D38"

  temporary_image <- tempfile(fileext = ".png")
  ragg::agg_png(temporary_image, width = figure_width_mm, height = 100,
                units = "mm", res = 300)
  measurement_device <- grDevices::dev.cur()
  on.exit({
    grDevices::dev.off(measurement_device)
    unlink(temporary_image)
  }, add = TRUE)

  em <- grid::convertHeight(grid::unit(font_size, "pt"), "mm", valueOnly = TRUE)

  text_grob <- function(label) {
    grid::textGrob(label, gp = grid::gpar(
      fontsize = font_size, lineheight = line_height,
      fontfamily = "sans", col = ink))
  }
  text_width <- function(label) {
    grid::convertWidth(grid::grobWidth(text_grob(label)), "mm", valueOnly = TRUE)
  }
  text_height <- function(label) {
    grid::convertHeight(grid::grobHeight(text_grob(label)), "mm", valueOnly = TRUE)
  }

  label_budget <- 0.58 * (column_width - 4 * em)
  if (label_budget <= 2 * em) stop("Increase figure width.")

  wrap_text <- function(label, width_fun, budget) {
    words <- strsplit(stringr::str_squish(label), " ", fixed = TRUE)[[1]]
    lines <- character()
    current <- ""
    for (word in words) {
      candidate <- if (nzchar(current)) paste(current, word) else word
      if (width_fun(candidate) <= budget) {
        current <- candidate
        next
      }
      if (nzchar(current)) lines <- c(lines, current)
      current <- ""
      while (nchar(word) > 1 && width_fun(word) > budget) {
        lo <- 1L
        hi <- nchar(word)
        while (lo < hi) {
          mid <- ceiling((lo + hi) / 2)
          if (width_fun(substr(word, 1L, mid)) <= budget) lo <- mid else hi <- mid - 1L
        }
        lines <- c(lines, substr(word, 1L, lo))
        word <- substring(word, lo + 1L)
      }
      current <- word
    }
    paste(c(lines, current), collapse = "\n")
  }

  wrap_label <- function(label) wrap_text(label, text_width, label_budget)

  title_text_width <- function(label) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(
      label, gp = grid::gpar(
        fontsize = font_size + 1.5, lineheight = line_height,
        fontfamily = "sans", col = ink))), "mm", valueOnly = TRUE)
  }

  groups <- lapply(groups, function(d) {
    d$label <- vapply(d$term, wrap_label, character(1), USE.NAMES = FALSE)
    d$text_height <- vapply(d$label, text_height, numeric(1))
    d
  })

  tallest_label <- max(vapply(groups, function(d) max(d$text_height), numeric(1)))
  row_height <- max(2 * em, tallest_label + 0.5 * em)
  bar_height <- min(1.75 * em, 0.62 * row_height)
  factor_gap <- 2.5 * em

  reds <- RColorBrewer::brewer.pal(9, "Reds")[3:8]
  color_limits <- c(1, 3)  # fixed -log10(FDR) scale

  panel_theme <- theme_classic(base_size = font_size + 1.5, base_family = "sans") +
    theme(
      text = element_text(color = ink),
      plot.title = element_text(size = font_size + 1.5, lineheight = 1.25,
                                margin = margin(b = 1 * em, unit = "mm")),
      plot.title.position = "panel",
      axis.text.y = element_text(size = font_size, lineheight = line_height,
                                 vjust = 0.5, margin = margin(r = 0.6 * em, unit = "mm")),
      axis.text.x = element_text(size = font_size + 1,
                                 margin = margin(t = 0.6 * em, unit = "mm")),
      axis.line = element_line(color = "black", linewidth = 0.25),
      axis.ticks = element_line(color = "#77717C", linewidth = 0.22),
      axis.ticks.length = grid::unit(0.5 * em, "mm"),
      legend.position = "none",
      plot.margin = margin(0.6 * em, 1.2 * em, 0, 1.2 * em, unit = "mm")
    )

  make_panel <- function(d) {
    panel_height <- nrow(d) * row_height
    d$position <- (rev(seq_len(nrow(d))) - 0.5) * row_height

    plot <- ggplot(d) +
      geom_rect(aes(xmin = 0, xmax = fold_enrichment,
                    ymin = position - bar_height / 2,
                    ymax = position + bar_height / 2, fill = score)) +
      scale_x_continuous(limits = c(0, max(d$fold_enrichment) * 1.1),
                         breaks = scales::breaks_pretty(n = 3),
                         expand = expansion(mult = 0)) +
      scale_y_continuous(breaks = d$position, labels = d$label,
                         limits = c(0, panel_height), expand = expansion(mult = 0)) +
      scale_fill_gradientn(colors = reds, limits = color_limits,
                           oob = scales::squish, breaks = c(1, 2, 3),
                           name = expression(-log[10](FDR))) +
      labs(title = d$panel_title[1], x = NULL, y = NULL) +
      panel_theme

    drawing <- ggplotGrob(plot)
    panel_row <- drawing$layout$t[drawing$layout$name == "panel"]
    drawing$heights[panel_row] <- grid::unit(panel_height, "mm")

    list(plot = plot, drawing = drawing,
         height = grid::convertHeight(sum(drawing$heights), "mm", valueOnly = TRUE))
  }

  panels <- lapply(groups, make_panel)
  widths <- do.call(grid::unit.pmax, lapply(panels, function(p) p$drawing$widths))
  for (i in seq_along(panels)) panels[[i]]$drawing$widths <- widths

  drawing_widths <- panels[[1]]$drawing$widths
  fixed_width_mm <- sum(grid::convertWidth(
    drawing_widths[grid::unitType(drawing_widths) != "null"], "mm", valueOnly = TRUE))
  title_budget <- column_width - fixed_width_mm - 0.5 * em
  if (title_budget <= 2 * em) stop("Increase figure width for panel titles.")

  groups <- lapply(groups, function(d) {
    d$panel_title <- paste(wrap_text(d$panel_title[1], title_text_width, title_budget),
                           collapse = "\n")
    d
  })

  panels <- lapply(groups, make_panel)
  widths <- do.call(grid::unit.pmax, lapply(panels, function(p) p$drawing$widths))
  for (i in seq_along(panels)) panels[[i]]$drawing$widths <- widths
  heights <- vapply(panels, function(panel) panel$height, numeric(1))

  legend_box <- get_legend(panels[[1]]$plot + theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.key.width = grid::unit(3 * em, "mm"),
    legend.key.height = grid::unit(em, "mm"),
    legend.title = element_text(size = font_size + 0.5),
    legend.text = element_text(size = font_size),
    legend.margin = margin(0, 0, 0, 0)))
  guide_id <- which(legend_box$layout$name == "guides")[1]
  legend <- legend_box$grobs[[guide_id]]
  bar_cell <- legend$layout[legend$layout$name == "bar", , drop = FALSE]
  if (nrow(bar_cell) != 1L) stop("Could not locate the color strip in the legend.")

  title_idx <- which(grepl("title", legend$layout$name))
  if (length(title_idx) == 1L) {
    legend$layout$t[title_idx] <- bar_cell$t
    legend$layout$b[title_idx] <- bar_cell$b
  }

  legend <- gtable::gtable_add_rows(legend, grid::unit(2.6 * em, "mm"), pos = 0)
  legend <- gtable::gtable_add_grob(
    legend,
    grid::textGrob("Fold enrichment", x = 0.5,
                   gp = grid::gpar(fontsize = font_size + 2.5,
                                   fontfamily = "sans", col = ink)),
    t = 1, l = bar_cell$l, r = bar_cell$r, clip = "off",
    name = "fold-enrichment-title")
  legend_height <- grid::convertHeight(sum(legend$heights), "mm", valueOnly = TRUE)

  # pack panels into two columns to balance total height (tallest first),
  # then stack them tightly within each column so short panels don't get
  # padded to a shared row height.
  panel_gap <- 1.5 * em
  columns <- vector("list", n_columns)
  col_heights <- numeric(n_columns)
  for (i in order(-heights)) {
    col <- which.min(col_heights)
    if (length(columns[[col]])) col_heights[col] <- col_heights[col] + panel_gap
    columns[[col]] <- c(columns[[col]], i)
    col_heights[col] <- col_heights[col] + heights[i]
  }

  panel_col <- integer(n_panels)
  panel_top <- numeric(n_panels)
  for (col in seq_len(n_columns)) {
    cum <- 0
    for (i in columns[[col]]) {
      panel_col[i] <- col
      panel_top[i] <- cum
      cum <- cum + heights[i] + panel_gap
    }
  }

  figure_height <- max(col_heights)

  body <- ggdraw()
  for (i in seq_along(panels)) {
    body <- body + draw_grob(
      panels[[i]]$drawing,
      x = (panel_col[i] - 1) / n_columns,
      y = (figure_height - panel_top[i] - heights[i]) / figure_height,
      width = 1 / n_columns,
      height = heights[i] / figure_height)
  }

  legend_gap <- 1.5 * em
  figure <- plot_grid(
    body, NULL, ggdraw() + draw_grob(legend),
    ncol = 1,
    rel_heights = c(figure_height, legend_gap, legend_height)
  )

  list(plot = figure, width = figure_width_mm,
       height = figure_height + legend_gap + legend_height)
}

result <- plot_factor_go(enrichment)
ggsave("figures/6_neuronal_factor_enrichment.png", result$plot,
       width = result$width, height = result$height, units = "mm",
       dpi = 300, device = ragg::agg_png, bg = "white", limitsize = FALSE)
