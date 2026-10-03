#!/usr/bin/env Rscript
# Replogle PerturbVI fit: theme enrichment and marker-gene heatmaps -> merged figure.
#
# Run: Rscript merged_replogle.R
# Inputs (results/): BW.csv
# Inputs (input/): replogle_factors_enrichment.csv, replogle_theme_membership.csv,
#                  top25_annotation.csv, K562_essential_downstream_gene.tsv
# Outputs (figures/): merged_replogle.png

library(tidyverse)
library(patchwork)
library(cowplot)

pdf(NULL)

RESULTS <- "results"
INPUT <- "input"
FIGURES <- "figures"
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)

# Layout tuning knobs ------------------------------------------------------

font_size          <- 7           # pt, theme-enrichment panels
line_height        <- 1.2
ink                <- "#342D38"
column_width       <- 99          # mm, one theme-enrichment column
fig3_height        <- 6.05        # in, marker-gene heatmap row
left_margin        <- .04
right_margin       <- .08
heatmap_gap        <- .20
legend_left_inset  <- .22
effect_gap         <- .4
bottom_margin      <- .03
heatmap_legend_gap <- .14
panel_gap_mm       <- 2


# Panel D markers: ribosome
ribosome_perturbations <- c(
  "SNRPF",
  "SNRPD3",
  "SNRPD2",
  "SNRPE",
  "SNRPD1",
  "UTP6",
  "IMP4",
  "UTP3",
  "RRP12",
  "NOL6",
  "PDCD11",
  "DDX47",
  "TSR2",
  "RPS7",
  "RPS24",
  "RPS19",
  "RPS6",
  "RPS8"
)
ribosome_annot <- list(
  "rRNA metabolic process" = c(
    "METTL5",
    "DKC1",
    "PA2G4",
    "EBNA1BP2",
    "EXOSC8",
    "LYAR",
    "GTF3A",
    "DIMT1",
    "NOL7",
    "NOP56"
  ),
  "protein-RNA complex organization" = c("PRKDC", "PRMT5"),
  "translational initiation" = c("EIF3E", "EIF3F", "EIF3L", "EIF4B", "RPL13A", "RPS5"),
  "cytoplasmic translation" = c(
    "RPL35A",
    "RPL10A",
    "RPL26",
    "RPL7",
    "RPL14",
    "RPL24",
    "RPLP0",
    "RPL23A",
    "RPL6",
    "RPL5",
    "RPL10",
    "RPL3",
    "RPS16",
    "RPS13",
    "RPS21",
    "RPS15",
    "RPS28",
    "RPS27",
    "RPSA",
    "RPS19",
    "RPS24",
    "RPS14"
  )
)

# Panel C markers: cell cycle
cell_cycle_perturbations <- c(
  "TAF1",
  "TAF10",
  "MED1",
  "TAF2",
  "RRM1",
  "POLE",
  "ERCC3",
  "INTS7",
  "E4F1",
  "NLE1",
  "EIF4E",
  "RPL24",
  "CDC5L",
  "RPS6",
  "GSPT1",
  "ATF5",
  "CHMP3",
  "PHB2"
)
cell_cycle_annot <- list(
  "kinetochore organization" = c("CENPF", "CENPE", "NDC80", "CENPA", "DLGAP5", "NUF2", "SMC4"),
  "chromosome condensation" = c("TOP2A", "PLK1", "CDK1", "NUSAP1", "NCAPG"),
  "cell cycle G2/M phase transition" = c("CCNB1", "AURKA", "KIF14", "CCNA2", "CCNB2", "AURKB"),
  "meiotic cell cycle" = c("ASPM", "CDC20", "CKS2", "SGO2", "KNL1", "TTK", "NEK2"),
  "cell cycle G1/S phase transition" = c("GTSE1"),
  "positive regulation of cell cycle" = c("CDCA8", "KIF23", "BUB1", "BIRC5", "KIF20B", "RACGAP1"),
  "negative regulation of cell cycle" = c("CCNF", "SPC25", "BUB1B"),
  "mitotic cell cycle phase transition" = c("UBE2S", "TACC3", "CKS1B"),
  "regulation of mitotic cell cycle" = c("MKI67", "CDCA2")
)


theme_keywords <- list(
  "cell_cycle" = c(
    "cell[- ]cycle",
    "mitotic",
    "chromosome condensation",
    "kinetochore"
  ),
  "mitochondria" = c(
    "mitochond",
    "energy derivation",
    "electron transport",
    "proton.*transport",
    "NADH",
    "cytochrome"
  ),
  "ribosome" = c(
    "ribosom",
    "\\brRNA\\b",
    "\\btranslation(al)?\\b",
    "ribonucleoprotein",
    "protein[- ]RNA complex",
    "\\bncRNA\\b"
  )
)

theme_order <- c(
  "cell_cycle",
  "mitochondria",
  "ribosome",
  "signaling",
  "transport",
  "metabolism",
  "immune",
  "stress",
  "death"
)

# ---------------------------------------------------------------------------
# Panel A: sum of squared effects (p)
# ---------------------------------------------------------------------------

bw <- read_csv(file.path(RESULTS, "BW.csv")) |>
  rename(perturbation_id = 1) |>
  column_to_rownames("perturbation_id")

top <- bw |>
  rownames_to_column("perturbation_id") |>
  mutate(ssq = rowSums(across(-perturbation_id, ~ .x^2))) |>
  select(perturbation_id, ssq) |>
  arrange(desc(ssq)) |>
  slice_head(n = 25)

CATEGORIES <- read_csv(file.path(INPUT, "top25_annotation.csv"))
assignment <- setNames(CATEGORIES$category, CATEGORIES$perturbation_id)

title_case <- function(label) {
  sapply(strsplit(label, " "), function(w) {
    paste(ifelse(w == toupper(w), w, str_to_title(w)), collapse = " ")
  })
}

top <- top |>
  mutate(
    category = assignment[perturbation_id],
    perturbation_id = factor(perturbation_id, levels = perturbation_id),
    x = seq_len(n())
  )

marker_colors <- setNames(RColorBrewer::brewer.pal(5, "Set1"), unique(top$category))

top <- top |>
  mutate(category = factor(category, levels = names(marker_colors)))

bar_width <- 0.42
x_spacing <- 0.62

p <- ggplot(top, aes(x = x * x_spacing, y = ssq, fill = category)) +
  geom_col(width = bar_width, linewidth = 0) +
  scale_x_continuous(
    breaks = top$x * x_spacing,
    labels = top$perturbation_id,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_fill_manual(
    values = marker_colors,
    name = "Marker annotation",
    breaks = names(marker_colors),
    labels = title_case
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(y = "Sum of squared effects") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 8.5, base_family = "sans") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = ink),
    axis.text.y = element_text(size = 7, color = ink),
    axis.title.y = element_text(size = 8.5, color = ink),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.title = element_text(size = 8.5, color = ink, margin = margin(t = 6, b = 4)),
    legend.text = element_text(size = 8, color = ink),
    legend.key.size = unit(0.38, "cm"),
    legend.key.spacing.y = unit(0.07, "cm"),
    legend.key.spacing.x = unit(0.25, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(1, 1, 1, 1)
  )

# ---------------------------------------------------------------------------
# Panel B: theme enrichment panels (p2)
# ---------------------------------------------------------------------------

input_file      <- file.path(INPUT, "replogle_factors_enrichment.csv")
membership_file <- file.path(INPUT, "replogle_theme_membership.csv")
membership <- read.csv(membership_file, stringsAsFactors = FALSE)
theme_ids <- setNames(strsplit(membership$term_ids, ";", fixed = TRUE), membership$theme_key)
theme_labels <- c(
  setNames(membership$theme_label, membership$theme_key),
  cell_cycle   = "Cell Cycle",
  mitochondria = "Mitochondria",
  ribosome     = "Ribosome"
)

plot_theme_enrichment <- function(theme_keywords,
                                  input_file,
                                  factors = NULL,
                                  theme_ids = NULL,
                                  theme_labels = NULL) {
  if (length(factors)) {
    factors <- lapply(factors, function(ids) {
      if (is.numeric(ids)) ids <- paste0("factor_", as.integer(ids))
      unique(ids)
    })
  }
  
  enrichment <- read.csv(input_file, stringsAsFactors = FALSE)
  enrichment <- enrichment[is.finite(enrichment$fdr) &
                             enrichment$fdr >= 0 & enrichment$fdr < 0.05 &
                             is.finite(enrichment$fold_enrichment) &
                             enrichment$fold_enrichment > 1, , drop = FALSE]
  
  selected_terms <- do.call(rbind, lapply(theme_order, function(theme) {
    if (theme %in% names(theme_ids)) {
      keep <- enrichment$term_id %in% theme_ids[[theme]]
    } else {
      words <- theme_keywords[[theme]]
      keep <- rep(FALSE, nrow(enrichment))
      if (length(words)) {
        pattern <- regex(paste(words, collapse = "|"), ignore_case = TRUE)
        keep <- !is.na(enrichment$term) &
          str_detect(enrichment$term, pattern)
      }
    }
    if (theme %in% names(factors)) {
      keep <- keep & enrichment$factor %in% factors[[theme]]
    }
    d <- enrichment[keep, , drop = FALSE]
    d$theme <- rep(theme, nrow(d))
    d
  }))
  
  enrichment <- selected_terms
  enrichment$score <- -log10(pmax(enrichment$fdr, .Machine$double.xmin))
  
  groups <- list()
  for (theme in theme_order) {
    d <- enrichment[enrichment$theme == theme, , drop = FALSE]
    if (!nrow(d))
      next
    
    if (theme %in% names(factors)) {
      ids <- factors[[theme]]
    } else {
      pairs <- unique(d[c("factor", "term_id")])
      counts <- table(pairs$factor)
      numbers <- as.integer(sub("^factor_", "", names(counts)))
      ids <- names(counts)[order(-as.integer(counts), numbers, names(counts))[1]]
    }
    
    for (i in seq_along(ids)) {
      panel <- d[d$factor == ids[i], , drop = FALSE]
      if (!nrow(panel))
        next
      panel <- panel[order(-panel$fold_enrichment,
                           panel$fdr,
                           panel$p_value,
                           panel$term_id), , drop = FALSE]
      panel <- panel[!duplicated(panel$term_id), , drop = FALSE]
      panel <- utils::head(panel, 5)
      number <- as.integer(sub("^factor_", "", ids[i])) + 1L
      panel$panel_title <- if (theme %in% names(theme_labels)) theme_labels[[theme]] else theme
      panel$factor_label <- paste0("(Factor ", number, ")")
      panel$factor_position <- i
      groups[[paste(theme, ids[i], sep = "::")]] <- panel
    }
  }
  
  make_figure <- function(groups) {
    themes <- vapply(groups, function(d)
      d$theme[1], character(1))
    theme_index <- match(themes, unique(themes))
    theme_block <- (theme_index - 1L) %/% 3L + 1L
    factor_rows <- vapply(groups, function(d)
      d$factor_position[1], integer(1))
    block_depths <- vapply(seq_len(max(theme_block)), function(block) {
      max(factor_rows[theme_block == block])
    }, integer(1))
    block_offsets <- c(0L, head(cumsum(block_depths), -1L))
    panel_rows <- factor_rows + block_offsets[theme_block]
    panel_rows <- match(panel_rows, sort(unique(panel_rows)))
    panel_columns <- (theme_index - 1L) %% 3L + 1L
    n_columns <- min(3L, length(unique(themes)))
    n_rows <- max(panel_rows)
    figure_width_mm <- column_width * n_columns
    
    temporary_image <- tempfile(fileext = ".png")
    ragg::agg_png(
      temporary_image,
      width = figure_width_mm,
      height = 100,
      units = "mm",
      res = 300
    )
    measurement_device <- grDevices::dev.cur()
    on.exit({
      grDevices::dev.off(measurement_device)
      unlink(temporary_image)
    }, add = TRUE)
    
    em <- grid::convertHeight(grid::unit(font_size, "pt"), "mm", valueOnly = TRUE)
    
    text_grob <- function(label) {
      grid::textGrob(
        label,
        gp = grid::gpar(
          fontsize = font_size,
          lineheight = line_height,
          fontfamily = "sans",
          col = ink
        )
      )
    }
    text_width <- function(label) {
      grid::convertWidth(grid::grobWidth(text_grob(label)), "mm", valueOnly = TRUE)
    }
    text_height <- function(label) {
      grid::convertHeight(grid::grobHeight(text_grob(label)), "mm", valueOnly = TRUE)
    }
    
    label_budget <- 0.58 * (column_width - 4 * em)
    
    wrap_text <- function(label, width_fun, budget) {
      words <- strsplit(str_squish(label), " ", fixed = TRUE)[[1]]
      lines <- character()
      current <- ""
      
      for (word in words) {
        candidate <- if (nzchar(current))
          paste(current, word)
        else
          word
        if (width_fun(candidate) <= budget) {
          current <- candidate
          next
        }
        
        if (nzchar(current))
          lines <- c(lines, current)
        current <- ""
        
        while (nchar(word) > 1 && width_fun(word) > budget) {
          lo <- 1L
          hi <- nchar(word)
          while (lo < hi) {
            mid <- ceiling((lo + hi) / 2)
            if (width_fun(substr(word, 1L, mid)) <= budget) {
              lo <- mid
            } else {
              hi <- mid - 1L
            }
          }
          lines <- c(lines, substr(word, 1L, lo))
          word <- substring(word, lo + 1L)
        }
        current <- word
      }
      paste(c(lines, current), collapse = "\n")
    }
    
    wrap_label <- function(label) {
      wrap_text(label, text_width, label_budget)
    }
    
    title_text_width <- function(label) {
      grid::convertWidth(grid::grobWidth(grid::textGrob(
        label,
        gp = grid::gpar(
          fontsize = font_size + 1.5,
          lineheight = line_height,
          fontfamily = "sans",
          col = ink
        )
      )), "mm", valueOnly = TRUE)
    }
    
    groups <- lapply(groups, function(d) {
      labels <- d$term
      d$label <- vapply(labels, wrap_label, character(1), USE.NAMES = FALSE)
      d$text_height <- vapply(d$label, text_height, numeric(1))
      d
    })
    
    tallest_label <- max(vapply(groups, function(d)
      max(d$text_height), numeric(1)))
    row_height <- max(2 * em, tallest_label + 0.5 * em)
    bar_height <- min(1.75 * em, 0.62 * row_height)
    factor_gap <- 2.5 * em
    
    reds <- RColorBrewer::brewer.pal(9, "Reds")[3:7]
    color_limits <- range(unlist(lapply(groups, function(d)
      d$score)), finite = TRUE)
    
    panel_theme <- theme_classic(base_size = font_size + 1.5, base_family = "sans") +
      theme(
        text = element_text(color = ink),
        plot.title = element_text(
          size = font_size + 1.5,
          lineheight = 1.25,
          margin = margin(b = 1 * em, unit = "mm")
        ),
        plot.title.position = "panel",
        axis.text.y = element_text(
          size = font_size,
          lineheight = line_height,
          vjust = 0.5,
          margin = margin(r = 0.6 * em, unit = "mm")
        ),
        axis.text.x = element_text(
          size = font_size + 1,
          margin = margin(t = 0.6 * em, unit = "mm")
        ),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(color = "#77717C", linewidth = 0.22),
        axis.ticks.length = grid::unit(0.5 * em, "mm"),
        legend.position = "none",
        plot.margin = margin(
          0.6 * em, 1.2 * em, 0, 1.2 * em,
          unit = "mm"
        )
      )
    
    make_panel <- function(d) {
      panel_height <- nrow(d) * row_height
      d$position <- (rev(seq_len(nrow(d))) - 0.5) * row_height
      
      plot <- ggplot(d) +
        geom_rect(
          aes(
            xmin = 0,
            xmax = fold_enrichment,
            ymin = position - bar_height / 2,
            ymax = position + bar_height / 2,
            fill = score
          )
        ) +
        scale_x_continuous(
          limits = c(0, max(d$fold_enrichment) * 1.1),
          breaks = scales::breaks_pretty(n = 3),
          expand = expansion(mult = 0)
        ) +
        scale_y_continuous(
          breaks = d$position,
          labels = d$label,
          limits = c(0, panel_height),
          expand = expansion(mult = 0)
        ) +
        scale_fill_gradientn(
          colors = reds,
          limits = color_limits,
          breaks = scales::breaks_pretty(n = 4),
          name = expression(-log[10](FDR))
        ) +
        labs(title = d$panel_title[1],
             x = NULL,
             y = NULL) +
        panel_theme
      
      drawing <- ggplotGrob(plot)
      panel_row <- drawing$layout$t[drawing$layout$name == "panel"]
      drawing$heights[panel_row] <- grid::unit(panel_height, "mm")
      
      list(
        plot = plot,
        drawing = drawing,
        height = grid::convertHeight(sum(drawing$heights), "mm", valueOnly = TRUE)
      )
    }
    
    panels <- lapply(groups, make_panel)
    widths <- do.call(grid::unit.pmax,
                      lapply(panels, function(p)
                        p$drawing$widths))
    for (i in seq_along(panels))
      panels[[i]]$drawing$widths <- widths
    
    drawing_widths <- panels[[1]]$drawing$widths
    fixed_width_mm <- sum(grid::convertWidth(drawing_widths[grid::unitType(drawing_widths) != "null"], "mm", valueOnly = TRUE))
    title_budget <- column_width - fixed_width_mm - 0.5 * em
    
    groups <- lapply(groups, function(d) {
      label_lines <- strsplit(wrap_text(d$panel_title[1], title_text_width, title_budget),
                              "\n",
                              fixed = TRUE)[[1]]
      candidate <- paste0(label_lines[length(label_lines)], " ", d$factor_label[1])
      if (title_text_width(candidate) <= title_budget) {
        label_lines[length(label_lines)] <- candidate
      } else {
        label_lines <- c(label_lines, d$factor_label[1])
      }
      d$panel_title <- paste(label_lines, collapse = "\n")
      d
    })
    
    panels <- lapply(groups, make_panel)
    widths <- do.call(grid::unit.pmax,
                      lapply(panels, function(p)
                        p$drawing$widths))
    for (i in seq_along(panels))
      panels[[i]]$drawing$widths <- widths
    heights <- vapply(panels, function(panel)
      panel$height, numeric(1))
    
    legend_box <- get_legend(
      panels[[1]]$plot + theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = grid::unit(4 * em, "mm"),
        legend.key.height = grid::unit(em, "mm"),
        legend.title = element_text(size = font_size + 0.5),
        legend.text = element_text(size = font_size),
        legend.margin = margin(0, 0, 0, 0)
      )
    )
    guide_id <- which(legend_box$layout$name == "guides")[1]
    legend <- legend_box$grobs[[guide_id]]
    bar_cell <- legend$layout[legend$layout$name == "bar", , drop = FALSE]
    if (nrow(bar_cell) != 1L)
      stop("Could not locate the color strip in the legend.")
    
    # Vertically center -log10(FDR) with the colorbar
    title_idx <- which(grepl("title", legend$layout$name))
    
    if (length(title_idx) == 1L) {
      legend$layout$t[title_idx] <- bar_cell$t
      legend$layout$b[title_idx] <- bar_cell$b
    }
    
    legend <- gtable::gtable_add_rows(legend, grid::unit(2.6 * em, "mm"), pos = 0)
    legend <- gtable::gtable_add_grob(
      legend,
      grid::textGrob(
        "Fold enrichment",
        x = 0.5,
        gp = grid::gpar(
          fontsize = font_size + 2.5,
          fontfamily = "sans",
          col = ink
        )
      ),
      t = 1,
      l = bar_cell$l,
      r = bar_cell$r,
      clip = "off",
      name = "fold-enrichment-title"
    )
    
    legend_height <- grid::convertHeight(sum(legend$heights), "mm", valueOnly = TRUE)
    grid_row_heights <- vapply(seq_len(n_rows), function(row) {
      max(heights[panel_rows == row])
    }, numeric(1))
    row_gap <- 1.5 * em
    row_offsets <- c(0, cumsum(head(grid_row_heights, -1L) + row_gap))
    figure_height <- sum(grid_row_heights) +
      (n_rows - 1L) * row_gap +
      factor_gap + legend_height
    figure <- ggdraw()
    
    for (i in seq_along(panels)) {
      figure <- figure + draw_grob(
        panels[[i]]$drawing,
        x = (panel_columns[i] - 1) / n_columns,
        y = (figure_height - row_offsets[panel_rows[i]] - heights[i]) / figure_height,
        width = 1 / n_columns,
        height = heights[i] / figure_height
      )
    }
    figure <- figure + draw_grob(
      legend, y = 0, height = legend_height / figure_height
    )
    
    list(
      plot = figure,
      width = figure_width_mm,
      height = figure_height,
      title_budget_mm = title_budget,
      panel_layout = data.frame(
        theme = themes,
        factor = vapply(groups, function(d)
          d$factor[1], character(1)),
        title = vapply(groups, function(d)
          d$panel_title[1], character(1)),
        row = panel_rows,
        column = panel_columns,
        n_terms = vapply(groups, nrow, integer(1)),
        row.names = NULL
      )
    )
  }
  
  result <- make_figure(groups)
  invisible(result)
}

result <- plot_theme_enrichment(
  input_file = input_file,
  theme_keywords = theme_keywords,
  theme_ids = theme_ids,
  theme_labels = theme_labels,
  factors = list(
    "cell_cycle" = 3,
    "mitochondria" = 15,
    "ribosome" = 1,
    "death" = 8
  )
)

p2 <- result$plot

writeLines(
  sprintf(
    "%-14s %-10s %2d terms",
    result$panel_layout$theme,
    result$panel_layout$factor,
    result$panel_layout$n_terms
  )
)

# ---------------------------------------------------------------------------
# Panels C / D: marker-gene heatmaps
# ---------------------------------------------------------------------------

gene_meta <- read_tsv(file.path(INPUT, "K562_essential_downstream_gene.tsv"))
by_symbol <- gene_meta |>
  distinct(gene_name, .keep_all = TRUE) |>
  select(gene_name, gene_id) |>
  deframe()

make_marker_annotations <- function(annot) {
  tibble(
    gene_name = unlist(annot, use.names = FALSE),
    annotation = rep(names(annot), lengths(annot))
  ) |>
    mutate(gene_ID = by_symbol[gene_name]) |>
    select(gene_ID, gene_name, annotation)
}

ribosome_df     <- make_marker_annotations(ribosome_annot)
cell_cycle_df   <- make_marker_annotations(cell_cycle_annot)

# Marker-annotation colours, named per category so the mapping does not depend
# on gene order. tab10 is matplotlib's default qualitative palette.
tab10 <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
           "#8C564B", "#E377C2", "#17BECF", "#BCBD22")

ribosome_colors     <- setNames(RColorBrewer::brewer.pal(length(ribosome_annot), "Set2"),
                                names(ribosome_annot))
cell_cycle_colors   <- setNames(tab10[seq_along(cell_cycle_annot)], names(cell_cycle_annot))

# ---------------------------------------------------------------------------
# Heatmap helpers
# ---------------------------------------------------------------------------

marker_annotation_key <- function(categories, palette, height) {
  ncol <- 2L
  chip_size <- .085
  text_gap <- .03
  column_gap <- .3
  row_step <- .16
  title_space <- .19
  
  rows <- (seq_along(categories) - 1L) %/% ncol
  cols <- (seq_along(categories) - 1L) %% ncol
  
  text_width <- function(x) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(
        x, gp = grid::gpar(fontsize = 8, fontfamily = "sans")
      )),
      "inches", valueOnly = TRUE
    )
  }
  
  # Dynamic width of each column from its longest rendered label.
  col_widths <- vapply(0:(ncol - 1L), function(j) {
    labels <- categories[cols == j]
    if (!length(labels)) return(0)
    max(vapply(labels, text_width, numeric(1)))
  }, numeric(1))
  
  # Dynamic column positions with a fixed physical gap.
  column_lefts <- numeric(ncol)
  for (j in 2:ncol) {
    column_lefts[j] <- column_lefts[j - 1L] +
      chip_size + text_gap + col_widths[j - 1L] + column_gap
  }
  
  content_width <- column_lefts[ncol] + chip_size + text_gap + col_widths[ncol]
  row_centers <- height - title_space - (rows * row_step + row_step / 2)
  
  p <- ggdraw() +
    draw_label(
      "Marker annotation", x = 0, y = .99,
      hjust = 0, vjust = 1, size = 8.5, color = ink
    )
  
  for (i in seq_along(categories)) {
    left <- column_lefts[cols[i] + 1L]
    y <- row_centers[i] / height
    label <- categories[i]
    
    p <- p +
      draw_grob(grid::rectGrob(
        x = grid::unit(left / content_width, "npc"),
        y = grid::unit(y, "npc"),
        width = grid::unit(chip_size / content_width, "npc"),
        height = grid::unit(chip_size / height, "npc"),
        just = c("left", "centre"),
        gp = grid::gpar(fill = palette[label], col = NA)
      )) +
      draw_grob(grid::textGrob(
        label,
        x = grid::unit((left + chip_size + text_gap) / content_width, "npc"),
        y = grid::unit(y, "npc") + grid::stringDescent(label) / 2,
        just = c("left", "centre"),
        gp = grid::gpar(
          fontsize = 8, fontfamily = "sans",
          col = ink, lineheight = 1
        )
      ))
  }
  
  list(plot = p, width = content_width)
}

effect_size_key <- function(limit, height) {
  colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  palette <- scales::gradient_n_pal(colors)(seq(0, 1, length.out = 257))
  
  bar_length <- 1.00
  bar_height <- .10
  bar_left <- .05
  key_width <- bar_length + .10
  
  title_top <- height - .025
  title_height <- .12
  title_bar_gap <- .06
  bar_top <- title_top - title_height - title_bar_gap
  bar_bottom <- bar_top - bar_height
  tick_y <- bar_bottom - .03
  
  # Nice symmetric ticks INSIDE the actual scale endpoints.
  pretty_vals <- pretty(c(-limit, limit), n = 5)
  pretty_vals <- pretty_vals[abs(pretty_vals) < limit]
  
  # Keep zero plus the outermost rounded value on either side.
  tick_max <- max(pretty_vals[pretty_vals > 0])
  ticks <- c(-tick_max, 0, tick_max)
  
  n <- length(palette)
  centers <- (seq_len(n) - .5) / n
  
  p <- ggdraw() +
    draw_label(
      "Effect size",
      x = bar_left / key_width, y = title_top / height,
      hjust = 0, vjust = 1, size = 8.5, color = ink
    ) +
    draw_grob(grid::rectGrob(
      x = grid::unit((bar_left + centers * bar_length) / key_width, "npc"),
      y = grid::unit((bar_bottom + bar_height / 2) / height, "npc"),
      width = grid::unit(bar_length / (n * key_width), "npc"),
      height = grid::unit(bar_height / height, "npc"),
      gp = grid::gpar(fill = palette, col = NA)
    ))
  
  q <- (ticks + limit) / (2 * limit)
  tick_depth <- bar_height * .36
  
  for (i in seq_along(ticks)) {
    tick_x <- (bar_left + q[i] * bar_length) / key_width
    
    # White ticks pointing inward from BOTH long edges.
    p <- p +
      draw_grob(grid::segmentsGrob(
        x0 = grid::unit(rep(tick_x, 2), "npc"),
        x1 = grid::unit(rep(tick_x, 2), "npc"),
        y0 = grid::unit(c(bar_bottom, bar_top) / height, "npc"),
        y1 = grid::unit(c(
          bar_bottom + tick_depth,
          bar_top - tick_depth
        ) / height, "npc"),
        gp = grid::gpar(col = "white", lwd = .65)
      )) +
      draw_label(
        format(ticks[i], trim = TRUE, scientific = FALSE),
        x = tick_x, y = tick_y / height,
        hjust = .5, vjust = 1, size = 7.5, color = ink
      )
  }
  
  list(plot = p, width = key_width)
}

annotation_key_height <- function(n_categories) {
  .21 + ceiling(n_categories / 2) * .16
}


# ---------------------------------------------------------------------------
# Heatmap builder
# ---------------------------------------------------------------------------

make_gene_heatmap <- function(bw, gene_df, perturbations, colors) {
  n_genes <- nrow(gene_df)
  n_pert <- length(perturbations)
  gene_df <- gene_df |> mutate(Row = n_genes - row_number() + 1)
  
  mat <- as.data.frame(bw[perturbations, gene_df$gene_ID, drop = FALSE]) |>
    rownames_to_column("perturbation_id") |>
    pivot_longer(-perturbation_id, names_to = "gene_ID", values_to = "effect") |>
    left_join(gene_df, by = "gene_ID") |>
    mutate(Column = match(perturbation_id, perturbations))
  
  categories <- unique(gene_df$annotation)
  
  palette <- colors[categories]
  limit <- max(abs(mat$effect))
  gene_df <- gene_df |> mutate(strip_color = palette[annotation])
  
  body <- ggplot(mat, aes(Column, Row)) +
    geom_tile(aes(fill = effect), width = 1, height = 1) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
      limits = c(-limit, limit), guide = "none"
    ) +
    
    # Thin annotation strip.
    geom_rect(
      data = gene_df,
      aes(
        xmin = n_pert + .63, xmax = n_pert + .86,
        ymin = Row - .5, ymax = Row + .5
      ),
      inherit.aes = FALSE, fill = gene_df$strip_color
    ) +
    
    annotate(
      "rect",
      xmin = .5, xmax = n_pert + .5,
      ymin = .5, ymax = n_genes + .5,
      fill = NA, color = "#555555", linewidth = .25
    ) +
    
    # Perturbation names on top; title on bottom.
    scale_x_continuous(
      limits = c(.5, n_pert + .86),
      breaks = seq_len(n_pert), labels = perturbations,
      position = "top", expand = c(0, 0),
      sec.axis = dup_axis(name = "Perturbations", breaks = NULL)
    ) +
    
    # Marker-gene names on right; title on left.
    scale_y_continuous(
      limits = c(.5, n_genes + .5),
      breaks = gene_df$Row, labels = gene_df$gene_name,
      position = "right", expand = c(0, 0),
      sec.axis = dup_axis(name = "Marker genes", breaks = NULL)
    ) +
    
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 7, base_family = "sans") +
    theme(
      axis.text.x.top = element_text(
        angle = 60, hjust = 0, vjust = 0,
        face = "italic", size = 7, color = ink,
        margin = margin(b = 4)
      ),
      axis.text.y.right = element_text(
        face = "italic", size = 7, color = ink,
        hjust = 0, margin = margin(l = 4)
      ),
      
      # Individual titles with a little breathing room from matrix.
      axis.title.x.bottom = element_text(
        size = 8.5, color = ink,
        margin = margin(t = 6)
      ),
      axis.title.y.left = element_text(
        size = 8.5, color = ink, angle = 90,
        margin = margin(r = 6)
      ),
      
      plot.margin = margin(t = 0, r = 4, b = 2, l = 0)
    )
  
  list(
    top = body,
    categories = categories,
    palette = palette,
    limit = limit
  )
}


p4 <- make_gene_heatmap(
  bw, cell_cycle_df, cell_cycle_perturbations, cell_cycle_colors
)

p5 <- make_gene_heatmap(
  bw, ribosome_df, ribosome_perturbations, ribosome_colors
)


place_panel <- function(canvas, plot, left, bottom, width, height,
                        canvas_width, canvas_height) {
  canvas + draw_plot(
    plot,
    x = left / canvas_width, y = bottom / canvas_height,
    width = width / canvas_width, height = height / canvas_height
  )
}

fig3_width <- result$width / 25.4

panel_w <- (fig3_width - left_margin - right_margin - heatmap_gap) / 2

x1 <- left_margin
x2 <- x1 + panel_w + heatmap_gap

legend_x1 <- x1 + legend_left_inset
legend_x2 <- x2 + legend_left_inset

ann1_h <- annotation_key_height(length(p4$categories))
ann2_h <- annotation_key_height(length(p5$categories))
legend_h <- max(ann1_h, ann2_h)

ann1 <- marker_annotation_key(
  p4$categories, p4$palette, legend_h
)

ann2 <- marker_annotation_key(
  p5$categories, p5$palette, legend_h
)

eff1 <- effect_size_key(
  p4$limit, legend_h
)

eff2 <- effect_size_key(
  p5$limit, legend_h
)

eff1_x <- legend_x1 + ann1$width + effect_gap
eff2_x <- legend_x2 + ann2$width + effect_gap

legend_bottom <- bottom_margin
legend_top <- legend_bottom + legend_h

heat_bottom <- legend_top + heatmap_legend_gap
heat_h <- fig3_height - heat_bottom

# each heatmap is its own patch so the enclosing patchwork can tag them
# separately (C, D). They are split through the middle of heatmap_gap so the
# two halves abut seamlessly.
split_x <- x2 - heatmap_gap / 2

fig3_left <- ggdraw()
fig3_left <- place_panel(
  fig3_left, p4$top,
  x1, heat_bottom, panel_w, heat_h,
  split_x, fig3_height
)
fig3_left <- place_panel(
  fig3_left, ann1$plot,
  legend_x1, legend_bottom, ann1$width, legend_h,
  split_x, fig3_height
)
fig3_left <- place_panel(
  fig3_left, eff1$plot,
  eff1_x, legend_bottom, eff1$width, legend_h,
  split_x, fig3_height
)

fig3_right <- ggdraw()
fig3_right <- place_panel(
  fig3_right, p5$top,
  x2 - split_x, heat_bottom, panel_w, heat_h,
  fig3_width - split_x, fig3_height
)
fig3_right <- place_panel(
  fig3_right, ann2$plot,
  legend_x2 - split_x, legend_bottom, ann2$width, legend_h,
  fig3_width - split_x, fig3_height
)
fig3_right <- place_panel(
  fig3_right, eff2$plot,
  eff2_x - split_x, legend_bottom, eff2$width, legend_h,
  fig3_width - split_x, fig3_height
)

fig3_height_mm <- fig3_height * 25.4

# Column where the two heatmap patches meet (the row spans columns 2-39).
split_col <- 1L + round(38 * split_x / fig3_width)

layout <- c(
  area(t = 1, l = 2, b = 1, r = 39),
  area(t = 2, l = 1, b = 2, r = 40),
  area(t = 3, l = 2, b = 3, r = 39),
  area(t = 4, l = 2, b = 4, r = split_col),
  area(t = 4, l = split_col + 1, b = 4, r = 39)
)

combined <- p +
  plot_spacer() +
  wrap_elements(full = p2) +
  wrap_elements(full = fig3_left) +
  wrap_elements(full = fig3_right) +
  plot_layout(
    design = layout,
    heights = c(48, panel_gap_mm, result$height, fig3_height_mm)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

ggsave(
  filename = file.path(FIGURES, "merged_replogle.png"),
  plot = combined,
  width = result$width,
  height = 48 + panel_gap_mm + result$height +
    panel_gap_mm + fig3_height_mm + 10,
  units = "mm",
  dpi = 450,
  device = ragg::agg_png
)

