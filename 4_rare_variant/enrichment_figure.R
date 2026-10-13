#!/usr/bin/env Rscript
# Figure 4: forest plot of the rare-variant ASD enrichment
#
# Input:  results/glm_asd_enrichment.csv
# Output: figures/rare_enrichment.png
#
# Run: Rscript enrichment_figure.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

dat <- read_csv("results/glm_asd_enrichment.csv", show_col_types = FALSE) |>
  mutate(
    sig_cat = if_else(q < 0.05, "q < 0.05", "q >= 0.05"),
    q_label = case_when(
      q < 0.001 ~ formatC(q, format = "e", digits = 1),
      q > 0.99 ~ ">0.99",
      q < 0.1 ~ formatC(q, format = "f", digits = 4),
      TRUE ~ formatC(q, format = "f", digits = 3)
    ),
    q_label = if_else(q < 0.05, paste0(q_label, "*"), q_label)
  )

# sort by OR
ord <- dat |> arrange(or_sd) |> pull(factor_label)
dat <- dat |> mutate(factor_label = factor(factor_label, levels = ord))

cols <- c("q >= 0.05" = "#5B6B73",
          "q < 0.05"  = "#D7263D")

# --- x-axis geometry --------------------------------------------------------
# Same construction as final_figure.R: the q-value column is drawn inside the
# panel, to the right of the longest confidence interval, so the panel has to be
# wide enough to hold that text. The right limit is derived from the geometry:
# the label text has a fixed physical width and the panel has a fixed physical
# width, so the fraction of the panel the labels occupy fixes how much
# log10-axis span is needed beyond the labels.
PANEL_IN_MARGINS <- 0.93  # y tick labels + left and right plot margins, inches
LABEL_IN         <- 0.47  # width of the longest q label ("7.4e-13*", size 3.0)
LABEL_GAP        <- 0.02  # blank strip after the labels, as a fraction of panel
PAD_CI           <- 1.06  # whisker headroom before the label column
PAD_LABEL        <- 1.03  # labels start just right of that headroom

# ggplot expands a continuous axis by 5% of its span on each side, so the
# visible span is 1.10 x the span of the limits we set. With lo, hi the extreme
# confidence limits and t the log10 position where the labels start, the
# requirement  visible_right - t >= 1.10 * (phi + gap) * span  solves to
#   span = (t - log10(lo)) / (1.05 - 1.10 * (phi + gap))
# because visible_right = right_limit + 0.05 * span = log10(lo) + 1.05 * span.
x_scale_limits <- function(width_in, min_lo, max_hi) {
  lo  <- min_lo * 0.95
  a   <- log10(lo)
  t   <- log10(max_hi * PAD_CI * PAD_LABEL)
  phi <- LABEL_IN / (width_in - PANEL_IN_MARGINS)
  d   <- (t - a) / (1.05 - 1.10 * (phi + LABEL_GAP))
  c(lo, 10^(a + d))
}

width  <- 6.00
height <- 5.80

lo     <- min(dat$or_lo_sd)
hi     <- max(dat$or_hi_sd)
x_text <- hi * PAD_CI
x_lim  <- x_scale_limits(width, lo, hi)

p <- ggplot(dat, aes(y = factor_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "grey55") +
  geom_linerange(aes(xmin = or_lo_sd, xmax = or_hi_sd, colour = sig_cat),
                 linewidth = 1.0, lineend = "round", show.legend = FALSE) +
  geom_point(aes(x = or_sd, colour = sig_cat), shape = 18, size = 3.2) +
  geom_text(aes(x = x_text * PAD_LABEL, label = q_label), hjust = 0, size = 3.0,
            colour = "#4A444E") +
  scale_colour_manual(
    values = cols,
    breaks = c("q >= 0.05", "q < 0.05"),
    labels = c("not significant", "significant"),
    name = NULL) +
  scale_x_log10(limits = x_lim, breaks = c(0.75, 1, 1.25, 1.5),
                labels = c("0.75", "1", "1.25", "1.5")) +
  labs(x = "OR", y = NULL) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    text = element_text(colour = "#1F1B23"),
    axis.text.y = element_text(size = 10, colour = "#2b2b2b"),
    axis.text.x = element_text(size = 9.5, colour = "#2b2b2b"),
    axis.title.x = element_text(size = 10.5, colour = "#1f1f1f",
                                margin = margin(t = 9)),
    axis.line = element_line(linewidth = 0.4, colour = "grey30"),
    axis.ticks = element_line(linewidth = 0.4, colour = "grey30"),
    legend.position = "bottom",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 9.5),
    legend.key.width = grid::unit(0.18, "cm"),
    legend.key.spacing = grid::unit(1.2, "mm"),
    legend.key.spacing.x = grid::unit(3, "mm"),
    legend.margin = margin(t = 2, b = 0),
    panel.grid = element_blank(),
    plot.margin = margin(4, 14, 4, 6)
  )

ggsave("figures/rare_enrichment.png", p, width = width, height = height, dpi = 600,
       bg = "white", device = ragg::agg_png)
