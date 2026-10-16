#!/usr/bin/env Rscript
# Figure 5: rare-variant ASD enrichment in PerturbVI factors.
# Left:  forest plot of odds ratios per 1-SD increase in PIP, all 20 factors.
# Right: PIP heatmap of the ASD genes that reach PIP > 0.95 in at least one of
#        the four enriched factors (14, 3, 11, 10); rows grouped by PIP pattern.
#
# The 0.95 threshold only picks the genes shown in the heatmap; the enrichment
# test on the left uses the continuous PIPs.
#
# Input:  data/PIP_W.csv,
#         data/top6k_genes.csv, data/hgnc_complete_set.tsv,
#         data/asd_rare_variant_genes.csv, results/rare_asd_enrichment.csv
# Output: figures/figure_5.png
#
# Run: Rscript plot_rare_asd_enrichment.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

# ---- forest plot -----------------------------------------------------------

dat <- read_csv("results/rare_asd_enrichment.csv",
                show_col_types = FALSE) |>
  mutate(
    sig_cat = if_else(q < 0.05, "q < 0.05", "q >= 0.05")
  )

# sort by OR
ord <- dat |> arrange(or_sd) |> pull(factor_label)
dat <- dat |> mutate(factor_label = factor(factor_label, levels = ord))

cols <- c("q >= 0.05" = "#5B6B73",
          "q < 0.05"  = "#D7263D")

S <- 0.76

p_forest <- ggplot(dat, aes(y = factor_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45 * S, colour = "grey55") +
  geom_linerange(aes(xmin = or_lo_sd, xmax = or_hi_sd, colour = sig_cat),
                 linewidth = 1.0 * S, lineend = "round", show.legend = FALSE) +
  geom_point(aes(x = or_sd, colour = sig_cat), shape = 18, size = 3.2 * S) +
  scale_colour_manual(
    values = cols,
    breaks = c("q >= 0.05", "q < 0.05"),
    labels = c("not significant", "BH-adjusted P < 0.05"),
    name = NULL) +
  scale_x_log10(breaks = c(0.75, 1, 1.25, 1.5),
                labels = c("0.75", "1", "1.25", "1.5")) +
  labs(x = "OR per 1-SD increase in PIP", y = NULL) +
  theme_classic(base_size = 11 * S, base_family = "sans") +
  theme(
    text = element_text(colour = "#1F1B23"),
    axis.text.y = element_text(size = 10 * S, colour = "#2b2b2b"),
    axis.text.x = element_text(size = 9.5 * S, colour = "#2b2b2b"),
    axis.title.x = element_text(size = 10.5 * S, colour = "#1f1f1f",
                                margin = margin(t = 7 * S)),
    axis.line = element_line(linewidth = 0.4 * S, colour = "grey30"),
    axis.ticks = element_line(linewidth = 0.4 * S, colour = "grey30"),
    legend.position = "bottom",
    legend.title = element_text(size = 9.5 * S),
    legend.text = element_text(size = 9.5 * S),
    legend.key.width = grid::unit(0.18 * S, "cm"),
    legend.key.spacing = grid::unit(1.2 * S, "mm"),
    legend.key.spacing.x = grid::unit(3 * S, "mm"),
    legend.margin = margin(t = 2 * S, b = 0),
    panel.grid = element_blank(),
    plot.margin = margin(3 * S, 7 * S, 3 * S, 5 * S)
  )

# ---- gene-level PIPs for the heatmap ---------------------------------------

prep_matrix <- function(path) {
  raw <- read.csv(path, row.names = 1, check.names = FALSE)
  m <- t(as.matrix(raw))
  colnames(m) <- paste0("w", sub("^factor_", "", colnames(m)))
  m
}

pip <- prep_matrix("data/PIP_W.csv")

id_map <- read.csv("data/top6k_genes.csv",
                   stringsAsFactors = FALSE)
symbols <- id_map$Name[match(rownames(pip), id_map$ID)]
keep <- !is.na(symbols) & symbols != ""
pip <- pip[keep, , drop = FALSE]
rownames(pip) <- symbols[keep]

hgnc <- read.delim("data/hgnc_complete_set.tsv",
                   check.names = FALSE,
                   quote = "", colClasses = "character")
apps <- trimws(hgnc$`Approved symbol`)
approved <- unique(apps[hgnc$Status == "Approved" & !is.na(apps) & apps != ""])
prev_sym <- strsplit(hgnc$`Previous symbols`, ",", fixed = TRUE)
alias_sym <- strsplit(hgnc$`Alias symbols`, ",", fixed = TRUE)
alias_long <- data.frame(
  alias = trimws(c(unlist(prev_sym), unlist(alias_sym))),
  app = c(rep(apps, lengths(prev_sym)), rep(apps, lengths(alias_sym))),
  stringsAsFactors = FALSE
)
alias_long <- unique(alias_long[!is.na(alias_long$alias) & alias_long$alias != "" &
                                  !alias_long$alias %in% approved & !is.na(alias_long$app), ])
cand_n <- table(alias_long$alias)
ambiguous <- names(cand_n)[cand_n > 1]
sym_map <- setNames(approved, approved)
unambig <- alias_long[alias_long$alias %in% names(cand_n)[cand_n == 1], ]
sym_map[unambig$alias] <- unambig$app

harmonize <- function(x) {
  x <- trimws(as.character(x))
  mapped <- sym_map[x]
  ifelse(is.na(mapped), x, mapped)
}

harmonize_table <- function(scores) {
  scores <- scores[!rownames(scores) %in% ambiguous, , drop = FALSE]
  df <- as.data.frame(scores, check.names = FALSE)
  df$approved <- harmonize(rownames(scores))
  agg <- aggregate(. ~ approved, data = df, FUN = max)
  rownames(agg) <- agg$approved
  agg$approved <- NULL
  as.matrix(agg)
}

pip <- harmonize_table(pip)

PERTURBED <- c("ADNP", "ARID1B", "ASH1L", "CHD2", "CHD8", "CTNND2", "DYRK1A",
               "HDAC5", "MECP2", "MYT1L", "POGZ", "PTEN", "RELN", "SETD5")
asd <- read.csv("data/asd_rare_variant_genes.csv",
                stringsAsFactors = FALSE)
asd <- asd[asd$ASD_253, ]
asd$approved <- harmonize(asd$Gene)
asd <- asd[!asd$Gene %in% ambiguous, ]

universe <- setdiff(rownames(pip), intersect(harmonize(PERTURBED), rownames(pip)))
asd_bg <- intersect(asd$approved, universe)

# ---- heatmap ---------------------------------------------------------------

FACTORS      <- c("w13", "w2", "w10", "w9")   # Factors 14, 3, 11, 10
FACTOR_NAMES <- c("Factor 14", "Factor 3", "Factor 11", "Factor 10")

mat <- pip[asd_bg, FACTORS, drop = FALSE]
sel <- rownames(mat)[rowSums(mat > 0.95) >= 1]  # display selection only
mat <- mat[sel, , drop = FALSE]
mat <- mat[hclust(dist(mat), method = "ward.D2")$order, , drop = FALSE]
cat(sprintf("heatmap: %d of %d ASD genes reach PIP > 0.95 in >= 1 enriched factor\n",
            nrow(mat), length(asd_bg)))

heat <- data.frame(
  gene   = factor(rep(rownames(mat), times = ncol(mat)), levels = rev(rownames(mat))),
  factor = factor(rep(FACTOR_NAMES, each = nrow(mat)), levels = FACTOR_NAMES),
  pip    = as.vector(mat)
)

p_heat <- ggplot(heat, aes(factor, gene, fill = pip)) +
  geom_tile(colour = "white", linewidth = 0.3 * S) +
  scale_fill_gradient(low = "#F0F0F0", high = "#D7263D", limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = "PIP") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right", expand = expansion(add = 0)) +
  annotate("text", x = 0.2, y = (nrow(mat) + 1) / 2, label = "Rare-variant ASD genes",
           angle = 90, size = 3.2 * S, colour = "#2b2b2b") +
  annotate("text", x = (ncol(mat) + 1) / 2, y = 0.1, label = "Factors",
           vjust = 1, size = 3.4 * S, colour = "#1f1f1f") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_colourbar(barwidth = grid::unit(2.4 * S, "cm"),
                                barheight = grid::unit(0.35 * S, "cm"),
                                title.position = "top", title.hjust = 0.5,
                                ticks.colour = "grey30", frame.colour = "grey30")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11 * S, base_family = "sans") +
  theme(
    text = element_text(colour = "#1F1B23"),
    axis.text.y.right = element_text(size = 8 * S, colour = "#2b2b2b", face = "italic"),
    axis.text.x.top = element_text(size = 8.6 * S, colour = "#2b2b2b"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9.5 * S),
    legend.text = element_text(size = 9 * S),
    legend.margin = margin(t = 2 * S, b = 0),
    panel.grid = element_blank(),
    plot.margin = margin(10 * S, 4 * S, 16 * S, 26 * S)
  )

# ---- combine ---------------------------------------------------------------

combined <- p_forest + p_heat +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13 * S, face = "bold", family = "sans"))

ggsave("figures/figure_5.png", combined, width = 17.4 / 2.54, height = 5.75,
       dpi = 600, bg = "white", device = ragg::agg_png)
