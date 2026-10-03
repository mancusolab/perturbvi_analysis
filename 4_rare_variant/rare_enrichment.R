#!/usr/bin/env Rscript
# Rare-variant ASD enrichment: LUHMES PerturbVI fit + external tables -> figure + result tables.
#
# Run: Rscript rare_enrichment.R
# Inputs (input/): PIP_W.csv, W.csv, top6k_genes.csv, asd_rare_variant_genes.csv, hgnc_complete_set.tsv
# Outputs (figures/): asd_rare_variant_enrichment.png, suppl_factor_effect_sizes.png
# Outputs (results/): suppl_factor_enrichment.csv, suppl_cluster_enrichment.csv, suppl_ranked_enrichment.csv

library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
library(ragg)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

PERTURBED <- c("ADNP", "ARID1B", "ASH1L", "CHD2", "CHD8", "CTNND2", "DYRK1A",
               "HDAC5", "MECP2", "MYT1L", "POGZ", "PTEN", "RELN", "SETD5")
FACTORS <- paste0("w", 0:19)
THRESHOLD <- 0.95

# 1. Prep: factor x gene -> gene x wN, ENSG -> symbol -----------------------------

prep_matrix <- function(path) {
  raw <- read.csv(path, row.names = 1, check.names = FALSE)
  m <- t(as.matrix(raw))
  colnames(m) <- paste0("w", sub("^factor_", "", colnames(m)))
  m
}

pip <- prep_matrix("input/PIP_W.csv")
W <- prep_matrix("input/W.csv")

id_map <- read.csv("input/top6k_genes.csv", stringsAsFactors = FALSE)
symbols <- id_map$Name[match(rownames(pip), id_map$ID)]
keep <- !is.na(symbols) & symbols != ""
pip <- pip[keep, , drop = FALSE]
W <- W[keep, , drop = FALSE]
rownames(pip) <- rownames(W) <- symbols[keep]

# 2. HGNC harmonization ----------------------------------------------------------

hgnc <- read.delim("input/hgnc_complete_set.tsv", check.names = FALSE, quote = "",
                   colClasses = "character")
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

# drop ambiguous, collapse aliases to one approved gene by max (matches binary membership)
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
W <- harmonize_table(abs(W))

# 3. ASD genes --------------------------------------------------------------------

asd <- read.csv("input/asd_rare_variant_genes.csv", stringsAsFactors = FALSE)
asd <- asd[asd$ASD_253, ]
asd$approved <- harmonize(asd$Gene)
asd <- asd[!asd$Gene %in% ambiguous, ]
asd$Cluster <- suppressWarnings(as.integer(asd$Cluster))

perturbed <- intersect(harmonize(PERTURBED), rownames(pip))
universe <- setdiff(rownames(pip), perturbed)
asd_bg <- intersect(asd$approved, universe)

fisher_enrich <- function(fset, gset_interest, background) {
  a <- length(intersect(fset, gset_interest))
  b <- length(setdiff(fset, gset_interest))
  cc <- length(setdiff(gset_interest, fset))
  d <- length(setdiff(background, union(fset, gset_interest)))
  tab <- matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE)
  ft_one <- fisher.test(tab, alternative = "greater")
  ft_two <- fisher.test(tab, alternative = "two.sided")
  or <- unname(ft_one$estimate)
  plot_or <- if (is.finite(or) && or > 0) or
             ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (cc + 0.5))
  list(a = a, b = b, c = cc, d = d, odds_ratio = or, plot_odds_ratio = plot_or,
       ci_low = ft_two$conf.int[1], ci_high = ft_two$conf.int[2], p = ft_one$p.value)
}

# 4. Per-factor binary enrichment -------------------------------------------------

factor_enrichment <- map_dfr(FACTORS, function(f) {
  fset <- intersect(rownames(pip)[pip[, f] >= THRESHOLD], universe)
  r <- fisher_enrich(fset, asd_bg, universe)
  tibble(
    factor = f,
    background_size = length(universe),
    n_asd_in_background = length(asd_bg),
    factor_size = length(intersect(fset, universe)),
    n_asd_overlap = r$a,
    overlap_genes = str_c(sort(intersect(fset, asd_bg)), collapse = ","),
    a = r$a, b = r$b, c = r$c, d = r$d,
    odds_ratio = r$odds_ratio, plot_odds_ratio = r$plot_odds_ratio,
    ci_low = r$ci_low, ci_high = r$ci_high, p = r$p
  )
})
factor_enrichment$q_bh_family <- p.adjust(factor_enrichment$p, "BH")

# 5. Per-factor x ASD cluster enrichment ------------------------------------------

asd_clu <- asd[!is.na(asd$Cluster) & asd$Cluster %in% 1:6, ]
asd_clu <- asd_clu[!duplicated(asd_clu$approved), ]
asd_clu <- asd_clu[asd_clu$approved %in% universe, ]

cluster_enrichment <- map_dfr(FACTORS, function(f) {
  fset <- intersect(rownames(pip)[pip[, f] >= THRESHOLD], universe)
  map_dfr(1:6, function(k) {
    cset <- asd_clu$approved[asd_clu$Cluster == k]
    r <- fisher_enrich(fset, cset, universe)
    tibble(
      factor = f, cluster = k,
      n_cluster_genes = length(cset), n_overlap = r$a,
      overlap_genes = str_c(sort(intersect(fset, cset)), collapse = ","),
      a = r$a, b = r$b, c = r$c, d = r$d,
      odds_ratio = r$odds_ratio, plot_odds_ratio = r$plot_odds_ratio,
      ci_low = r$ci_low, ci_high = r$ci_high, p = r$p
    )
  })
})
cluster_enrichment$q_bh <- p.adjust(cluster_enrichment$p, "BH")

# 6. Threshold-free ranked enrichment (Mann-Whitney AUC) --------------------------

ranked_enrichment <- map_dfr(FACTORS, function(f) {
  map_dfr(c("PIP", "abs_loading"), function(measure) {
    values <- if (measure == "PIP") pip[, f] else W[, f]
    values <- values[universe]
    is_asd <- universe %in% asd$approved
    v_asd <- values[is_asd]
    v_rest <- values[!is_asd]
    ranks <- rank(c(v_asd, v_rest))
    u <- sum(ranks[seq_along(v_asd)]) - length(v_asd) * (length(v_asd) + 1) / 2
    auc <- u / (length(v_asd) * length(v_rest))
    p <- wilcox.test(v_asd, v_rest, alternative = "greater", exact = FALSE)$p.value
    tibble(
      factor = f, measure = measure,
      n_asd = length(v_asd), n_rest = length(v_rest),
      median_asd = median(v_asd), median_rest = median(v_rest),
      mean_asd = mean(v_asd), mean_rest = mean(v_rest),
      AUC = auc, rank_biserial = 2 * auc - 1, p = p
    )
  })
})
ranked_enrichment$q_bh_family <- p.adjust(ranked_enrichment$p, "BH")

# 7. Figure -----------------------------------------------------------------------

binary <- factor_enrichment |>
  transmute(factor, binary_FDR = q_bh_family, odds_ratio)
ranked <- ranked_enrichment |> filter(measure == "PIP") |>
  transmute(factor, ranked_FDR = q_bh_family, AUC)

factor_labels <- setNames(paste("Factor", 1:20), FACTORS)
set1 <- brewer.pal(9, "Set1")
threshold <- -log10(0.05)

plot_df <- binary |>
  left_join(ranked, by = "factor") |>
  mutate(
    binary_logFDR = -log10(pmax(binary_FDR, 1e-300)),
    ranked_logFDR = -log10(pmax(ranked_FDR, 1e-300)),
    any_sig = binary_FDR < 0.05 | ranked_FDR < 0.05,
    factor_label = factor_labels[factor]
  )

p_concordance <- ggplot(plot_df, aes(binary_logFDR, ranked_logFDR)) +
  geom_vline(xintercept = threshold, linetype = "dashed", linewidth = 0.55, color = "grey55") +
  geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 0.55, color = "grey55") +
  geom_point(aes(color = any_sig), size = 3.2, alpha = 0.95) +
  geom_text_repel(
    data = filter(plot_df, any_sig), aes(label = factor_label),
    hjust = 0, vjust = 0, size = 3.1, family = "sans", color = "#222222",
    segment.color = "black", segment.size = 0.35, segment.alpha = 0.8,
    min.segment.length = 0, point.padding = 0.5, box.padding = 0.5,
    force = 1, max.overlaps = Inf, seed = 123
  ) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = set1[1])) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0.02, 0.10))) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Concordance of ASD enrichment across PerturbVI factors",
    subtitle = "Dashed lines indicate FDR = 0.05",
    x = expression("PIP > 0.95 enrichment, " * -log[10](FDR)),
    y = expression("PIP-ranked enrichment, " * -log[10](FDR))
  ) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.text = element_text(size = 9, color = "#333333"),
    axis.title = element_text(size = 10, color = "#222222"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, color = "#555555"),
    plot.margin = margin(4, 6, 4, 4),
    legend.position = "none"
  )

cluster_df <- cluster_enrichment |>
  mutate(
    factor = factor(factor, levels = rev(FACTORS)),
    cluster = factor(paste0("C", cluster), levels = paste0("C", 1:6)),
    log2_or = log2(plot_odds_ratio),
    significant = q_bh < 0.05,
    cell_label = paste0(n_overlap, "/", n_cluster_genes, if_else(significant, "*", ""))
  )
or_limit <- max(quantile(abs(cluster_df$log2_or), 0.95, na.rm = TRUE), 2)

p_cluster <- ggplot(cluster_df, aes(cluster, factor, fill = log2_or)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = cell_label, fontface = if_else(significant, "bold", "plain")),
            size = 2.8, family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-or_limit, or_limit), oob = scales::squish,
                       name = expression(log[2]("OR"))) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(labels = factor_labels, expand = c(0, 0)) +
  labs(
    title = "Rare-variant ASD cluster enrichment",
    subtitle = "Cells show overlapping / represented cluster genes; * FDR < 0.05",
    x = NULL, y = NULL,
    caption = paste0("C1  ASD    C2  ASD + modest DD     C3  ASD + SCZ\n",
                     "C4  ASD + DD     C5  ASD + DD + epilepsy     C6  ASD + DD + SCZ")
  ) +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold", color = "#222222"),
    axis.text.y = element_text(size = 8.5, color = "#333333"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, color = "#555555"),
    plot.caption = element_text(size = 7.5, color = "#444444", hjust = 0,
                                lineheight = 1.4, margin = margin(t = 7)),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.45, "in"),
    plot.margin = margin(4, 6, 4, 4)
  )

combined <- p_concordance / p_cluster +
  plot_layout(heights = c(0.85, 1.45)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(family = "sans", size = 14, face = "bold")))

ggsave("figures/asd_rare_variant_enrichment.png", combined, width = 7.2, height = 10,
       units = "in", dpi = 600, bg = "white")

# 7b. Supplementary effect-size figure -------------------------------------------

effect_theme <- theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.text = element_text(size = 9, color = "#333333"),
    axis.title = element_text(size = 10, color = "#222222"),
    plot.title = element_text(size = 10, face = "plain", hjust = 0, margin = margin(t = 8, b = 4)),
    plot.tag = element_text(family = "sans", size = 12, face = "plain"),
    plot.margin = margin(24, 8, 6, 8),
    legend.position = "none"
  )
effect_colors <- c(sig = set1[1], ns = "grey55")
effect_legend <- c(expression(italic(q) < 0.05), expression(italic(q) >= 0.05))

# Shared row order: strongest binary effect first.
ord <- factor_enrichment$factor[order(-factor_enrichment$odds_ratio)]
y_levels <- rev(factor_labels[ord])

forest_df <- factor_enrichment |>
  mutate(sig = if_else(q_bh_family < 0.05, "sig", "ns"),
         ylab = factor(factor_labels[factor], levels = y_levels))

p_odds <- ggplot(forest_df, aes(odds_ratio, ylab)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey60") +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high, color = sig), linewidth = 0.9) +
  geom_point(aes(color = sig), size = 2.8) +
  scale_x_log10(breaks = scales::pretty_breaks(n = 5)) +
  scale_color_manual(values = effect_colors, name = NULL, breaks = c("sig", "ns"),
                     labels = effect_legend) +
  labs(title = "Binary enrichment (loading PIP > 0.95)", tag = "A",
       x = "Odds ratio (95% CI)", y = NULL) +
  effect_theme

auc_df <- ranked_enrichment |>
  filter(measure == "PIP") |>
  mutate(sig = if_else(q_bh_family < 0.05, "sig", "ns"),
         ylab = factor(factor_labels[factor], levels = y_levels))

p_auc <- ggplot(auc_df, aes(AUC, ylab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.5, color = "grey60") +
  geom_point(aes(color = sig), size = 2.8) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(0.4, 0.72)) +
  scale_color_manual(values = effect_colors, name = NULL, breaks = c("sig", "ns"),
                     labels = effect_legend) +
  labs(title = "Threshold-free enrichment (PIP-ranked)", tag = "B",
       x = "AUC (Mann-Whitney U)", y = NULL) +
  effect_theme

shared_legend <- cowplot::get_legend(
  p_odds + theme(legend.position = "bottom", legend.direction = "horizontal",
                 legend.margin = margin(0, 0, 0, 0), legend.key.size = grid::unit(0.3, "cm"))
)

suppl_effect <- (p_odds + p_auc +
                   plot_layout(axes = "collect_y")) /
  wrap_elements(full = shared_legend) +
  plot_layout(heights = c(24, 1))

ggsave("figures/suppl_factor_effect_sizes.png", suppl_effect, width = 9.0, height = 5.6,
       units = "in", dpi = 600, bg = "white", device = ragg::agg_png)

# 8. Final result tables ----------------------------------------------------------

write.csv(factor_enrichment, "results/suppl_factor_enrichment.csv", row.names = FALSE)
write.csv(cluster_enrichment, "results/suppl_cluster_enrichment.csv", row.names = FALSE)
write.csv(ranked_enrichment, "results/suppl_ranked_enrichment.csv", row.names = FALSE)

cat("Done: asd_rare_variant_enrichment.png + 3 result tables\n")
