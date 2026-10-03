#!/usr/bin/env Rscript
# Rare-variant ASD enrichment: gene-constraint (LOEUF) supplementary analyses.
#
# Run: Rscript loeuf.R
# Inputs (input/): PIP_W.csv, top6k_genes.csv, hgnc_complete_set.tsv,
#                  gene_covariates.tsv, asd_rare_variant_genes.csv
# Outputs (results/): suppl_constraint_tests.csv, suppl_permutation_robustness.csv

library(tidyverse)

dir.create("results", showWarnings = FALSE)

PERTURBED <- c("ADNP", "ARID1B", "ASH1L", "CHD2", "CHD8", "CTNND2", "DYRK1A",
               "HDAC5", "MECP2", "MYT1L", "POGZ", "PTEN", "RELN", "SETD5")
FACTORS <- paste0("w", 0:19)
THRESHOLD <- 0.95
N_PERM <- 10000

# 1. Prep: factor x gene -> gene x wN, ENSG -> symbol -----------------------------

prep_matrix <- function(path) {
  raw <- read.csv(path, row.names = 1, check.names = FALSE)
  m <- t(as.matrix(raw))
  colnames(m) <- paste0("w", sub("^factor_", "", colnames(m)))
  m
}

pip <- prep_matrix("input/PIP_W.csv")

id_map <- read.csv("input/top6k_genes.csv", stringsAsFactors = FALSE)
symbols <- id_map$Name[match(rownames(pip), id_map$ID)]
keep <- !is.na(symbols) & symbols != ""
pip <- pip[keep, , drop = FALSE]
rownames(pip) <- symbols[keep]

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

# 3. ASD genes and LOEUF covariates ----------------------------------------------

asd <- read.csv("input/asd_rare_variant_genes.csv", stringsAsFactors = FALSE)
asd <- asd[asd$ASD_253, ]
asd$approved <- harmonize(asd$Gene)
asd <- asd[!asd$Gene %in% ambiguous, ]

cov <- read.delim("input/gene_covariates.tsv", check.names = FALSE, colClasses = "character")
cov <- cov[!cov$Gene %in% ambiguous, ]
cov$approved <- harmonize(cov$Gene)
cov$LOEUF <- suppressWarnings(as.numeric(cov$LOEUF))
cov <- aggregate(LOEUF ~ approved, data = cov, FUN = median)

perturbed <- intersect(harmonize(PERTURBED), rownames(pip))
universe <- setdiff(rownames(pip), perturbed)
uni <- data.frame(approved = universe) |> left_join(cov, by = "approved")
has_loeuf <- !is.na(uni$LOEUF)
asd_bg <- universe %in% intersect(asd$approved, universe)
cat(sprintf("background %d genes; LOEUF for %d (%.1f%%); median LOEUF %.3f\n",
            length(universe), sum(has_loeuf), 100 * mean(has_loeuf),
            median(uni$LOEUF[has_loeuf])))

# 4. Constraint checks (PIP vs LOEUF) --------------------------------------------

check1 <- map_dfr(FACTORS, function(f) {
  pipv <- pip[universe, f]
  mask <- has_loeuf & is.finite(pipv)
  ct <- suppressWarnings(cor.test(pipv[mask], uni$LOEUF[mask], method = "spearman"))
  tibble(check = "spearman_pip_loeuf", factor = f, n = sum(mask),
         rho = unname(ct$estimate), p = ct$p.value)
})
check1$q_bh <- p.adjust(check1$p, "BH")

check2 <- map_dfr(FACTORS, function(f) {
  fac <- universe[pip[universe, f] >= THRESHOLD]
  fac_loeuf <- uni$LOEUF[uni$approved %in% fac & has_loeuf]
  rest_loeuf <- uni$LOEUF[!uni$approved %in% fac & has_loeuf]
  pv <- wilcox.test(fac_loeuf, rest_loeuf, alternative = "two.sided", exact = FALSE)$p.value
  tibble(check = "loeuf_factor_vs_background", factor = f, n = length(fac_loeuf),
         rho = NA_real_, p = pv,
         median_factor = median(fac_loeuf), median_background = median(rest_loeuf))
})
check2$q_bh <- p.adjust(check2$p, "BH")

checks <- bind_rows(check1, check2) |>
  select(check, factor, n, rho, median_factor, median_background, p, q_bh)
write.csv(checks, "results/suppl_constraint_tests.csv", row.names = FALSE)

# 5. LOEUF-stratified permutation ------------------------------------------------

set.seed(0)
membership <- matrix(FALSE, nrow = length(universe), ncol = length(FACTORS))
for (j in seq_along(FACTORS)) {
  membership[, j] <- universe %in% rownames(pip)[pip[, FACTORS[j]] >= THRESHOLD]
}
storage.mode(membership) <- "integer"
asd_int <- as.integer(asd_bg)
observed <- as.vector(asd_int %*% membership)

strata <- rep(0L, length(universe))
breaks <- unique(quantile(uni$LOEUF[has_loeuf], probs = seq(0, 1, 0.1)))
strata[has_loeuf] <- as.integer(cut(uni$LOEUF[has_loeuf], breaks = breaks,
                                    include.lowest = TRUE, labels = FALSE))
strata_idx <- split(seq_along(universe), strata)

n_at_least <- sum_perm <- sumsq_perm <- numeric(length(FACTORS))
for (i in seq_len(N_PERM)) {
  permuted <- asd_int
  for (idx in strata_idx) permuted[idx] <- sample(permuted[idx])
  overlap <- as.vector(permuted %*% membership)
  n_at_least <- n_at_least + as.integer(overlap >= observed)
  sum_perm <- sum_perm + overlap
  sumsq_perm <- sumsq_perm + overlap^2
}
null_mean <- sum_perm / N_PERM
null_sd <- sqrt(pmax(sumsq_perm / N_PERM - null_mean^2, 0))

permutation <- tibble(
  factor = FACTORS,
  background_size = length(universe),
  loeuf_coverage = mean(has_loeuf),
  n_factor = colSums(membership),
  n_asd_in_background = sum(asd_bg),
  n_asd_overlap = observed,
  null_mean = null_mean,
  null_sd = null_sd,
  n_perm = N_PERM,
  perm_p = (n_at_least + 1) / (N_PERM + 1)
)
permutation$perm_q_bh <- p.adjust(permutation$perm_p, "BH")
write.csv(permutation, "results/suppl_permutation_robustness.csv", row.names = FALSE)

# 6. Supplementary figure: constraint effect -------------------------------------
# Everything below is commented out; the constraint results are reported as the
# two tables only. Uncomment to also write figures/suppl_constraint.png.
#
# library(RColorBrewer)
# library(patchwork)
# library(ragg)
# dir.create("figures", showWarnings = FALSE)
#
# factor_labels <- setNames(paste("Factor", 1:20), FACTORS)
# set1 <- brewer.pal(9, "Set1")
# effect_colors <- c(sig = set1[1], ns = "grey55")
# effect_legend <- c(expression(italic(q) < 0.05), expression(italic(q) >= 0.05))
#
# plot_df <- check1 |>
#   select(factor, rho, rho_sig = q_bh) |>
#   left_join(check2 |> select(factor, median_factor, constraint_sig = q_bh), by = "factor") |>
#   mutate(
#     rho_sig = if_else(rho_sig < 0.05, "sig", "ns"),
#     constraint_sig = if_else(constraint_sig < 0.05, "sig", "ns"),
#     ylab = factor(factor_labels[factor], levels = rev(factor_labels[order(median_factor)]))
#   )
#
# effect_theme <- theme_classic(base_size = 10, base_family = "sans") +
#   theme(
#     axis.text = element_text(size = 9, color = "#333333"),
#     axis.title = element_text(size = 10, color = "#222222"),
#     plot.title = element_text(size = 10, face = "plain", hjust = 0,
#                               margin = margin(t = 8, b = 4)),
#     plot.tag = element_text(family = "sans", size = 12, face = "plain"),
#     plot.margin = margin(24, 8, 6, 8),
#     legend.position = "none"
#   )
#
# p_rho <- ggplot(plot_df, aes(rho, ylab)) +
#   geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey60") +
#   geom_point(aes(color = rho_sig), size = 2.8) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
#   scale_color_manual(values = effect_colors, name = NULL, breaks = c("sig", "ns"),
#                      labels = effect_legend) +
#   labs(title = "PIP vs constraint (Spearman)", tag = "A",
#        x = expression("Spearman " * rho * " (PIP vs LOEUF)"), y = NULL) +
#   effect_theme
#
# p_constraint <- ggplot(plot_df, aes(median_factor, ylab)) +
#   geom_vline(xintercept = median(uni$LOEUF[has_loeuf]), linetype = "dashed",
#              linewidth = 0.5, color = "grey60") +
#   geom_point(aes(color = constraint_sig), size = 2.8) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
#   scale_color_manual(values = effect_colors, name = NULL, breaks = c("sig", "ns"),
#                      labels = effect_legend) +
#   labs(title = "Constraint of factor genes", tag = "B",
#        x = "Median LOEUF (PIP > 0.95)", y = NULL) +
#   effect_theme
#
# shared_legend <- cowplot::get_legend(
#   p_rho + theme(legend.position = "bottom", legend.direction = "horizontal",
#                 legend.margin = margin(0, 0, 0, 0), legend.key.size = grid::unit(0.3, "cm"))
# )
#
# suppl_constraint <- (p_rho + p_constraint + plot_layout(axes = "collect_y")) /
#   wrap_elements(full = shared_legend) +
#   plot_layout(heights = c(24, 1))
#
# ggsave("figures/suppl_constraint.png", suppl_constraint, width = 9.0, height = 5.6,
#        units = "in", dpi = 600, bg = "white", device = ragg::agg_png)

cat("Done: 2 result tables\n")
