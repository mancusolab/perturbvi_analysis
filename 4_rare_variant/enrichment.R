#!/usr/bin/env Rscript
# Rare-variant ASD enrichment in LUHMES PerturbVI factors.
# Model: logit(P(ASD)) = b0 + b1 * PIP
# Report PIP odds ratios per 1 SD, with two-sided 95% Wald confidence intervals.
# Test enrichment using one-sided (upper-tail) Wald z-tests.
# BH correction is applied across the 20 factors.
#
# Run: Rscript enrichment.R
# Inputs (input/): PIP_W.csv, top6k_genes.csv,
#                 hgnc_complete_set.tsv, asd_rare_variant_genes.csv
# Output (results/): glm_asd_enrichment.csv

library(tidyverse)

dir.create("results", showWarnings = FALSE)

PERTURBED <- c("ADNP", "ARID1B", "ASH1L", "CHD2", "CHD8", "CTNND2", "DYRK1A",
               "HDAC5", "MECP2", "MYT1L", "POGZ", "PTEN", "RELN", "SETD5")
FACTORS <- paste0("w", 0:19)
factor_labels <- setNames(paste("Factor", 1:20), FACTORS)

# ---- 1. Prep: factor x gene -> gene x wN, ENSG -> symbol -----------------------

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

# ---- 2. HGNC harmonization --------------------------------------------------

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

# ---- 3. ASD genes and background universe -----------------------------------

asd <- read.csv("input/asd_rare_variant_genes.csv", stringsAsFactors = FALSE)
asd <- asd[asd$ASD_253, ]
asd$approved <- harmonize(asd$Gene)
asd <- asd[!asd$Gene %in% ambiguous, ]

perturbed <- intersect(harmonize(PERTURBED), rownames(pip))
universe <- setdiff(rownames(pip), perturbed)
asd_bg <- universe %in% intersect(asd$approved, universe)

cat(sprintf("Background: %d genes; ASD genes: %d\n",
            length(universe), sum(asd_bg)))

# ---- 4. Per-factor logistic regressions ------------------------------------

extract_pip <- function(fit, sd_pip) {
  cf <- summary(fit)$coefficients
  beta <- cf["pipv", "Estimate"]
  se <- cf["pipv", "Std. Error"]
  z <- beta / se
  tibble(
    beta_pip = beta,
    se_pip = se,
    or_sd = exp(beta * sd_pip),
    or_lo_sd = exp((beta - 1.96 * se) * sd_pip),
    or_hi_sd = exp((beta + 1.96 * se) * sd_pip),
    p = pnorm(z, lower.tail = FALSE)
  )
}

res <- map_dfr(FACTORS, function(f) {
  pipv_all <- pip[universe, f]
  ok <- is.finite(pipv_all)
  pipv <- pipv_all[ok]
  y <- asd_bg[ok]
  sd_pip <- sd(pipv)

  fit_u <- glm(y ~ pipv, family = binomial)
  u <- extract_pip(fit_u, sd_pip)

  bind_cols(
    tibble(
      factor = f,
      factor_label = unname(factor_labels[f]),
      n_genes = sum(ok),
      n_asd = sum(y),
      sd_pip = sd_pip
    ),
    u
  )
})

res$q <- p.adjust(res$p, method = "BH")

write.csv(res, "results/glm_asd_enrichment.csv", row.names = FALSE)

cat("\n--- ASD enrichment: PIP odds ratio per 1 SD ---\n")
print(as.data.frame(res |>
  transmute(
    factor_label,
    n_genes,
    n_asd,
    or = round(or_sd, 3),
    ci_lower = round(or_lo_sd, 3),
    ci_upper = round(or_hi_sd, 3),
    p = signif(p, 3),
    q = signif(q, 3)
  )), row.names = FALSE)

cat("\nDone: results/glm_asd_enrichment.csv\n")
