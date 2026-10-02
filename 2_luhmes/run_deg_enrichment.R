#!/usr/bin/env Rscript
# DEG GO enrichment: GO BP ORA per perturbation for GSFA and PerturbVI
# (genes with LFSR < 0.05), using the frozen WebGestalt annotation.
#
# Run:
#   Rscript run_deg_enrichment.R
#
# Outputs:
#   input/deg_go_enrichment.csv - per-term associations (all tested terms)
#   input/deg_go_counts.csv     - per-perturbation significant-term counts

library(WebGestaltR)

pvi <- read.csv("results/LFSR_BW.csv", row.names = 1, check.names = FALSE)
background <- colnames(pvi)

gsfa <- readRDS("input/gsfa_fit.light.rds")$lfsr
stopifnot(setequal(background, rownames(gsfa)))

# drop Nontargeting / offset columns in the GSFA fit
perturbations <- sort(intersect(rownames(pvi), colnames(gsfa)))

run_ora <- function(degs, method, perturbation) {
  if (!length(degs)) return(NULL)
  ans <- WebGestaltR::WebGestaltR(
    enrichMethod = "ORA",
    organism = "hsapiens",
    interestGene = degs,
    interestGeneType = "ensembl_gene_id",
    referenceGene = background,
    referenceGeneType = "ensembl_gene_id",
    # for reproducibility reasons, we use the frozen annotation
    enrichDatabaseFile = "input/frozen_annotation/go_bp_noRedundant.gmt",
    enrichDatabaseDescriptionFile = "input/frozen_annotation/go_bp_noRedundant.des",
    enrichDatabaseType = "entrezgene",
    # instead, you can use the live annotation, 
    # enrichDatabase = "geneontology_Biological_Process_noRedundant",
    minNum = 10, maxNum = 500,
    fdrMethod = "BH", sigMethod = "fdr", fdrThr = 1,
    isOutput = FALSE
  )
  if (is.null(ans) || !nrow(ans)) return(NULL)
  data.frame(
    method = method,
    perturbation = perturbation,
    term_id = ans$geneSet,
    description = ans$description,
    fold_enrichment = ans$enrichmentRatio,
    fdr = ans$FDR,
    overlap = ans$overlap,
    p_value = ans$pValue,
    significant = ans$FDR < 0.05,
    stringsAsFactors = FALSE
  )
}

rows <- lapply(perturbations, function(target) {
  lapply(c("GSFA", "PerturbVI"), function(method) {
    if (method == "PerturbVI") {
      degs <- background[pvi[target, background] < 0.05]
    } else {
      degs <- background[gsfa[background, target] < 0.05]
    }
    run_ora(degs, method, target)
  })
})

enrichment <- do.call(rbind, unlist(rows, recursive = FALSE))
write.csv(enrichment, "input/deg_go_enrichment.csv", row.names = FALSE)

counts <- data.frame(
  perturbation = perturbations,
  GSFA = vapply(perturbations, function(t)
    sum(enrichment$method == "GSFA" & enrichment$perturbation == t &
          enrichment$significant), integer(1)),
  PerturbVI = vapply(perturbations, function(t)
    sum(enrichment$method == "PerturbVI" & enrichment$perturbation == t &
          enrichment$significant), integer(1)),
  check.names = FALSE
)
write.csv(counts, "input/deg_go_counts.csv", row.names = FALSE)
