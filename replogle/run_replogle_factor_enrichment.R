#!/usr/bin/env Rscript
# Recompute factor-level GO-BP enrichment for Replogle/K562.

suppressPackageStartupMessages(library(WebGestaltR))

pip_data <- read.csv("fit/PIP_W.csv", check.names = FALSE, stringsAsFactors = FALSE)
factor_ids <- as.character(pip_data[[1]])
pip <- as.matrix(pip_data[, -1, drop = FALSE])
storage.mode(pip) <- "numeric"
rownames(pip) <- factor_ids
background <- colnames(pip)

results <- lapply(seq_len(nrow(pip)), function(i) {
  factor_id <- factor_ids[[i]]
  foreground <- background[pip[i, ] > 0.95]
  if (!length(foreground)) return(NULL)

  ans <- WebGestaltR::WebGestaltR(
    enrichMethod = "ORA",
    organism = "hsapiens",
    enrichDatabase = "geneontology_Biological_Process_noRedundant",
    interestGene = foreground,
    interestGeneType = "ensembl_gene_id",
    referenceGene = background,
    referenceGeneType = "ensembl_gene_id",
    isOutput = FALSE
  )

  if (is.null(ans) || !nrow(ans)) return(NULL)

  data.frame(
    factor = factor_id,
    term_id = ans$geneSet,
    term = ans$description,
    term_genes = if ("size" %in% names(ans)) ans$size else NA_integer_,
    overlap = ans$overlap,
    fold_enrichment = ans$enrichmentRatio,
    p_value = ans$pValue,
    fdr = ans$FDR,
    overlap_genes = if ("overlapId" %in% names(ans)) ans$overlapId else NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

enrichment <- do.call(rbind, results)
write.csv(enrichment, "results/suppl_fig_3B_factors_enrichment.csv", row.names = FALSE, na = "")