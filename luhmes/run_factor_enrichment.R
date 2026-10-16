#!/usr/bin/env Rscript

# GO-BP enrichment by PerturbVI factor.

if (as.character(utils::packageVersion("WebGestaltR")) != "1.0.1") {
  stop("This analysis uses WebGestaltR 1.0.1.")
}

PIP_W <- read.csv("fit/PIP_W.csv", row.names = 1, check.names = FALSE)
background <- colnames(PIP_W)
foregrounds <- lapply(as.data.frame(t(PIP_W)), function(probability) {
  background[probability > 0.95]
})

enrichment_rows <- lapply(names(foregrounds), function(factor) {
  foreground <- foregrounds[[factor]]
  if (!length(foreground)) {
    return(NULL)
  }

  result <- WebGestaltR::WebGestaltR(
    enrichMethod = "ORA",
    organism = "hsapiens",
    enrichDatabase = "geneontology_Biological_Process_noRedundant",
    interestGene = foreground,
    interestGeneType = "ensembl_gene_id",
    referenceGene = background,
    referenceGeneType = "ensembl_gene_id",
    minNum = 10,
    maxNum = 500,
    fdrMethod = "BH",
    sigMethod = "fdr",
    fdrThr = 0.05,
    isOutput = FALSE
  )

  if (is.null(result) || !nrow(result)) {
    return(NULL)
  }

  data.frame(
    factor = factor,
    term_id = result$geneSet,
    description = result$description,
    term_size = result$size,
    overlap_count = result$overlap,
    expected_overlap = result$expect,
    fold_enrichment = result$enrichmentRatio,
    p_value = result$pValue,
    fdr = result$FDR,
    overlap_ids = result$overlapId,
    stringsAsFactors = FALSE
  )
})

enrichment <- do.call(rbind, enrichment_rows)
write.csv(enrichment, "results/suppl_fig_2F_perturbvi.csv", row.names = FALSE)
