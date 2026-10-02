#!/usr/bin/env Rscript
# LUHMES factor GO enrichment.
#
# Computes GO Biological Process ORA (WebGestaltR) for each factor's
# PIP_W > 0.95 gene set, with BH correction, and saves the combined table.
#
# Run:
#   Rscript run_factor_enrichment.R
#
# Outputs: input/factor_enrichment.csv

library(WebGestaltR)

# 1. select factor genes with loading PIP > 0.95.
PIP_W <- read.csv("results/PIP_W.csv", row.names = 1, check.names = FALSE)

background <- colnames(PIP_W)
foregrounds <- lapply(as.data.frame(t(PIP_W)), function(probability) {
  background[probability > 0.95]
})

# 2. run GO Biological Process ORA with BH correction per factor.
results <- lapply(foregrounds, function(foreground) {
  if (!length(foreground)) return(NULL)
  WebGestaltR::WebGestaltR(
    enrichMethod = "ORA",
    organism = "hsapiens",
    interestGene = foreground,
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
    fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05,
    isOutput = FALSE
  )
})

# 3. save the combined enrichment table.
enrichment <- do.call(rbind, lapply(names(results), function(id) {
  ans <- results[[id]]
  if (is.null(ans) || !nrow(ans)) return(NULL)
  data.frame(group = id, term_id = ans$geneSet, description = ans$description,
             fold_enrichment = ans$enrichmentRatio, fdr = ans$FDR,
             overlap = ans$overlap, p_value = ans$pValue)
}))
write.csv(enrichment, "input/factor_enrichment.csv", row.names = FALSE)
