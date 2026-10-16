#!/usr/bin/env Rscript

# GO-BP enrichment by perturbation for GSFA and PerturbVI.

if (as.character(utils::packageVersion("WebGestaltR")) != "1.0.1") {
  stop("This analysis uses WebGestaltR 1.0.1.")
}

pvi <- as.matrix(read.csv("fit/LFSR_BW.csv", row.names = 1, check.names = FALSE))
gsfa_fit <- readRDS("fit/gsfa_fit.light.rds")
gsfa <- as.matrix(gsfa_fit$lfsr)
background <- colnames(pvi)
gsfa <- gsfa[background, , drop = FALSE]
perturbations <- sort(intersect(rownames(pvi), colnames(gsfa)))

database <- "geneontology_Biological_Process_noRedundant"
all_results <- list()

for (target in perturbations) {
  for (method in c("GSFA", "PerturbVI")) {
    degs <- if (method == "GSFA") {
      background[gsfa[background, target] < 0.05]
    } else {
      background[pvi[target, background] < 0.05]
    }

    if (!length(degs)) next

    message(target, " / ", method, ": ", length(degs), " DEGs")
    out <- WebGestaltR::WebGestaltR(
      enrichMethod = "ORA",
      organism = "hsapiens",
      enrichDatabase = database,
      interestGene = degs,
      interestGeneType = "ensembl_gene_id",
      referenceGene = background,
      referenceGeneType = "ensembl_gene_id",
      isOutput = FALSE
    )

    if (is.null(out) || !nrow(out)) next

    all_results[[length(all_results) + 1L]] <- data.frame(
      method = method,
      perturbation = target,
      term_id = out$geneSet,
      description = out$description,
      term_size = out$size,
      overlap_count = out$overlap,
      expected_overlap = out$expect,
      fold_enrichment = out$enrichmentRatio,
      p_value = out$pValue,
      fdr = out$FDR,
      overlap_ids = out$overlapId,
      significant = TRUE,
      stringsAsFactors = FALSE
    )
  }
}

enrichment <- do.call(rbind, all_results)
write.csv(enrichment, "results/go_bp_enrichment_results.csv", row.names = FALSE)
write.csv(enrichment[enrichment$method == "GSFA", -1],
  "results/suppl_fig_2D_gsfa.csv",
  row.names = FALSE
)
write.csv(enrichment[enrichment$method == "PerturbVI", -1],
  "results/suppl_fig_2D_perturbvi.csv",
  row.names = FALSE
)

counts <- data.frame(
  perturbation = perturbations,
  GSFA = vapply(perturbations, function(target) {
    sum(enrichment$method == "GSFA" & enrichment$perturbation == target)
  }, integer(1)),
  PerturbVI = vapply(perturbations, function(target) {
    sum(enrichment$method == "PerturbVI" & enrichment$perturbation == target)
  }, integer(1))
)
write.csv(counts, "results/go_bp_significant_term_counts.csv", row.names = FALSE)

writeLines(c(
  paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("WebGestaltR: ", utils::packageVersion("WebGestaltR")),
  paste0("Enrichment: ORA; ", database, "."),
  "WebGestaltR default settings; BH-adjusted FDR < 0.05.",
  "Foreground: genes with LFSR < 0.05 for each method and perturbation.",
  paste0("Shared reference background: ", length(background), " modeled genes."),
  paste0("Perturbations analyzed: ", length(perturbations), ".")
), "results/analysis_metadata.txt")
