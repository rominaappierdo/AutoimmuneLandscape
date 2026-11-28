## =====================================================================
## Grafici dei top pathway nella meta-analisi ssGSEA maschi vs femmine (solo pazienti)
## - NON rilancia nessuna meta-analisi: legge i *_complete_results.tsv e gli oggetti *_metaintegrator_object.qs
##   gia' prodotti dagli script perform_meta_analysis_gsva_sex_*.R
## - filtro: effectSizeFDR < 0.05 e |effectSize| > 0.4 (nessun filtro su numStudies)
## - scrive SOLO in: /mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_top_plots/output_pathways/
##   (non tocca le cartelle dei risultati; non crea PDF vuoti)
## =====================================================================
library(MetaIntegrator)   # forestPlot

fdr_cutoff <- 0.05
es_cutoff  <- 0.4

root     <- "/mnt/data3/appierdo/AutoimmuneLandscape"
out_root <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_top_plots/output_pathways"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_top_plots")
source("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_top_plots/functions_plot_top_features.R")

contexts  <- c("BLOOD", "TISSUE")
gene_sets <- c("hallmarks", "kegg", "wikipathways")

summary_list <- list()
for (ctx in contexts) {
  for (gs in gene_sets) {
    base <- file.path(root, "Analysis", paste0(ctx, "_sex_differences"), "meta_analysis_pathways_activity", gs)
    if (!dir.exists(base)) next
    for (disease in list.dirs(base, full.names = FALSE, recursive = FALSE)) {
      in_dir <- file.path(base, disease)
      for (rf in list.files(in_dir, pattern = "_complete_results\\.tsv$", full.names = TRUE)) {
        prefix  <- sub("_complete_results\\.tsv$", "", basename(rf))
        qs_file <- file.path(in_dir, paste0(prefix, "_metaintegrator_object.qs"))
        out_dir <- file.path(out_root, ctx, gs, disease)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        message("### ", ctx, " | ", gs, " | ", disease, " | ", prefix)
        info <- plot_top_features(rf, qs_file, file.path(out_dir, prefix), fdr_cutoff, es_cutoff)

        summary_list[[length(summary_list) + 1]] <- data.frame(
          context = ctx, gene_set = gs, disease = disease, comparison = prefix,
          n_pathways = info$n_total, n_pass_filter = info$n_selected,
          n_higher_in_male = info$n_male, n_higher_in_female = info$n_female,
          n_plotted = info$n_plotted, pdf = info$pdf, stringsAsFactors = FALSE)
      }
    }
  }
}

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  write.table(summary_df, file.path(out_root, sprintf("summary_top_pathways_FDR%s_ES%s.tsv", fdr_cutoff, es_cutoff)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(summary_df[, c("context", "gene_set", "disease", "n_pathways", "n_pass_filter", "n_higher_in_male", "n_higher_in_female", "n_plotted")])
} else {
  message("Nessuna cartella di risultati trovata sotto ", root)
}

print("###############DONE! (top pathways M vs F)")
