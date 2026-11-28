## =====================================================================
## Funzioni per i grafici "top features" (forest plot) delle meta-analisi M vs F
## Filtro: SOLO effectSizeFDR < fdr_cutoff  e  |effectSize| > es_cutoff
##         (nessun filtro su numStudies, a differenza di forestplot_top_genes della pipeline)
## effectSize > 0 = piu' alto nei MASCHI;  effectSize < 0 = piu' alto nelle FEMMINE
## =====================================================================

select_top_features <- function(results_tsv, fdr_cutoff = 0.05, es_cutoff = 0.4) {
  # read.delim (e non read.table): i nomi dei pathway possono contenere apostrofi
  res <- read.delim(results_tsv, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- !is.na(res$effectSize) & !is.na(res$effectSizeFDR) &
          res$effectSizeFDR < fdr_cutoff & abs(res$effectSize) > es_cutoff
  sel <- res[keep, , drop = FALSE]
  sel <- sel[order(sel$effectSize, decreasing = TRUE), , drop = FALSE]
  sel$higher_in <- ifelse(sel$effectSize > 0, "male", "female")
  list(n_total = nrow(res), selected = sel)
}

# Scrive: <out_prefix>_top_features_FDR<..>_ES<..>.tsv  (sempre, se c'e' almeno una feature)
#         <out_prefix>_top_features_FDR<..>_ES<..>.pdf  (solo se c'e' almeno un grafico: mai PDF vuoti)
plot_top_features <- function(results_tsv, qs_file, out_prefix, fdr_cutoff = 0.05, es_cutoff = 0.4) {
  tag <- sprintf("FDR%s_ES%s", fdr_cutoff, es_cutoff)
  out <- select_top_features(results_tsv, fdr_cutoff, es_cutoff)
  sel <- out$selected
  info <- list(n_total = out$n_total, n_selected = nrow(sel),
               n_male = sum(sel$higher_in == "male"), n_female = sum(sel$higher_in == "female"),
               n_plotted = 0L, pdf = NA_character_)

  if (nrow(sel) == 0) {
    message("   nessuna feature passa il filtro: nessun file scritto")
    return(info)
  }

  write.table(sel[, c("gene", "effectSize", "effectSizeStandardError", "effectSizeFDR", "numStudies", "higher_in")],
              paste0(out_prefix, "_top_features_", tag, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (!file.exists(qs_file)) {
    warning("oggetto .qs non trovato, scritta solo la tabella (nessun PDF): ", qs_file)
    return(info)
  }

  meta.res <- qs::qread(qs_file)
  pdf_file <- paste0(out_prefix, "_top_features_", tag, ".pdf")
  pdf(pdf_file)
  n_ok <- 0L
  for (g in sel$gene) {
    ok <- tryCatch({ MetaIntegrator::forestPlot(meta.res, g); TRUE },
                   error = function(e) { message("   forestPlot fallito per ", g, ": ", conditionMessage(e)); FALSE })
    if (ok) n_ok <- n_ok + 1L
  }
  dev.off()

  if (n_ok == 0L) {
    file.remove(pdf_file)     # niente PDF senza pagine
  } else {
    info$pdf <- pdf_file
  }
  info$n_plotted <- n_ok
  info
}
