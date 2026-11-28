## =====================================================================
## Differenze di sesso: PAZIENTI vs SANI  (interazione sesso x malattia sui pathway ssGSEA)
## Per ogni (contesto, gene set, malattia) con risultati in ENTRAMBE le analisi:
##   ES_pazienti = effect size M-vs-F nei pazienti   (Analysis/<ctx>_sex_differences/meta_analysis_pathways_activity/...)
##   ES_sani     = effect size M-vs-F nei sani       (Analysis/sex_differences_healthy_controls/<ctx>/...)
##   diff = ES_pazienti - ES_sani ;  z = diff / sqrt(SE_pazienti^2 + SE_sani^2) ;  p a due code ; FDR BH per (contesto, gene set, malattia)
## Le due stime vengono da campioni diversi (pazienti vs sani), quindi sono indipendenti.
## Classi (per |ES| > es_cutoff e FDR < 0.05 nell'analisi corrispondente):
##   disease_specific : differenza significativa (FDR_int < 0.05) e |ES_pazienti| > es_cutoff
##   shared_sex_effect: significativo (FDR<0.05, |ES|>es_cutoff) sia nei pazienti sia nei sani, stessa direzione, interazione non significativa
##   opposite         : significativo in entrambi con segno opposto
##   healthy_only     : significativo solo nei sani
##   patient_only_nonsig_interaction: significativo solo nei pazienti ma interazione NON significativa (da non sovrainterpretare)
## effect size > 0 = piu' alto nei maschi.
## Scrive SOLO in .../Analysis/sex_differences_healthy_controls/output_comparison/
## =====================================================================
library(ggplot2)

es_cutoff <- 0.4
root <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis"
out  <- file.path(root, "sex_differences_healthy_controls", "output_comparison")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
setwd(out)

read_one <- function(dir) {
  f <- list.files(dir, pattern = "_complete_results\\.tsv$", full.names = TRUE)
  if (length(f) != 1) return(NULL)
  r <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  r[!is.na(r$effectSize) & !is.na(r$effectSizeStandardError) & r$effectSizeStandardError > 0,
    c("gene", "effectSize", "effectSizeStandardError", "effectSizeFDR", "numStudies")]
}

res <- list(); cors <- list()
for (ctx in c("BLOOD", "TISSUE")) {
  for (gs in c("hallmarks", "kegg", "wikipathways")) {
    pbase <- file.path(root, paste0(ctx, "_sex_differences"), "meta_analysis_pathways_activity", gs)
    hbase <- file.path(root, "sex_differences_healthy_controls", ctx, gs)
    if (!dir.exists(pbase) || !dir.exists(hbase)) next
    for (dis in intersect(list.dirs(pbase, full.names = FALSE, recursive = FALSE),
                          list.dirs(hbase, full.names = FALSE, recursive = FALSE))) {
      p <- read_one(file.path(pbase, dis)); h <- read_one(file.path(hbase, dis))
      if (is.null(p) || is.null(h)) next
      m <- merge(p, h, by = "gene", suffixes = c("_pat", "_healthy"))
      if (nrow(m) < 10) next
      m$diff <- m$effectSize_pat - m$effectSize_healthy
      m$z    <- m$diff / sqrt(m$effectSizeStandardError_pat^2 + m$effectSizeStandardError_healthy^2)
      m$p_int   <- 2 * pnorm(-abs(m$z))
      m$FDR_int <- p.adjust(m$p_int, "BH")
      sp <- m$effectSizeFDR_pat < 0.05 & abs(m$effectSize_pat) > es_cutoff
      sh <- m$effectSizeFDR_healthy < 0.05 & abs(m$effectSize_healthy) > es_cutoff
      si <- m$FDR_int < 0.05
      m$class <- "none"
      m$class[sp & !sh & si]                                  <- "disease_specific"
      m$class[sp & !sh & !si]                                 <- "patient_only_nonsig_interaction"
      m$class[!sp & sh]                                       <- "healthy_only"
      m$class[sp & sh & sign(m$effectSize_pat) == sign(m$effectSize_healthy)] <- "shared_sex_effect"
      m$class[sp & sh & sign(m$effectSize_pat) != sign(m$effectSize_healthy)] <- "opposite"
      m$class[sp & sh & si & sign(m$effectSize_pat) == sign(m$effectSize_healthy)] <- "disease_specific"   # stessa direzione ma molto piu' forte/debole
      m$context <- ctx; m$gene_set <- gs; m$disease <- dis
      res[[length(res) + 1]] <- m
      cors[[length(cors) + 1]] <- data.frame(context = ctx, gene_set = gs, disease = dis, n_pathways = nrow(m),
        spearman_pat_vs_healthy = cor(m$effectSize_pat, m$effectSize_healthy, method = "spearman"),
        n_sig_pat = sum(sp), n_sig_healthy = sum(sh), n_sig_interaction = sum(si),
        n_disease_specific = sum(m$class == "disease_specific"), n_shared = sum(m$class == "shared_sex_effect"),
        n_opposite = sum(m$class == "opposite"), n_healthy_only = sum(m$class == "healthy_only"))
      message(ctx, " | ", gs, " | ", dis, ": ", nrow(m), " pathway")
    }
  }
}
if (length(res) == 0) stop("Nessuna coppia pazienti/sani trovata: lancia prima gli script dei sani.")
all <- do.call(rbind, res); summ <- do.call(rbind, cors)
colnames(all) <- sub("^gene$", "pathway", colnames(all))
all <- all[, c("context", "gene_set", "disease", "pathway", "effectSize_pat", "effectSizeFDR_pat", "numStudies_pat",
               "effectSize_healthy", "effectSizeFDR_healthy", "numStudies_healthy", "diff", "z", "p_int", "FDR_int", "class")]
write.table(all, "sex_effect_patients_vs_healthy_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all[all$class %in% c("disease_specific", "opposite"), ], "sex_effect_disease_specific_candidates.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ, "summary_patients_vs_healthy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

## scatter ES pazienti vs ES sani, un pannello per analisi
all$an <- paste(all$context, all$disease, all$gene_set, sep = " | ")
sc <- ggplot(all, aes(effectSize_healthy, effectSize_pat, color = class)) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_point(size = 0.8, alpha = 0.7) + facet_wrap(~ an, scales = "free", ncol = 4) +
  scale_color_manual(values = c(none = "grey80", disease_specific = "#D55E00", shared_sex_effect = "#0072B2",
                                opposite = "#CC79A7", healthy_only = "#009E73", patient_only_nonsig_interaction = "#E69F00")) +
  labs(x = "effect size M vs F nei SANI", y = "effect size M vs F nei PAZIENTI",
       title = "Differenze di sesso: pazienti vs sani (tratteggiata = uguali)") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
ggsave("scatter_patients_vs_healthy.pdf", sc, width = 14, height = max(6, 2.2 * ceiling(length(unique(all$an)) / 4)), limitsize = FALSE)
print("###############DONE! (patients vs healthy sex effects)")
