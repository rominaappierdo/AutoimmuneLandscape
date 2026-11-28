## =====================================================================
## UNICA heatmap (sangue + tessuto) dei pathway differenti tra maschi e femmine nei pazienti
## - legge i *_complete_results.tsv della meta-analisi ssGSEA M vs F (hallmarks, kegg, wikipathways)
## - tiene i pathway con |effectSize| > es_cutoff (e FDR < fdr_cutoff, se use_fdr = TRUE) in ALMENO una analisi
## - colonne = contesto_malattia (BLOOD_/TISSUE_), righe = pathway, un pannello per gene set
## - colore = direzione: ROSSO = piu' alto nei MASCHI (effectSize > 0), BLU = piu' alto nelle FEMMINE (< 0)
## - le celle che non passano il filtro in quella analisi restano grigio chiaro; le celle con * hanno FDR < 0.05
## - scrive SOLO in .../Analysis/sex_differences_top_plots/output_heatmap/
## =====================================================================
library(ggplot2)

es_cutoff  <- 0.4
fdr_cutoff <- 0.05
use_fdr    <- TRUE      # FALSE = solo |effectSize| > es_cutoff, senza filtro FDR

root     <- "/mnt/data3/appierdo/AutoimmuneLandscape"
out_root <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_top_plots/output_heatmap"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
setwd(out_root)

gene_sets <- c(hallmarks = "HALLMARKS", kegg = "KEGG", wikipathways = "WIKIPATHWAYS")

## ---- 1) leggo tutti i risultati ------------------------------------------------
all_res <- list()
for (ctx in c("BLOOD", "TISSUE")) {
  for (gs in names(gene_sets)) {
    base <- file.path(root, "Analysis", paste0(ctx, "_sex_differences"), "meta_analysis_pathways_activity", gs)
    if (!dir.exists(base)) next
    for (disease in list.dirs(base, full.names = FALSE, recursive = FALSE)) {
      for (rf in list.files(file.path(base, disease), pattern = "_complete_results\\.tsv$", full.names = TRUE)) {
        r <- read.delim(rf, stringsAsFactors = FALSE, check.names = FALSE)
        if (nrow(r) == 0) next
        all_res[[length(all_res) + 1]] <- data.frame(
          analysis = paste0(ctx, "_", disease), context = ctx, disease = disease,
          gene_set = gene_sets[[gs]], pathway = r$gene,
          effectSize = r$effectSize, FDR = r$effectSizeFDR, numStudies = r$numStudies,
          stringsAsFactors = FALSE)
      }
    }
  }
}
stopifnot(length(all_res) > 0)
d <- do.call(rbind, all_res)
d <- d[!is.na(d$effectSize) & !is.na(d$FDR), ]

## ---- 2) filtro ------------------------------------------------------------------
d$pass <- abs(d$effectSize) > es_cutoff & (!use_fdr | d$FDR < fdr_cutoff)
keep   <- unique(d[d$pass, c("gene_set", "pathway")])
if (nrow(keep) == 0) stop("Nessun pathway passa il filtro.")
d <- merge(d, keep, by = c("gene_set", "pathway"))
d$value <- ifelse(d$pass, d$effectSize, NA_real_)
d$star  <- ifelse(d$pass & d$FDR < 0.05, "*", "")

## ---- 3) ordine di righe e colonne -----------------------------------------------
d$pathway_label <- gsub("_", " ", sub("^(HALLMARK|KEGG|WP)_", "", d$pathway))
d$pathway_label <- ifelse(nchar(d$pathway_label) > 60, paste0(substr(d$pathway_label, 1, 57), "..."), d$pathway_label)
d$row_id <- paste(d$gene_set, d$pathway_label, sep = "||")
# righe: dal piu' "maschile" al piu' "femminile" (media dei valori che passano) dentro ogni gene set
rowmean <- tapply(d$value, d$row_id, function(x) mean(x, na.rm = TRUE))
d$row_id <- factor(d$row_id, levels = names(sort(rowmean, decreasing = FALSE)))
d$analysis <- factor(d$analysis, levels = sort(unique(d$analysis)))   # BLOOD_* poi TISSUE_*
d$gene_set <- factor(d$gene_set, levels = unname(gene_sets))

## ---- 4) tabella della matrice (pathway x analisi) -------------------------------
mat <- reshape(d[, c("gene_set", "pathway", "analysis", "value")], idvar = c("gene_set", "pathway"),
               timevar = "analysis", direction = "wide")
colnames(mat) <- sub("^value\\.", "", colnames(mat))
write.table(mat, "heatmap_matrix_pathways_ES0.4.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")

## ---- 5) heatmap -----------------------------------------------------------------
lim <- max(abs(d$value), na.rm = TRUE)
p <- ggplot(d, aes(x = analysis, y = row_id, fill = value)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 3, vjust = 0.75) +
  facet_grid(gene_set ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), na.value = "grey93",
                       name = "Effect size\n(rosso = maschi,\nblu = femmine)") +
  labs(x = NULL, y = NULL,
       title = "Differenze di attivita' dei pathway tra maschi e femmine (solo pazienti)",
       subtitle = sprintf("|effect size| > %s%s; grigio = non passa il filtro; * = FDR < 0.05",
                          es_cutoff, if (use_fdr) sprintf(" e FDR < %s", fdr_cutoff) else "")) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"))

n_rows <- length(unique(d$row_id)); n_cols <- length(unique(d$analysis))
w <- max(7, 4 + 0.45 * n_cols)
h <- max(5, 2 + 0.17 * n_rows)
ggsave("heatmap_pathways_sex_differences_ES0.4.pdf", p, width = w, height = h, limitsize = FALSE)
ggsave("heatmap_pathways_sex_differences_ES0.4.png", p, width = w, height = h, dpi = 200, limitsize = FALSE)

message(n_rows, " pathway, ", n_cols, " analisi -> ", out_root)
print("###############DONE! (heatmap pathways M vs F)")
