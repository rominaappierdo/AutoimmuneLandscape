## =====================================================================
## LA FIRMA DI MALATTIA E' DIVERSA TRA MASCHI E FEMMINE?   (interazione sesso x malattia)
## Input (prodotti dagli script perform_meta_analysis_sexstrat_*.R): per ogni malattia, la meta-analisi malati-vs-sani
##   fatta solo nei maschi e solo nelle femmine (geni, Hallmarks, KEGG, WikiPathways; sangue e tessuto).
## Livello 1 - per malattia e per feature:   diff = ES_M - ES_F ;  z = diff / sqrt(SE_M^2 + SE_F^2) ;  FDR (BH) dentro ogni analisi.
##   (le due stime vengono da persone diverse: sono indipendenti)
##   Classi: male_specific / female_specific  (significativo in un sesso, |ES| > es_cutoff, e interazione FDR < 0.05),
##           shared (significativo in entrambi, stesso segno), opposite (significativo in entrambi, segno opposto),
##           one_sex_only_no_interaction (significativo in un sesso ma NON dimostrabile come diverso: non e' un risultato),
##           none.
## Livello 2 - pan-autoimmune: la differenza (ES_M - ES_F) di ogni malattia e' un'unita'; random effects DerSimonian-Laird con
##   Hartung-Knapp (stesso metodo della meta-analisi pan-autoimmune precedente), p, FDR, I2, coerenza di segno, leave-one-out.
##   Si poolano anche separatamente la firma nei maschi e nelle femmine.
## effect size > 0 = piu' alto nei MALATI (dentro quel sesso). diff > 0 = la firma e' piu' forte/positiva nei maschi.
## Scrive SOLO in: Analysis/sex_specific_disease_signature/output/
## =====================================================================
library(ggplot2)

SS  <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_specific_disease_signature"
out <- file.path(SS, "output"); dir.create(out, recursive = TRUE, showWarnings = FALSE); setwd(out)

min_studies <- 2; es_cutoff <- 0.4; fdr_cutoff <- 0.05
min_k <- 3; pan_es_cutoff <- 0.15; same_sign_frac <- 0.8
gene_sets <- c("hallmarks", "kegg", "wikipathways")
SEX_LINKED <- c("XIST","TSIX","JPX","FTX","KDM6A","KDM5C","DDX3X","EIF2S3","EIF1AX","ZFX","PRKX","USP9X","UBA1","SMC1A","STS","RPS4X","SRY",
                "RPS4Y1","RPS4Y2","ZFY","AMELY","TBL1Y","PCDH11Y","TGIF2LY","USP9Y","DDX3Y","UTY","TMSB4Y","NLGN4Y","KDM5D","EIF1AY","PRKY",
                "TXLNGY","LINC00278","ZFY-AS1","TTTY15","TTTY14","SHROOM2","ZRSR2","CD99","SLC25A6","ASMTL","GTPBP6","PLCXD1","SHOX","ZBED1","IL9R")

read_res <- function(f) {
  r <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  r <- r[!is.na(r$effectSize) & !is.na(r$effectSizeStandardError) & r$effectSizeStandardError > 0 &
         !is.na(r$numStudies) & r$numStudies >= min_studies, c("gene", "effectSize", "effectSizeStandardError", "effectSizeFDR", "numStudies")]
  r[!duplicated(r$gene), ]
}
find_res <- function(base, tag) {
  f <- list.files(base, pattern = paste0(tag, ".*complete_results\\.tsv$"), full.names = TRUE)
  f <- f[!grepl("backup", basename(f))]
  if (length(f) == 1) f else NULL
}

# random effects DerSimonian-Laird + Hartung-Knapp, vettorializzato (identico alla meta-analisi pan-autoimmune)
pool_re <- function(ES, SE) {
  nas <- is.na(ES) | is.na(SE); ES0 <- ES; ES0[nas] <- 0
  W <- 1 / SE^2; W[nas] <- 0; k <- rowSums(!nas); sw <- rowSums(W)
  fe <- rowSums(W * ES0) / sw; Q <- rowSums(W * (ES0 - fe)^2); C <- sw - rowSums(W^2) / sw
  tau2 <- pmax(0, (Q - (k - 1)) / C); Ws <- 1 / (SE^2 + tau2); Ws[nas] <- 0
  mu <- rowSums(Ws * ES0) / rowSums(Ws); se <- sqrt(1 / rowSums(Ws))
  qhk <- rowSums(Ws * (ES0 - mu)^2) / pmax(k - 1, 1); se_hk <- sqrt(pmax(1, qhk) / rowSums(Ws))
  out <- data.frame(k = k, ES = mu, SE = se, SE_HK = se_hk, p = 2 * pnorm(-abs(mu / se)),
                    p_HK = 2 * pt(-abs(mu / se_hk), pmax(k - 1, 1)), tau2 = tau2,
                    I2 = ifelse(Q > 0, pmax(0, (Q - (k - 1)) / Q), 0),
                    n_pos = rowSums(ES > 0, na.rm = TRUE), n_neg = rowSums(ES < 0, na.rm = TRUE), row.names = rownames(ES))
  out[k < min_k, c("ES", "SE", "SE_HK", "p", "p_HK", "tau2", "I2")] <- NA
  out
}
loo_p <- function(ES, SE) {
  if (ncol(ES) <= min_k) return(rep(NA_real_, nrow(ES)))
  m <- sapply(seq_len(ncol(ES)), function(j) { r <- pool_re(ES[, -j, drop = FALSE], SE[, -j, drop = FALSE]); r$p })
  suppressWarnings(apply(m, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)))
}

## ---------- Livello 1: per malattia -------------------------------------------------------
l1 <- list(); summ <- list()
set_dirs <- function(sex, ctx, ft) {            # ritorna lista disease -> cartella
  base <- if (ft == "genes") file.path(SS, sex, paste0(ctx, "_genes")) else file.path(SS, sex, paste0(ctx, "_pathways"), ft)
  if (!dir.exists(base)) return(list())
  ds <- list.dirs(base, full.names = FALSE, recursive = FALSE)
  setNames(file.path(base, ds), ds)
}
for (ctx in c("BLOOD", "TISSUE")) for (ft in c("genes", gene_sets)) {
  dm <- set_dirs("males", ctx, ft); df <- set_dirs("females", ctx, ft)
  for (dis in intersect(names(dm), names(df))) {
    fm <- find_res(dm[[dis]], "ONLYMALES"); ff <- find_res(df[[dis]], "ONLYFEMALES")
    if (is.null(fm) || is.null(ff)) { message("   ", ctx, "/", ft, "/", dis, ": mancano i risultati di un sesso, salto"); next }
    m <- merge(read_res(fm), read_res(ff), by = "gene", suffixes = c("_M", "_F"))
    if (nrow(m) < 10) next
    m$diff <- m$effectSize_M - m$effectSize_F
    m$z_int <- m$diff / sqrt(m$effectSizeStandardError_M^2 + m$effectSizeStandardError_F^2)
    m$p_int <- 2 * pnorm(-abs(m$z_int)); m$FDR_int <- p.adjust(m$p_int, "BH")
    sM <- m$effectSizeFDR_M < fdr_cutoff & abs(m$effectSize_M) > es_cutoff
    sF <- m$effectSizeFDR_F < fdr_cutoff & abs(m$effectSize_F) > es_cutoff
    si <- m$FDR_int < fdr_cutoff
    m$class <- "none"
    m$class[(sM | sF) & !si] <- "one_sex_only_no_interaction"
    m$class[sM & sF & sign(m$effectSize_M) == sign(m$effectSize_F)] <- "shared"
    m$class[sM & sF & sign(m$effectSize_M) != sign(m$effectSize_F)] <- "opposite"
    m$class[sM & !sF & si] <- "male_specific"
    m$class[sF & !sM & si] <- "female_specific"
    m$context <- ctx; m$feature_type <- ft; m$disease <- dis; m$sex_linked <- m$gene %in% SEX_LINKED
    l1[[length(l1) + 1]] <- m
    summ[[length(summ) + 1]] <- data.frame(context = ctx, feature_type = ft, disease = dis, n_features = nrow(m),
      spearman_ES_M_vs_F = cor(m$effectSize_M, m$effectSize_F, method = "spearman"),
      n_sig_M = sum(sM), n_sig_F = sum(sF), n_shared = sum(m$class == "shared"), n_opposite = sum(m$class == "opposite"),
      n_male_specific = sum(m$class == "male_specific"), n_female_specific = sum(m$class == "female_specific"),
      n_interaction_fdr05 = sum(si), n_one_sex_only_no_interaction = sum(m$class == "one_sex_only_no_interaction"))
    message(ctx, " | ", ft, " | ", dis, ": ", nrow(m), " feature")
  }
}
if (length(l1) == 0) stop("Nessuna coppia maschi/femmine trovata: lancia prima gli script perform_meta_analysis_sexstrat_*.R")
cols <- c("context", "feature_type", "disease", "gene", "sex_linked", "effectSize_M", "effectSizeStandardError_M", "effectSizeFDR_M", "numStudies_M",
          "effectSize_F", "effectSizeStandardError_F", "effectSizeFDR_F", "numStudies_F", "diff", "z_int", "p_int", "FDR_int", "class")
all1 <- do.call(rbind, lapply(l1, function(x) x[, cols])); colnames(all1)[colnames(all1) == "gene"] <- "feature"
summ <- do.call(rbind, summ)
write.table(summ, "summary_sex_specific_by_disease.tsv", sep = "\t", quote = FALSE, row.names = FALSE); print(summ)
write.table(all1[all1$class %in% c("male_specific", "female_specific", "opposite"), ], "candidates_sex_specific_by_disease.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
for (ft in unique(all1$feature_type)) for (ctx in unique(all1$context)) {
  x <- all1[all1$feature_type == ft & all1$context == ctx, ]
  if (nrow(x) > 0) write.table(x[x$class != "none", ], sprintf("sex_specific_by_disease_%s_%s_nonnull.tsv", ctx, ft), sep = "\t", quote = FALSE, row.names = FALSE)
}

# scatter: firma nei maschi vs firma nelle femmine (per i geni solo i non-"none" per leggerezza)
sc <- all1[all1$feature_type != "genes" | all1$class != "none", ]
sc$an <- paste(sc$context, sc$disease, sc$feature_type, sep = " | ")
p <- ggplot(sc, aes(effectSize_F, effectSize_M, color = class)) +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") + geom_abline(linetype = 2, color = "grey60") +
  geom_point(size = 0.7, alpha = 0.7) + facet_wrap(~ an, scales = "free", ncol = 4) +
  scale_color_manual(values = c(none = "grey80", shared = "#0072B2", opposite = "#CC79A7", male_specific = "#D55E00",
                                female_specific = "#009E73", one_sex_only_no_interaction = "#E69F00")) +
  labs(x = "firma di malattia nelle FEMMINE (ES malati vs sani)", y = "firma di malattia nei MASCHI", title = "Firma di malattia: maschi vs femmine (tratteggiata = uguali)") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
ggsave("scatter_disease_signature_M_vs_F.pdf", p, width = 14, height = max(6, 2.3 * ceiling(length(unique(sc$an)) / 4)), limitsize = FALSE)

## ---------- Livello 2: pan-autoimmune ---------------------------------------------------------
forest <- function(long, feats, file, title) {
  if (length(feats) == 0) return(invisible())
  long <- long[long$feature %in% feats, ]; long$feature <- factor(long$feature, levels = feats)
  long$lo <- long$es - 1.96 * long$se; long$hi <- long$es + 1.96 * long$se
  long$unit <- factor(long$unit, levels = rev(unique(long$unit)))
  pdf(file, width = 11, height = 8.5)
  for (s in seq(1, length(feats), by = 12)) {
    sel <- feats[s:min(length(feats), s + 11)]
    print(ggplot(long[long$feature %in% sel, ], aes(es, unit, color = pooled, shape = pooled)) + geom_vline(xintercept = 0, color = "grey70") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) + geom_point(size = 2) +
      scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "#B2182B"), guide = "none") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
      facet_wrap(~ feature, scales = "free_x", ncol = 3, labeller = labeller(feature = function(x) substr(gsub("_", " ", x), 1, 45))) +
      labs(x = "differenza M - F della firma di malattia (>0 = piu' forte nei maschi)", y = NULL, title = title) + theme_bw(base_size = 8))
  }
  dev.off()
}
pan_all <- list()
for (ctx in unique(all1$context)) for (ft in unique(all1$feature_type)) {
  x <- all1[all1$context == ctx & all1$feature_type == ft, ]
  ds <- unique(x$disease); if (length(ds) < min_k) { message("   pan ", ctx, " ", ft, ": solo ", length(ds), " malattie, salto"); next }
  feats <- unique(x$feature)
  mk <- function(v) { m <- matrix(NA_real_, length(feats), length(ds), dimnames = list(feats, ds))
                      for (d in ds) { y <- x[x$disease == d, ]; m[y$feature, d] <- y[[v]] }; m }
  ESd <- mk("diff"); SEd <- sqrt(mk("effectSizeStandardError_M")^2 + mk("effectSizeStandardError_F")^2)
  r <- pool_re(ESd, SEd)
  pm <- pool_re(mk("effectSize_M"), mk("effectSizeStandardError_M")); pf <- pool_re(mk("effectSize_F"), mk("effectSizeStandardError_F"))
  r$p_LOO_max <- loo_p(ESd, SEd)
  r$FDR_HK <- p.adjust(r$p_HK, "BH"); r$same_sign_frac <- pmax(r$n_pos, r$n_neg) / r$k
  r$ES_pooled_M <- pm$ES; r$ES_pooled_F <- pf$ES
  r$feature <- rownames(r); r$context <- ctx; r$feature_type <- ft; r$sex_linked <- r$feature %in% SEX_LINKED
  r$direction <- ifelse(r$ES > 0, "stronger_in_males", "stronger_in_females")
  r$robust <- !is.na(r$FDR_HK) & r$FDR_HK < fdr_cutoff & abs(r$ES) > pan_es_cutoff & r$same_sign_frac >= same_sign_frac
  r <- cbind(r, setNames(as.data.frame(ESd), paste0("diff_", colnames(ESd)))); r <- r[order(r$p_HK), ]
  write.table(r, sprintf("pan_sex_by_disease_%s_%s.tsv", ctx, ft), sep = "\t", quote = FALSE, row.names = FALSE)
  pan_all[[paste(ctx, ft)]] <- r
  top <- r$feature[r$robust & !r$sex_linked]; top <- top[seq_len(min(length(top), 60))]
  long <- do.call(rbind, lapply(top, function(f) {
    d <- data.frame(feature = f, unit = colnames(ESd), es = ESd[f, ], se = SEd[f, ], pooled = FALSE); d <- d[!is.na(d$es), ]
    rbind(d, data.frame(feature = f, unit = "POOLED", es = r[f, "ES"], se = r[f, "SE_HK"], pooled = TRUE)) }))
  if (length(top) > 0) forest(long, top, sprintf("forest_pan_sex_by_disease_%s_%s.pdf", ctx, ft), paste(ctx, ft, "- interazione sesso x malattia"))
}
if (length(pan_all) > 0) {
  ps <- do.call(rbind, lapply(pan_all, function(r) data.frame(context = r$context[1], feature_type = r$feature_type[1], n_diseases = max(r$k),
          n_pooled = sum(!is.na(r$p_HK)), n_robust = sum(r$robust), n_robust_non_sex_linked = sum(r$robust & !r$sex_linked),
          n_stronger_in_males = sum(r$robust & r$direction == "stronger_in_males" & !r$sex_linked),
          n_stronger_in_females = sum(r$robust & r$direction == "stronger_in_females" & !r$sex_linked))))
  write.table(ps, "summary_pan_sex_by_disease.tsv", sep = "\t", quote = FALSE, row.names = FALSE); print(ps)
}
print("###############DONE! (sex-specific disease signature)")
