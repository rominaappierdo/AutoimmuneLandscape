## =====================================================================
## META-ANALISI PAN-AUTOIMMUNE DELLE DIFFERENZE DI SESSO (M vs F nei pazienti)
## Secondo livello di meta-analisi: ogni MALATTIA e' una unita'. Si combinano gli effect size gia' calcolati
## DENTRO ogni dataset (e poi meta-analizzati per malattia): nessun campione di GSE diversi viene mai messo
## nella stessa matrice, quindi i batch effect tra dataset non entrano nel confronto.
##
## Input (gia' prodotti, solo lettura):
##   geni     : Analysis/<ctx>_sex_differences/meta_analysis_gene_expression/<malattia>/*male_vs_female*complete_results.tsv
##   pathway  : Analysis/<ctx>_sex_differences/meta_analysis_pathways_activity/<hallmarks|kegg|wikipathways>/<malattia>/*complete_results.tsv
## Metodo per ogni feature (gene o pathway) e per ogni contesto (BLOOD, TISSUE separati):
##   - unita' = malattie con numStudies >= min_studies_per_disease (default 2) e SE valido
##   - random effects DerSimonian-Laird tra malattie (tau2, I2, Q) ; p-value primario con aggiustamento Hartung-Knapp (t, k-1 gdl,
##     varianza non inferiore a 1: conservativo, adatto a pochi k) ; FDR di Benjamini-Hochberg per (contesto, tipo di feature)
##   - leave-one-disease-out: p massimo (z) rimuovendo ciascuna malattia
##   - coerenza di direzione: quante malattie vanno nello stesso verso
##   - flag "robust": k >= min_k, FDR_HK < 0.05, |ES| > es_cutoff e >= 80% delle malattie con stesso segno
## effect size > 0 = piu' alto nei MASCHI.
## Controllo di qualita' (controllo positivo): geni Y e XIST devono risultare fortemente M>F / F>M in tutte le malattie.
## Opzionale: se esistono i risultati dei SANI (cartella sex_differences_healthy_controls) calcola lo stesso pooled nei sani
##   e la differenza pazienti - sani per i pathway.
## Scrive SOLO in: Analysis/sex_differences_pan_autoimmune/output/
## =====================================================================
library(ggplot2)

root      <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis"
out_root  <- file.path(root, "sex_differences_pan_autoimmune", "output")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
setwd(out_root)

min_studies_per_disease <- 2     # studi minimi dentro ogni malattia
min_k       <- 3                 # malattie minime per poolare una feature
es_cutoff   <- 0.2               # soglia sul pooled effect size per il flag "robust" (il pooled e' una media: e' piu' basso dei singoli)
fdr_cutoff  <- 0.05
same_sign_frac <- 0.8
gene_sets   <- c("hallmarks", "kegg", "wikipathways")

SEX_LINKED <- c("XIST","TSIX","JPX","FTX","KDM6A","KDM5C","DDX3X","EIF2S3","EIF1AX","ZFX","PRKX","USP9X","UBA1","SMC1A","STS","RPS4X",
                "SRY","RPS4Y1","RPS4Y2","ZFY","AMELY","TBL1Y","PCDH11Y","TGIF2LY","USP9Y","DDX3Y","UTY","TMSB4Y","NLGN4Y","KDM5D","EIF1AY",
                "PRKY","TXLNGY","LINC00278","ZFY-AS1","TTTY15","TTTY14","BCORP1","CYorf15A","CYorf15B","SHROOM2","ZRSR2","CD99","SLC25A6","ASMTL","GTPBP6","PLCXD1","SHOX","ZBED1","IL9R","WASH6P")
QC_GENES <- c("XIST","RPS4Y1","DDX3Y","KDM5D","UTY","USP9Y","EIF1AY","ZFY","TXLNGY")

## ---------- funzioni -----------------------------------------------------------------
read_res <- function(f) {
  r <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  r <- r[!is.na(r$effectSize) & !is.na(r$effectSizeStandardError) & r$effectSizeStandardError > 0 &
           !is.na(r$numStudies) & r$numStudies >= min_studies_per_disease, c("gene", "effectSize", "effectSizeStandardError", "numStudies")]
  r[!duplicated(r$gene), ]
}

# matrici feature x malattia
build_mats <- function(files) {            # files: named vector disease -> path
  L <- lapply(files, read_res)
  feats <- unique(unlist(lapply(L, `[[`, "gene")))
  ES <- SE <- matrix(NA_real_, length(feats), length(L), dimnames = list(feats, names(L)))
  for (d in names(L)) { i <- match(L[[d]]$gene, feats); ES[i, d] <- L[[d]]$effectSize; SE[i, d] <- L[[d]]$effectSizeStandardError }
  list(ES = ES, SE = SE)
}

# random effects DerSimonian-Laird + Hartung-Knapp, vettorializzato sulle righe
pool_re <- function(ES, SE) {
  nas <- is.na(ES) | is.na(SE)
  ES0 <- ES; ES0[nas] <- 0
  W <- 1 / SE^2; W[nas] <- 0
  k  <- rowSums(!nas)
  sw <- rowSums(W); fe <- rowSums(W * ES0) / sw
  Q  <- rowSums(W * (ES0 - fe)^2)
  C  <- sw - rowSums(W^2) / sw
  tau2 <- pmax(0, (Q - (k - 1)) / C)
  Ws <- 1 / (SE^2 + tau2); Ws[nas] <- 0
  mu <- rowSums(Ws * ES0) / rowSums(Ws)
  se <- sqrt(1 / rowSums(Ws))
  z  <- mu / se; p <- 2 * pnorm(-abs(z))
  qhk <- rowSums(Ws * (ES0 - mu)^2) / pmax(k - 1, 1)
  se_hk <- sqrt(pmax(1, qhk) / rowSums(Ws))
  p_hk  <- 2 * pt(-abs(mu / se_hk), pmax(k - 1, 1))
  I2 <- ifelse(Q > 0, pmax(0, (Q - (k - 1)) / Q), 0)
  out <- data.frame(k = k, ES = mu, SE = se, SE_HK = se_hk, z = z, p = p, p_HK = p_hk,
                    tau2 = tau2, I2 = I2, Q = Q, pQ = pchisq(Q, pmax(k - 1, 1), lower.tail = FALSE),
                    n_pos = rowSums(ES > 0, na.rm = TRUE), n_neg = rowSums(ES < 0, na.rm = TRUE),
                    row.names = rownames(ES))
  out[k < min_k, c("ES", "SE", "SE_HK", "z", "p", "p_HK", "tau2", "I2", "Q", "pQ")] <- NA
  out
}

pan_analysis <- function(m, label) {
  if (ncol(m$ES) < min_k) { message("   ", label, ": meno di ", min_k, " malattie, salto"); return(NULL) }
  res <- pool_re(m$ES, m$SE)
  # leave-one-disease-out: p massimo (z) rimuovendo ciascuna malattia
  loo_p <- matrix(NA_real_, nrow(res), ncol(m$ES))
  if (ncol(m$ES) > min_k) {            # con k = min_k il LOO lascerebbe < min_k malattie
    for (j in seq_len(ncol(m$ES))) {
      sub <- suppressWarnings(pool_re_k(m$ES[, -j, drop = FALSE], m$SE[, -j, drop = FALSE]))
      loo_p[, j] <- sub
    }
  }
  res$p_LOO_max <- suppressWarnings(apply(loo_p, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)))
  res$FDR_HK <- p.adjust(res$p_HK, "BH")
  res$FDR    <- p.adjust(res$p, "BH")
  res$same_sign_frac <- pmax(res$n_pos, res$n_neg) / res$k
  res$direction <- ifelse(res$ES > 0, "male", "female")
  res$robust <- !is.na(res$FDR_HK) & res$FDR_HK < fdr_cutoff & abs(res$ES) > es_cutoff & res$same_sign_frac >= same_sign_frac
  res$feature <- rownames(res)
  res$diseases <- apply(!is.na(m$ES), 1, function(x) paste(colnames(m$ES)[x], collapse = ","))
  cbind(res, setNames(as.data.frame(m$ES), paste0("ES_", colnames(m$ES))))
}
# p (z) del random effects senza HK, per il LOO (richiede k >= 2 nella matrice ridotta; qui sempre >= min_k)
pool_re_k <- function(ES, SE) {
  nas <- is.na(ES) | is.na(SE); ES0 <- ES; ES0[nas] <- 0
  W <- 1 / SE^2; W[nas] <- 0; k <- rowSums(!nas); sw <- rowSums(W)
  fe <- rowSums(W * ES0) / sw; Q <- rowSums(W * (ES0 - fe)^2); C <- sw - rowSums(W^2) / sw
  tau2 <- pmax(0, (Q - (k - 1)) / C); Ws <- 1 / (SE^2 + tau2); Ws[nas] <- 0
  mu <- rowSums(Ws * ES0) / rowSums(Ws); p <- 2 * pnorm(-abs(mu / sqrt(1 / rowSums(Ws))))
  p[k < 2] <- NA; p
}

forest_pages <- function(m, res, feats, file, title, per_page = 12) {
  if (length(feats) == 0) return(invisible())
  long <- do.call(rbind, lapply(feats, function(f) {
    d <- data.frame(feature = f, unit = colnames(m$ES), es = m$ES[f, ], se = m$SE[f, ], pooled = FALSE)
    d <- d[!is.na(d$es), ]
    rbind(d, data.frame(feature = f, unit = "POOLED", es = res[f, "ES"], se = res[f, "SE_HK"], pooled = TRUE))
  }))
  long$lo <- long$es - 1.96 * long$se; long$hi <- long$es + 1.96 * long$se
  long$feature <- factor(long$feature, levels = feats)
  long$unit <- factor(long$unit, levels = rev(c(colnames(m$ES), "POOLED")))
  pdf(file, width = 11, height = 8.5)
  for (s in seq(1, length(feats), by = per_page)) {
    sel <- feats[s:min(length(feats), s + per_page - 1)]
    p <- ggplot(long[long$feature %in% sel, ], aes(es, unit, color = pooled, shape = pooled)) +
      geom_vline(xintercept = 0, color = "grey70") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) + geom_point(size = 2) +
      scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "#B2182B"), guide = "none") +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
      facet_wrap(~ feature, scales = "free_x", ncol = 3, labeller = labeller(feature = function(x) substr(gsub("_", " ", x), 1, 45))) +
      labs(x = "effect size M vs F (>0 = maschi)", y = NULL, title = title) + theme_bw(base_size = 8)
    print(p)
  }
  dev.off()
}

## ---------- 1) analisi per contesto e tipo di feature ---------------------------------
all_res <- list(); qc <- list()
for (ctx in c("BLOOD", "TISSUE")) {
  ctx_dir <- file.path(root, paste0(ctx, "_sex_differences"))
  sets <- list(genes = file.path(ctx_dir, "meta_analysis_gene_expression"))
  for (gs in gene_sets) sets[[gs]] <- file.path(ctx_dir, "meta_analysis_pathways_activity", gs)
  for (ft in names(sets)) {
    base <- sets[[ft]]
    if (!dir.exists(base)) next
    files <- c()
    for (dis in list.dirs(base, full.names = FALSE, recursive = FALSE)) {
      f <- list.files(file.path(base, dis), pattern = "male_vs_female.*complete_results\\.tsv$", full.names = TRUE)
      f <- f[!grepl("backup|healthy", basename(f))]
      if (length(f) == 1) files[dis] <- f else if (length(f) > 1) message("   ", ctx, "/", ft, "/", dis, ": piu' file, salto")
    }
    if (length(files) < min_k) next
    message("### ", ctx, " | ", ft, " | ", length(files), " malattie: ", paste(names(files), collapse = ","))
    m <- build_mats(files)
    res <- pan_analysis(m, paste(ctx, ft))
    if (is.null(res)) next
    res$context <- ctx; res$feature_type <- ft
    res$sex_linked <- res$feature %in% SEX_LINKED
    res <- res[order(res$p_HK), ]
    all_res[[paste(ctx, ft)]] <- res
    write.table(res, sprintf("pan_autoimmune_%s_%s.tsv", ctx, ft), sep = "\t", quote = FALSE, row.names = FALSE)
    if (ft == "genes") {
      q <- res[res$feature %in% QC_GENES, ]
      # formato lungo (le malattie cambiano tra sangue e tessuto, quindi non si puo' fare rbind in formato largo)
      qc[[ctx]] <- do.call(rbind, lapply(grep("^ES_", colnames(q), value = TRUE), function(cn)
        data.frame(context = ctx, gene = q$feature, disease = sub("^ES_", "", cn), ES = q[[cn]], stringsAsFactors = FALSE)))
    }
    top <- rownames(res)[which(res$robust)]
    top <- top[seq_len(min(length(top), 60))]
    if (ft != "genes" && length(top) > 0)
      forest_pages(m, res, top, sprintf("forest_pan_autoimmune_%s_%s.pdf", ctx, ft), paste(ctx, ft, "- pan-autoimmune M vs F"))
    if (ft == "genes") {      # per i geni solo quelli NON sessuo-legati
      top_ns <- rownames(res)[which(res$robust & !res$sex_linked)]
      top_ns <- top_ns[seq_len(min(length(top_ns), 36))]
      forest_pages(m, res, top_ns, sprintf("forest_pan_autoimmune_%s_genes_non_sex_linked.pdf", ctx), paste(ctx, "geni (non sessuo-legati) - pan-autoimmune M vs F"))
    }
  }
}
if (length(all_res) == 0) stop("Nessun risultato trovato.")

## ---------- 2) riepilogo e QC ---------------------------------------------------------
summ <- do.call(rbind, lapply(all_res, function(r) data.frame(
  context = r$context[1], feature_type = r$feature_type[1], n_diseases = max(r$k),
  n_features_pooled = sum(!is.na(r$p_HK)), n_robust = sum(r$robust),
  n_robust_male = sum(r$robust & r$direction == "male"), n_robust_female = sum(r$robust & r$direction == "female"),
  n_robust_non_sex_linked = sum(r$robust & !r$sex_linked))))
write.table(summ, "summary_pan_autoimmune.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
if (length(qc) > 0) {
  qc_all <- do.call(rbind, qc)
  qc_all$expected <- ifelse(qc_all$gene == "XIST", "female (ES<0)", "male (ES>0)")
  qc_all$ok <- ifelse(is.na(qc_all$ES), NA, ifelse(qc_all$gene == "XIST", qc_all$ES < 0, qc_all$ES > 0))
  qc_all <- qc_all[!is.na(qc_all$ES), ]
  write.table(qc_all, "QC_sex_genes_positive_controls.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  message("QC: i geni Y (RPS4Y1, DDX3Y, KDM5D, ...) devono essere > 0 e XIST < 0 in TUTTE le malattie (guarda QC_sex_genes_positive_controls.tsv)")
}

## ---------- 3) pazienti vs sani (opzionale, solo pathway) -----------------------------
hroot <- file.path(root, "sex_differences_healthy_controls")
if (dir.exists(hroot)) {
  cmp <- list()
  for (ctx in c("BLOOD", "TISSUE")) for (gs in gene_sets) {
    key <- paste(ctx, gs); hb <- file.path(hroot, ctx, gs)
    if (!dir.exists(hb) || is.null(all_res[[key]])) next
    files <- c()
    for (dis in list.dirs(hb, full.names = FALSE, recursive = FALSE)) {
      f <- list.files(file.path(hb, dis), pattern = "complete_results\\.tsv$", full.names = TRUE)
      if (length(f) == 1) files[dis] <- f
    }
    if (length(files) < min_k) { message("   sani ", key, ": meno di ", min_k, " malattie, salto"); next }
    mh <- build_mats(files); rh <- pool_re(mh$ES, mh$SE)
    rp <- all_res[[key]]; common <- intersect(rownames(rp)[!is.na(rp$ES)], rownames(rh)[!is.na(rh$ES)])
    d <- data.frame(context = ctx, gene_set = gs, pathway = common,
                    ES_patients = rp[common, "ES"], SE_patients = rp[common, "SE_HK"], k_patients = rp[common, "k"],
                    ES_healthy = rh[common, "ES"], SE_healthy = rh[common, "SE_HK"], k_healthy = rh[common, "k"])
    d$diff <- d$ES_patients - d$ES_healthy
    d$z_int <- d$diff / sqrt(d$SE_patients^2 + d$SE_healthy^2); d$p_int <- 2 * pnorm(-abs(d$z_int)); d$FDR_int <- p.adjust(d$p_int, "BH")
    cmp[[key]] <- d[order(d$p_int), ]
  }
  if (length(cmp) > 0) write.table(do.call(rbind, cmp), "pan_autoimmune_patients_vs_healthy_pathways.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}
print("###############DONE! (pan-autoimmune sex differences)")
