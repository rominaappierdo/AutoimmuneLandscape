setwd("/mnt/data3/appierdo/PhD-project/Analysis/TISSUE/meta_analysis_cell_proportions/ipf")

df <- data.table::fread("LUNG_ipf_vs_healthy_tissue_complete_results.tsv")
df <- as.data.frame(df[, c("effectSize", "gene")])
rownames(df) <- df$gene
df$gene <- NULL

ComplexHeatmap::Heatmap(as.matrix(df))

metaobj <-qs::qread("LUNG_ipf_vs_healthy_tissue_metaintegrator_object.qs")

MetaIntegrator::forestPlot(metaObject = metaobj, geneName = "plasmacytoid_dendritic_cell")
