setwd("/mnt/data3/appierdo/PhD-project/Analysis/TISSUE/meta_analysis_pathways_activity/hallmarks/ipf")

metaboj <- qs::qread("LUNG_gsva_ipf_vs_healthy_tissue_metaintegrator_object.qs")
#MetaIntegrator::forestPlot(metaObject = metaboj, geneName = )

tsv <- read.csv("LUNG_gsva_ipf_vs_healthy_tissue_complete_results.tsv", sep = '\t', row.names = 15)
tsv <- tsv["effectSize"]
ComplexHeatmap::Heatmap(as.matrix(tsv))
