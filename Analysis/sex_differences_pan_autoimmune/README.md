# Meta-analisi pan-autoimmune delle differenze M vs F
Posizione: /mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_pan_autoimmune/
Lancio (R): source("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_pan_autoimmune/pan_autoimmune_sex_meta_analysis.R")
Non rilancia nessuna meta-analisi e non serve il dataset: legge i *_complete_results.tsv gia' prodotti. Scrive solo in output/.
Output: pan_autoimmune_<BLOOD|TISSUE>_<genes|hallmarks|kegg|wikipathways>.tsv (tutte le feature, colonna robust), forest_pan_autoimmune_*.pdf,
summary_pan_autoimmune.tsv, QC_sex_genes_positive_controls.tsv, e (se hai gia' i sani) pan_autoimmune_patients_vs_healthy_pathways.tsv.
Rationale: vedi la risposta in chat.
