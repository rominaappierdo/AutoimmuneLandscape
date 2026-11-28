## =====================================================================
## FIRMA DI MALATTIA STRATIFICATA PER SESSO  [TISSUE - psoriasis - pathways]
## Derivato da: Analysis/TISSUE_sex_differences/meta_analysis_pathways_activity/perform_meta_analysis_gsva_sex_psoriasis.R
## malati vs sani, ESEGUITO SEPARATAMENTE nei soli MASCHI e nelle sole FEMMINE (stessi filtri GSE/tessuto/timepoint/etichette dello script originale)
## - dentro ogni GSE si confrontano solo persone dello stesso sesso (nessun confronto tra dataset diversi -> nessun batch effect incrociato)
## - un GSE entra in un sesso solo con >= 2 casi e >= 2 sani di quel sesso
## - effect size > 0 = piu' alto nei malati che nei sani (dentro quel sesso)
## - il confronto M vs F della firma (interazione sesso x malattia) si fa con compare_sex_specific_disease_signature.R
## Risultati: Analysis/sex_specific_disease_signature/<males|females>/TISSUE_pathways/...   (cartella NUOVA, non tocca nulla di esistente)
## =====================================================================
# Percorsi assoluti (radice del progetto: /mnt/data3/appierdo/AutoimmuneLandscape)
SS_ROOT <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_specific_disease_signature"
gene_sets <- list(
  HALLMARKS    = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/hallmarks_list_GENE_SYMBOL.RData",
                      wd   = paste0(SS_ROOT, "/%s/TISSUE_pathways/hallmarks/psoriasis")),
  KEGG         = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/KEGG_list_GENE_SYMBOL.RData",
                      wd   = paste0(SS_ROOT, "/%s/TISSUE_pathways/kegg/psoriasis")),
  WIKIPATHWAYS = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData",
                      wd   = paste0(SS_ROOT, "/%s/TISSUE_pathways/wikipathways/psoriasis"))
)
for (sd in c("males", "females")) {
  for (g in gene_sets) dir.create(sprintf(g$wd, sd), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(SS_ROOT, sd, "gene_set_variation_analysis"), recursive = TRUE, showWarnings = FALSE)
}
setwd(sprintf(gene_sets$HALLMARKS$wd, "males"))
source("../../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../../metapipeline-master_gsva/")
library(GSVA)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

data<-read.csv("../../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/psoriasis_metadata.tsv", sep="\t", header=T)

#remove checked, single cell and cultured and sorted cells
setDT(data) #turn dataset into data.table object
data <- data[, check:treatment] #subset dataset until treatment
data[data==""] <- NA
data <- data[!grepl("x",data$check),]
data[,check:=NULL] #drop check column
#drop cultured cells e sscrnaseq
data <- data[!grepl("cultured",data$cell_type),]
data <- data[!grepl("single",data$sample_source),]
data <- data[is.na(data$cell_type),]
data<- data[data$sample_type=="tissue"]

names(data)<-c("GSE", "geo_accession", "subject.id","sample.type", "sample.tissue","sample_source","age", "gender", "cell_type", "timepoint.id","treatment")

directory_downloaded_rnaseq<- "../../../../../../Data_Collection/downloaded_datasets/psoriasis/downloaded_RNAseq/"
rnaseq_series <- "../../../../../../Data_Collection/gse_and_gsm_lists/psoriasis/rnaseq_series_psoriasis.txt"
microarray_series<-readLines("../../../../../../Data_Collection/gse_and_gsm_lists/psoriasis/microarray_series_psoriasis.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-series[!series %in% c("GSE160932", "GSE16161")]
series<-unique(series[!grepl("GPL", series)])



#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]


comp<-list(c("psoriasis", "psoriasis (dryness syndrome)", "psoriasis (erythrodermic)", "psoriasis (generalised pustular)", 
             "psoriasis (heat syndrom)", "psoriasis (inverse)", "psoriasis (mild)", "psoriasis (palmoplantar pustular)", "psoriasis (palmoplantar)", 
             "psoriasis (scalp)", "psoriasis (severe)", "psoriasis (stasis syndrome)"), c("healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("SKIN_pso_vs_healthy")

tissues<-c("abdomen skin (lesional)", "arm skin (lesional)", "back skin (lesional)", "buttock skin (lesional)",
           "center of plaque skin (lesional)", "dermis", "dermis (lesional)","elbow skin (lesional)", "epidermis", "epidermis (lesional)", 
           "leg skin (lesional)", "lower back skin (lesional)", "palm skin", "scalp skin", "scalp skin (lesional)",
           "skin", "skin ", "kin (chronic lesion)", "skin (epidermis+dermis)","skin (lesional)", "skin (limb)",
           "skin (perilesional)", "skin (psoriasis)", "skin (trunk)", "skin epidermis", "skin epidermis (lesional)",
           "sole skin", "thigh skin (lesional)", "tissue (lesional)","upper arm skin (lesional)" )


data<-data[data$sample_source %in% tissues]



## =====================================================================
## PREPARAZIONE: malati + sani, con sesso valido; poi si esegue la pipeline separatamente per M e F
## =====================================================================
data$gender <- toupper(trimws(as.character(data$gender)))
data$gender[data$gender == "MALE"]   <- "M"
data$gender[data$gender == "FEMALE"] <- "F"
data <- data[data$gender %in% c("M", "F"), ]      # scarta sesso mancante o non valido
data <- data[data$sample.type %in% c(comp$cases, comp$controls), ]   # malati e sani, con le etichette dello script originale
base_name <- names(confronti)[1]                   # nome originale del confronto (<malattia>_vs_healthy_)
message("Campioni per GSE, gruppo e sesso (servono >= 2 casi e >= 2 sani per sesso in ogni GSE):")
print(table(data$GSE, ifelse(data$sample.type %in% comp$cases, "case", "control"), data$gender))
## =====================================================================


for (sx in c("M", "F")) {
  sex_dir <- ifelse(sx == "M", "males", "females")
  data_sx <- data[data$gender == sx, ]
  confronti <- list(comp[c("cases", "controls")])
  names(confronti) <- paste0(base_name, ifelse(sx == "M", "ONLYMALES", "ONLYFEMALES"))
  names(confronti) <- sub("^((GUT|SKIN|SYNOVIUM)_)?", "\\1gsva_", names(confronti))   # stesso schema dei nomi dei pathway esistenti
  message("######## FIRMA DI MALATTIA - SOLO ", ifelse(sx == "M", "MASCHI", "FEMMINE"), " (", sum(data_sx$sample.type %in% comp$cases), " casi, ", sum(data_sx$sample.type %in% comp$controls), " sani)")
for (gs in names(gene_sets)) {
  wd_now <- sprintf(gene_sets[[gs]]$wd, sex_dir); dir.create(wd_now, recursive = TRUE, showWarnings = FALSE)
  setwd(wd_now)
  pathways_list <- readRDS(gene_sets[[gs]]$file)
  message("########## ssGSEA firma di malattia stratificata per sesso - gene set: ", gs)

  results <- run_pipeline(data_sx, confronti,sample.tissue="tissue", timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      specify_tissue_subtypes = TRUE,
                      list_tissue_subtypes = list(tissues, tissues),
                          pathways = pathways_list,
                          gs = gs
  )
}

}   # fine ciclo sesso
print("###############DONE! (sex-stratified disease signature, pathways)")