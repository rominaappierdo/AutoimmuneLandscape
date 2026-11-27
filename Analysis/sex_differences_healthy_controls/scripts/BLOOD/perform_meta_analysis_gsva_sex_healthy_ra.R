## =====================================================================
## ssGSEA + META-ANALISI SUI PATHWAY: MASCHI vs FEMMINE NEI SOLI SANI (CONTROLLI)   [BLOOD - ra]
## Derivato da: Analysis/BLOOD_sex_differences/meta_analysis_pathways_activity/perform_meta_analysis_gsva_sex_ra.R  (stesso script dei pazienti, ma al posto dei pazienti tiene i CONTROLLI SANI)
## - stessi filtri su GSE / tessuto / timepoint / sorgente, stessa pipeline metapipeline-master_gsva (ssGSEA + scale per GSE)
## - sample.type = sesso: casi = "M", controlli = "F"  =>  effect size > 0 = attivita' piu' alta nei MASCHI SANI
## - serve per distinguere le differenze di sesso "normali" da quelle specifiche della malattia (vedi compare_patients_vs_healthy_sex.R)
## Risultati in: Analysis/sex_differences_healthy_controls/BLOOD/<gene set>/ra/   (cartella NUOVA, non tocca nulla di esistente)
## =====================================================================
# Percorsi assoluti (radice del progetto: /mnt/data3/appierdo/AutoimmuneLandscape)
gene_sets <- list(
  HALLMARKS    = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/hallmarks_list_GENE_SYMBOL.RData",
                      wd   = "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_healthy_controls/BLOOD/hallmarks/ra"),
  KEGG         = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/KEGG_list_GENE_SYMBOL.RData",
                      wd   = "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_healthy_controls/BLOOD/kegg/ra"),
  WIKIPATHWAYS = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData",
                      wd   = "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_healthy_controls/BLOOD/wikipathways/ra")
)
for (g in gene_sets) dir.create(g$wd, recursive = TRUE, showWarnings = FALSE)
dir.create("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_healthy_controls/gene_set_variation_analysis", recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_differences_healthy_controls/BLOOD/hallmarks/ra")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)


options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/rheumatoid_arthritis_metadata.tsv", sep="\t", header=T)
setDT(data) #turn dataset into data.table object
data <- data[, check:treatment] #subset dataset until treatment column
data[data==""] <- NA
data <- data[!grepl("x",data$check),]
data[,check:=NULL] #drop check column
#drop cultured cells e sscrnaseq
data <- data[!grepl("cultured",data$cell_type),]
data <- data[!grepl("single",data$sample_source),]
data <- data[is.na(data$cell_type),]

names(data)<-c("GSE", "geo_accession", "subject.id","sample.type", "sample.tissue","sample_source","age", "gender", "cell_type", "timepoint.id","treatment")
directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/ra/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/ra/rnaseq_series_ra.txt"
microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/ra/microarray_series_ra.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])

#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]


comp<-list(c("rheumatoid arthritis", "rheumatoid arthritis", "rheumatoid arthritis (DMARD-IR)", "rheumatoid arthritis (TNF-IR)",
             "rheumatoid arthritis (early)", "rheumatoid arthritis (end-stage)", "rheumatoid arthritis (high)", "rheumatoid arthritis (late stage)",
             "rheumatoid arthritis (long-standing)", "rheumatoid arthritis (low)", "rheumatoid arthritis (moderate)", 
             "rheumatoid arthritis (newly diagnosed)", "rheumatoid arthritis (remission)", "rheumatoid arthritis (severe)"), c("healthy", "control, healthy"))


names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("ra_vs_healthy_")

data[data == "blood"] <- "whole blood"




## =====================================================================
## CONFRONTO MASCHI vs FEMMINE NEI SOLI PAZIENTI  (derivato dallo script originale)
## - tiene solo i pazienti, con le STESSE etichette di malattia dello script originale
## - sample.type diventa il sesso: casi = "M", controlli = "F"
## - effect size > 0  =>  espressione piu' alta nei maschi
## =====================================================================
new_name <- sub("vs_healthy", "healthy_male_vs_female", names(confronti)[1])
if (identical(new_name, names(confronti)[1])) new_name <- paste0(new_name, "_healthy_male_vs_female")

patient_labels <- comp$controls   # SANI al posto dei pazienti
data <- data[data$sample.type %in% patient_labels, ]   # solo controlli sani

data$gender <- toupper(trimws(as.character(data$gender)))
data$gender[data$gender == "MALE"]   <- "M"
data$gender[data$gender == "FEMALE"] <- "F"
data <- data[data$gender %in% c("M", "F"), ]      # scarta sesso mancante o non valido

data$sample.type <- data$gender
comp_sex <- list(cases = "M", controls = "F")
confronti <- list(comp_sex)
names(confronti) <- new_name
# stesso schema dei nomi degli script dei pathway esistenti (gsva_<malattia>_..., oppure <TESSUTO>_gsva_<malattia>_...)
names(confronti) <- sub("^((GUT|SKIN|SYNOVIUM)_)?", "\\1gsva_", names(confronti))

message("Sani per GSE e sesso (servono almeno 2 M e 2 F per ogni GSE):")
print(table(data$GSE, data$gender))
## =====================================================================


for (gs in names(gene_sets)) {
  dir.create(gene_sets[[gs]]$wd, recursive = TRUE, showWarnings = FALSE)
  setwd(gene_sets[[gs]]$wd)
  pathways_list <- readRDS(gene_sets[[gs]]$file)
  message("########## ssGSEA maschi vs femmine SANI - gene set: ", gs)

  results <- run_pipeline(data, confronti,sample.tissue=c("whole blood"), timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                          pathways = pathways_list,
                          gs = gs
  )
}

print("###############DONE! (ssGSEA sex differences, HEALTHY)")