## =====================================================================
## FIRMA DI MALATTIA STRATIFICATA PER SESSO  [TISSUE - uc - genes]
## Derivato da: Analysis/TISSUE_sex_differences/meta_analysis_gene_expression/uc/perform_meta_analysis_sex_GUT.R
## malati vs sani, ESEGUITO SEPARATAMENTE nei soli MASCHI e nelle sole FEMMINE (stessi filtri GSE/tessuto/timepoint/etichette dello script originale)
## - dentro ogni GSE si confrontano solo persone dello stesso sesso (nessun confronto tra dataset diversi -> nessun batch effect incrociato)
## - un GSE entra in un sesso solo con >= 2 casi e >= 2 sani di quel sesso
## - effect size > 0 = piu' alto nei malati che nei sani (dentro quel sesso)
## - il confronto M vs F della firma (interazione sesso x malattia) si fa con compare_sex_specific_disease_signature.R
## Risultati: Analysis/sex_specific_disease_signature/<males|females>/TISSUE_genes/...   (cartella NUOVA, non tocca nulla di esistente)
## =====================================================================
SS_ROOT <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_specific_disease_signature"
for (sd in c("males", "females")) dir.create(file.path(SS_ROOT, sd, "TISSUE_genes", "uc"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(SS_ROOT, "males", "TISSUE_genes", "uc"))
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master/")
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/ulcerative_colitis_metadata.tsv", sep="\t", header=T)

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

directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/uc/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/uc/rnaseq_series_uc.txt"
microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/uc/microarray_series_uc.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-series[!series %in% c("GSE1710")]
series<-append(series, "GSE87466")#microarray series that cannot be downloaded from getGEOData
series<-unique(series[!grepl("GPL", series)])

#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]


comp<-list(c("ulcerative colitis", "ulcerative colitis (active)", "ulcerative colitis (extensive)", "ulcerative colitis (inactive)",
             "ulcerative colitis (left-sided)", "ulcerative colitis (limited)", "ulcerative colitis (quiescent)", "ulcerative colitis (remission)"
), c("healthy", "control, healthy", "non-IBD control, healthy"))

names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("GUT_uc_vs_healthy")


tissues<-c("intestinal mucosa", "intestinal mucosa (inflamed)", "intestinum",
           
           "ascending (right-sided) colon", "ascending colon", "colon", "colon (affected)", "colon (inflamed)", 
           "colon (whole mucosa)", "colon epithelium", "colon mucosa", "colonic mucosa", "colonic mucosa (inflamed)", 
           "descending (left-sided) colon", "descending colon", "descending colon (inflamed)", "descending colon (inflammed)",
           "left colon", "right colon", "right colon (inflamed)", "sigmoid colon", "sigmoid colon (inflamed)", "transverse colon",
           "tranverse colon (inflamed)",
           
           
           "ileum", "ileum slightly inflamed", "terminal ileum")

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
  message("######## FIRMA DI MALATTIA - SOLO ", ifelse(sx == "M", "MASCHI", "FEMMINE"), " (", sum(data_sx$sample.type %in% comp$cases), " casi, ", sum(data_sx$sample.type %in% comp$controls), " sani)")
  wd_now <- file.path(SS_ROOT, sex_dir, "TISSUE_genes", "uc"); dir.create(wd_now, recursive = TRUE, showWarnings = FALSE); setwd(wd_now)
  results <- run_pipeline(data_sx, confronti,sample.tissue="tissue", timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      specify_tissue_subtypes = TRUE,
                      list_tissue_subtypes = list(tissues, tissues)
  )
}

print("###############DONE! (sex-stratified disease signature, genes)")