## =====================================================================
## FIRMA DI MALATTIA STRATIFICATA PER SESSO  [BLOOD - lupus - genes]
## Derivato da: Analysis/BLOOD_sex_differences/meta_analysis_gene_expression/lupus/perform_meta_analysis_sex.R
## malati vs sani, ESEGUITO SEPARATAMENTE nei soli MASCHI e nelle sole FEMMINE (stessi filtri GSE/tessuto/timepoint/etichette dello script originale)
## - dentro ogni GSE si confrontano solo persone dello stesso sesso (nessun confronto tra dataset diversi -> nessun batch effect incrociato)
## - un GSE entra in un sesso solo con >= 2 casi e >= 2 sani di quel sesso
## - effect size > 0 = piu' alto nei malati che nei sani (dentro quel sesso)
## - il confronto M vs F della firma (interazione sesso x malattia) si fa con compare_sex_specific_disease_signature.R
## Risultati: Analysis/sex_specific_disease_signature/<males|females>/BLOOD_genes/...   (cartella NUOVA, non tocca nulla di esistente)
## =====================================================================
SS_ROOT <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_specific_disease_signature"
for (sd in c("males", "females")) dir.create(file.path(SS_ROOT, sd, "BLOOD_genes", "lupus"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(SS_ROOT, "males", "BLOOD_genes", "lupus"))
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master/")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv", sep="\t", header=T)
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

directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/lupus/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/lupus/rnaseq_series_lupus.txt"

microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/lupus/microarray_series_lupus.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])
series<-append(series,"GSE11908")

#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]

comp<-list(c("systemic lupus erythematosus", "systemic lupus erythematosus (active)", "systemic lupus erythematosus (childhood-onset)",
             "systemic lupus erythematosus (less active)", "systemic lupus erythematosus (pediatric)", "systemic lupus erythematosus (quiescent)"
), c("healthy", "healthy (ANA-)", "healthy (ANA+)", "control/healthy"))

#comp<-list(c("cutaneous lupus", "cutaneous lupus (acute)", "cutaneous lupus (chronic)", "cutaneous lupus (subacute)",
#            "cutaneous lupus erythematosus (chronic)", "subacute cutaneous lupus" ), c("healthy"))

#comp<-list(c("lupus nephritis", "lupus nephritis (class 1)", "lupus nephritis (class 2)", "lupus nephritis (class 3), 
#lupus nephritis (class 3+5)", "lupus nephritis (class 4)", "lupus nephritis (class 4+5)", "lupus nephritis (class 5)",
#          "lupus nephritis (mixed 3/4+5)"), c("healthy"))

#comp<-list(c("discoid lupus", "discoid lupus erythematosus", "discoid lupus erythematosus (chronic)"), c("healthy"))

names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("lupus_vs_healthy")

data[data == "blood"] <- "whole blood"


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
  wd_now <- file.path(SS_ROOT, sex_dir, "BLOOD_genes", "lupus"); dir.create(wd_now, recursive = TRUE, showWarnings = FALSE); setwd(wd_now)
  results <- run_pipeline(data_sx, confronti,sample.tissue=c("whole blood"), timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series
  )
}

print("###############DONE! (sex-stratified disease signature, genes)")