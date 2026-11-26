setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE_sex_differences/meta_analysis_gene_expression/uc/")
source("../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../metapipeline-master/")
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

data<-read.csv("../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/ulcerative_colitis_metadata.tsv", sep="\t", header=T)

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

directory_downloaded_rnaseq<- "../../../../Data_Collection/downloaded_datasets/uc/downloaded_RNAseq/"
rnaseq_series <- "../../../../Data_Collection/gse_and_gsm_lists/uc/rnaseq_series_uc.txt"
microarray_series<-readLines("../../../../Data_Collection/gse_and_gsm_lists/uc/microarray_series_uc.txt")

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
## CONFRONTO MASCHI vs FEMMINE NEI SOLI PAZIENTI  (derivato dallo script originale)
## - tiene solo i pazienti, con le STESSE etichette di malattia dello script originale
## - sample.type diventa il sesso: casi = "M", controlli = "F"
## - effect size > 0  =>  espressione piu' alta nei maschi
## =====================================================================
new_name <- sub("vs_healthy", "male_vs_female", names(confronti)[1])
if (identical(new_name, names(confronti)[1])) new_name <- paste0(new_name, "_male_vs_female")

patient_labels <- comp$cases
data <- data[data$sample.type %in% patient_labels, ]

data$gender <- toupper(trimws(as.character(data$gender)))
data$gender[data$gender == "MALE"]   <- "M"
data$gender[data$gender == "FEMALE"] <- "F"
data <- data[data$gender %in% c("M", "F"), ]      # scarta sesso mancante o non valido

data$sample.type <- data$gender
comp_sex <- list(cases = "M", controls = "F")
confronti <- list(comp_sex)
names(confronti) <- new_name

message("Pazienti per GSE e sesso (servono almeno 2 M e 2 F per ogni GSE):")
print(table(data$GSE, data$gender))
## =====================================================================


results<-run_pipeline(data, confronti,sample.tissue="tissue", timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      specify_tissue_subtypes = TRUE,
                      list_tissue_subtypes = list(tissues, tissues)
)

print("##########DONE!")