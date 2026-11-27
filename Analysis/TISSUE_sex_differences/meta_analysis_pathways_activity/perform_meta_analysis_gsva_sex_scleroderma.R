## =====================================================================
## ssGSEA + META-ANALISI SUI PATHWAY: MASCHI vs FEMMINE NEI SOLI PAZIENTI
## Derivato da: Analysis/TISSUE_sex_differences/meta_analysis_gene_expression/scleroderma/perform_meta_analysis_sex.R
## - stessi pazienti, filtri su GSE/tessuto/timepoint e confronto M vs F dello script sull'espressione genica
## - per ogni GSE: ssGSEA (GSVA, method = 'ssgsea') + scale() sui campioni selezionati, poi meta-analisi dei punteggi
##   con metapipeline-master_gsva (la stessa pipeline dei pathway "vs healthy")
## - cicla su 3 gene set: HALLMARKS, KEGG, WIKIPATHWAYS
## - effect size > 0  =>  attivita' del pathway piu' alta nei MASCHI
## Risultati: Analysis/TISSUE_sex_differences/meta_analysis_pathways_activity/<gene set>/scleroderma/
## Punteggi per campione: Analysis/TISSUE_sex_differences/gene_set_variation_analysis/  (cartella NUOVA: non tocca quella originale)
## =====================================================================
# Percorsi assoluti (radice del progetto: /mnt/data3/appierdo/AutoimmuneLandscape)
gene_sets <- list(
  HALLMARKS    = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/hallmarks_list_GENE_SYMBOL.RData",
                      wd   = "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE_sex_differences/meta_analysis_pathways_activity/hallmarks/scleroderma"),
  KEGG         = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/KEGG_list_GENE_SYMBOL.RData",
                      wd   = "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE_sex_differences/meta_analysis_pathways_activity/kegg/scleroderma"),
  WIKIPATHWAYS = list(file = "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData",
                      wd   = "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE_sex_differences/meta_analysis_pathways_activity/wikipathways/scleroderma")
)
for (g in gene_sets) dir.create(g$wd, recursive = TRUE, showWarnings = FALSE)
dir.create("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE_sex_differences/gene_set_variation_analysis", recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE_sex_differences/meta_analysis_pathways_activity/hallmarks/scleroderma")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/scleroderma_metadata.tsv", sep="\t", header=T)

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

directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/scleroderma/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/scleroderma/rnaseq_series_ss.txt"
microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/scleroderma/microarray_series_ss.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])



#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]


comp<-list(c("systemic sclerosis", "systemic sclerosis (diffuse cutaneous)", "systemic sclerosis (diffuse)", 
             "systemic sclerosis (early)", "systemic sclerosis (limited cutaneous)", "systemic sclerosis (limited)", 
             "systemic sclerosis (localized)", "systemic sclerosis (non-cutaneous)", "systemic sclerosis sine scleroderma"), c("healthy", "control, healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("SKIN_scleroderma_vs_healthy")

#####start getgeodata modified for not rnaseq

get_geo_data <- function(sel.datasets, dest.dir = "../../../../../geo_data/", from.cache = c(), all.samples.metadata = NULL,
                         aggregate.subject.data = FALSE, path_to_rnaseq_datasets, rnaseq_series_file, pathways) {
  
  to.parse <- sel.datasets[!sel.datasets %in% from.cache]
  from.cache <- intersect(sel.datasets, from.cache)
  data <- list()
  if(length(to.parse) > 0)
    data <- MetaIntegrator::getGEOData(to.parse, destdir = dest.dir)
  if(length(from.cache) > 0)
    for(x in from.cache) {
      message(sprintf("Loading dataset %s from cache", x))
      df <- readRDS(file.path(dest.dir, sprintf("%s.rds", x)))
      data$originalData[[x]] <- df
    }
  if(aggregate.subject.data) {
    if(is.null(all.samples.metadata))
      stop("You must provide all.samples.metadata if aggregate.subject.data is TRUE")
    for(i in 1:length(data$originalData))
      data$originalData[[i]] <- aggregate_subject_data(data$originalData[[i]], all.samples.metadata)
  }
  matrices<-lapply(data$originalData, function(x) map_probes_to_gene_symbols_from_metaobj(x$expr, x$keys))
  matrices<-unname(matrices)
  
  counter<-1
  while (counter<=length(matrices)){
    data$originalData[[counter]]$expr <- matrices[[counter]] 
    counter<-counter+1
  }
  
  counter<-1
  while (counter<=length(matrices)){
    gem<-data$originalData[[counter]]
    
    #check for log scale
    if (.GEM_log_check(gem$expr) == FALSE) { #se non sono in log scale li trasforma in log scale
      min_value <- min(gem$expr, na.rm = T)
      
      # There is possibility of negative values in processed non-normalized 
      # data from the background subtraction algorithm. For e.g., GSE103842 (Illumina data)
      # if minimum value is negative, add small offset to all values to make them all positive
      if (min_value < 0) {
        gem$expr <- gem$expr - min_value + 1
      }
      
      gem$expr <- log2(gem$expr)
      gem$exp_comment = "Data was not in log scale originally"
      
      data$originalData[[counter]]$expr<-log2(gem$expr)
      data$originalData[[counter]]$exp_comment="Data was not in log scale originally"
    }  
    counter<-counter+1
  }
  
  log_check<-sapply(data$originalData, function(x) .GEM_log_check(x$expr))
  if (any(!log_check)==TRUE){
    stop("there are data not in log scale")
  } else{
    print("all data are in log scale")
  }
  
  #compute pas metanalysis
  pathway_activy_scores <- list()
  counter <- 0
  for (x in data$originalData){
    counter <- counter +1
    print(x$formattedName)
    out_gsva <- gsva(as.matrix(x$expr), gset.idx.list = pathways, min.sz=5, parallel.sz=8, method = 'ssgsea', ssgsea.norm = F)
    out_gsva <- scale(out_gsva)
    pathway_activy_scores[[counter]] <- out_gsva
  }
  for (i in seq(1, length(data$originalData))){
    data$originalData[[i]]$expr<-pathway_activy_scores[[i]]
    lista_pathway<- rownames(pathway_activy_scores[[i]])
    names(lista_pathway)<- rownames(pathway_activy_scores[[i]])
    data$originalData[[i]]$keys<-lista_pathway
  }
  return(data)
}
####end getgeodata modified for not rnaseq


tissues<-c("forearm skin","skin")

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
# stesso schema dei nomi degli script dei pathway esistenti (gsva_<malattia>_..., oppure <TESSUTO>_gsva_<malattia>_...)
names(confronti) <- sub("^((GUT|SKIN|SYNOVIUM)_)?", "\\1gsva_", names(confronti))

message("Pazienti per GSE e sesso (servono almeno 2 M e 2 F per ogni GSE):")
print(table(data$GSE, data$gender))
## =====================================================================


for (gs in names(gene_sets)) {
  setwd(gene_sets[[gs]]$wd)
  pathways_list <- readRDS(gene_sets[[gs]]$file)
  message("########## ssGSEA maschi vs femmine - gene set: ", gs)

  results <- run_pipeline(data, confronti,sample.tissue="tissue", timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
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

print("##########DONE!")