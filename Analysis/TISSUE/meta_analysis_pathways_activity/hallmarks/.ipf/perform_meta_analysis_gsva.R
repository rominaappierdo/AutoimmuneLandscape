setwd("/mnt/data3/appierdo/PhD-project/Analysis/TISSUE/meta_analysis_pathways_activity/hallmarks/ipf/")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)
library(RMySQL)
source("../../../../../utils_functions/modify_metaintegrator_to_get_study_level_effectSizes_tissue.R")
assignInNamespace(x = '.combine.effect.sizes', .combine.effect.sizes, ns = "MetaIntegrator")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

pathways_list<-readRDS("../../../../../gene_sets/hallmarks_list_GENE_SYMBOL.RData")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/idiopathic_pulmonary_fibrosis_metadata.tsv", sep="\t", header=T)

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

#change names to use pakg to perform any comparison
names(data)<-c("GSE", "geo_accession", "subject.id","sample.type", "sample.tissue","sample_source","age", "gender", "cell_type", "timepoint.id","treatment")

directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/ipf/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/ipf/rnaseq_series_ipf.txt"
microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/ipf/microarray_series_ipf.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])
#"GSE52463" da sistemare e inserire
series<-series[!series %in% c("GSE2052")]

#AGGIUNGERE A series LE SERIES RNASEQ

series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]


comp<-list(c("idiopathic pulmonary fibrosis"), c("healthy", "control, healthy", "other, control, healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("LUNG_gsva_ipf_vs_healthy")


tissues<-c("central lung", "fibrotic peripheral area (lung)", "lung", "lung (affected)", "lung, alveolar septae",
           "lung, fibroblast foci", "lung, lower lobe", "lung, upper lobe", "peripheral lung")



data<-data[data$sample_source %in% tissues]


results<-run_pipeline(data, confronti,sample.tissue="tissue", timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      pathways = pathways_list,
                      gs = 'HALLMARKS',
                      specify_tissue_subtypes = TRUE,
                      list_tissue_subtypes = list(tissues, tissues)
)

print("##########DONE!")