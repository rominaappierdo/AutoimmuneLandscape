setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/meta_analysis_pathways_activity/wikipathways/sjogren")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

pathways_list<-readRDS("../../../../../gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/sjogren_syndrome_metadata.tsv", sep="\t", header=T)
setDT(data) #turn dataset into data.table object
data <- data[, check:treatment] #subset dataset until treatment
data[data==""] <- NA
data <- data[!grepl("x",data$check),]
data[,check:=NULL] #drop check column
#drop cultured cells e sscrnaseq
data <- data[!grepl("cultured",data$cell_type),]
data <- data[!grepl("single",data$sample_source),]
data <- data[is.na(data$cell_type),]

names(data)<-c("GSE", "geo_accession", "subject.id","sample.type", "sample.tissue","sample_source","age", "gender", "cell_type", "timepoint.id","treatment")

directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/sjogren/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/sjogren/rnaseq_series_sjogren.txt"
microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/sjogren/microarray_series_sjogren.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])

#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]

comp<-list(c("sjogren's syndrome", "sjogren's syndrome (advanced)", "sjogren's syndrome (early)", "sjogren's syndrome (keratoconjunctivitis sicca)", 
             "sjogren's syndrome (moderate)", "sjogren's syndrome (primary)" ), c("control, sjogren's syndrome, healthy", "control, healthy", "non-sjogren's sicca", "healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("gsva_sjogren_vs_healthy")

data[data == "blood"] <- "whole blood"


results<-run_pipeline(data, confronti,sample.tissue=c("whole blood"), timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      pathways = pathways_list,
                      gs = 'WIKIPATHWAYS'
)

print("######DONE!")