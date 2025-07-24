setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/meta_analysis_pathways_activity/hallmarks/jia")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

pathways_list<-readRDS("../../../../../gene_sets/hallmarks_list_GENE_SYMBOL.RData")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/JIA_metadata.tsv", sep="\t", header=T)
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
directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/jia/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/jia/rnaseq_series_jia.txt"
microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/jia/microarray_series_jia.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])

#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]

comp<-list(c("juvenile idiopathic arthritis", "juvenile idiopathic arthritis (active)", "juvenile idiopathic arthritis (childhood-onset)", 
             "juvenile idiopathic arthritis (oligoarticular)", "juvenile idiopathic arthritis (pauciarticular)", "juvenile idiopathic arthritis (polyarticular, active)", 
             "juvenile idiopathic arthritis (polyarticular, remission)", "juvenile idiopathic arthritis (polyarticular)", "juvenile idiopathic arthritis (remission)", 
             "juvenile idiopathic arthritis (RF- polyarticular)", "juvenile idiopathic arthritis (systemic, active)", "juvenile idiopathic arthritis (systemic, inactive)",
             "juvenile idiopathic arthritis (systemic)", "new-onset juvenile idiopathic arthritis (active)", "new-onset juvenile idiopathic arthritis (systemic)", "
systemic-onset juvenile idiopathic arthritis"), c("healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("gsva_jia_vs_healthy")

data[data == "blood"] <- "whole blood"

results<-run_pipeline(data, confronti,sample.tissue=c("whole blood"), timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      pathways = pathways_list,
                      gs = 'HALLMARKS'
)

print("#########DONE!!!")