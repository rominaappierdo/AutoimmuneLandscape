setwd("Analysis/BLOOD/meta_analysis_cell_proportions/uc/")
source("../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../metapipeline-master_deconvolution/")
source("../../../../utils_functions/fix_metaintegrator_immunoStates/function.R")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

data<-read.csv("../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/ulcerative_colitis_metadata.tsv", sep="\t", header=T)
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

directory_downloaded_rnaseq<- "../../../../Data_Collection/downloaded_datasets/uc/downloaded_RNAseq/"
rnaseq_series <- "../../../../Data_Collection/gse_and_gsm_lists/uc/rnaseq_series_uc.txt"
microarray_series<-readLines("../../../../Data_Collection/gse_and_gsm_lists/uc/microarray_series_uc.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])

#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]

comp<-list(c("ulcerative colitis", "ulcerative colitis (active)", "ulcerative colitis (extensive)", "ulcerative colitis (inactive)",
             "ulcerative colitis (left-sided)", "ulcerative colitis (limited)", "ulcerative colitis (quiescent)", "ulcerative colitis (remission)"
), c("healthy", "control, healthy", "non-IBD control, healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("cell_prop_uc_vs_healthy")


data[data == "blood"] <- "whole blood"



results<-run_pipeline(data, confronti,sample.tissue=c("whole blood"), timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series
)

print("######DONE!")