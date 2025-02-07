setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/TFs/meta_analysis_tf_activity/multiple_sclerosis/")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva_tf/")
library(GSVA)

source("../../../../../utils_functions/modify_metaintegrator_to_get_study_level_effectSizes.R")
assignInNamespace(x = '.combine.effect.sizes', .combine.effect.sizes, ns = "MetaIntegrator")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

pathways_list<-readRDS("../../../../../gene_sets/collecTRI_regulons_dorothea_08122023.rds")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/multiple_sclerosis_metadata.tsv", sep="\t", header=T)
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


directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/multiple_sclerosis/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/multiple_sclerosis/rnaseq_series_ms.txt"

microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/multiple_sclerosis/microarray_series_ms.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])


#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]


comp<-list(c("multiple sclerosis", "multiple sclerosis (clinically isolated syndrome)", "multiple sclerosis (non-recovered)", 
             "multiple sclerosis (primary progressive)", "multiple sclerosis (progressive relapsing)", "multiple sclerosis (radiologically isolated syndrome)", 
             "multiple sclerosis (recovered)", "multiple sclerosis (relapsing-remitting active)", "multiple sclerosis (relapsing-remitting stable)", 
             "multiple sclerosis (relapsing-remitting)", "multiple sclerosis (secondary progressive)", "multuple sclerosis (definite)"), c("healthy", "control, healthy"))


names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("tf_gsva_ms_vs_healthy")


data[data == "blood"] <- "whole blood"

results<-run_pipeline(data, confronti,sample.tissue=c("whole blood"), timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      pathways = pathways_list
)

print("#################DONE!")
