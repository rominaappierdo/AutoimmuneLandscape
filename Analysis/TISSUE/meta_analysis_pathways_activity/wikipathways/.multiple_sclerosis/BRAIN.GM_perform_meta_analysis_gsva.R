setwd("/mnt/data3/appierdo/PhD-project/Analysis/TISSUE/meta_analysis_pathways_activity/kegg/multiple_sclerosis/")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)


options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

pathways_list<-readRDS("../../../../../gene_sets/KEGG_list_GENE_SYMBOL.RData")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/multiple_sclerosis_metadata.tsv", sep="\t", header=T)

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


directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/multiple_sclerosis/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/multiple_sclerosis/rnaseq_series_ms.txt"

microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/multiple_sclerosis/microarray_series_ms.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])


#AGGIUNGERE A series LE SERIES RNASEQ

series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]

comp<-list(c("multiple sclerosis", "multiple sclerosis (clinically isolated syndrome)", "multiple sclerosis (non-recovered)", 
             "multiple sclerosis (primary progressive)", "multiple sclerosis (progressive relapsing)", "multiple sclerosis (radiologically isolated syndrome)", 
             "multiple sclerosis (recovered)", "multiple sclerosis (relapsing-remitting active)", "multiple sclerosis (relapsing-remitting stable)", 
             "multiple sclerosis (relapsing-remitting)", "multiple sclerosis (secondary progressive)", "multuple sclerosis (definite)"), c("healthy", "control, healthy"))


names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("BRAIN.GM_gsva_ms_vs_healthy")


tissues<-c("brain", "brain (active plaque)", "brain (chronic plaque)", "brain (chronice active plaque)", "grey matter  (lesional), frontal cortex",
           "grey matter  (lesional), parietal cortex", "grey matter (lesional), motor cortex", "grey matter in brodmann area 22",
           "grey matter in brodmann area 23", "grey matter in brodmann area 24", "grey matter in brodmann area 25", "grey matter in brodmann area 26",
           "grey matter in brodmann area 27", "grey matter in brodmann area 28", "grey matter in brodmann area 29", "grey matter in brodmann area 30",
           "grey matter in brodmann area 31", "grey matter in brodmann area 32", "grey matter in brodmann area 33", "grey matter in brodmann area 34", "grey matter in brodmann area 35",
           "grey matter in brodmann area 36", "grey matter in brodmann area 37", "grey matter in brodmann area 38", "grey matter in brodmann area 39",
           "grey matter in brodmann area 40", "grey matter in brodmann area 41", "grey matter, frontal cortex", "grey matter, motor cortex", "grey matter, parietal cortex",
           "leukocortex, grey matter (lesional)", "subpial grey matter lesions from the frontal gyri")

data<-data[data$sample_source %in% tissues]


results<-run_pipeline(data, confronti,sample.tissue="tissue", timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      pathways = pathways_list,
                      gs = 'KEGG'
)

print("#################DONE!")
