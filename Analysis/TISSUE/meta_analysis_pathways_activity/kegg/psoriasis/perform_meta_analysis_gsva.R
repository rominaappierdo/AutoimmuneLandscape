setwd("Analysis/TISSUE/meta_analysis_pathways_activity/kegg/psoriasis/")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)
library(RMySQL)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

pathways_list<-readRDS("../../../../../gene_sets/KEGG_list_GENE_SYMBOL.RData")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/psoriasis_metadata.tsv", sep="\t", header=T)

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

directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/psoriasis/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/psoriasis/rnaseq_series_psoriasis.txt"
microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/psoriasis/microarray_series_psoriasis.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])
series<-series[!series %in% c("GSE160932", "GSE16161")]

#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]
comp<-list(c("psoriasis", "psoriasis (dryness syndrome)", "psoriasis (erythrodermic)", "psoriasis (generalised pustular)", 
             "psoriasis (heat syndrom)", "psoriasis (inverse)", "psoriasis (mild)", "psoriasis (palmoplantar pustular)", "psoriasis (palmoplantar)", 
             "psoriasis (scalp)", "psoriasis (severe)", "psoriasis (stasis syndrome)"), c("healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("SKIN_gsva_pso_vs_healthy")


tissues<-c("abdomen skin (lesional)", "arm skin (lesional)", "back skin (lesional)", "buttock skin (lesional)",
           "center of plaque skin (lesional)", "dermis", "dermis (lesional)","elbow skin (lesional)", "epidermis", "epidermis (lesional)", 
           "leg skin (lesional)", "lower back skin (lesional)", "palm skin", "scalp skin", "scalp skin (lesional)",
           "skin", "skin ", "kin (chronic lesion)", "skin (epidermis+dermis)","skin (lesional)", "skin (limb)",
           "skin (perilesional)", "skin (psoriasis)", "skin (trunk)", "skin epidermis", "skin epidermis (lesional)",
           "sole skin", "thigh skin (lesional)", "tissue (lesional)","upper arm skin (lesional)" )


data<-data[data$sample_source %in% tissues]


results<-run_pipeline(data, confronti,sample.tissue="tissue", timepoint.id=c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),
                      genes.of.interest = NULL,
                      dest.dir = "../../../../../geo_data/",
                      from.cache = c(),
                      aggregate.subject.data = FALSE,
                      path_to_rnaseq_datasets = directory_downloaded_rnaseq,
                      rnaseq_series_file = rnaseq_series,
                      pathways = pathways_list,
                      gs = 'KEGG',
                      specify_tissue_subtypes = TRUE,
                      list_tissue_subtypes = list(tissues, tissues)
)
print("#########DONE!")