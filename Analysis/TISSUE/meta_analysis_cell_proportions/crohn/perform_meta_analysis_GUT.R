setwd("/mnt/data3/appierdo/PhD-project/Analysis/TISSUE/meta_analysis_cell_proportions/crohn/")
source("../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../metapipeline-master_deconvolution/")
source("../../../../utils_functions/fix_metaintegrator_immunoStates/function.R")
#assignInNamespace(x = 'iSdeconvolution', iSdeconvolution, ns = "MetaIntegrator")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

source("../../../../utils_functions/modify_metaintegrator_to_get_study_level_effectSizes_tissue.R")
assignInNamespace(x = '.combine.effect.sizes', .combine.effect.sizes, ns = "MetaIntegrator")


data<-read.csv("../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/crohn_metadata.tsv", sep="\t", header=T)

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


directory_downloaded_rnaseq<- "../../../../Data_Collection/downloaded_datasets/crohn/downloaded_RNAseq/"
rnaseq_series <- "../../../../Data_Collection/gse_and_gsm_lists/crohn/rnaseq_series_crohn.txt"
microarray_series<-readLines("../../../../Data_Collection/gse_and_gsm_lists/crohn/microarray_series_crohn.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])
series<-series[!series %in% c("GSE1710", "GSE83982", "GSE77387") ]


#Add to "series" the RNAseq series
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]


comp<-list(c("crohn's disease", "crohn's disease (active)", "crohn's disease (colon-only)", "crohn's disease (colonic)",
             "crohn's disease (fibrotic)", "crohn's disease (ileal)", "crohn's disease (ileo-colonic)", "crohn's disease (inactive)", 
             "crohn's disease (late stage)", "crohn's disease (newly diagnosed)", "crohn's disease (newly-diagnosed)", "crohn's disease (notfibrotic)",
             "crohn's disease (post-operative recurrent)", "crohn's disease (post-operative)" ,"crohn's disease (UACL)"), c("healthy", "control, healthy", "non-IBD control, healthy"))
names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("GUT_crohn_vs_healthy")


tissues<-c("ileal margin", "ileal mucosa", "ileum", "ileum ", "ileum (active)", "ileum (inflamed)", "ileum (involved)",
           "ileum lower", "ileum middle", "ileum upper", "neoterminal ileum", "neoterminal ileum (inflamed)", "terminal ileal mucosa", 
           "terminal ileum", "terminal ileum (inflamed)",
           
           "colon", "colon (affected)", "colon epithelium", "colonic mucosa", "colonic mucosa (inflamed)","left colon", "right colon",
           "right colon (inflamed)", "sigmoid colon", "sigmoid colon ", "sigmoid colon (inflamed)", "sigmoid colon (involved)",
           "transverse colon", "transverse colon (inflamed)", "transverse colon (involved)",
           
           "intestinal mucosa", "intestinal mucosa (inflamed margin)", "intestinal tissue", "intestinal tissue (involved)", 
           "intestine", "intestinum", "small intestine",
           
           "ascending (right-sided) colon", "ascending colon", "ascending colon ", "ascending colon (inflamed)", "ascending colon (involved)",
           
           "descending (left-sided) colon","descending colon", "descending colon (inflamed)", "descending colon (involved)",
           
           "rectum")


data<-data[data$sample_source %in% tissues]


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


print("#########DONE!")