setwd("/mnt/data3/appierdo/PhD-project/Analysis/TISSUE/meta_analysis_pathways_activity/kegg/lupus_subtypes_kidney/")
source("../../../../../utils_functions/functions_to_process_RNAseq_Data.R")
devtools::load_all("../../../../../metapipeline-master_gsva/")
library(GSVA)
library(RMySQL)
# source("../../../../../utils_functions/modify_metaintegrator_to_get_study_level_effectSizes_tissue.R")
# assignInNamespace(x = '.combine.effect.sizes', .combine.effect.sizes, ns = "MetaIntegrator")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

pathways_list<-readRDS("../../../../../gene_sets/KEGG_list_GENE_SYMBOL.RData")

data<-read.csv("../../../../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv", sep="\t", header=T)

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

directory_downloaded_rnaseq<- "../../../../../Data_Collection/downloaded_datasets/lupus/downloaded_RNAseq/"
rnaseq_series <- "../../../../../Data_Collection/gse_and_gsm_lists/lupus/rnaseq_series_lupus.txt"

microarray_series<-readLines("../../../../../Data_Collection/gse_and_gsm_lists/lupus/microarray_series_lupus.txt")

series <- unname(unlist(sapply(microarray_series, function(x) strsplit(x, "_"))))
series<-unique(series[!grepl("GPL", series)])
series<-append(series,"GSE11908")

#AGGIUNGERE A series LE SERIES RNASEQ
series<-append(series, readLines(rnaseq_series))
data<-data[data$GSE %in% series,]

comp<-list(c("systemic lupus erythematosus", "systemic lupus erythematosus (active)", "systemic lupus erythematosus (childhood-onset)",
             "systemic lupus erythematosus (less active)", "systemic lupus erythematosus (pediatric)", "systemic lupus erythematosus (quiescent)",
             "lupus nephritis", "lupus nephritis (class 1)", "lupus nephritis (class 2)", "lupus nephritis (class 3),
lupus nephritis (class 3+5)", "lupus nephritis (class 4)", "lupus nephritis (class 4+5)", "lupus nephritis (class 5)",
             "lupus nephritis (mixed 3/4+5)"), c("healthy", "healthy (ANA-)", "healthy (ANA+)", "control/healthy"))


names(comp)<-c("cases", "controls")
confronti<- list(comp)
names(confronti)<-c("KIDNEY.GLOMERULUS_gsva_lupus_vs_healthy")


#######start getmetanalisi modified
get_meta_analysis <- function(sel.samples, dest.dir = "../../../../../geo_data", from.cache = c(), comparisons, all.samples.metadata = NULL,
                              aggregate.subject.data = FALSE, path_to_rnaseq_datasets, rnaseq_series_file, pathways, gs, specify_tissue_subtypes) {
  data <- get_geo_data(names(sel.samples), dest.dir = dest.dir, from.cache = from.cache,
                       all.samples.metadata = all.samples.metadata, aggregate.subject.data = aggregate.subject.data,
                       path_to_rnaseq_datasets, rnaseq_series_file, pathways)
  
  print(names(data$originalData))
  data$originalData$GSE104948_GPL22945<-NULL
  names(data$originalData)<-c("GSE104948","GSE113342","GSE32591","GSE101437","GSE110999",
                              "GSE112087","GSE139350","GSE149050","GSE153781","GSE162828","GSE167923","GSE169080",
                              "GSE175759","GSE175839","GSE80183","GSE97264")
  
  data$originalData<-data$originalData[names(data$originalData) %in% names(sel.samples)] #modification
  print(names(data$originalData))
  
  for(gse in names(data$originalData)) {
    obj <- data$originalData[[gse]]
    cases <- intersect(colnames(obj$expr), unlist(sel.samples[[gse]]$cases))
    controls <- intersect(colnames(obj$expr), unlist(sel.samples[[gse]]$controls))
    obj <- MetaIntegrator::subsetOriginalData(obj, union(cases, controls))
    v <- setNames(rep(NA, nrow(obj$pheno)), row.names(obj$pheno))
    v[cases] <- 1
    v[controls] <- 0
    obj$class <- v
    data$originalData[[gse]] <- obj
  }
  
  for (pas in data$originalData) {
          pas_prop_matrix<-pas$expr
          gse<-pas$formattedName
          disease<-strsplit(strsplit(directory_downloaded_rnaseq, "downloaded_datasets/")[[1]][[2]], "/")[[1]][1]
          
          if (specify_tissue_subtypes == TRUE) {
            disease <- paste0(disease, '.', strsplit(names(comparisons), '_')[[1]][1], sep='')
          }
          
          file_name<-paste0("../../../gene_set_variation_analysis/", gs, "_", disease, "_",gse, "_pathway_activity_scores_matrix.tsv", sep='')
          write.table(pas_prop_matrix, file=file_name, quote = F, col.names = NA, sep='\t')
  }
  
  meta.res <- runMetaAnalysis(data, runLeaveOneOutAnalysis = F, maxCores = 1)
  return(meta.res)
}
#######end getmetanalisi


tissues<-c("kidney","kidney (glomerulus)")

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
                      list_tissue_subtypes = list(c("kidney (glomerulus)"), c("kidney (glomerulus)"))
                      
)

print("###############DONE!")