#importing libraries
library(tidyverse)
library(data.table)
library(MetaIntegrator)
library(dplyr)
library(rlang)
library(tidyr)
requireNamespace(package = "RMySQL", quietly=TRUE)
#####################################


counts_to_tpm <- function(counts, median.transcript.length) {
    TPM <- counts %>%
        data.frame() %>%
        rownames_to_column('Gene_symbol') %>% #adds the row_names column at the first index position of our data frame
        inner_join(median.transcript.length) %>%
        relocate(all_of(c('length')), .after='Gene_symbol') %>%
        column_to_rownames('Gene_symbol') %>%
        mutate(
            across(
                -c(length),
                ~(.x+1)/length #qua fa +1 per non far venire 0 credo
            )
        ) %>%
        select(-length)
    
    TPM <- t( t(TPM)*1e6 / colSums(TPM) )
    return(TPM)
}


replace_zero_with_min <- function(m) {
    x <- min(m[m > 0], na.rm = T)
    m[m == 0] <- x
    return(m)
}

fpkm_to_tpm <- function(fpkm) {                                                                                                                          
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))                                                                                                           
}

check_if_can_become_rownames<-function(dataframe, column_name){
   dataframe<-as.data.frame(dataframe)
   if (length(dataframe[, (column_name)])== length(unique(dataframe[, (column_name)])) ) {
       print("no duplicated gene symbols")
       rownames(dataframe)<-dataframe[,column_name]
       dataframe[,(column_name)]<-NULL
       
       
   } else {
       print("duplicated gene symbols")
       print("deleting duplicated gene symbols")
       dataframe<-as.data.frame(dataframe)
       
       
    
       dataframe<-dataframe %>% add_count(across({{column_name}}), name = 'occurrence')
       
       dataframe<-dataframe[dataframe[,"occurrence"]==1,]
        
       dataframe["occurrence"]<-NULL

       rownames(dataframe)<-dataframe$column_name
       dataframe$column_name<-NULL
       print("deleted duplicated gene symbols")
       
       rownames(dataframe)<-dataframe[,column_name]
       dataframe[,(column_name)]<-NULL
       
       
   }
    
return(dataframe)
}





remove_all_zeros_rows<- function(data){

      data<-data[rowSums(data[])>0,]
      #data<-data[apply(data[,-1], 1, function(x) !all(x==0)),]
      return(data)

}
                       
                       
order_by_col_names<-function(data){

    data<-data[ , order(names(data))]
    return(data)

}       
                       
sort_by_col<-function(d, colname){
    d<-as.data.frame(d)
    d <- d[order(d[,colname]),] #sort by column
    return(d)
}
                       
merge_samples_from_multiple_files<- function(path, patt=".txt", colname="symbol") {
    
    for (file in list.files(path, pattern=patt)[1]){
        gsm<-strsplit(file, "_")[[1]][1]
        data<-read.csv(file, sep="\t")
        setDT(data)
        colnames(data)<-c(colname, gsm)
    }

    for (file in list.files(path, pattern=patt)[2:length(list.files(path, pattern=patt))]){
        gsm<-strsplit(file, "_")[[1]][1]
        single<-read.csv(file, sep="\t")
        colnames(single)<-c(colname, gsm)
        setDT(single)
        data[[gsm]]<-single[[gsm]]
    
    }
    return(data)                                
                                }
                       
      
remove_version_from_ensid<-function(data){    #devono essere rownames
    
        ids_without_version<-unname(sapply(rownames(data), function(x) strsplit(x, "\\.")[[1]][1]))
        data<-(cbind(data,as.data.frame(ids_without_version)))
        return(data)
  }     
                                           
                                           
                                           
                                           
                                           


                      
get_metadata_for_a_series<-function(gse, path){

    metadata<-fread(path)
    #metadata<-read.csv(path, sep="\t")
    setDT(metadata)
    metadata<-metadata[metadata$series_id==gse]
    return(metadata)

}
                      
                      
                         
################################################################################################################################## 
                                 
.ucsc_refseq_table <- function(){
  return(ucsc_refseq_table)
}    


.ensembl_entrez_connector <- function(){
  
  #Establish connection to EnsEMBL
  con <- DBI::dbConnect(DBI::dbDriver("MySQL"),
                        username= 'anonymous',
                        host    = 'ensembldb.ensembl.org',
                        dbname  = 'homo_sapiens_core_76_38',
                        password= '',
                        port    = 3306)
  
  #build the parts for the query string
  query  <- "select xref.dbprimary_acc,
  xref.display_label 
  from external_db 
  join xref 
  where external_db.db_name = 'EntrezGene' AND 
  xref.external_db_id = external_db.external_db_id"
  
  #run query and get table
  ens_table           <- DBI::dbGetQuery(con,query)
  colnames(ens_table) <- c("EntrezID","GeneName")
  
  #Close connection
  DBI::dbDisconnect(con)
  
  #return table
  return(ens_table)
}


.ensembl_ensgID_table <- function(){
  #return table from data
  return(ens_ensgID_table)
}


                         
new_mapping_ids <- function(expr_data, colname, var="ens"){
    
  expr_data<-as.data.frame(expr_data)
  
  symbol_pos <- NULL;
  comment    <- NULL;
  keys       <- NULL;

    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case3-> genes are encoded in ENSG ID [use Human for this]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
      if(var=="ens"){
        
        #convert genes into gene names
        ens_table   <- .ensembl_ensgID_table();
        probe2table <- match(expr_data[, colname],ens_table$stable_id);
        
        #store data
        keys        <- ens_table$display_label[probe2table];
        #names(keys) <- as.character(GEO_fData[[1]]);
        
        comment     <- "Human ENSG format";
      }
    
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case4-> genes are encoded in REFSEQ format [use Human for this]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
     
      if(var=="refseq"){
        
        ucsc_table  <- .ucsc_refseq_table()
        
        #if REFSEQ IDs have _at at the end add it to the table
        if(length(grep("_at",expr_data[, colname]))>0){
          ucsc_table$name <- gsub("$","_at",ucsc_table$name)
        }
        
        probe2table <- match(gsub("\\..$","",expr_data[, colname]),
                             ucsc_table$name)
        
        #store data
        keys        <- ucsc_table$name2[probe2table]
        #names(keys) <- as.character(GEO_fData[[1]])
        
        comment     <- "Human REFSEQ format"
      }
    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case6-> genes are encoded in EntrezID format [use Human for this]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      if(var=="entrez"){
        #get table
        entrezid_table <- .ensembl_entrez_connector()
        
        probe2table    <- match(expr_data$symbol,entrezid_table$EntrezID)
        
        #store data
        keys        <- entrezid_table$GeneName[probe2table]
        #names(keys) <- as.character(GEO_fData[[1]])
        
        comment     <- "Human EntrezID format"
      }
   
  
  #remove probes that do not map back to genes
  
    #genes      <- .ensembl_ensgID_table()
    #fake_genes <- which(!(unlist(sapply(keys,function(i)strsplit(i,",")[[1]][1])) 
    #                      %in% 
    #                        genes$display_label))
    #keys[fake_genes] <- NA
  
 expr_data$mapped_symbol<-keys
  #return output
  return(expr_data)
}         
     
                                 
                                 
                                 
                                 
                                 
                                 
conte <- function(data, gse, plt='NULL'){

    mtl<-read.csv("/mnt//data3/appierdo/new_median_transcript_length.tsv",sep="\t", header=T)
    #devo per forza rimuovere i duplicati
    mtl<-check_if_can_become_rownames(mtl, "Gene_symbol")
    mtl$Gene_symbol<-rownames(mtl)
    rownames(mtl)<-NULL
    mtl<-mtl[,c("length", "Gene_symbol")]
    
    tpm_data<-(counts_to_tpm(data, mtl))

    expression_comment="Data in log scale"

    if (max(tpm_data,na.rm=T) > 50){
        tpm_data<- replace_zero_with_min(tpm_data)
        tpm_data<- log2(tpm_data)
        expression_comment = "Data was not in log scale originally"
    }


    plot(density(unlist((tpm_data))))
    
    class_vector<-rep(0, ncol(tpm_data))
    names(class_vector)<-colnames(tpm_data)

    keys<-rownames(tpm_data)
    names(keys)<-rownames(tpm_data)

    pheno <- data.frame(sample.id = colnames(tpm_data), gse = gse)                                                                                       
    row.names(pheno) <- colnames(tpm_data)
    
    res<-list(                                                                                                                                         
        name = gse,                                                                                                                               
        load_comment = "Loaded with GEOquery",
        class=class_vector,   
        keys = keys,                                                                                                               
        key_comment = "Gene Symbol",                                                                                                                                                
        pheno = pheno,
        platform = plt,
        expr = as.matrix(tpm_data),
        exp_comment = expression_comment,
        formattedName = gse                                                                                                                                         
    )  

    if (checkDataObject(res, "Dataset")==TRUE){
        print("ok")
        return(res)
        }

}

                                 
                                 
fpkm<-function(data, gse, plt= 'NULL'){
    expression_comment="Data in log scale"
    tpm_data<- apply(data, 2, fpkm_to_tpm) #per colonna

    if (max(tpm_data,na.rm=T) > 50){
        tpm_data<- replace_zero_with_min(tpm_data)
        tpm_data<- log2(tpm_data)
        expression_comment = "Data was not in log scale originally"
    }


    plot(density(unlist((tpm_data))))

    class_vector<-rep(0, ncol(tpm_data))
    names(class_vector)<-colnames(tpm_data)

    pheno <- data.frame(sample.id = colnames(tpm_data), gse = gse)                                                                                       
    row.names(pheno) <- colnames(tpm_data) 

    keys<-rownames(tpm_data)
    names(keys)<-rownames(tpm_data)


    res<-list(                                                                                                                 
        name = gse,    
        load_comment = "Loaded with GEOquery",
        class=class_vector,   
        keys = keys, 
        key_comment = "Gene symbol",                              
        pheno = pheno,
        platform = plt,
        expr = as.matrix(tpm_data),
        exp_comment = expression_comment,
        formattedName = gse                                                                                                                                         
    )  

    if (checkDataObject(res, "Dataset")==TRUE){
            print("ok")
            return(res)
            }
}
                                 
                                 
                          
                                 
create_metaobject_from_rnaseq<- function(path_series, downloaded_RNAseq_directory){
    series<-readLines(path_series)
    list_series<-c()
    for (i in series){
        path<-paste0(downloaded_RNAseq_directory,i,"/", i, ".rds",sep="")
        obj<-readRDS(path)
        list_series[[i]]<-obj
     }

    metaobj<-list()
    metaobj$originalData<-list_series
    if (checkDataObject(metaobj, "Meta", "Pre-Analysis")==TRUE){
        return(metaobj)
        
    }else{
        stop("errore")
    }
    
}
                                 
############################################################################################################################################################################################################################################################################################################################################################################################                                 
                                 
                                 
mapping_ids<-function(data, name_of_mapping_column="ensembl_id", name_of_data_column){
    data<-as.data.frame(data)
    #mapping ids
    mapping_genes<-read.csv("/mnt/data3/appierdo/gene_id_conversion_table.tsv",header=T, sep="\t")
    mapping_genes<-mapping_genes[,c(name_of_mapping_column, "Gene_name")]
    #remove duplicated rows
    mapping_genes<-mapping_genes[!duplicated(mapping_genes), ]

    #remove non mapped ids
    mapping_genes<-mapping_genes[!mapping_genes$Gene_name=="",]
    mapping_genes<-mapping_genes[!is.na(mapping_genes$Gene_name),]

    conversion_vector<-mapping_genes$Gene_name
    names(conversion_vector)<-mapping_genes[,name_of_mapping_column]

    
    v<-unname((lapply((data)[,name_of_data_column], function(x) unname(conversion_vector[x]))))
    v<-as.data.frame(do.call(rbind, v))
    data["Gene_symbol"]<-v
                      
    data<-check_if_can_become_rownames(data, "Gene_symbol")
    data[,name_of_data_column]<-NULL
    return(data)
    
}


#####merge multiple rows for same gene performing a mean

mean_expression_gene_rowise <- function(dataset, colname) {
  
  dataset<- dataset %>%
    group_by(get(colname)) %>%
    summarise_if(is.numeric, mean) %>%
    rename("Symbol" = "get(colname)") %>%
    as.data.frame()
  
  rownames(dataset) <- dataset$Symbol
  dataset$Symbol<-NULL
  
  return(dataset)
  
}



#####merge multiple rows for same gene performing a sum

sum_expression_gene_rowise <- function(dataset, colname) {
  
  dataset<- dataset %>%
    group_by(get(colname)) %>%
    summarise_if(is.numeric, sum) %>%
    dplyr::rename("Symbol" = "get(colname)") %>%
    as.data.frame()
  
  rownames(dataset) <- dataset$Symbol
  dataset$Symbol<-NULL
  
  return(dataset)
  
}

#####old

# expand_gene_symbols<- function(df){
#   
#   df<-df[!is.na(rownames(df)), ]
#   keys <- as.character(rownames(df))
#   skey <- strsplit( keys, split = ",")
#   df <- df[rep(1:nrow(df), sapply(skey, length)), ]
#   rownames(df) <- unlist(skey)
#   df<-do.call(rbind, by(df, row.names(df), FUN = colMeans))
#   return(df)
#   
# }

expand_gene_symbols<- function(df, colname){
  k <- as.character(df[[colname]])
  skey <- strsplit(k, split = ",")
  df <- df[rep(1:nrow(df), sapply(skey, length)), ]
  df[[colname]] <- unlist(skey)
  return(df)
}


rename_single_column <- function(df, old_col_name, new_col_name){
  
  colnames(df)[colnames(df) == old_col_name] <- new_col_name
  return(as.data.frame(df))
  
}


wide_to_long <- function(dataframe, name_row_column, name_column_column, name_value_column) {
  dataframe <- as.data.frame(dataframe)
  # Add row names as a separate column
  out_long <- dataframe %>%
    rownames_to_column(var = name_row_column)
  
  # Reshape the data to long format
  out_long <- out_long %>%
    pivot_longer(cols = -name_row_column, 
                 names_to = name_column_column, 
                 values_to = name_value_column)
  
  return(out_long)
}


remove_rows_with_less_than_5_reads_in_more_than_half_samples <- function(data, n) {
  if ((n %% 2) == 0) {
    half_samples<- (n/2) +1
  } else {
    half_samples<- (n/2) +0.5
  }
  
  data$counts_filter <- matrixStats::rowCounts(data > 5)
  data <- data[data$counts_filter >= half_samples, ]
  data$counts_filter<-NULL
  data<-order_by_col_names(data)
  return(data)
}

create_metaobject <- function(expr_data, studyID, cases_regex, metadata, metadata_sample_column_name, metadata_dis_state_column_name) {
  
  obj <- list()
  obj$originalData <- list()
  
  obj$originalData[[studyID]] <- list()  
  obj$originalData[[studyID]]$expr <- as.matrix(expr_data)
  obj$originalData[[studyID]]$keys <- rownames(expr_data)
  names(obj$originalData[[studyID]]$keys) <- rownames(expr_data)
  obj$originalData[[studyID]]$class <- ifelse(grepl(cases_regex, metadata[, metadata_dis_state_column_name]), 1, 0)
  names(obj$originalData[[studyID]]$class) <- metadata[, metadata_sample_column_name]
  obj$originalData[[studyID]]$formattedName <- studyID
  pheno_df <- metadata
  rownames(pheno_df) <- pheno_df[, metadata_sample_column_name]
  obj$originalData[[studyID]]$pheno <- pheno_df
  
  if (!checkDataObject(obj, "Meta", "Pre-Analysis")) {
    stop("invalid metaobject")
  }
  
  return(obj)
}


remove_rows_with_more_than_N_zeros <- function(DF, n) {
  DF <- DF[rowSums(DF == 0) <= n, ]
  return(DF)
}