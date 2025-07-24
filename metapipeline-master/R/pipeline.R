
aggregate_subject_data <- function(data, all.samples.metadata) {
    metadata.cols <- c("sample.type", "sample.tissue", "geo_accession", "subject.id", "timepoint.id")
    metadata <- all.samples.metadata[all.samples.metadata$GSE == data$name, metadata.cols]

    if(all(apply(metadata, 2, function(x) {all(!is.na(x))}))) {
        grouping.cols <- setdiff(metadata.cols, "geo_accession")
        if(anyDuplicated(metadata[, grouping.cols]) != 0) {
            message("Aggregating duplicate samples for dataset ", data$name)

            expr <- data$expr
            m <- as.data.frame(t(expr))

            m$geo_accession <- row.names(m)
            m <- merge(m, metadata[, metadata.cols], by = "geo_accession")

            m <- plyr::ddply(m, grouping.cols, function(x) {
                geo.accession <- x[1, "geo_accession"]
                x <- x[, !names(x) %in% metadata.cols]
                ret <- data.frame(value = colMeans(as.matrix(x)))
                ret$geo_accession <- geo.accession
                ret$gene <- row.names(ret)
                return(ret)
            })

            m <- reshape2::acast(m, geo_accession ~ gene, value.var = "value")
            data$keys <- data$keys[colnames(m)]
            data$pheno <- data$pheno[row.names(m),]
        }

    } else {
        message("Skipping subject aggreation because metadata is not available for all samples in dataset ", data$name)
    }
    return(data)
}



#path_to_rnaseq_datasets sarà diverso per ogni malattia

get_geo_data <- function(sel.datasets, dest.dir = "/mnt/data3/appierdo/geo_data/", from.cache = c(), all.samples.metadata = NULL,
                         aggregate.subject.data = FALSE, path_to_rnaseq_datasets, rnaseq_series_file) {

  to.parse <- sel.datasets[!sel.datasets %in% from.cache]
  from.cache <- intersect(sel.datasets, from.cache)

  #devo togliere le series che non sono microarray, perchè con getgeodata non le posso scaricare correttamente

  to.parse<-to.parse[!to.parse %in% readLines(rnaseq_series_file)]

  data <- list()

  if(length(to.parse) > 0)
    data <- MetaIntegrator::getGEOData(to.parse, destdir = dest.dir)

  if(length(from.cache) > 0)
    for(x in from.cache) {
      message(sprintf("Loading dataset %s from cache", x))
      df <- readRDS(file.path(dest.dir, sprintf("%s.rds", x)))
      data$originalData[[x]] <- df
    }


  ################rnaseq part
  rnaseq_data<-create_metaobject_from_rnaseq(rnaseq_series_file, path_to_rnaseq_datasets)
  data$originalData<-append(data$originalData, rnaseq_data$originalData)
  ###########

  if(aggregate.subject.data) {
    if(is.null(all.samples.metadata))
      stop("You must provide all.samples.metadata if aggregate.subject.data is TRUE")
    for(i in 1:length(data$originalData))
      data$originalData[[i]] <- aggregate_subject_data(data$originalData[[i]], all.samples.metadata)
  }

  return(data)
}
# ####end getgeodata modified for rnaseq


select_samples_character_column <- function(all.samples.metadata, cases, controls, sample.tissue = NULL, timepoint.id = NULL, specify_tissue_subtypes, list_tissue_subtypes) {
    ret <- list()

    for(gse in unique(all.samples.metadata$GSE)) {
        temp <- all.samples.metadata[all.samples.metadata$GSE == gse, ]
        if(!is.null(sample.tissue))
            temp <- temp[temp$sample.tissue %in% sample.tissue, ]
        if(!is.null(timepoint.id))
            temp <- temp[temp$timepoint.id %in% timepoint.id, ]
        
        if (specify_tissue_subtypes == TRUE) {
          tissue_subtypes_case <- list_tissue_subtypes[[1]] #first position cases
          tissue_subtypes_control <- list_tissue_subtypes[[2]] #second position control
          
          temp_cases <- temp[temp$sample.type %in% cases, ]
          temp_cases <- temp_cases[temp_cases$sample_source %in% tissue_subtypes_case,]
          
          temp_control <- temp[temp$sample.type %in% controls, ]
          temp_control <- temp_control[temp_control$sample_source %in% tissue_subtypes_control, ]
          
          # print(gse)
          # print((unique(temp_cases[, "sample_source"])))
          # print((unique(temp_control[, "sample_source"])))
          
          case.samples <- temp_cases[, "geo_accession"]
          control.samples <- temp_control[, "geo_accession"]
          
        } else {
          case.samples <- temp[temp$sample.type %in% cases, "geo_accession"]
          control.samples <- temp[temp$sample.type %in% controls, "geo_accession"]
        }
        
        if(length(case.samples) > 1 && length(control.samples) > 1)
            ret <- c(ret, setNames(list(list(cases = case.samples, controls = control.samples)), gse))
    }
    
    
    return(ret)
}


select_samples_numeric_column <- function(all.samples.metadata, col.name, mid.point = NULL, sample.tissue = NULL, timepoint.id = NULL) {
    all.samples.metadata <- all.samples.metadata[!is.na(all.samples.metadata[, col.name]), ]
    ret <- list()
    if(is.null(mid.point))
        mid.point <- median(all.samples.metadata[, col.name], na.rm = T)

    message(sprintf("%s midpoint: %d", col.name, mid.point))

    for(gse in unique(all.samples.metadata$GSE)) {
        temp <- all.samples.metadata[all.samples.metadata$GSE == gse, ]
        if(!is.null(sample.tissue))
            temp <- temp[temp$sample.tissue %in% sample.tissue, ]
        if(!is.null(timepoint.id))
            temp <- temp[temp$timepoint.id %in% timepoint.id, ]

        case.samples <- temp[temp[, col.name] >= mid.point, "geo_accession"]
        control.samples <- temp[temp[, col.name] < mid.point, "geo_accession"]
        if(length(case.samples) > 1 && length(control.samples) > 1)
            ret <- c(ret, setNames(list(list(cases = case.samples, controls = control.samples)), gse))
    }
    return(ret)
}



#' @export
select_samples <- function(all.samples.metadata, comparison, sample.tissue = NULL, timepoint.id = NULL, specify_tissue_subtypes, list_tissue_subtypes) {
    ret <- NULL
    if("col.name" %in% names(comparison)) {
        ret <- select_samples_numeric_column(all.samples.metadata, col.name = comparison$col.name, mid.point = comparison$mid.point,
                                                     sample.tissue = sample.tissue, timepoint.id = timepoint.id)
    } else {
        ret <- select_samples_character_column(all.samples.metadata, comparison$cases, comparison$controls, sample.tissue = sample.tissue, timepoint.id = timepoint.id, specify_tissue_subtypes, list_tissue_subtypes)
    }

    return(ret)
}



#gestire quando ci sono più di una gpl

#' @export
get_meta_analysis <- function(sel.samples, dest.dir = "./geo_data", from.cache = c(), all.samples.metadata = NULL,
                              aggregate.subject.data = FALSE, path_to_rnaseq_datasets, rnaseq_series_file) {
    data <- get_geo_data(names(sel.samples), dest.dir = dest.dir, from.cache = from.cache,
                         all.samples.metadata = all.samples.metadata, aggregate.subject.data = aggregate.subject.data,
                         path_to_rnaseq_datasets = path_to_rnaseq_datasets, rnaseq_series_file = rnaseq_series_file)

    print(names(data$originalData))
    data$originalData<-data$originalData[names(data$originalData) %in% names(sel.samples)] #modification
    print(names(data$originalData))

    for(gse in names(data$originalData)) {
        obj <- data$originalData[[gse]]
        # Some samples may have been removed during data processing because of problems
        cases <- intersect(colnames(obj$expr), sel.samples[[gse]]$cases)
        controls <- intersect(colnames(obj$expr), sel.samples[[gse]]$controls)

        obj <- MetaIntegrator::subsetOriginalData(obj, union(cases, controls))
        v <- setNames(rep(NA, nrow(obj$pheno)), row.names(obj$pheno))
        v[cases] <- 1
        v[controls] <- 0
        obj$class <- v
        data$originalData[[gse]] <- obj
    }

    meta.res <- runMetaAnalysis(data, runLeaveOneOutAnalysis = F, maxCores = 1)
    return(meta.res)
}



#' @export
run_pipeline <- function(
    all.samples.metadata,
    comparisons,
    sample.tissue,
    timepoint.id,
    genes.of.interest = NULL,
    dest.dir = "/mnt/data3/appierdo/PhD-project/geo_data",
    from.cache = c(),
    aggregate.subject.data = FALSE,
    path_to_rnaseq_datasets,
    rnaseq_series_file,
    specify_tissue_subtypes = FALSE,
    list_tissue_subtypes = NULL
) {
    comparison.table <- list()

    for(st in sample.tissue) {
        for(comparison in names(comparisons)) {
            sel.samples <- select_samples(all.samples.metadata, comparisons[[comparison]], sample.tissue = st, timepoint.id = timepoint.id, specify_tissue_subtypes, list_tissue_subtypes)

            if(length(sel.samples) > 0) {
                meta.res <- get_meta_analysis(sel.samples, dest.dir = dest.dir, from.cache = from.cache,
                                              all.samples.metadata = all.samples.metadata, aggregate.subject.data = aggregate.subject.data,
                                              path_to_rnaseq_datasets = path_to_rnaseq_datasets, rnaseq_series_file = rnaseq_series_file)

                sample.tissue.string <- paste(st, collapse = "-")


                pdf(sprintf("%s_%s_top_genes_plot.pdf", comparison, sample.tissue.string))
                forestplot_top_genes(meta.res)
                dev.off()

                df <- meta.res$metaAnalysis$pooledResults
                df$gene <- row.names(df)


                if(!is.null(genes.of.interest)) {
                    plot_genes_of_interest(meta.res, genes.of.interest, sprintf("%s_%s", comparison, sample.tissue.string))
                    write.table(df[df$gene %in% unlist(genes.of.interest), ], sprintf("%s_%s_genes_of_interest_results.tsv", comparison, sample.tissue.string), col.names = T, row.names = F, sep = "\t")
                }

                write.table(df, sprintf("%s_%s_complete_results.tsv", comparison, sample.tissue.string), col.names = T, row.names = F, sep = "\t")
                comparison.table <- c(comparison.table, list(
                    list(
                        comparison = comparison,
                        sample.tissue = sample.tissue.string,
                        datasets = paste(names(sel.samples), collapse = ","),
                        num.datasets = length(sel.samples),
                        num.samples = length(unlist(sel.samples))
                    )
                ))
                qs::qsave(meta.res, sprintf("%s_%s_metaintegrator_object.qs", comparison, sample.tissue.string))

            }
        }
    }

    df <- do.call(rbind, lapply(comparison.table, data.frame))
    return(invisible(df))

}

