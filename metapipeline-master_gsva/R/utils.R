#' @export
get_datasets_summary <- function(tab, sel.columns = c("sample.type", "sample.tissue", "timepoint.id")) {
    tab <- tab[, c("GSE", sel.columns)]

    #tab$sample.type <- toupper(tab$sample.type)

    m <- plyr::ddply(tab, ~GSE, function(x) {
        ret <- data.frame(lapply(x, function(v) {paste(sort(unique(v)), collapse = ";")}))
        ret$total.samples <- nrow(x)
        return(ret)
    })

    return(m)
}

#' @export
merge_geo_datasets <- function(a, b, out.name = NULL) {
    ret <- list()

    ret$name <- paste(a$name, b$name, sep = "_")
    ret$load_comment <- a$load_comment
    ret$class <- c(a$class, b$class)

    common.keys <- intersect(names(a$keys), names(b$keys))
    if(any(a$keys[common.keys] != b$keys[common.keys], na.rm = TRUE))
        stop("Keys mismatch, cannot merge dataset")
    ret$keys <- c(a$keys, b$keys[setdiff(names(b$keys), names(a$keys))])
    ret$key_comment <- a$key_comment
    ret$pheno <- plyr::rbind.fill(a$pheno, b$pheno)
    row.names(ret$pheno) <- c(row.names(a$pheno), row.names(b$pheno))
    ret$platform <- a$platform

    if(length(intersect(colnames(a$expr), colnames(b$expr))))
        stop("Duplicated samples, cannot merge dataset")

    expr.a <- data.frame(a$expr, probe.id = row.names(a$expr))
    expr.b <- data.frame(b$expr, probe.id = row.names(b$expr))

    expr.merge <- plyr::join(expr.a, expr.b, by = "probe.id")
    row.names(expr.merge) <- expr.merge$probe.id
    expr.merge$probe.id <- NULL
    ret$expr <- as.matrix(expr.merge)

    if(!(all(row.names(ret$expr) %in% names(ret$keys))))
        stop("Missing keys, cannot merge datasets")

    ret$keys <- ret$keys[row.names(ret$expr)]

    ret$exp_comment <- a$exp_comment
    ret$formattedName <- paste(a$formattedName, b$formattedName, sep = "_")

    if(!is.null(out.name))
        ret$name <- ret$formattedName <- out.name

    return(ret)
}



#' Get metadata for a GSE entry
#' @param dest.dir The location where the metadata is to be saved
#' @param GSE The GSE accession for which metadata should be retrieved
#' @param geo.cache.dir The directory where GEO data shuold be downloaded
#'
#' @return The metadata is written in a file called \code{GSE}_meta_all.csv in
#'   \code{dest.dir} and is also returned invisbly
#'
#'
#'@export
get_gse_metadata <- function(GSE, dest.dir = "./", geo.cache.dir = tempdir()) {
    dir.create(dest.dir, recursive = TRUE, showWarnings = FALSE)

    gse1 <- GEOquery::getGEO(GSE, GSEMatrix = FALSE, parseCharacteristics = TRUE, destdir = geo.cache.dir)
    gpls <- names(GEOquery::GPLList(gse1))
    l <- length(gpls)
    gse2 <-GEOquery::getGEO(GSE, GSEMatrix = TRUE, parseCharacteristics = TRUE, destdir = geo.cache.dir)
    df <- NULL
    for (i in 1:l){
        gpl <- (Biobase::pData(Biobase::phenoData(gse2[[i]])))
        df <- dplyr::bind_rows(df, gpl)
    }

    file <- file.path(dest.dir, sprintf("%s_meta_all.csv", GSE))
    write.csv(df, file)
    return(invisible(df))
}


#' @export
get_overall_design <- function(accession, geo.cache.dir = tempdir()) {
    exp.data <- NA
    x <- tryCatch(
        GEOquery::getGEO(accession, getGPL = FALSE, destdir = geo.cache.dir),
        error = function(cond) {NA}
    )

    if(!is.na(x) && class(x) != "GDS") {
        exp.data <- experimentData(x[[1]])
        exp.data <- exp.data@other$overall_design
    }
    if(is.null(exp.data))
        exp.data <- NA
    return(exp.data)
}





#' @export
bind_all_metadata <- function(dir) {
    files.list <- list.files(pattern = "^GSE.*", path = dir, full.names = T)
    all.data <- lapply(files.list, read.table, row.names = 1, sep = ",", header = T, colClasses = "character")
    names(all.data) <- sapply(strsplit(basename(files.list), "_"), "[[", 1)
    all.data <- dplyr::bind_rows(all.data, .id = "GSE")
    return(all.data)
}




