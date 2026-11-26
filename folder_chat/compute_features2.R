
compute_features2 <- function(expr, module, deconvolution = TRUE, ncore = 8) {
  
  # Input: genes x GSM
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  
  if (is.null(rownames(expr)) || is.null(colnames(expr)))
    stop("expr must have gene symbols as rownames and GSM as colnames")
  
  if (anyDuplicated(rownames(expr)))
    stop("Duplicated gene symbols in expr")
  
  if (anyDuplicated(colnames(expr)))
    stop("Duplicated GSM identifiers in expr")
  
  if (ncol(expr) < 2)
    stop("At least two GSM are required for standardization")
  
  module <- unique(module)
  
  # Standardize each feature across samples
  zscore <- function(x) {
    x <- as.matrix(x)
    z <- t(scale(t(x)))
    z[!is.finite(z)] <- NA_real_
    z
  }
  
  # Run ssGSEA and standardize each gene-set score
  run_ssgsea <- function(sets, suffix = "", invert = FALSE) {
    
    if (length(sets) == 0)
      return(NULL)
    
    # Keep gene sets with at least 5 genes in the expression matrix
    sets <- lapply(sets, function(x) {
      intersect(unique(x), rownames(expr))
    })
    
    sets <- sets[lengths(sets) >= 5]
    
    if (length(sets) == 0)
      return(NULL)
    
    scores <- gsva(
      expr,
      gset.idx.list = sets,
      min.sz = 5,
      parallel.sz = ncore,
      method = "ssgsea",
      ssgsea.norm = FALSE,
      verbose = FALSE
    )
    
    scores <- zscore(scores)
    
    if (invert)
      scores <- -scores
    
    rownames(scores) <- paste0(rownames(scores), suffix)
    
    scores
  }
  
  features <- list()
  
  # 1. CYTOKINES
  cyto_genes <- rownames(expr)[
    paste0(rownames(expr), "_cyto") %in% module
  ]
  
  if (length(cyto_genes) > 0) {
    cyto <- zscore(expr[cyto_genes, , drop = FALSE])
    rownames(cyto) <- paste0(rownames(cyto), "_cyto")
    features[["cyto"]] <- cyto
  }
  
  # 2. CELL PROPORTIONS
  if (deconvolution) {
    
    cell_prop <- as.data.frame(compute_cell_proportions(expr))
    
    # Use the rn column as GSM identifiers
    if ("rn" %in% colnames(cell_prop)) {
      rownames(cell_prop) <- as.character(cell_prop$rn)
      cell_prop$rn <- NULL
    }
    
    # Check GSM correspondence
    if (!all(colnames(expr) %in% rownames(cell_prop))) {
      stop("Cell proportion output does not match GSM identifiers")
    }
    
    # Keep selected cell types and preserve GSM order
    cell_prop <- cell_prop[
      colnames(expr),
      intersect(colnames(cell_prop), module),
      drop = FALSE
    ]
    
    if (ncol(cell_prop) > 0) {
      features[["cells"]] <- t(as.matrix(cell_prop))
    }
  }
  
  # 3. HALLMARK
  hall <- readRDS(
    "/mnt/data3/appierdo/PhD-project/gene_sets/hallmarks_list_GENE_SYMBOL.RData"
  )
  
  features[["hallmark"]] <- run_ssgsea(
    hall[names(hall) %in% module]
  )
  
  # 4. KEGG
  kegg <- readRDS(
    "/mnt/data3/appierdo/PhD-project/gene_sets/KEGG_list_GENE_SYMBOL.RData"
  )
  
  features[["kegg"]] <- run_ssgsea(
    kegg[names(kegg) %in% module]
  )
  
  # 5. WIKIPATHWAYS
  wp <- readRDS(
    "/mnt/data3/appierdo/PhD-project/gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData"
  )
  
  features[["wp"]] <- run_ssgsea(
    wp[names(wp) %in% module]
  )
  
  # 6. TRANSCRIPTION FACTORS
  tf <- readRDS(
    "/mnt/data3/appierdo/PhD-project/gene_sets/collecTRI_regulons_dorothea_08122023.rds"
  )
  
  features[["tf"]] <- run_ssgsea(
    tf[paste0(names(tf), "_tf") %in% module],
    suffix = "_tf"
  )
  
  # 7. miRNA
  mirna <- readRDS(
    "/mnt/data3/appierdo/PhD-project/gene_sets/miRTarBase_regulons.rds"
  )
  
  features[["mirna"]] <- run_ssgsea(
    mirna[names(mirna) %in% module],
    invert = TRUE
  )
  
  # Combine feature matrices
  features <- Filter(Negate(is.null), features)
  
  if (length(features) == 0) {
    return(matrix(
      numeric(0),
      nrow = ncol(expr),
      ncol = 0,
      dimnames = list(colnames(expr), character(0))
    ))
  }
  
  result <- do.call(rbind, features)
  
  # Remove features that cannot be calculated
  result <- result[
    rowSums(is.finite(result)) > 0,
    ,
    drop = FALSE
  ]
  
  # Return GSM x features
  result <- t(result)
  
  if (anyDuplicated(colnames(result)))
    stop("Duplicated feature names in output")
  
  return(result)
}
