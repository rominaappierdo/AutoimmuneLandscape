setwd("/mnt/data3/appierdo/AutoimmuneLandscape/")
setwd("Analysis/TISSUE/meta_analysis_gene_expression")

files <- list.files(".", recursive = TRUE, pattern = "comple")

# Read one file -> named numeric vector: names = gene, values = effectSize
read_effect_vec <- function(file) {
  df <- read.csv(file, sep = "\t")
  df <- df[df$numStudies >= 2, ]
  df <- df[df$effectSizeFDR < 0.05, ]
  setNames(df$effectSize, df$gene)
}

# Disease label = top-level directory (before first "/")
disease_from_path <- function(path) sub("/.*", "", path)

# Load vectors and labels
vecs <- lapply(files, read_effect_vec)
labels <- vapply(files, disease_from_path, character(1))

# Make labels unique if a disease appears multiple times
labels_unique <- make.unique(labels, sep = "_")

names(vecs) <- labels_unique
n <- length(vecs)

# Prepare correlation matrix
cor_mat <- matrix(NA_real_, n, n, dimnames = list(labels_unique, labels_unique))
diag(cor_mat) <- 1

# Fill upper triangle with pairwise-overlap Spearman correlations
for (i in 1:(n-1)) {
  vi <- vecs[[i]]
  for (j in (i+1):n) {
    vj <- vecs[[j]]
    g <- intersect(names(vi), names(vj))
    if (length(g) >= 3) {
      # Spearman on the overlap only
      r <- suppressWarnings(cor(vi[g], vj[g], method = "spearman"))
      cor_mat[i, j] <- r
      cor_mat[j, i] <- r
    }
  }
}

# Optional: nicer printing
print(round(cor_mat, 3))

#write.table(round(cor_mat, 3), "pairwise_spearman_results.tsv", col.names = NA, row.names = T, quote = F, sep = '\t')

# If you want the original disease names (possibly duplicated) shown instead:
# dimnames(cor_mat) <- list(labels, labels); print(round(cor_mat, 3))
