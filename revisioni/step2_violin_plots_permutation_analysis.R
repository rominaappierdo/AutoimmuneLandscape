setwd("/mnt/data3/appierdo/AutoimmuneLandscape/revisioni")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
library(singscore)
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")
#source("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/module_validation_figs/psoriasis_GSE117239/Untitled.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/function_compute_all_features_for_module.R")

gse<-'GSE78193'


#INTERSECTION GENES AMONG IFN PATHWAYS

ifn_genes <- c("CXCL10", "EIF2AK2", "IFIT2", "IRF1", "IRF2", "IRF9", "ISG15", "PSMB9",  
               "STAT2", "TAP1" )


########## set module
module_name <- "mod5"
module <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod5/mod5.txt")
module_interpretation <- 'IFN-g module'

########## get expression data
data <- getGEOData(gse)
expr <- data$originalData$GSE78193_GPL6848$expr
map_probes <- (as.data.frame(data$originalData$GSE78193_GPL6848$keys))
colnames(map_probes) <- 'keys'
expr <- merge(expr, map_probes, by = 'row.names')
expr$Row.names <- NULL
expr <- na.omit(expr)
expr <- mean_expression_gene_rowise(expr, colname = 'keys')
expr <- order_by_col_names(expr)

########## get metadata
meta <- fread("/mnt/data3/appierdo/PhD-project/Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv")
meta <- meta[meta$series_id==gse]
meta <- meta[meta$cell_type=='cultured sample', ]
meta <- sort_by_col(meta, colname = 'geo_accession')

if (!identical(meta$geo_accession, colnames(expr))) {
  stop("different samples from expression and metadata")
}

rank_expr <- rankGenes(expr)


# ============================================================
# 2. CONVERT RANKS INTO PERCENTILES
# ============================================================

percentage_rank_df <- apply(
  rank_expr,
  MARGIN = 2,
  function(x) {
    x / max(x, na.rm = TRUE) * 100
  }
)

# Interpretation:
#
# percentile rank ~100 = among the highest expressed genes
# percentile rank ~50  = around the middle
# percentile rank ~0   = among the lowest expressed genes

# ============================================================
# 3. KEEP ONLY THE SHARED IFN GENES PRESENT IN THE DATASET
# ============================================================

ifn_genes_present <- intersect(
  ifn_genes,
  rownames(percentage_rank_df)
)

cat(
  "Shared IFN genes present:",
  length(ifn_genes_present),
  "out of",
  length(ifn_genes),
  "\n"
)

print(ifn_genes_present)


# ============================================================
# 4. OBSERVED MEAN PERCENTILE RANK
# ============================================================

sub_rank_expr <- percentage_rank_df[
  ifn_genes_present,
  ,
  drop = FALSE
]

observed <- colMeans(
  sub_rank_expr,
  na.rm = TRUE
)

# observed contains one value per sample.
#
# Example:
# Sample1 = 82
#
# means that, on average, the shared IFN genes lie around
# the 82nd percentile of expression in that sample.



# ============================================================
# 5. PERMUTATION TEST
# ============================================================

set.seed(1234)

n_perm <- 10000
n_genes <- length(ifn_genes_present)

# All measured genes can potentially be selected
background_genes <- rownames(percentage_rank_df)

# Matrix:
# rows    = permutations
# columns = samples

null_distribution <- matrix(
  NA_real_,
  nrow = n_perm,
  ncol = ncol(percentage_rank_df)
)

colnames(null_distribution) <- colnames(percentage_rank_df)


# ------------------------------------------------------------
# Randomization
# ------------------------------------------------------------

for(i in seq_len(n_perm)) {
  
  random_genes <- sample(
    background_genes,
    size = n_genes,
    replace = FALSE
  )
  
  null_distribution[i, ] <- colMeans(
    percentage_rank_df[
      random_genes,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
}


# ============================================================
# 6. EMPIRICAL P-VALUE
# ============================================================

# Alternative hypothesis:
#
# The shared IFN genes are MORE HIGHLY RANKED
# than random genes.

p_values <- sapply(
  seq_along(observed),
  function(j) {
    
    (
      sum(
        null_distribution[, j] >= observed[j]
      ) + 1
    ) /
      (n_perm + 1)
  }
)

names(p_values) <- names(observed)


# Multiple-testing correction across samples

padj <- p.adjust(
  p_values,
  method = "BH"
)


# ============================================================
# 7. EFFECT SIZE: Z-SCORE
# ============================================================

null_mean <- colMeans(
  null_distribution,
  na.rm = TRUE
)

null_sd <- apply(
  null_distribution,
  2,
  sd,
  na.rm = TRUE
)

z_score <- (
  observed - null_mean
) / null_sd


# ============================================================
# 8. RESULT TABLE
# ============================================================

results <- data.frame(
  sample = names(observed),
  observed_rank = as.numeric(observed),
  expected_random_rank = as.numeric(null_mean),
  difference = as.numeric(
    observed - null_mean
  ),
  z_score = as.numeric(z_score),
  pvalue = as.numeric(p_values),
  padj = as.numeric(padj)
)

results <- results %>%
  arrange(pvalue)

print(results)


# ============================================================
# 9. PREPARE DATA FOR PLOT
# ============================================================

null_plot_df <- data.frame(
  sample = rep(
    colnames(null_distribution),
    each = nrow(null_distribution)
  ),
  mean_percentile_rank = as.vector(
    null_distributions
  )
)

observed_plot_df <- data.frame(
  sample = names(observed),
  mean_percentile_rank = as.numeric(observed)
)


# ============================================================
# 10. VIOLIN PLOT
# ============================================================

#metadata
submeta <- meta[, c('geo_accession', 'timepoint', 'treatment')]
dim(null_plot_df)
colnames(submeta) <- c('sample', 'timepoint', 'treatment')


null_plot_df <- merge(null_plot_df, submeta, by = 'sample')
dim(null_plot_df)

observed_plot_df <- merge(observed_plot_df, submeta, by = 'sample')

null_plot_df$timepoint <- factor(
  null_plot_df$timepoint,
  levels = c("t0", "24h", "48h")
)

observed_plot_df$timepoint <- factor(
  observed_plot_df$timepoint,
  levels = c("t0", "24h", "48h")
)
ggplot(
  null_plot_df,
  aes(
    x = sample,
    y = mean_percentile_rank,
    fill = treatment
  )
) +
  geom_violin(
    trim = FALSE
  ) +
  geom_point(
    data = observed_plot_df,
    aes(
      x = sample,
      y = mean_percentile_rank
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean percentile rank",
    title = "Rank enrichment of shared IFN genes",
    subtitle = paste0(
      "Points = observed ",
      n_genes,
      "-gene set; violins = random gene sets"
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  facet_wrap(
    ~ timepoint,
    scales = "free_x",
    nrow = 1
  ) 


# ============================================================
# 11. Z-SCORE PLOT
# ============================================================

results_plot <- results %>%
  mutate(
    sample = factor(
      sample,
      levels = sample[order(z_score)]
    )
  )

ggplot(
  results_plot,
  aes(
    x = sample,
    y = z_score
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Permutation z-score",
    title = "Enrichment of shared IFN genes in expression ranks"
  )


# ============================================================
# GENE-LEVEL DATA
# ============================================================

gene_level_df <- as.data.frame(
  t(
    percentage_rank_df[
      ifn_genes_present,
      ,
      drop = FALSE
    ]
  )
)

gene_level_df$sample <- rownames(gene_level_df)

gene_level_long <- gene_level_df %>%
  pivot_longer(
    cols = all_of(ifn_genes_present),
    names_to = "gene",
    values_to = "percentile_rank"
  )


# ============================================================
# PLOT
# Distribution across samples for each of the 10 genes
# ============================================================

ggplot(
  gene_level_long,
  aes(
    x = gene,
    y = percentile_rank
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.6,
    size = 2
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed"
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Expression percentile rank",
    title = "Expression ranks of shared IFN genes"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

gene_order <- gene_level_long %>%
  group_by(gene) %>%
  summarise(
    mean_rank = mean(percentile_rank, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_rank))

gene_level_long$gene <- factor(
  gene_level_long$gene,
  levels = gene_order$gene
)

ggplot(
  gene_level_long,
  aes(
    x = gene,
    y = percentile_rank
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.5,
    size = 2
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed"
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Expression percentile rank",
    title = "Shared IFN genes ranked by expression percentile"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )



#USANDO COME BACKGROUND GENI DEI GENE SETS

# Load gene sets
hallmarks <- readRDS(
  "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/hallmarks_list_GENE_SYMBOL.RData"
)

wiki <- readRDS(
  "/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData"
)

# Extract gene sets
alpha <- unique(
  hallmarks[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]
)

gamma_hallmark <- unique(
  hallmarks[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]
)

gamma_wiki <- unique(
  wiki[["WP_TYPE_II_INTERFERON_SIGNALING_IFNG"]]
)


background_ifn <- setdiff(
  unique(c(alpha, gamma_hallmark, gamma_wiki)),
  ifn_genes
)

library(singscore)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# 1. Rank genes
# ============================================================

rank_expr <- rankGenes(expr)

# Convert ranks to percentile
percentage_rank_df <- apply(
  rank_expr,
  2,
  function(x) x / max(x, na.rm = TRUE) * 100
)


# ============================================================
# 2. Shared 10 genes
# ============================================================

ifn_genes <- Reduce(
  intersect,
  list(alpha, gamma_hallmark, gamma_wiki)
)

ifn_genes_present <- intersect(
  ifn_genes,
  rownames(percentage_rank_df)
)

length(ifn_genes_present)
print(ifn_genes_present)


# ============================================================
# 3. IFN background genes
#    all genes from the 3 gene sets,
#    excluding the 10 shared genes
# ============================================================

background_ifn <- unique(
  c(alpha, gamma_hallmark, gamma_wiki)
)

background_ifn <- setdiff(
  background_ifn,
  ifn_genes_present
)

# Keep only genes actually measured in expr
background_ifn <- intersect(
  background_ifn,
  rownames(percentage_rank_df)
)

cat("Shared genes:", length(ifn_genes_present), "\n")
cat("IFN background genes:", length(background_ifn), "\n")


# ============================================================
# 4. Observed mean percentile rank of the 10 shared genes
# ============================================================

observed <- colMeans(
  percentage_rank_df[
    ifn_genes_present,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)


# ============================================================
# 5. Permutation test
#    Randomly select the same number of genes
#    from the IFN gene-set background
# ============================================================

set.seed(1234)

n_perm <- 10000
n_genes <- length(ifn_genes_present)

null_distribution <- matrix(
  NA_real_,
  nrow = n_perm,
  ncol = ncol(percentage_rank_df)
)

colnames(null_distribution) <- colnames(percentage_rank_df)


for(i in seq_len(n_perm)) {
  
  random_genes <- sample(
    background_ifn,
    size = n_genes,
    replace = FALSE
  )
  
  null_distribution[i, ] <- colMeans(
    percentage_rank_df[
      random_genes,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
}


# ============================================================
# 6. Empirical p-values
#
# H1:
# shared 10 genes have higher percentile ranks
# than other IFN genes
# ============================================================

p_values <- sapply(
  seq_along(observed),
  function(j) {
    
    (
      sum(
        null_distribution[, j] >= observed[j]
      ) + 1
    ) /
      (n_perm + 1)
    
  }
)

names(p_values) <- names(observed)

padj <- p.adjust(
  p_values,
  method = "BH"
)


# ============================================================
# 7. Effect size
# ============================================================

null_mean <- colMeans(
  null_distribution,
  na.rm = TRUE
)

null_sd <- apply(
  null_distribution,
  2,
  sd,
  na.rm = TRUE
)

z_score <- (
  observed - null_mean
) / null_sd


# ============================================================
# 8. Results table
# ============================================================

results <- data.frame(
  sample = names(observed),
  observed_rank = observed,
  expected_IFN_background = null_mean,
  difference = observed - null_mean,
  z_score = z_score,
  pvalue = p_values,
  padj = padj
)

results <- results %>%
  arrange(pvalue)

print(results)


# ============================================================
# 9. Plot:
# observed 10 shared genes vs random IFN genes
# ============================================================

null_plot_df <- data.frame(
  sample = rep(
    colnames(null_distribution),
    each = nrow(null_distribution)
  ),
  mean_percentile_rank = as.vector(
    null_distribution
  )
)

observed_plot_df <- data.frame(
  sample = names(observed),
  mean_percentile_rank = observed
)

#metadata
submeta <- meta[, c('geo_accession', 'timepoint', 'treatment')]
dim(null_plot_df)
colnames(submeta) <- c('sample', 'timepoint', 'treatment')


null_plot_df <- merge(null_plot_df, submeta, by = 'sample')
dim(null_plot_df)

observed_plot_df <- merge(observed_plot_df, submeta, by = 'sample')

null_plot_df$timepoint <- factor(
  null_plot_df$timepoint,
  levels = c("t0", "24h", "48h")
)

observed_plot_df$timepoint <- factor(
  observed_plot_df$timepoint,
  levels = c("t0", "24h", "48h")
)
ggplot(
  null_plot_df,
  aes(
    x = sample,
    y = mean_percentile_rank,
    fill = treatment
  )
) +
  geom_violin(
    trim = FALSE
  ) +
  geom_point(
    data = observed_plot_df,
    aes(
      x = sample,
      y = mean_percentile_rank
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_classic() +
  labs(
     x = NULL,
     y = "Mean percentile rank",
     title = "Rank enrichment of shared IFN genes",
     subtitle = "Points = 10 shared genes; violins = random genes from IFN gene sets"
    ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  facet_wrap(
    ~ timepoint,
    scales = "free_x",
    nrow = 1
  ) 


# ============================================================
# TABLE WITH PERMUTATION P-VALUES
# ============================================================

results_table <- data.frame(
  Sample = names(observed),
  Observed_mean_rank = as.numeric(observed),
  Expected_mean_rank = as.numeric(null_mean),
  Difference = as.numeric(observed - null_mean),
  Z_score = as.numeric(z_score),
  P_value = as.numeric(p_values),
  FDR = as.numeric(padj)
)

# Round values for visualization
results_table_print <- results_table %>%
  mutate(
    Observed_mean_rank = round(Observed_mean_rank, 2),
    Expected_mean_rank = round(Expected_mean_rank, 2),
    Difference = round(Difference, 2),
    Z_score = round(Z_score, 2),
    P_value = signif(P_value, 3),
    FDR = signif(FDR, 3)
  )

print(results_table_print)


####

library(singscore)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# ============================================================
# 1. RANK ALL GENES WITHIN EACH SAMPLE
# ============================================================

rank_expr <- rankGenes(expr)

percentage_rank_df <- apply(
  rank_expr,
  2,
  function(x) x / max(x, na.rm = TRUE) * 100
)


# ============================================================
# 2. GENES SHARED BY ALL THREE GENE SETS
# ============================================================

shared_genes <- Reduce(
  intersect,
  list(
    alpha,
    gamma_hallmark,
    gamma_wiki
  )
)

shared_genes
length(shared_genes)


# ============================================================
# 3. FUNCTION TO RANK GENES WITHIN EACH GENE SET
# ============================================================

rank_geneset <- function(
    geneset,
    geneset_name,
    percentage_rank_df,
    shared_genes
) {
  
  genes_present <- intersect(
    geneset,
    rownames(percentage_rank_df)
  )
  
  x <- percentage_rank_df[
    genes_present,
    ,
    drop = FALSE
  ]
  
  result <- data.frame(
    gene = rownames(x),
    
    mean_percentile_rank = rowMeans(
      x,
      na.rm = TRUE
    ),
    
    median_percentile_rank = apply(
      x,
      1,
      median,
      na.rm = TRUE
    ),
    
    sd_percentile_rank = apply(
      x,
      1,
      sd,
      na.rm = TRUE
    )
  )
  
  result$shared <- ifelse(
    result$gene %in% shared_genes,
    "Shared",
    "Not shared"
  )
  
  result$geneset <- geneset_name
  
  result <- result %>%
    arrange(desc(mean_percentile_rank)) %>%
    mutate(
      rank_within_geneset = row_number(),
      relative_position = rank_within_geneset / n() * 100
    )
  
  return(result)
}


# ============================================================
# 4. APPLY TO YOUR THREE GENE SETS
# ============================================================

alpha_rank <- rank_geneset(
  geneset = alpha,
  geneset_name = "IFN-alpha Hallmark",
  percentage_rank_df = percentage_rank_df,
  shared_genes = shared_genes
)

gamma_hallmark_rank <- rank_geneset(
  geneset = gamma_hallmark,
  geneset_name = "IFN-gamma Hallmark",
  percentage_rank_df = percentage_rank_df,
  shared_genes = shared_genes
)

gamma_wiki_rank <- rank_geneset(
  geneset = gamma_wiki,
  geneset_name = "IFN-gamma WikiPathways",
  percentage_rank_df = percentage_rank_df,
  shared_genes = shared_genes
)


# ============================================================
# 5. COMBINE RESULTS
# ============================================================

all_ranked <- bind_rows(
  alpha_rank,
  gamma_hallmark_rank,
  gamma_wiki_rank
)


# ============================================================
# 6. TABLE WITH ONLY THE 10 SHARED GENES
# ============================================================

shared_positions <- all_ranked %>%
  filter(shared == "Shared") %>%
  select(
    geneset,
    gene,
    mean_percentile_rank,
    median_percentile_rank,
    sd_percentile_rank,
    rank_within_geneset,
    relative_position
  ) %>%
  arrange(
    geneset,
    rank_within_geneset
  )

print(shared_positions)


# ============================================================
# 7. PLOT ALL GENES
# ============================================================

ggplot(
  all_ranked,
  aes(
    x = rank_within_geneset,
    y = mean_percentile_rank
  )
) +
  geom_point(
    aes(color = shared),
    size = 2
  ) +
  geom_text_repel(
    data = all_ranked %>%
      filter(shared == "Shared"),
    aes(label = gene),
    color = "red",
    size = 3,
    max.overlaps = Inf
  ) +
  facet_wrap(
    ~ geneset,
    scales = "free_x"
  ) +
  scale_color_manual(
    values = c(
      "Not shared" = "grey70",
      "Shared" = "red"
    ),
    labels = c(
      "Not shared" = "Other genes",
      "Shared" = "Shared by all three"
    )
  ) +
  theme_classic() +
  labs(
    x = "Gene rank within gene set",
    y = "Mean expression percentile rank",
    color = NULL,
    title = "Rank-based prioritization of IFN pathway genes"
  )
