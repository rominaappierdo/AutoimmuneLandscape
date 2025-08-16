setwd("Analysis/BLOOD/meta_analysis_gene_expression")

ra <- read.csv("ra/ra_vs_healthy__whole blood_complete_results.tsv", sep = '\t')
ra <- ra[ra$effectSizeFDR < 0.05, ]
ra <- ra[ra$numStudies >= 2, ]

sle <- read.csv("lupus/lupus_vs_healthy_whole blood_complete_results.tsv", sep = '\t')
sle <- sle[sle$effectSizeFDR < 0.05, ]
sle <- sle[sle$numStudies >= 2, ]

ssc <- read.csv("scleroderma/scleroderma_vs_healthy_whole blood_complete_results.tsv", sep = '\t')
ssc <- ssc[ssc$effectSizeFDR < 0.05, ]
ssc <- ssc[ssc$numStudies >= 2, ]

ss <- read.csv("sjogren/sjogren_vs_healthy__whole blood_complete_results.tsv", sep = '\t')
ss <- ss[ss$effectSizeFDR < 0.05, ]
ss <- ss[ss$numStudies >= 2, ]


merged<- merge(ra, sle , by = 'gene')
cor(merged$effectSize.x, merged$effectSize.y)
cor.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')



merged<- merge(ra, ssc , by = 'gene')
cor(merged$effectSize.x, merged$effectSize.y)
cor.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')


merged<- merge(ra, ss , by = 'gene')
cor(merged$effectSize.x, merged$effectSize.y)
cor.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')

merged<- merge(ss, sle , by = 'gene')
cor(merged$effectSize.x, merged$effectSize.y)
cor.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')



merged<- merge(ssc, sle , by = 'gene')
cor(merged$effectSize.x, merged$effectSize.y)
cor.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')


merged<- merge(ss, ssc , by = 'gene')
cor(merged$effectSize.x, merged$effectSize.y)
cor.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')







# Load and filter function
load_and_filter <- function(path) {
  df <- read.csv(path, sep = '\t')
  df <- df[df$effectSizeFDR < 0.05 & df$numStudies >= 2, ]
  return(df)
}

# Load datasets
ra  <- load_and_filter("ra/ra_vs_healthy__whole blood_complete_results.tsv")
sle <- load_and_filter("lupus/lupus_vs_healthy_whole blood_complete_results.tsv")
ssc <- load_and_filter("scleroderma/scleroderma_vs_healthy_whole blood_complete_results.tsv")
ss  <- load_and_filter("sjogren/sjogren_vs_healthy__whole blood_complete_results.tsv")

datasets <- list(RA = ra, SLE = sle, SSC = ssc, SS = ss)

# Function to compute correlation between two diseases
correlation_value <- function(df1, df2) {
  merged <- merge(df1, df2, by = "gene")
  if (nrow(merged) == 0) return(NA)  # no overlap
  cor(merged$effectSize.x, merged$effectSize.y, method = "spearman")
}

# Build correlation matrix
diseases <- names(datasets)
corr_mat <- matrix(NA, nrow = length(diseases), ncol = length(diseases),
                   dimnames = list(diseases, diseases))

for (i in diseases) {
  for (j in diseases) {
    if (i == j) {
      corr_mat[i, j] <- 1
    } else {
      corr_mat[i, j] <- correlation_value(datasets[[i]], datasets[[j]])
    }
  }
}

print(corr_mat)

# Nice formatting
library(knitr)
kable(corr_mat, digits = 3, caption = "Spearman correlations of effect sizes")

# Round to 2 digits
corr_mat_round <- round(corr_mat, 2)

write.table(corr_mat_round,
            file = "RA_SLE_SS_SSc_correlations.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA) 
