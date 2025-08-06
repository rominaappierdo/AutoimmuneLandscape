setwd("Supplementary/supp_figures")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
source("../../utils_functions/function_compute_all_features_for_module.R")
source("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")
set.seed(1234)

gse<-'GSE78193'

########## set module
module_name <- "mod5"
module <- readLines("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod5/mod5.txt")
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
meta <- fread("../..//Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv")
meta <- meta[meta$series_id==gse]
meta <- meta[meta$cell_type=='cultured sample', ]
meta <- sort_by_col(meta, colname = 'geo_accession')
meta <- meta[meta$timepoint == '48h',]

expr <- expr[, colnames(expr) %in% meta$geo_accession]

if (!identical(meta$geo_accession, colnames(expr))) {
  stop("different samples from expression and metadata")
}

#real one

module <- mod5
########## compute features
total_df <- compute_features(expr = expr, module = module, deconvolution = F)
total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable

score_df <- as.data.frame(colSums(total_df))
colnames(score_df) <- "module"
score_df$Condition <- paste0(meta$timepoint, "_", meta$treatment)

stim_df <- score_df[score_df$Condition == "48h_IFN-gamma", ]
unstim_df <- score_df[score_df$Condition == "48h_unstimulated", ]

n1 <- nrow(stim_df)
n2 <- nrow(unstim_df)

m1 <- mean(stim_df$module)
m2 <- mean(unstim_df$module)
diff <- m1 - m2

#compute pooled standard deviation
sd1  <- stats::sd(stim_df$module)
sd2  <- stats::sd(unstim_df$module)
sp   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) )

#correct for small sample size [from Cohen's D to Hedge's G]
cf   <- 1 - 3/( 4*(n1 + n2) - 9 ) # small sample size correction
g_real    <- cf * diff/sp
se.g_real <- sqrt( (n1+n2)/(n1*n2) + 0.5*g_real^2 /(n1+n2-3.94) )

#random one
#module has 49 features


dist_random_es <- c()
while (length(dist_random_es) < 1000){
  print(paste0("number: ", length(dist_random_es)+1, sep = ''))
  module <- c(mod1, mod2, mod3, mod4, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
  module <- sample(module, size = 49)
  ########## compute features 
  total_df <- compute_features(expr = expr, module = module, deconvolution = F)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  score_df$Condition <- paste0(meta$timepoint, "_", meta$treatment)
  
  stim_df <- score_df[score_df$Condition == "48h_IFN-gamma", ]
  
  unstim_df <- score_df[score_df$Condition == "48h_unstimulated", ]
  
  n1 <- nrow(stim_df)
  n2 <- nrow(unstim_df)
  
  m1 <- mean(stim_df$module)
  m2 <- mean(unstim_df$module)
  diff <- m1 - m2
  
  #compute pooled standard deviation
  sd1  <- stats::sd(stim_df$module)
  sd2  <- stats::sd(unstim_df$module)
  sp   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) )
  
  #correct for small sample size [from Cohen's D to Hedge's G]
  cf   <- 1 - 3/( 4*(n1 + n2) - 9 ) # small sample size correction
  g    <- cf * diff/sp
  se.g <- sqrt( (n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94) )
  
  if (is.na(g)) {
    next
  } else {
    dist_random_es <- append(dist_random_es, g)
  }
  
}


# plot(density(dist_random_es))
# abline(v = g_real, col = "red", lwd = 2, lty = 2)  # Vertical line
# 

# Compute the sample mean and standard deviation
sample_mean <- mean(dist_random_es)
sample_sd <- sd(dist_random_es)

# Calculate the z-score
z_score <- (g_real - sample_mean) / sample_sd

# Calculate the p-value (two-tailed)
p_value <- 2 * (1 - pnorm(abs(z_score)))

df <- data.frame('density' = dist_random_es)


p <- ggplot(df, aes(x = density)) +
  geom_density(fill = "orange", alpha = 0.3) +
  theme_classic() +
  geom_vline(xintercept = g_real, color = "red", linetype = "dashed", linewidth = 0.8) +  
  labs(title = "Density Plot - Z-Score Comparison. -  48h", x = "Hedges' g", y = "Density") +
  annotate("text", #x = g_real + 0.5, y = 0.1, 
           x = Inf, y = Inf,hjust = 1, vjust = 1,
           label = paste("p =", signif(p_value, digits = 3), 
                         "\nz =", signif(z_score, digits = 3)), 
           color = "red", size = 5
           ) 

jpeg(file = 'SuppFig10_48h.png', res = 400, width = 3000, height =2000)
p
dev.off()
