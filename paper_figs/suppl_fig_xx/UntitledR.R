setwd("paper_figs/suppl_fig_xx")

files <- list.files("../../Analysis/", pattern = 'complete_res', recursive = T, full.names = T)
files <- files[grepl('gene', files)]

diseases_blood <- c('crohn_vs_healthy_w')
diseases_tissue <- c('GUT_crohn')

for (i in 1:length(diseases_blood)) {
  fb <- files[grepl(diseases_blood[i], files)]
  ft <- files[grepl(diseases_tissue[i], files)]
  
  b <- read.csv(fb, sep = '\t')
  b <- b[b$numStudies>=2, ]
  b <- b[b$effectSizeFDR < 0.05, ]
  t <- read.csv(ft, sep = '\t')
  t <- t[t$numStudies>=2, ]
  t <- t[t$effectSizeFDR < 0.05, ]
  
  merged <- merge(b, t, by = 'gene')
  
  cor(merged$effectSize.x, merged$effectSize.y)
  
  stop()
}