setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/meta_analysis_gene_expression/lupus")
library(qs)
a <- qread(("lupus_vs_healthy_whole blood_metaintegrator_object.qs"))

l <- a$originalData
for (gse in names(l)) {
  print(gse)
  print("cases")
  print(sum(l[[gse]]$class))
  print("controls")
  print(sum(l[[gse]]$class == 0))
  
}