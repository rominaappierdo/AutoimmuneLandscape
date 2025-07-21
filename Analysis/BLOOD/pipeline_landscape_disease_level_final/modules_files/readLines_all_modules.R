mod1 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod1/mod1.txt")
mod2 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod2/mod2.txt")
mod3 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod3/mod3.txt")
mod4 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod4/mod4.txt")
mod5 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod5/mod5.txt")
mod6 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod6/mod6.txt")
mod7 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod7/mod7.txt")
mod8 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod8/mod8.txt")
mod9 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod9/mod9.txt")
mod10 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod10/mod10.txt")
mod11 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod11/mod11.txt")
mod12 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod12/mod12.txt")
mod13 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod13/mod13.txt")
mod14 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod14/mod14.txt")
mod15 <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod15/mod15.txt")

# 
# library(tidyverse)
# 
# # 1. Get module names
# module_names <- ls(pattern = "^mod\\d+$")
# 
# # 2. Build the data frame
# df_modules <- module_names %>%
#   set_names() %>%                     # sets names to the same vector (auto-naming)
#   map(~ get(.x)) %>%                  # retrieve each object
#   enframe(name = "module", value = "gene") %>%  # make long table
#   unnest(gene) %>%                    # one gene per row
#   group_by(module) %>%
#   mutate(row = row_number()) %>%      # row ID within each module
#   pivot_wider(names_from = module, values_from = gene) %>%  # wide format
#   select(order(as.numeric(gsub("\\D", "", names(.)))))      # optional: sort columns mod1-mod15
# 
# df_modules$row <- NULL
# colnames(df_modules) <- paste0("Module_", c(1:15))
# df_modules[is.na(df_modules)] <- ""
# 
# write.table(df_modules, file = 'SupplementaryFile87.tsv', row.names = F, col.names = T, sep = '\t', quote = F)