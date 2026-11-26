setwd("/mnt/data3/appierdo/AutoimmuneLandscape/revisioni")

library(ggplot2)
library(ggVennDiagram)

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


# --------------------------------------------------
# VENN DIAGRAM
# --------------------------------------------------

gene_lists <- list(
  "Type I IFN-alpha\n(Hallmark)" = alpha,
  "Type II IFN-gamma\n(Hallmark)" = gamma_hallmark,
  "Type II IFN-gamma\n(WikiPathways)" = gamma_wiki
)

ggVennDiagram(gene_lists)



all_three <- Reduce(
  intersect,
  list(
    alpha,
    gamma_hallmark,
    gamma_wiki
  )
)


p <- ggVennDiagram(gene_lists) +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen"
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(
      t = 30,
      r = 50,
      b = 30,
      l = 80
    )
  )

p
