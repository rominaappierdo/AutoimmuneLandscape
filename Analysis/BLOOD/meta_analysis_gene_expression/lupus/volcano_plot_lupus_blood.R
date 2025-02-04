setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/meta_analysis_gene_expression/lupus")
library(EnhancedVolcano)
library(ggplot2)
library(MetaIntegrator)
data <- read.csv("lupus_vs_healthy_whole blood_complete_results.tsv", sep = '\t')
data <- data[data$numStudies >= 2,]
data <- data[data$effectSizeFDR < 0.05,]
data <- data[abs(data$effectSize) > 0.8,]
sum(data$effectSize>0)
sum(data$effectSize<0)

stop()
#data <- data[data$heterogeneityPval > 0.05,]


jpeg(filename = "Gene_Expr_volcanoplot.jpeg", res = 700, width = 17000, height = 13000)
EnhancedVolcano(data,
                lab = data$gene,
                #lab = NA,
                x = 'effectSize',
                y = 'effectSizeFDR',
                title = 'Lupus vs Control (Blood)',
                selectLab=data[data$effectSizeFDR<0.05 & abs(data$effectSize) > 1.3,]$gene,
                #selectLab=data[data$gene == 'IFI27',],
                
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 4.0,
                labSize = 10.0,
                col=c('black', 'black', 'black', 'red3'),
                legendPosition = 'top',
                legendLabSize = 16,
                legendIconSize = 8.0,
                drawConnectors = T,
                colAlpha = 1,
                #widthConnectors = 0.5,
                boxedLabels=T,
                subtitle = "Differential gene expression") + 
                ylim(0, 50) +
                xlim(-3.5, 3.5)+
                xlab('pooled_effectSize') + 
                ylab('pvalueFDR') 

dev.off()
stop()

data <- read.csv("metaanalyses_results/results_GeneExpr_meta.tsv", sep = '\t')
data <- data[data$numStudies >=2,]
data <- data[data$effectSizeFDR < 0.05,]
#data <- data[data$heterogeneityPval > 0.05,]
data <- data[abs(data$effectSize) > 1,]
write.table(data$X, file = "DEG_genes_metaanalysis.txt", sep ='\n', quote = F, col.names = F, row.names = F)


write.table(data$X[data$effectSize>0], file = 'upreg_genes.txt', quote = F, col.names = F, row.names = F, sep = '\n')
write.table(data$X[data$effectSize<0], file = 'down_genes.txt', quote = F, col.names = F, row.names = F, sep = '\n')