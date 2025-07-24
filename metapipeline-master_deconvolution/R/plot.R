plot_gene_expression <- function(data, sel.gene) {
    df <- list()

    df <- lapply(names(data$originalData), function(dataset) {
        w <- which(data$originalData[[dataset]]$keys == sel.gene)

        if(any(w)) {
            gx <- data$originalData[[dataset]]$expr[w,, drop = F]
            # colMeans is necessary in case there are multiple probes that map to the same gene
            data.frame(gx = as.vector(colMeans(gx)), class = factor(data$originalData[[dataset]]$class[colnames(gx)]))
        }
        else
            NULL
    })
    names(df) <- names(data$originalData)
    df <- bind_rows(df, .id = "dataset")
    df$group <- "controls"
    df[df$class == 1, "group"] <- "cases"


    (p <- ggplot2::ggplot(ggplot2::aes(x = group, y = gx), data = df)
        + ggplot2::geom_boxplot()
        + ggplot2::geom_jitter(ggplot2::aes(color = group))
        + ggplot2::facet_wrap(~dataset, scales = "free_y")
        + ggplot2::labs(title = sel.gene)
    )
}


plot_genes_of_interest <- function(meta.res, genes.of.interest, filename.prefix) {
    for(gene.group in names(genes.of.interest)) {
        sel.genes <- genes.of.interest[[gene.group]]
        pdf(sprintf("%s_%s_genes_plot.pdf", filename.prefix, gene.group), width = 11)
        par(mfrow = c(2, ceiling(length(sel.genes) / 2)))
        for(g in sel.genes)
            if(g %in% row.names(meta.res$metaAnalysis$pooledResults))
                MetaIntegrator::forestPlot(meta.res, g)

        dev.off()
    }
}

forestplot_top_genes <- function(meta.res) {
    x <- meta.res$metaAnalysis$pooledResults
    x <- x[x$effectSizeFDR < 0.05 & abs(x$effectSize) > 1, ]

    if(!(all(x$numStudies == 1)))
        x <- x[x$numStudies > 1, ]

    x <- x[order(x$effectSize, decreasing = T), ]

    for(g in row.names(x))
        MetaIntegrator::forestPlot(meta.res, g)
}
