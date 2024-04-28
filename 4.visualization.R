##R version 4.3.2

library(scMINER)
library(NetBID2)
library(ggplot2)

load('/mnt/sda/Public/Project/B-ALL/Bdev/activity/acs_sc')
pData(acs_sc) <- filter(pData(acs_sc), pData(acs_sc)[,'scminer'] != 'ActivatedB')
pData(acs_sc)[,'scminer'] <- droplevels(pData(acs_sc)[,'scminer'])
load('SJARACNe/Input.eset')
genes_of_interest <-c("FLT3", "EBF1", "MS4A1","BCL2", "BCL2L1", "MCL1")
genes_of_interest_network <-c("FLT3.SIG", "EBF1.TF", "MS4A1.SIG","BCL2.SIG", "BCL2L1.SIG", "MCL1.SIG")

# 1.feature_heatmap
feature_heatmap <- function (input_eset, target, feature = "geneSymbol", group_name = "label", 
    name = "log2Exp", save_plot = TRUE, width = 4, height = 8, 
    cluster_rows = FALSE, colors = rev(colorRampPalette(brewer.pal(10, 
        "RdYlBu"))(256)), plot_name = "GeneHeatmap.pdf", ...) 
{
    input <- exprs(input_eset)[,rownames(pData(input_eset))]
    gn <- intersect(target, fData(input_eset)[, feature])
    indx <- match(gn, fData(input_eset)[, feature])
    exp <- exprs(input_eset)[indx, ]
    rownames(exp) <- gn
    lab <- pData(input_eset)[, group_name]
    names(lab) <- sampleNames(input_eset)
    ranks <- names(sort(lab, decreasing = FALSE))
    exp.ordered <- as.matrix(exp[, ranks])
    lab.ordered <- lab[ranks]
    df <- data.frame(scMINER = lab.ordered)
    n <- length(unique(lab.ordered))
    ncols <- (scales::hue_pal())(n)
    names(ncols) <- unique(lab.ordered)
    myanndf = HeatmapAnnotation(df = df, col = list(scMINER = ncols))
    mycolors = colors
    hmp <- Heatmap(exp.ordered, col = mycolors, name = name, 
        show_row_names = TRUE, show_column_names = FALSE, cluster_rows = cluster_rows, 
        cluster_columns = FALSE, top_annotation = myanndf, width = width, height = height, ...)
    if (save_plot) {
        pdf(file = plot_name, width = width, height = height)
        ComplexHeatmap::draw(hmp)
        dev.off()
    }
    return(hmp)
}

feature_heatmap(input_eset = input_eset, target = genes_of_interest, group_name = "cellName",save_plot = TRUE, width = 8, height = 4, name = "log2Exp", plot_name='/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/input_heatmap.pdf')
feature_heatmap(input_eset = acs_sc, target = genes_of_interest_network, feature="ID", group_name = "scminer",save_plot = TRUE, width = 8, height = 4, name = "activity", plot_name='/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/acs_sc_heatmap.pdf')

# 2.feature_vlnplot
feature_vlnplot <- function (input_eset, target = NULL, feature = "geneSymbol", 
    group_by = "celltype", ylabel = "Expression", color_by = NULL, 
    colors = NULL, ncol = 3, stat = "median", boxplot = FALSE, 
    title.size = 5) 
{
    if (!group_by %in% colnames(pData(input_eset))) 
        stop("Please check your group_by information!", "\n")
    if (!feature %in% colnames(fData(input_eset))) 
        stop("Please check your feature information!", "\n")
    input <- exprs(input_eset)[,rownames(pData(input_eset))]
    gn <- intersect(target, fData(input_eset)[, feature])
    if (length(indx) == 0) 
        stop("No target feature found in data!", "\n")
    label <- as.factor(pData(input_eset)[, group_by])
    if (is.null(color_by)) 
        color_by = group_by
    condition <- as.factor(pData(input_eset)[, color_by])
    if (length(target) != 1) {
        target_values <- t(as.matrix(input[indx, ]))
        colnames(target_values) <- gn
        df <- data.frame(target_values, cluster = label, condition = condition)
    }
    else {
        target_values <- input[indx, ]
        df <- data.frame(target_values, cluster = label, condition = condition)
        colnames(df)[1] <- gn
    }
    df_melt <- reshape2::melt(df, id.vars = c("cluster", "condition"))
    p <- ggplot(df_melt, aes(x = cluster, y = value, fill = condition)) + 
        theme_classic() + geom_violin(trim = TRUE, scale = "width", 
        na.rm = TRUE, size = 0.4, width = 0.5)
    if (!is.null(stat)) {
        if (stat == "median") 
            p <- p + stat_summary(fun.y = median, geom = "point", 
                size = 1.2, color = "black", position = position_dodge(width = 1))
        else if (stat == "mean") 
            p <- p + stat_summary(fun.y = mean, geom = "point", 
                size = 1.2, color = "black", position = position_dodge(width = 1))
        else cat("Stat not supported, please check your spelling.", 
            "\n")
    }
    if (boxplot) 
        p <- p + geom_boxplot(fill = "white", width = 0.1, size = 0.1, 
            outlier.size = 0.001, show.legend = FALSE, na.rm = TRUE)
    p <- p + facet_wrap(~variable, scales = "free", ncol = ncol) + 
        labs(x = group_by, y = ylabel) + theme(axis.text.x = element_text(size = 10), 
        plot.title = element_text(size = title.size, face = "bold"), 
        strip.background = element_rect(fill = "#FFFFFF"))
    if (!is.null(colors)) 
        p <- p + scale_fill_manual(values = colors)
    if (ylabel == "Activity") {
        p <- p + geom_boxplot(width = 0.2, size = 0.1, outlier.size = 0.001, 
            show.legend = FALSE, na.rm = TRUE)
    }
    return(p)
}
feature_vlnplot(input_eset = input_eset, target = genes_of_interest, 
                feature = "geneSymbol", group_by = "cellName", ylabel = "log2Exp", ncol = 1)
ggsave('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/input_Violin.pdf')

feature_vlnplot(input_eset = acs_sc, target = genes_of_interest_network, 
                feature = "ID", group_by = "scminer", ylabel = "activity", ncol = 1)
ggsave('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/acs_sc_Violin.pdf')

# 3.feature_highlighting: UMAP scatter plot
load('./B_sub.rdata')
pData(input_eset) <- merge(x = pData(input_eset), y = as.data.frame(B_sub@reductions$umap@cell.embeddings), by = 0, all.x = TRUE)
rownames(pData(input_eset)) <- pData(input_eset)[['Row.names']]
pData(input_eset) <- pData(input_eset)[,-1]

pData(acs_sc) <- merge(x = pData(acs_sc), y = as.data.frame(B_sub@reductions$umap@cell.embeddings), by = 0, all.x = TRUE)
rownames(pData(acs_sc)) <- pData(acs_sc)[['Row.names']]
pData(acs_sc) <- pData(acs_sc)[,-1]

feature_highlighting <- function (input_eset, target = NULL, feature = "geneSymbol", 
    x = "X", y = "Y", wrap_by = NULL, ylabel = "Expression", 
    pct.size = 0.8, title.size = 15, ncol = 4, alpha = 0.8, colors = colorRampPalette(c("#E3E3E3", 
        "#BCA2FC", "#4900FE"), interpolate = "linear")(8)) 
{
    input <- as.matrix(exprs(input_eset))
    gn <- intersect(target, fData(input_eset)[, feature])
    indx <- match(gn, fData(input_eset)[, feature])
    id.vars <- c(x, y, wrap_by)
    projection <- pData(input_eset)[colnames(input), id.vars]
    if (length(indx) != 1) {
        target_values <- t(as.matrix(input[indx, ]))
        colnames(target_values) <- gn
        proj_target <- cbind(projection, target_values)
        proj_target_melt <- reshape2::melt(proj_target, id.vars = id.vars)
        p <- ggplot(proj_target_melt, aes_string(x, y)) + theme_classic() + 
            facet_wrap(c("variable", wrap_by), scales = "free", 
                ncol = ncol)
        labs(title = "")
    }
    else {
        target_values <- input[indx, ]
        proj_target <- cbind(projection, target = target_values)
        proj_target_melt <- reshape2::melt(proj_target, id.vars = id.vars)
        p <- ggplot(proj_target_melt, aes_string(x, y)) + theme_classic() + 
            labs(title = target, scales = "free")
        if (!is.null(wrap_by)) 
            p <- p + facet_wrap(c(wrap_by), scales = "free", 
                ncol = ncol)
    }
    p <- p + geom_point(aes(colour = value), size = pct.size, 
        alpha = alpha) + scale_colour_gradientn(colors = colors) + 
        theme(plot.title = element_text(size = title.size, face = "bold"), 
            axis.title = element_text(size = 10), legend.title = element_text(size = 10)) + 
        labs(color = ylabel)
    return(p)
}

feature_highlighting(input_eset = input_eset, target = genes_of_interest, 
                     feature = "gene_short_name", ylabel = "log2Exp", x = "UMAP_1", y = "UMAP_2", pct.size = 0.5, ncol=2)
ggsave('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/input_feature_highlighting.pdf')

feature_highlighting(input_eset = acs_sc, target = genes_of_interest_network, 
                     feature = "ID", ylabel = "activity", x = "UMAP_1", y = "UMAP_2", pct.size = 0.5,ncol=2)
ggsave('/mnt/sda/Public/Project/B-ALL/grant/scMINER/fig/acs_sc_feature_highlighting.pdf')
