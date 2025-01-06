library(Seurat)
library(tidyseurat)
library(patchwork)
library(SeuratExtend)
tar_load(sc_opt)
genes <- c("FIBCD1", "HAPLN3", "CDYL2", "FHOD3", "SIT1")
# 获取基因表达数据
library(furrr)
plan(multisession, workers = 5)
options(future.globals.maxSize = 40 * 1024^3)
process_gene <- function(gene) {
    library(tidyseurat)
    gene_expression <- FetchData(sc_opt, vars = gene)

    # 定义表达量大于0的细胞
    cells.highlight <- rownames(gene_expression)[gene_expression[[gene]] > 0]
    p1 <- DimPlot2(sc_opt,
        features = c(gene, "cell_type_opt"), label = TRUE,
        # 高亮表达量大于0的细胞
        cells.highlight = cells.highlight,
        nrow = 1
    )
    p2 <- VlnPlot2(sc_opt, features = c(gene), group.by = "cell_type_opt")
    toplot <- CalcStats(sc_opt, features = gene, method = "zscore", order = "p", group.by = "cell_type_opt")
    p3 <- Heatmap(toplot, lab_fill = "zscore")
    toplot <- CalcStats(
        sc_opt %>%
            filter(cell_type_opt == "RCC"),
        features = gene, method = "mean", order = "p", group.by = "batch"
    )
    p4 <- Heatmap(toplot, lab_fill = "mean")
    p <- (p1 / p2) / p3 / p4 + plot_layout(heights = c(1, 1, 0.2, 0.2))

    ggsave(paste0("result/207.gene/cell_type_opt_", gene, ".png"), p, width = 20, height = 14)

    # see percentage of AGER expression in each cell type
    high_cells <- cells.highlight
    sc_opt <- sc_opt %>%
        mutate(AGER_exp = ifelse(cell %in% high_cells, paste0(gene, "_exp"), paste0(gene, "_not_exp")))

    p <- ClusterDistrBar(origin = sc_opt$cell_type_opt, cluster = sc_opt$AGER_exp)
    ggsave(paste0("result/207.gene/exp_", gene, ".png"), p, width = 10, height = 6)
}

future_map(genes, process_gene)

# AGER_dist <- function(sc_opt, dataset_one) {
#     sc_opt <- sc_opt %>%
#         filter(dataset == dataset_one)
#     p <- ClusterDistrBar(origin = sc_opt$cell_type_opt, cluster = sc_opt$AGER_exp) +
#         ggtitle(dataset_one)
#     p
# }
# p_list <- map(unique(sc_opt$batch), ~ AGER_dist(sc_opt, .x))
# ggsave(paste0("result/207.gene/AGER_exp_dist_", gene, ".png"), patchwork::wrap_plots(p_list, ncol = 3), width = 20, height = 20)
