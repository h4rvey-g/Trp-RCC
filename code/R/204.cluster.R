run_clusters <- function(sc_anno) {
    sc_cluster <- sc_anno %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 0.1, algorithm = 4, method = "igraph") %>%
        FindClusters(resolution = 0.2, algorithm = 4, method = "igraph")
    sc_cluster
}

cluster_tree <- function(sc_cluster) {
    p1 <- DimPlot(sc_cluster, group.by = "cell_type", label = TRUE, reduction = "umap")
    p2 <- DimPlot(sc_cluster, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap")
    ggsave("result/204.cluster/cluster_res_0.1.png", p1 + p2, width = 20, height = 10)

    p2 <- DimPlot(sc_cluster, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap")
    ggsave("result/204.cluster/cluster_res_0.2.png", p1 + p2, width = 20, height = 10)

    # p <- clustree(sc_cluster, prefix = "RNA_snn_res.")
    # ggsave("result/204.cluster/clustree.png", p, width = 20, height = 15)
    # p <- clustree_overlay(sc_cluster,
    #     prefix = "RNA_snn_res.", exprs = "data",
    #     x_value = "pca_1", y_value = "pca_2", red_dim = "pca"
    # )
    # ggsave("result/204.cluster/clustree_overlay_umap.png", p, width = 10, height = 10)
    "result/204.cluster"
}

optimize_clusters <- function(sc_cluster) {
    # marker_genes <- c(
    #     "CD19", "CD20", "MS4A1", "CD79A", "CD79B", # B cell
    #     "CD3D", "CD3E", "CD3G", "CD3", "TRBC2", # T cell
    #     "NKG7", "KLRD1", "KLRF1", "GNLY", "NCR1", # NK cell
    #     "PECAM1", "CD34", "CD31", "KDR", # endothelial cell
    #     "MZB1", "IGHA1", "IGHG1", # plasma cell
    #     "CD14", "CD16", "LYZ", # monocyte
    #     "CD68", "C1QA", # macrophage
    #     "CD1C", "CLEC4C", "IRF8", # dendritic cell
    #     "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R", # Granulocyte
    #     "TPSB2", "TPSAB1", "MS4A2", # Mast cell
    #     "ACTA2", "COL1A1", "COL1A2", "TAGLN", # Fibroblast
    #     "PDZK1IP1", "LRP2", "ALDOB", # proximal tubule cell
    #     "SLC4A4", "SLC5A12", "SLC22A19", # kidney cell
    #     "CA9", "NDUFA4L2", "VEGFA" # RCC cell
    # )
    # p <- DotPlot(sc_cluster, features = marker_genes, group.by = "RNA_snn_res.0.2", idents = c(1, 37, 6))
    # p <- p + theme(plot.background = element_rect(fill = "white")) +
    #     # make y axis in the order of "B-cell", "T-cell", "NK", "EC", "Plasma", "Myeloid", "pDC", "Mast", "Fibro", "RCC", "Epi_PT", "Epi_non-PT"
    #     # scale_y_discrete(limits = c("B-cell", "T-cell", "NK", "EC", "Plasma", "Myeloid", "pDC", "Mast", "Fibro", "RCC", "Epi_PT", "Epi_non-PT")) +
    #     # x axis label angle 45
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # ggsave("result/204.cluster/T_or_NK.png", p, width = 30, height = 10)
    # p <- DotPlot(sc_cluster, features = marker_genes, group.by = "RNA_snn_res.0.2", idents = c(12, 40, 23))
    # p <- p + theme(plot.background = element_rect(fill = "white")) +
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # ggsave("result/204.cluster/T_or_RCC.png", p, width = 30, height = 10)
    # p <- DotPlot(sc_cluster, features = marker_genes, group.by = "RNA_snn_res.0.2", idents = c(11, 28, 25, 36, 33, 26, 40))
    # p <- p + theme(plot.background = element_rect(fill = "white")) +
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # ggsave("result/204.cluster/TBD.png", p, width = 30, height = 10)
    # p <- DotPlot(sc_cluster, features = marker_genes, group.by = "RNA_snn_res.0.2", idents = c(21, 22, 23, 29, 31, 32, 34, 14, 38, 19, 8, 17, 27, 35, 15))
    # p <- p + theme(plot.background = element_rect(fill = "white")) +
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # ggsave("result/204.cluster/RCC.png", p, width = 30, height = 10)

    # # use FindMarkers to find marker genes for cluster 10, 11, 13, 26, 28, 33, 36, 37, 40
    # cluster_markers <- map(c(10, 11, 13, 26, 28, 33, 36, 37, 40), ~ FindMarkers(sc_cluster, ident.1 = .x))
    # cluster_markers <- setNames(cluster_markers, c(10, 11, 13, 26, 28, 33, 36, 37, 40))
    # cluster_markers_df <- bind_rows(cluster_markers, .id = "cluster")
    # # filter by p_val_adj < 0.05, arrange by desc abs(avg_logFC), then get top 10 genes of each cluster
    # cluster_markers_df_top <- cluster_markers_df %>%
    #     as_tibble(rownames = "gene") %>%
    #     # remove the "...\d+" part in gene
    #     mutate(gene = str_remove(gene, "\\.\\.\\.\\d+")) %>%
    #     # str_remove("\\.\\d+")) %>%
    #     filter(p_val_adj < 0.05) %>%
    #     arrange(desc(avg_log2FC)) %>%
    #     group_by(cluster) %>%
    #     slice_max(avg_log2FC, n = 20) %>%
    #     ungroup() %>%
    #     mutate(gene_fix = transId(gene, transTo = "symbol", keepNA = TRUE, unique = TRUE) %>% pull(symbol))
    # # write to xlsx
    # writexl::write_xlsx(cluster_markers_df_top, "data/204.cluster/cluster_markers.xlsx")
    # # see featureplot of cluster 11
    # p <- FeaturePlot(sc_cluster,
    #     features = c("PECAM1", "CD14", "LYZ", "CD68", "C1QA", "S100A9", "CSF3R"),
    #     cells = WhichCells(sc_cluster, ident = 11)
    # )
    # ggsave("result/204.cluster/cluster_11.png", p, width = 10, height = 10)
    
    # myeloid_marker <- read_tsv("data/204.cluster/myeloid_markers.tsv")
    # # group by each first_type of myeloid_marker, then run DotPlot on gene
    # p <- myeloid_marker %>%
    #     group_by(first_type) %>%
    #     group_map(~ {
    #         plot <- DotPlot(sc_cluster, features = .x$marker, group.by = "RNA_snn_res.0.2")
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #     }) %>%
    #     wrap_plots(., ncol = 1)
    # ggsave("result/204.cluster/myeloid_marker_all.png", p, width = 10, height = 40)
    # # group by second_type of myeloid_marker, then run DotPlot on gene
    # p <- myeloid_marker %>%
    #     filter(second_type != "None") %>%
    #     group_by(second_type) %>%
    #     group_map(~ {
    #         plot <- DotPlot(sc_cluster, features = .x$marker, group.by = "RNA_snn_res.0.2")
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #     }) %>%
    #     wrap_plots(., ncol = 1)
    # ggsave("result/204.cluster/myeloid_marker_second.png", p, width = 10, height = 40)
    # # group by third_type of myeloid_marker, then run DotPlot on gene
    # p <- myeloid_marker %>%
    #     filter(third_type != "None") %>%
    #     group_by(third_type) %>%
    #     group_map(~ {
    #         plot <- DotPlot(sc_cluster, features = .x$marker, group.by = "RNA_snn_res.0.2")
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #     }) %>%
    #     wrap_plots(., ncol = 1)
    # ggsave("result/204.cluster/myeloid_marker_third.png", p, width = 10, height = 40)
    # # draw FeaturePlot of myeloid_marker second_type, on cluster 3, 28, 11
    # p <- myeloid_marker %>%
    #     filter(second_type != "None") %>%
    #     group_by(second_type) %>%
    #     group_walk(~ {
    #         plot <- FeaturePlot(sc_cluster, features = .x$marker, cells = WhichCells(sc_cluster, ident = c(3, 28, 11)))
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #         ggsave(paste0("result/204.cluster/myeloid_marker_second_featureplot_", .y, ".png"), plot, width = 10, height = 10)
    #     })
    
    # kidney_marker <- read_tsv("data/204.cluster/kidney_markers.tsv")
    # # group by each first_type of kidney_marker, then run DotPlot on gene
    # p <- kidney_marker %>%
    #     group_by(first_type) %>%
    #     group_map(~ {
    #         plot <- DotPlot(sc_cluster, features = .x$marker, group.by = "RNA_snn_res.0.2")
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #     }) %>%
    #     wrap_plots(., ncol = 1)
    # ggsave("result/204.cluster/kidney_marker_all.png", p, width = 10, height = 40)
    # # group by second_type of kidney_marker, then run DotPlot on gene
    # p <- kidney_marker %>%
    #     filter(second_type != "None") %>%
    #     group_by(second_type) %>%
    #     group_map(~ {
    #         plot <- DotPlot(sc_cluster, features = .x$marker, group.by = "RNA_snn_res.0.2")
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #             # remove legend
    #             theme(legend.position = "none")
    #     }) %>%
    #     wrap_plots(., ncol = 3)
    # ggsave("result/204.cluster/kidney_marker_second.png", p, width = 30, height = 40)

    # lymphoid_marker <- read_tsv("data/204.cluster/lymphoid_markers.tsv")
    # # group by each first_type of lymphoid_marker, then run DotPlot on gene
    # p <- lymphoid_marker %>%
    #     group_by(first_type) %>%
    #     group_map(~ {
    #         plot <- DotPlot(sc_cluster, features = .x$marker, group.by = "RNA_snn_res.0.2")
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #     }) %>%
    #     wrap_plots(., ncol = 1)
    # ggsave("result/204.cluster/lymphoid_marker_all.png", p, width = 10, height = 40)
    # # group by second_type of lymphoid_marker, then run DotPlot on gene
    # p <- lymphoid_marker %>%
    #     filter(second_type != "None") %>%
    #     group_by(second_type) %>%
    #     group_map(~ {
    #         plot <- DotPlot(sc_cluster, features = .x$marker, group.by = "RNA_snn_res.0.2")
    #         plot + ggtitle(.y) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #     }) %>%
    #     wrap_plots(., ncol = 1)
    # ggsave("result/204.cluster/lymphoid_marker_second.png", p, width = 10, height = 40)
    "result/204.cluster/final_annotation.tsv"
}

get_final_annotation <- function(sc_cluster, final_annotation) {
    final_annotation_df <- read_tsv(final_annotation) %>%
        # convert cluster to factor
        mutate(cluster = as.factor(cluster))
    # left join to sc_anno, by cluster
    sc_opt <- sc_cluster %>%
        left_join(final_annotation_df, by = c("RNA_snn_res.0.2" = "cluster")) %>%
        # if cell_type_opt is "celltypist", then use cell_type
        mutate(cell_type_opt = if_else(cell_type_opt == "cell_typist", cell_type, cell_type_opt)) %>%
        # convert Epi_PT and Epi_non-PT to Nephron
        mutate(cell_type_opt = if_else(cell_type_opt %in% c("Epi_PT", "Epi_non-PT"), "Nephron", cell_type_opt)) %>%
        # convert Myeloid to Monocytic lineage
        mutate(cell_type_opt = if_else(cell_type_opt == "Myeloid", "Monocytic lineage", cell_type_opt))
    # drop type == "normal" but cell_type_opt == "RCC"
    sc_opt <- sc_opt %>%
        filter(!(type == "normal" & cell_type_opt == "RCC"))
    # plot DimPlot with cell_type and cell_type_opt
    p <- DimPlot(sc_opt, group.by = c("cell_type", "cell_type_opt"), label = TRUE, reduction = "umap")
    ggsave("result/202.annotation/cluster_opt.png", p, width = 20, height = 10)
    # save cell of cell_type_opt == "T-cell" to data/204.cluster/T_cell.tsv
    sc_opt %>%
        filter(cell_type_opt == "T-cell") %>%
        as_tibble() %>%
        mutate(cell_name = cell) %>%
        select(cell_name) %>%
        write_tsv("data/204.cluster/T_cell_names.tsv")
    sc_opt
}
