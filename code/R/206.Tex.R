get_T_monocle <- function(sc_T_anno) {
    cds <- as.cell_data_set(sc_T_anno)
    cds <- cluster_cells(cds, resolution = 0.4)
    # change clusters(cds) to colData(cds)$seurat_clusters, with names attr
    seurat_clusters <- colData(cds)$T_cell_sub_cluster
    names(seurat_clusters) <- rownames(colData(cds))
    cds@clusters@listData$UMAP$clusters <- seurat_clusters
    cds@int_colData@listData$reducedDims$UMAP <- sc_T_anno@reductions$umap_T_cell@cell.embeddings
    p1 <- plot_cells(cds,
        reduction_method = "UMAP", color_cells_by = "cell_type_dtl", show_trajectory_graph = FALSE, label_cell_groups = FALSE,
        label_groups_by_cluster = FALSE
    )
    p2 <- plot_cells(cds,
        reduction_method = "UMAP", color_cells_by = "partition", show_trajectory_graph = FALSE, label_cell_groups = FALSE,
        label_groups_by_cluster = FALSE
    )
    patchwork::wrap_plots(p1, p2) %>%
        ggsave(filename = "result/206.Tex/cells_partition.png", width = 13)

    cds_trajectory <- learn_graph(cds, use_partition = FALSE, learn_graph_control = list(nn.cores = 30))

    # order cds, root is M0
    cds_trajectory <- order_cells(cds_trajectory,
        root_cells = colnames(cds_trajectory[, colData(cds_trajectory)$cell_type_dtl %in% c("CD8+T_EFF")])
    )
    p <- plot_cells(cds_trajectory, color_cells_by = "cell_type_dtl", show_trajectory_graph = TRUE, label_cell_groups = FALSE)
    ggsave(p, filename = "result/206.Tex/trajectory_cell_types.png", width = 7)
    p <- plot_cells(cds_trajectory, color_cells_by = "cell_type_dtl", show_trajectory_graph = TRUE, label_roots = FALSE)
    ggsave(p, filename = "result/206.Tex/trajectory_cell_types_2.png", width = 7)
    p <- plot_cells(cds_trajectory,
        color_cells_by = "pseudotime",
        group_cells_by = "cell_type_dtl",
        label_cell_groups = FALSE,
        label_groups_by_cluster = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        label_roots = FALSE,
        trajectory_graph_color = "grey60"
    )
    ggsave(p, filename = "result/206.Tex/trajectory_pseudotime.png", width = 7)
    # p <- ggplot()
    # cds_trajectory
    "results/07.monocle/trajectory_cell_types.png"
}

DEG_Tex <- function(sc_T_anno, sc_opt, Tex_score) {
    sc_opt <- sc_opt %>%
        left_join(sc_T_anno %>% as_tibble() %>% select(cell, cell_type_dtl), by = "cell")
    # fill sc_opt cell_type_dtl NA with values in cell_type
    sc_opt <- sc_opt %>%
        mutate(cell_type_dtl = case_when(
            is.na(cell_type_dtl) ~ cell_type,
            TRUE ~ cell_type_dtl
        ))
    sc_opt <- sc_opt %>%
        left_join(Tex_score, by = c("dataset" = "dataset"))
    sc_opt <- sc_opt %>%
        mutate(Tex_group = case_when(
            Tex_value > 0 ~ "Tex_high",
            Tex_value <= 0 ~ "Tex_low"
        ))
    sc_opt <- sc_opt %>%
        filter(cell_type_dtl == "RCC")
    pseudo_sc <- AggregateExpression(sc_opt, return.seurat = TRUE, group.by = c("dataset", "Tex_group"))
    Tex_marker <- pseudo_sc %>%
        FindMarkers(
            group.by = "Tex_group",
            ident.1 = "Tex-high",
            ident.2 = "Tex-low",
            min.pct = 0.25,
            logfc.threshold = 0.25,
            test.use = "DESeq2"
        )
    Tex_marker <- Tex_marker %>%
        rownames_to_column("gene") %>%
        as_tibble()
    Tex_marker_2 <- sc_opt %>%
        FindMarkers(
            group.by = "Tex_group",
            ident.1 = "Tex_high",
            ident.2 = "Tex_low",
            min.pct = 0.25,
            logfc.threshold = 0.25,
            reduction = "scvi"
        )
    Tex_marker_2 <- Tex_marker_2 %>%
        rownames_to_column("gene") %>%
        as_tibble()
    # find genes not in Tex_marker but in Tex_marker_2
    Tex_marker_2 <- Tex_marker_2 %>%
        anti_join(Tex_marker, by = "gene")
    Tex_marker <- bind_rows(Tex_marker, Tex_marker_2)

    Tex_marker <- Tex_marker %>%
        arrange(desc(abs(avg_log2FC))) %>%
        filter(abs(avg_log2FC) > 1 & p_val < 0.05)
    write_tsv(Tex_marker, "result/206.Tex/Tex_marker.tsv")
    Tex_marker
}

Tex_gsea <- function(Tex_marker, msigdb) {
    gene_up <- Tex_marker %>%
        filter(avg_log2FC > 0) %>%
        pull(gene) %>%
        genekitr::transId(transTo = "symbol", unique = TRUE)
    gene_down <- Tex_marker %>%
        filter(avg_log2FC < 0) %>%
        pull(gene) %>%
        genekitr::transId(transTo = "symbol", unique = TRUE)
    gene_all <- c(gene_up, gene_down)
    gs <- geneset::getMsigdb(org = "human")
    gsea_up <- genGSEA(gene_up, gs)
    gsea_down <- genGSEA(gene_down, gs)
    gsea_all <- genGSEA(gene_all, gs)
    p <- plotGSEA(gsea_up, plot_type = "ridge", stats_metric = "p.adjust", show_pathway = nrow(gsea_up$gsea_df))
    ggsave("result/108.direct_enrich/gsea_up.png", p)
    p <- plotGSEA(gsea_up, plot_type = "bar", colour = c("navyblue", "orange"))
    ggsave("result/108.direct_enrich/gsea_up_bar.png", p)
    write_tsv(gsea_up$gsea_df, "data/108.direct_enrich/gsea_up.tsv")
    gs$geneset <- msigdb
    gsea_up <- genGSEA(gene_up, gs)
}
