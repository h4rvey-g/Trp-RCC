get_sub_cluster <- function(sc_opt, res) {
    sc_T_cell <- sc_opt %>%
        filter(cell_type_opt == "T-cell") %>%
        FindClusters(
            resolution = res,
            method = "igraph",
            algorithm = 4,
            cluster.name = "T_cell_sub_cluster"
        )
    # Run FindAllMarkers
    T_sub_markers <- FindAllMarkers(
        object = sc_T_cell,
        only.pos = TRUE,
        min.pct = 0.25,
        verbose = FALSE
    )
    # filter top 5 features of each sub cluster
    top_markers <- T_sub_markers %>%
        filter(p_val_adj < 0.05) %>%
        group_by(cluster) %>%
        top_n(5, avg_log2FC) %>%
        ungroup() %>%
        pull(gene) %>%
        unique()
    # DotPlot
    p <- DotPlot(
        object = sc_T_cell,
        features = top_markers,
        assay = "RNA",
        dot.scale = 4
    ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white"),
        )
    ggsave(
        filename = paste0("result/205.sub_cluster/T_cell_sub_cluster_", res, ".png"),
        plot = p,
        width = 20,
        height = 10
    )
    top_markers_10 <- T_sub_markers %>%
        filter(p_val_adj < 0.05) %>%
        group_by(cluster) %>%
        top_n(10, avg_log2FC) %>%
        ungroup()
    write_tsv(top_markers_10, paste0("result/205.sub_cluster/T_cell_sub_cluster_", res, ".tsv"))
    paste0("result/205.sub_cluster/T_cell_sub_cluster_", res, ".png")
}

get_scGate_T <- function(sc_opt) {
    sc_T_cell <- sc_opt %>%
        filter(cell_type_opt == "T-cell") %>%
        RunUMAP(
            dims = 1:30,
            reduction = "scvi",
            n.neighbors = 30,
            reduction.name = "umap_T_cell"
        ) %>%
        FindClusters(
            resolution = 0.08,
            method = "igraph",
            algorithm = 4,
            cluster.name = "T_cell_sub_cluster"
        )
    p <- DimPlot(sc_T_cell,
        group.by = c("T_cell_sub_cluster", "cell_type_opt"), label = TRUE,
        reduction = "umap_T_cell"
    ) +
        # no legend
        theme(legend.position = "none")
    ggsave("result/205.sub_cluster/T_cell_all.png", p, width = 20, height = 10)
    # models.DB <- scGate::get_scGateDB()
    # use qs to save models.DB
    # qs::qsave(models.DB, file = "data/205.sub_cluster/models.DB.qs")
    # load
    models.DB <- qs::qread("data/205.sub_cluster/models.DB.qs")
    models.list <- c(models.DB$human$CD8_TIL, models.DB$human$CD4_TIL)
    sc_T_cell <- scGate(sc_T_cell, models.list, ncores = 10, reduction = "umap_T_cell")
    p1 <- DimPlot(sc_T_cell, group.by = "T_cell_sub_cluster", label = TRUE, reduction = "umap_T_cell") +
        # no legend
        theme(legend.position = "none")
    p2 <- DimPlot(sc_T_cell, group.by = "scGate_multi", label = TRUE, reduction = "umap_T_cell")
        # no legend
        # theme(legend.position = "none")
    p <- p1 + p2
    ggsave(
        filename = "result/205.sub_cluster/T_cell_sub_cluster_scGate.png",
        plot = p,
        width = 20,
        height = 10
    )
    sc_T_cell
}
