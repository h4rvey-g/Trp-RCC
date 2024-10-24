run_clusters <- function(sc_anno) {
    sc_cluster <- sc_anno %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 0.1, algorithm = 4, method = "igraph")
    sc_cluster
}

cluster_tree <- function(sc_cluster) {
    p1 <- DimPlot(sc_cluster, group.by = "cell_type", label = TRUE, reduction = "umap")
    p2 <- DimPlot(sc_cluster, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap")
    ggsave("result/204.cluster/cluster_res_0.1.png", p1+ p2, width = 20, height = 10)

    p <- clustree(sc_cluster, prefix = "RNA_snn_res.")
    ggsave("result/204.cluster/clustree.png", p, width = 20, height = 15)
    p <- clustree_overlay(sc_cluster,
        prefix = "RNA_snn_res.", exprs = "data",
        x_value = "pca_1", y_value = "pca_2", red_dim = "pca"
    )
    ggsave("result/204.cluster/clustree_overlay_umap.png", p, width = 10, height = 10)
}