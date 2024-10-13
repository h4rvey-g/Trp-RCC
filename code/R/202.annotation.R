run_annotation_train <- function() {
    # run code/Python/202.annotation.ipynb
    # use_python("/home/vscode/.pyenv/shims/python", required = TRUE)
    # py_run_file("code/Python/202.annotation.py")
    "data/202.annotation/predicted_labels.csv"
}

combine_annotation <- function(sc, predicted_labels_path) {
    predicted_labels <- readr::read_csv(predicted_labels_path) %>%
        # rename first column to cell
        rename(cell_name = "...1") %>%
        rename(cell_type = majority_voting) %>% 
        dplyr::select("cell_name", "cell_type")
    sc_anno <- left_join(sc, predicted_labels, by = c(".cell" = "cell_name"))
    sc_anno <- FindNeighbors(sc_anno, dims = 1:10, reduction = "scvi")
    sc_anno <- FindClusters(sc_anno, resolution = 0.3)
    sc_anno <- RunUMAP(sc_anno, dims = 1:10, reduction = "scvi")
    p <- DimPlot(sc_anno, reduction = "umap", group.by = "cell_type", label = TRUE)
    dir.create("result/202.annotation", showWarnings = FALSE)
    ggsave("result/202.annotation/umap_annotation_v1.png", p, width = 10, height = 10)
    sc_anno
}
