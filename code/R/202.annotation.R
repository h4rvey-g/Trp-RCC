run_annotation_train <- function() {
    # run code/Python/202.annotation.ipynb
    # use_python("/home/vscode/.pyenv/shims/python", required = TRUE)
    # py_run_file("code/Python/202.annotation.py")
    # list all data/202.annotation/predicted_labels_*.csv files
    list.files("data/202.annotation", pattern = "predicted_labels_.*.csv", full.names = TRUE)
}

combine_annotation <- function(sc, predicted_labels_path) {
    # read all files in predicted_labels_path
    predicted_labels_list <- lapply(predicted_labels_path, readr::read_csv) %>%
        # rename first column to cell
        lapply(function(df) {
            df %>%
                rename(cell_name = "...1") %>%
                rename(cell_type = majority_voting) %>%
                dplyr::select("cell_name", "cell_type")
        })
    predicted_labels <- do.call(rbind, predicted_labels_list)
    sc_anno <- left_join(sc, predicted_labels, by = c(".cell" = "cell_name"))
    sc_anno <- FindNeighbors(sc_anno, dims = 1:10, reduction = "scvi")
    sc_anno <- FindClusters(sc_anno, resolution = 0.3)
    sc_anno <- RunUMAP(sc_anno, dims = 1:10, reduction = "scvi")
    p <- DimPlot(sc_anno, reduction = "umap", group.by = "cell_type", label = TRUE)
    dir.create("result/202.annotation", showWarnings = FALSE)
    ggsave("result/202.annotation/umap_annotation_v1.png", p, width = 10, height = 10)
    sc_anno
}
