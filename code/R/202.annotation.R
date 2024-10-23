run_annotation_train <- function() {
    # run code/Python/202.annotation.ipynb
    # use_python("/home/vscode/.pyenv/shims/python", required = TRUE)
    # py_run_file("code/Python/202.annotation.py")
    "data/202.annotation/predicted_labels.csv"
}

combine_annotation <- function(sc, predicted_labels_path) {
    predicted_labels <- readr::read_csv(predicted_labels_path) %>%
        # rename first column to cell
        rename(cell_name_label = "...1") %>%
        rename(cell_type = majority_voting) %>% 
        dplyr::select("cell_name_label", "cell_type")
    sc_anno <- left_join(sc, predicted_labels, by = c(".cell" = "cell_name_label"))
    # for rows with NA in cell_type, fill with value in predict
    sc_anno$cell_type[is.na(sc_anno$cell_type)] <- sc_anno$predict[is.na(sc_anno$cell_type)]
    # convert tumor to RCC
    sc_anno <- sc_anno %>%
        mutate(cell_type = if_else(cell_type == "tumor", "RCC", cell_type))
    sc_anno <- FindNeighbors(sc_anno, dims = 1:10, reduction = "scvi")
    sc_anno <- FindClusters(sc_anno, resolution = 0.3)
    sc_anno <- RunUMAP(sc_anno, dims = 1:10, reduction = "scvi")
    p <- DimPlot(sc_anno, reduction = "umap", group.by = "cell_type", label = TRUE)
    dir.create("result/202.annotation", showWarnings = FALSE)
    ggsave("result/202.annotation/umap_annotation.png", p, width = 10, height = 10)
    sc_anno
}
