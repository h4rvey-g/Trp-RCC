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

test_annotation <- function(sc_anno) {
    marker_genes <- c(
        "CD19", "CD20", "MS4A1", "CD79A", "CD79B", # B cell
        "CD3D", "CD3E", "CD3G", "CD3", "TRBC2", # T cell
        "NKG7", "KLRD1", "KLRF1", "GNLY", "NCR1", # NK cell
        "PECAM1", "CD34", "CD31", "KDR", # endothelial cell
        "MZB1", "IGHA1", "IGHG1", # plasma cell
        "CD14", "CD16", "LYZ", # monocyte
        "CD68", "C1QA", # macrophage
        "CD1C", "CLEC4C", "IRF8", # dendritic cell
        "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R", # Granulocyte
        "TPSB2", "TPSAB1", "MS4A2", # Mast cell
        "ACTA2", "COL1A1", "COL1A2", "TAGLN", # Fibroblast
        "SLC4A4", "SLC5A12", "SLC22A19", "LRP2", "ALDOB" # kidney cell
    )
    p <- DotPlot(sc_anno, features = marker_genes, group.by = "cell_type")
    p <- p + theme(plot.background = element_rect(fill = "white")) +
        # make y axis in the order of "B-cell", "T-cell", "NK", "EC", "Plasma", "Myeloid", "pDC", "Mast", "Fibro", "RCC", "Epi_PT", "Epi_non-PT"
        scale_y_discrete(limits = c("B-cell", "T-cell", "NK", "EC", "Plasma", "Myeloid", "pDC", "Mast", "Fibro", "RCC", "Epi_PT", "Epi_non-PT")) +
        # x axis label angle 45
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("result/202.annotation/marker_dotplot.png", p, width = 30, height = 10)
    "result/202.annotation/marker_dotplot.png"
}
