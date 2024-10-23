prepare_infer_cnv <- function(sc) {
    all_datasets <- sc %>%
        pull(dataset) %>%
        unique()
    raw_list <- map(all_datasets, function(one_dataset) {
        sc_filtered <- sc %>%
            dplyr::filter(dataset == one_dataset)
        sc_filtered[["RNA"]]$counts
    })
    names(raw_list) <- all_datasets
    # # if detected respective qs files in data/203.inferCNV, remove relative datasets in raw_list
    # if (length(list.files("data/203.inferCNV", pattern = ".qs")) > 0) {
    #     qs_files <- list.files("data/203.inferCNV", pattern = ".qs")
    #     datasets <- gsub(".qs", "", qs_files)
    #     raw_list <- raw_list[!names(raw_list) %in% datasets]
    # }
    raw_list
    # save all raw_count_matrix as tsv in data/203.inferCNV/raw_count_matrix
    # for (i in 1:length(all_datasets)) {
    #     write.table(raw_list[[i]], file = paste0("data/203.inferCNV/raw_count_matrix/", all_datasets[i], ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    # }
    # "data/203.inferCNV/raw_count_matrix"
}

infer_cnv <- function() {
    # # unlist raw_matrix and name_matrix
    # raw_matrix <- raw_matrix[[1]]
    # name_matrix <- name_matrix[[1]]
    # res <- pipelineCNA(raw_matrix, par_cores = 10)
    # qs::qsave(res, file = paste0("data/203.inferCNV/", name_matrix, ".qs"))
    # bash code/Shell/201.inferCNV.sh
    "data/203.inferCNV/result"
}

post_cnv <- function(cnv_results, sc) {
    # use read_tsv to read all tsv files in cnv_results
    cnv_results_list <- list.files(cnv_results, pattern = ".csv", full.names = TRUE) %>%
        map(readr::read_csv)
    # names are file names without extension
    names(cnv_results_list) <- list.files(cnv_results, pattern = ".csv") %>%
        str_remove(".csv")
    # bind_rows all dataframes in cnv_results_list, add new column "dataset" to each dataframe
    cnv_results_df <- dplyr::bind_rows(cnv_results_list) %>%
        # convert predict column, 0 to normal, 1 to tumor
        mutate(predict = if_else(predict == 0, "normal", "tumor")) %>%
        # remove the part before the first dot in sample
        mutate(sample = sub("^[^.]*\\.", "", sample)) %>%
        # remove the end .1 or -1 in sample
        mutate(sample = sub("\\.1$", "", sample)) %>%
        mutate(sample = sub("-1$", "", sample)) %>%
        distinct()
    sc <- sc %>%
        # create col cell_name based on .cell, remove the end .1 or -1
        mutate(cell_name = sub("\\.1$", "", .cell)) %>%
        mutate(cell_name = sub("-1$", "", cell_name))

    # left_join cnv_results_df with sc by "dataset"
    sc_cnv <- left_join(sc, cnv_results_df, by = c("cell_name" = "sample")) %>%
        # remove predict is NA
        filter(!is.na(predict))
    # p <- DimPlot(sc_cnv, reduction = "umap", group.by = "predict", label = TRUE)
    # ggsave("result/203.inferCNV/umap_cnv.png", p, width = 10, height = 10)

    sc_cnv_normal <- sc_cnv %>%
        dplyr::filter(predict == "normal")
    sc_anndata <- AnnData(
        X = sc_cnv_normal[["RNA"]]$counts %>% t(),
        obs = sc_cnv_normal[[]],
        obsm = list(
            pca = as.matrix(Embeddings(sc_cnv_normal, reduction = "pca"))
        )
    )
    write_h5ad(sc_anndata, "data/201.load_sc/sc_pre.h5ad")
    sc_cnv
}
