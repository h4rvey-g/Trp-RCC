get_sub_cluster <- function(sc_opt) {
    sc_T_sub_cluster <- sc_opt %>%
        filter(cell_type_opt == "T-cell") %>%
        RunUMAP(
            dims = 1:30,
            reduction = "scvi",
            n.neighbors = 30,
            reduction.name = "umap_T_cell"
        ) %>%
        FindNeighbors(
            dims = 1:30,
            reduction = "scvi",
            n.neighbors = 30
        ) %>%
        FindClusters(
            resolution = 0.4,
            method = "igraph",
            algorithm = 4,
            cluster.name = "T_cell_sub_cluster"
        )
    sc_T_sub_cluster
}

anno_Tex <- function(sc_T_sub_cluster, predicted_labels_T) {
    predicted_labels_T <- read_csv("data/205.sub_cluster/predicted_labels_T.csv")
    sc_T_cell <- sc_T_sub_cluster %>%
        left_join(predicted_labels_T, by = c("cell" = "cell_name"))
    sc_T_cell <- sc_T_cell %>%
        # mutate cell_type_dtl = majority_voting, if contain CD8+T_EFF, then change to CD8+T_EFF
        # if contain CD4+T_Act, then change to CD4+T_Act
        # if contain CD8+T_Cycling, then change to CD8+T_Cycling
        # if contain CD4+T_EX, then change to CD4+T_EX
        # if contain CD8+T_preEX, then change to CD8+T_preEX
        mutate(
            cell_type_dtl = case_when(
                str_detect(majority_voting, "CD8\\+T_EFF") ~ "CD8+T_EFF",
                str_detect(majority_voting, "CD4\\+T_Act") ~ "CD4+T_Act",
                str_detect(majority_voting, "CD8\\+T_Cycling") ~ "CD8+T_Cycling",
                str_detect(majority_voting, "CD8\\+T_EX") ~ "CD8+T_EX",
                str_detect(majority_voting, "CD8\\+T_preEX") ~ "CD8+T_EX",
                str_detect(majority_voting, "gdT") ~ "gdT",
                TRUE ~ majority_voting
            )
        )

    sc_T_cell_test <- sc_T_cell %>%
        filter(!is.na(cell_type_dtl))
    p <- ClusterDistrBar(
        origin = sc_T_cell_test$T_cell_sub_cluster %>% as.character() %>% paste0("_c"),
        cluster = sc_T_cell_test$cell_type_dtl
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("result/205.sub_cluster/T_cell_all_bar.png", p, width = 20, height = 10)
    p <- DimPlot2(
        sc_T_cell_test,
        cells.highlight = sc_T_cell_test %>% filter(cell_type_dtl == "CD8+T_EX") %>% pull(cell),
        label = TRUE,
        reduction = "umap_T_cell"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("result/205.sub_cluster/T_cell_CD8+T_EX.png", p, width = 10, height = 10)

    cluster_anno_mtx <- ClusterDistrBar(
        origin = sc_T_cell_test$T_cell_sub_cluster %>% as.character() %>% paste0("_c"),
        cluster = sc_T_cell_test$cell_type_dtl,
        plot = FALSE
    ) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("T_cell_sub_cluster") %>%
        pivot_longer(cols = -T_cell_sub_cluster, names_to = "cell_type_dtl", values_to = "percentage")

    # select the most frequent cell_type_dtl and value for each T_cell_sub_cluster
    cluster_anno_mtx <- cluster_anno_mtx %>%
        group_by(T_cell_sub_cluster) %>%
        filter(percentage == max(percentage)) %>%
        ungroup() %>%
        arrange(desc(percentage))
    write_tsv(cluster_anno_mtx, "data/205.sub_cluster/T_cell_sub_cluster_anno.tsv")
    # sc_T_cell <- sc_T_cell %>%
    #     # change cluster 12,7,16,2  to CD8+T_EX, rest to NA
    #     mutate(
    #         cell_type_dtl = case_when(
    #             T_cell_sub_cluster == 12 | T_cell_sub_cluster == 7 | T_cell_sub_cluster == 16 | T_cell_sub_cluster == 2 ~ "CD8+T_EX",
    #             TRUE ~ NA_character_
    #         )
    #     )
    sc_T_cell
}

get_module_score_T <- function(sc_T_anno) {
    cellmatch_T <- list(
        "T_CD4" = "CD4",
        "T_CD8" = c("CD8A", "CD8B"),
        "Tex" = c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT"),
        "Tn_CD4" = c("CCR7", "PTPRC", "CD4"),
        "Tn_CD8" = c("PTPRC", "CD69", "CCR7", "CD8A", "SELL", "LEF1", "TCF7"),
        "Trm_CD4" = c("CD69", "ITGAE", "CXCR6", "CD4"),
        "Trm_CD8" = c("CD69", "ITGAE", "CXCR6", "CD8A"),
        "Treg" = c("FOXP3", "CD4"),
        "Tem" = c("GZMK", "NKG7", "CD8A", "IFNG"),
        "Teff" = c("CX3CR1", "KLRF1", "FGFBP2", "FCGR3A"),
        "MAIT" = c("SLC4A10", "RORC")
    )
    module_score <- sc_T_anno %>%
        AddModuleScore(
            features = cellmatch_T,
            name = paste0("module_"),
            ctrl = 5,
            search = TRUE
        )
    meta <- module_score[[]]
    # rename module_1 to module_11 to paste0 module_ names(cellmatch_T)
    new_name <- paste0("module_", names(cellmatch_T))
    meta <- meta %>%
        rename_at(vars(starts_with("module_")), ~new_name)
    module_score[[]] <- meta
    toplot <- CalcStats(
        module_score,
        features = new_name,
        group.by = "T_cell_sub_cluster"
    )
    p <- Heatmap(
        toplot,
        lab_fill = "zscore"
    )
    ggsave("result/205.sub_cluster/module_score_T.png", p, width = 20, height = 10)
    toplot <- CalcStats(
        module_score,
        features = new_name,
        group.by = "dataset"
    )
    p <- Heatmap(
        toplot,
        lab_fill = "zscore"
    )
    ggsave("result/205.sub_cluster/module_score_T_dataset.png", p, width = 40, height = 10)
    Tex_score <- toplot %>%
        rownames_to_column("module") %>%
        pivot_longer(cols = -module, names_to = "dataset", values_to = "Tex_value") %>%
        filter(module == "module_Tex")
    Tex_score
}

Tex_in_RCC <- function(sc_T_anno, sc_opt) {
    cellmatch_T <- list(
        "Tex" = c("HAVCR2", "LAG3", "PDCD1", "TIGIT")
    )
    sc_T_anno <- sc_T_anno %>%
        filter(batch == "RCC")
    p <- DotPlot2(
        sc_T_anno,
        features = cellmatch_T,
        group.by = "dataset"
    )
    ggsave("result/206.Tex/Tex_in_RCC_dotplot.png", p, width = 10, height = 10)
    df_zscore <- CalcStats(
        sc_T_anno,
        features = cellmatch_T$Tex,
        method = "zscore",
        order = "p",
        group.by = "dataset"
    ) %>%
        as.data.frame() %>%
        rownames_to_column("Tex_genes") %>%
        pivot_longer(cols = -Tex_genes, names_to = "dataset", values_to = "zscore")
    df_proportions <- feature_percent(
        sc_T_anno,
        feature = cellmatch_T$Tex,
        group.by = "dataset"
    ) %>%
        as.data.frame() %>%
        rownames_to_column("Tex_genes") %>%
        pivot_longer(cols = -Tex_genes, names_to = "dataset", values_to = "proportion")
    df <- df_zscore %>%
        left_join(df_proportions, by = c("dataset", "Tex_genes"))
    # for each in Tex_genes, identify the dataset to positive or negative, based on if zsore and proportion > mean zscore and proportion in this Tex_genes
    df_2 <- df %>%
        group_by(Tex_genes) %>%
        mutate(
            zscore_mean = mean(zscore),
            proportion_mean = mean(proportion),
            zscore_positive = zscore > zscore_mean,
            proportion_positive = proportion > proportion_mean,
            positive = zscore_positive & proportion_positive
        ) %>%
        ungroup() %>%
        select(Tex_genes, dataset, positive)
    # filter positive, convert to list, each element is a Tex_genes, containing vector of dataset
    positive_Tex <- df_2 %>%
        filter(positive) %>%
        group_by(Tex_genes) %>%
        summarize(dataset = list(unique(dataset))) %>%
        deframe()
    p <- ggVennDiagram(
        positive_Tex
    )
    ggsave("result/206.Tex/Tex_in_RCC_venn.png", p, width = 10, height = 10)
    # intersect positive == FALSE in all Tex_genes
    negative_Tex <- df_2 %>%
        filter(!positive) %>%
        group_by(Tex_genes) %>%
        summarize(dataset = list(unique(dataset))) %>%
        deframe()
    p <- ggVennDiagram(
        negative_Tex
    )
    ggsave("result/206.Tex/Tex_in_RCC_negative_venn.png", p, width = 10, height = 10)
    # get the dataset in full intersection of positive_Tex and negative_Tex
    positive_dataset <- Reduce(intersect, positive_Tex)
    negative_dataset <- Reduce(intersect, negative_Tex)

    # filter positive_dataset and negative_dataset in sc_T_anno
    sc_opt_RCC <- sc_opt %>%
        filter(dataset %in% c(positive_dataset, negative_dataset)) %>%
        filter(cell_type_opt == "RCC") %>%
        # add a new column to indicate if the dataset is positive or negative
        mutate(
            Tex_status = case_when(
                dataset %in% positive_dataset ~ "positive",
                dataset %in% negative_dataset ~ "negative"
            )
        )
    sc_pseudo <- AggregateExpression(sc_opt_RCC, assays = "RNA", return.seurat = T, group.by = c("dataset", "Tex_status"))
    sc_pseudo$dataset_Tex <- paste(sc_pseudo$dataset, sc_pseudo$Tex_status, sep = "_")
    Idents(sc_pseudo) <- "Tex_status"
    DEG_RCC <- FindMarkers(sc_pseudo,
        ident.1 = "positive",
        ident.2 = "negative",
        test.use = "DESeq2",
        min.cells.group = 2
    ) %>%
        Add_Pct_Diff() %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        arrange(desc(abs(avg_log2FC)))
    write_tsv(DEG_RCC, "result/206.Tex/DEG_RCC.tsv")
}
