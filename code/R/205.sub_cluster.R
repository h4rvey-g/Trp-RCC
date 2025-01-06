get_sub_cluster <- function(sc_opt, res) {
    sc_T_cell <- sc_opt %>%
        filter(cell_type_opt == "T-cell") %>%
        FindClusters(
            resolution = res,
            method = "igraph",
            algorithm = 4,
            cluster.name = "T_cell_sub_cluster"
        )
    # # Run FindAllMarkers
    T_sub_markers <- FindAllMarkers(
        object = sc_T_cell,
        only.pos = TRUE,
        min.pct = 0.25,
        verbose = FALSE
    )
    # filter top 5 features of each sub cluster
    top_markers <- T_sub_markers %>%
        filter(p_val_adj < 0.05) %>%
        filter(cluster %in% cluster_TBD) %>%
        group_by(cluster) %>%
        top_n(10, avg_log2FC) %>%
        ungroup() %>%
        arrange(cluster, desc(avg_log2FC)) %>%
        pull(gene) %>%
        unique()
    # DotPlot
    p <- DotPlot(
        object = sc_T_cell,
        features = top_markers,
        assay = "RNA",
        dot.scale = 4,
        idents = cluster_TBD,
        group.by = "T_cell_sub_cluster"
    ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white"),
        ) +
        # order y axis by numeric
        scale_y_discrete(limits = cluster_TBD %>% as.numeric() %>% sort() %>% as.character()) +
        # draw a grey vertical line after every 10th tick
        geom_vline(xintercept = seq(10, length(cluster_TBD)*10, by = 10), color = "grey")

    ggsave(
        filename = paste0("result/205.sub_cluster/T_cell_sub_cluster_", res, ".png"),
        plot = p,
        width = 20,
        height = 10
    )
    top_markers_10 <- T_sub_markers %>%
        filter(p_val_adj < 0.05) %>%
        group_by(cluster) %>%
        top_n(20, avg_log2FC) %>%
        ungroup() %>%
        arrange(cluster, desc(avg_log2FC)) %>%
        select(cluster, gene) %>%
        filter(cluster %in% cluster_TBD) %>%
        # make each cluster is a row, and gene is arranged after cluster
        reshape2::dcast(cluster ~ data.table::rowid(cluster), value.var = "gene") 
    write_csv(top_markers_10, paste0("result/205.sub_cluster/T_cell_sub_cluster_", res, ".csv"))
    # paste0("result/205.sub_cluster/T_cell_sub_cluster_", res, ".png")
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

    sub_annotations <- read_csv("data/205.sub_cluster/predicted_labels_T.csv")
    sc_T_cell <- sc_T_cell %>%
        left_join(sub_annotations, by = c("cell" = "cell_name"))
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
    p <- DimPlot2(sc_T_cell,
        features = c("T_cell_sub_cluster", "cell_type_dtl"), label = TRUE,
        reduction = "umap_T_cell", pt.size = 1
    ) +
        # no legend
        theme(legend.position = "none", plot.background = element_rect(fill = "white"))
    ggsave("result/205.sub_cluster/T_cell_all.png", p, width = 30, height = 10)
    sc_T_cell_test <- sc_T_cell %>%
        filter(!is.na(cell_type_dtl))
    p <- ClusterDistrBar(
        origin = sc_T_cell_test$T_cell_sub_cluster %>% as.character() %>% paste0("_c"),
        cluster = sc_T_cell_test$cell_type_dtl
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("result/205.sub_cluster/T_cell_all_bar.png", p, width = 20, height = 10)

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
    # left_join(
    #     # calculate cell number for each T_cell_sub_cluster
    #     sc_T_cell_test %>%
    #         as_tibble() %>%
    #         count(T_cell_sub_cluster) %>%
    #         rename(n = n) %>%
    #         # add _c to T_cell_sub_cluster
    #         mutate(T_cell_sub_cluster = T_cell_sub_cluster %>% paste0("_c")),
    #     by = "T_cell_sub_cluster"
    # )
    write_tsv(cluster_anno_mtx, "data/205.sub_cluster/T_cell_sub_cluster_anno.tsv")
    # pull cluster with percentage < 0.8, remove _c
    cluster_TBD <- cluster_anno_mtx %>%
        filter(percentage < 80 | cell_type_dtl == "Low quality") %>%
        pull(T_cell_sub_cluster) %>%
        str_remove("_c")

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
        "PDZK1IP1", "LRP2", "ALDOB", # proximal tubule cell
        "SLC4A4", "SLC5A12", "SLC22A19", # kidney cell
        "CA9", "NDUFA4L2", "VEGFA" # RCC cell
    ) %>%
        genekitr::transId(transTo = "symbol", keepNA = FALSE) %>%
        pull(symbol) %>%
        unique()
    stats1 <- CalcStats(sc_T_cell,
        features = cellmatch_T_vector,
        method = "zscore", order = "p",
        cells = WhichCells(sc_T_cell, idents = cluster_TBD),
        group.by = "T_cell_sub_cluster"
    ) %>%
        rownames_to_column("gene") %>%
        mutate(gene = factor(gene, levels = cellmatch_T_vector)) %>%
        arrange(gene) %>%
        column_to_rownames("gene")
    p1 <- Heatmap(stats1,
        lab_fill = "zscore"
    )
    stats2 <- CalcStats(sc_T_cell,
        features = marker_genes,
        method = "zscore", order = "p",
        cells = WhichCells(sc_T_cell, idents = cluster_TBD),
        group.by = "T_cell_sub_cluster"
    ) %>%
        rownames_to_column("gene") %>%
        mutate(gene = factor(gene, levels = marker_genes)) %>%
        arrange(gene) %>%
        column_to_rownames("gene")
    p2 <- Heatmap(stats2, lab_fill = "zscore")
    p <- plot_grid(p1, p2, ncol = 1)
    ggsave("result/205.sub_cluster/T_cell_TBD_heatmap.png", p, width = 20, height = 20)
    p1 <- DotPlot(
        object = sc_T_cell,
        features = cellmatch_T_vector,
        assay = "RNA",
        dot.scale = 4,
        idents = cluster_TBD,
        group.by = "T_cell_sub_cluster"
    ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white"),
        ) +
        # order y axis by numeric
        scale_y_discrete(limits = cluster_TBD %>% as.numeric() %>% sort() %>% as.character())
    p2 <- DotPlot(
        object = sc_T_cell,
        features = marker_genes,
        assay = "RNA",
        dot.scale = 4,
        idents = cluster_TBD,
        group.by = "T_cell_sub_cluster"
    ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white"),
        ) +
        # order y axis by numeric
        scale_y_discrete(limits = cluster_TBD %>% as.numeric() %>% sort() %>% as.character()) +
        # draw a grey vertical line after every 10th tick
        geom_vline(xintercept = seq(10, length(cluster_TBD) * 10, by = 10), color = "grey")
    p <- plot_grid(p1, p2, ncol = 2)
    ggsave("result/205.sub_cluster/T_cell_TBD_dotplot.png", p, width = 20, height = 20)


    sc_T_cell_mtx <- rev_gene(sc_T_cell[["RNA"]]@data,
        data_type = "data",
        species = "Human",
        geneinfo = geneinfo
    )
    # scCATCH_T <- createscCATCH(data = sc_T_cell_mtx, cluster = sc_T_cell$T_cell_sub_cluster %>%
    #     as.character()) %>%
    #     findmarkergene(., species = "Human", marker = cellmatch, tissue = "Kidney")
    # scCATCH_T <- findcelltype(scCATCH_T)

    cellmatch_T <- list(
        "T_CD4" = "CD4",
        "T_CD8" = c("CD8A", "CD8B"),
        "Tex" = c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT"),
        "Tn_CD4" = c("CCR7", "PTPRC", "CD4"),
        "Tn_CD8" = c("PTPRC", "CD69", "CCR7", "CD8A"),
        "Trm_CD4" = c("CD69", "CD103", "CXCR6", "CD4"),
        "Trm_CD8" = c("CD69", "CD103", "CXCR6", "CD8A"),
        "Treg" = c("FOXP3"),
        "Tem" = c("GZMK", "NKG7", "CD8A", "IFNG")
    )
    # use transId to convert gene to Symbol
    cellmatch_T_trans <- lapply(cellmatch_T, function(x) {
        genekitr::transId(x, transTo = "symbol", keepNA = TRUE)$symbol
    })
    cellmatch_T_vector <- unlist(cellmatch_T_trans) %>% unique()
    # Run each element in cellmatch_T on sc_T_cell, use VlnPlot2(sc_T_cell, features = genes, ncol = 1), save to paste0("result/205.sub_cluster/", celltype, ".png")
    for (celltype in names(cellmatch_T)) {
        genes <- cellmatch_T_trans[[celltype]]
        p <- VlnPlot2(sc_T_cell, features = genes, ncol = 2, group.by = "T_cell_sub_cluster", pt = FALSE)
        ggsave(
            filename = paste0("result/205.sub_cluster/", celltype, ".png"),
            plot = p,
            width = 14,
            height = 7 * ceiling(length(genes) / 2)
        )
        p <- FeaturePlot3.grid(sc_T_cell, features = genes, color = "ryb", reduction = "umap_T_cell") +
            theme(plot.background = element_rect(fill = "white"))
        ggsave(
            filename = paste0("result/205.sub_cluster/", celltype, "_umap.png"),
            plot = p,
            width = 14,
            height = 7,
            device = grDevices::png
        )
    }
    # convert cellmatch_T to data.frame, add celltype column
    cellmatch_T_df <- cellmatch_T %>%
        enframe(name = "celltype", value = "gene") %>%
        unnest(cols = gene) %>%
        mutate(celltype = factor(celltype, levels = names(cellmatch_T))) %>%
        # add column species tissue cancer               condition   subtype1   subtype2 subtype3
        mutate(
            species = "Human", tissue = "Kidney", cancer = "Normal", condition = "Normal cell", subtype1 = NA,
            subtype2 = NA, subtype3 = NA, pmid = "098333"
        )
    cellmatch_T_df <- read_tsv("data/205.sub_cluster/T_markers.tsv") %>%
        mutate(
            species = "Human", tissue = "Kidney", cancer = "Normal", condition = "Normal cell", subtype1 = NA,
            subtype2 = NA, subtype3 = NA, pmid = "098333"
        ) %>%
        # move the celltype to subtype3, then change celltype to CD4_T or CD8_T, based on the beginning of celltype
        mutate(
            subtype3 = celltype,
            celltype = case_when(
                str_detect(celltype, "^CD4") ~ "CD4_T",
                str_detect(celltype, "^CD8") ~ "CD8_T",
                TRUE ~ celltype
            )
        )
    # add a row gene = "CD4" for all celltype begin with "CD4"
    cellmatch_T_df <- cellmatch_T_df %>%
        filter(celltype != "CD4_T") %>%
        bind_rows(cellmatch_T_df %>%
            filter(celltype == "CD4_T") %>%
            mutate(gene = "CD4"))
    # add a row gene = "CD8" for all celltype begin with "CD8"
    cellmatch_T_df <- cellmatch_T_df %>%
        filter(celltype != "CD8_T") %>%
        bind_rows(cellmatch_T_df %>%
            filter(celltype == "CD8_T") %>%
            mutate(gene = "CD8"))

    # use genekitr::transId to convert gene to Symbol
    transId_result <- genekitr::transId(cellmatch_T_df$gene, transTo = "symbol", keepNA = TRUE)
    cellmatch_T_df <- cellmatch_T_df %>%
        left_join(transId_result, by = c("gene" = "input_id")) %>%
        # replace gene with symbol
        mutate(gene = symbol)
    scCATCH_T <- createscCATCH(data = sc_T_cell_mtx, cluster = sc_T_cell$T_cell_sub_cluster %>%
        as.character()) %>%
        findmarkergene(., marker = cellmatch_T_df, use_method = "2", if_use_custom_marker = TRUE) %>%
        findcelltype()
    scCATCH_T@celltype %>%
        as_tibble() %>%
        print(n = Inf)
    cellmatch_T_db <- cellmatch %>%
        as_tibble() %>%
        filter(celltype == "T Cell") %>%
        filter(condition == "Normal cell") %>%
        mutate(tissue = "Kidney") %>%
        filter(subtype3 == "Exhausted")
}
