uni_enrich <- function(cox_res, data_dds, msigdb) {
    gene_up_combine <- data_dds %>%
        as_tibble() %>%
        dplyr::select(.feature, log2FoldChange) %>%
        filter(log2FoldChange > 0) %>%
        distinct() %>%
        inner_join(cox_res %>% filter(coef > 0),
            by = c(".feature" = "gene")
        ) %>%
        # mutate(ENTREZID = bitr(.feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)$ENTREZID) %>%
        mutate(ENTREZID = transId(.feature, transTo = "entrez", keepNA = TRUE, unique = TRUE)$entrezid) %>%
        filter(!is.na(ENTREZID)) %>%
        arrange(desc(log2FoldChange))
    gene_up <- gene_up_combine$log2FoldChange %>%
        setNames(gene_up_combine$ENTREZID)

    gene_down_combine <- data_dds %>%
        as_tibble() %>%
        dplyr::select(.feature, log2FoldChange) %>%
        filter(log2FoldChange < 0) %>%
        distinct() %>%
        inner_join(cox_res %>% filter(coef < 0),
            by = c(".feature" = "gene")
        ) %>%
        # mutate(ENTREZID = bitr(.feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)$ENTREZID) %>%
        mutate(ENTREZID = transId(.feature, transTo = "entrez", keepNA = TRUE, unique = TRUE)$entrezid) %>%
        filter(!is.na(ENTREZID)) %>%
        mutate(log2FoldChange = abs(log2FoldChange)) %>%
        arrange(desc(log2FoldChange))
    gene_down <- gene_down_combine$log2FoldChange %>%
        setNames(gene_down_combine$ENTREZID)

    gs <- geneset::getMsigdb(org = "human")
    gs$geneset <- msigdb
    gsea_up <- genGSEA(gene_up, gs)
    p <- plotGSEA(gsea_up, plot_type = "ridge", stats_metric = "p.adjust", show_pathway = nrow(gsea_up$gsea_df))
    ggsave("result/108.direct_enrich/gsea_up.png", p)
    p <- plotGSEA(gsea_up, plot_type = "bar", colour = c("navyblue", "orange"))
    ggsave("result/108.direct_enrich/gsea_up_bar.png", p)
    write_tsv(gsea_up$gsea_df, "data/108.direct_enrich/gsea_up.tsv")

    gsea_down <- tryCatch(
        genGSEA(gene_down, gs),
        error = function(e) {
            message("Error in genGSEA for gene_down: ", e$message)
            return(NULL)
        }
    )

    if (!is.null(gsea_down)) {
        p <- plotGSEA(gsea_down, plot_type = "ridge", stats_metric = "p.adjust", show_pathway = nrow(gsea_down$gsea_df))
        ggsave("result/108.direct_enrich/gsea_down.png", p)
        p <- plotGSEA(gsea_down, plot_type = "bar", colour = c("navyblue", "orange"))
        ggsave("result/108.direct_enrich/gsea_down_bar.png", p)
        write_tsv(gsea_down$gsea_df, "data/108.direct_enrich/gsea_down.tsv")
    }

    gsea_all <- genGSEA(c(gene_up, gene_down %>%
        # revert gene_down, and make it negative
        rev() %>% `*`(-1)), gs)
    p <- plotGSEA(gsea_all, plot_type = "ridge", stats_metric = "p.adjust", show_pathway = nrow(gsea_all$gsea_df))
    ggsave("result/108.direct_enrich/gsea_all.png", p)
    p <- plotGSEA(gsea_all, plot_type = "bar", colour = c("navyblue", "orange"))
    ggsave("result/108.direct_enrich/gsea_all_bar.png", p)
    write_tsv(gsea_all$gsea_df, "data/108.direct_enrich/gsea_all.tsv")

    list(gsea_up = gsea_up, gsea_down = gsea_down, gsea_all = gsea_all)
}

plot_uni_enrich <- function(uni_EA_res, data_dds, cox_res) {
    gsea_up <- uni_EA_res$gsea_up$gsea_df
    gsea_down <- uni_EA_res$gsea_down$gsea_df
    gsea_all <- uni_EA_res$gsea_all$gsea_df
    data_dds <- data_dds %>%
        as_tibble() %>%
        select(.feature, log2FoldChange, padj) %>%
        distinct() %>%
        inner_join(cox_res, by = c(".feature" = "gene"))

    gsea_up <- gsea_up %>%
        mutate(gene_name = map(geneID, ~ {
            ids <- str_split(.x, "/")[[1]] # Split by "/" and keep as a list
            transId(ids, transTo = "symbol", keepNA = TRUE, unique = TRUE) %>%
                pull(symbol)
        })) %>%
        as_tibble()
    gsea_all <- gsea_all %>%
        mutate(gene_name = map(geneID, ~ {
            ids <- str_split(.x, "/")[[1]] # Split by "/" and keep as a list
            transId(ids, transTo = "symbol", keepNA = TRUE, unique = TRUE) %>%
                pull(symbol)
        })) %>%
        as_tibble()

    # traverse each row of gsea_up, and get the gene_name, use data in data_dds, draw volcano plot. Use purrr
    draw_volcano <- function(gsea_df, prefix) {
        pathway <- gsea_df$ID %>% tolower()
        genes <- gsea_df$gene_name %>% unlist()
        # filter genes in data_dds
        dds_res <- data_dds %>%
            dplyr::filter(.feature %in% genes)
        p <- dds_res %>%
            ggplot(aes(
                x = log2FoldChange, y = coef, color = log2FoldChange > 0
            )) +
            geom_point(size = 3) +
            geom_label_repel(aes(label = .feature), box.padding = 0.5, max.overlaps = 15) +
            scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            theme_minimal() +
            labs(
                title = pathway,
                x = "log2FoldChange", y = "coef"
            ) +
            theme_classic() +
            theme(panel.background = element_rect(fill = "white")) +
            theme(legend.position = "none")
        ggsave(paste0("result/108.direct_enrich/", prefix, "/volcano_", pathway, ".png"), p,
            width = 10, height = 10
        )
    }
    gsea_up %>%
        rowwise() %>%
        group_map(~ draw_volcano(.x, "up"))
    # gsea_down %>%
    #     rowwise() %>%
    #     group_map(~ draw_volcano(.x, "down"))
    gsea_all %>%
        rowwise() %>%
        group_map(~ draw_volcano(.x, "all"))

    "result/108.direct_enrich"
}

select_intersect_genes <- function(uni_EA_res, data_dds, cox_res) {
    gsea_up <- uni_EA_res$gsea_up$gsea_df
    gsea_down <- uni_EA_res$gsea_down$gsea_df
    gsea_all <- uni_EA_res$gsea_all$gsea_df
    data_dds <- data_dds %>%
        as_tibble() %>%
        select(.feature, log2FoldChange, padj) %>%
        distinct() %>%
        inner_join(cox_res, by = c(".feature" = "gene")) %>%
        select(.feature, log2FoldChange, coef)

    gsea_up <- gsea_up %>%
        mutate(gene_name = map(geneID, ~ {
            ids <- str_split(.x, "/")[[1]] # Split by "/" and keep as a list
            transId(ids, transTo = "symbol", keepNA = TRUE, unique = TRUE) %>%
                pull(symbol)
        })) %>%
        as_tibble()
    gsea_all <- gsea_all %>%
        mutate(gene_name = map(geneID, ~ {
            ids <- str_split(.x, "/")[[1]] # Split by "/" and keep as a list
            transId(ids, transTo = "symbol", keepNA = TRUE, unique = TRUE) %>%
                pull(symbol)
        })) %>%
        as_tibble()
    gsea_all_intersect <- gsea_all %>%
        unnest_longer(gene_name) %>%
        group_by(gene_name) %>%
        summarise(
            pathways = list(ID)
        ) %>%
        # make everything in pathways tolower
        mutate(pathways = map(pathways, ~ .x %>% tolower())) %>%
        mutate(
            pathway_count = lengths(pathways)
        ) %>%
        filter(pathway_count > 1) %>%
        arrange(desc(pathway_count))
    write_tsv(
        gsea_all_intersect %>%
            unnest_wider(pathways, names_sep = "_") %>%
            relocate(pathway_count, .after = gene_name) %>%
            left_join(data_dds, by = c("gene_name" = ".feature")) %>%
            relocate(log2FoldChange, coef, .after = gene_name),
        "result/108.direct_enrich/gene_pathway_count.tsv"
    )
    genes_to_plot <- gsea_all_intersect %>%
        unnest_longer(pathways) %>%
        group_by(pathways) %>%
        summarise(
            genes = list(gene_name)
        )
    write_tsv(genes_to_plot %>% unnest_wider(genes, names_sep = "_"), "result/108.direct_enrich/pathway_genes.tsv")
    # draw upset plot on gsea_all
    p <- ggplot(genes_to_plot, aes(genes)) +
        geom_bar() +
        scale_x_upset()
    ggsave("result/108.direct_enrich/upset.png", p, height = 20)

    "result/108.direct_enrich"
}
