run_risk_DEG <- function(data_filt, data_risk_score) {
    data_dds <- data_filt %>%
        left_join(data_risk_score, by = c(".sample" = ".sample")) %>%
        # add risk column, if risk_score is above median, set as high, otherwise low
        mutate(risk = ifelse(risk_score > median(risk_score, na.rm = TRUE), "high", "low")) %>%
        test_differential_abundance(
            .formula = ~ 0 + risk,
            .abundance = counts_scaled_adjusted,
            scaling_method = "none",
            contrasts = c("riskhigh - risklow"),
            action = "add",
            method = "limma_voom"
        )
    data_dds <- data_dds %>%
        as_tibble() %>%
        dplyr::mutate(log2FoldChange = `logFC___riskhigh - risklow`) %>%
        dplyr::mutate(pvalue = `P.Value___riskhigh - risklow`) %>%
        dplyr::mutate(padj = `adj.P.Val___riskhigh - risklow`) %>%
        dplyr::select(.sample, .feature, counts_scaled_adjusted, risk, risk_score, log2FoldChange, pvalue, padj) %>%
        filter(padj < 0.05 & abs(log2FoldChange) > 1)
    write_tsv(
        data_dds %>%
            dplyr::select(.feature, log2FoldChange, pvalue, padj) %>% distinct(),
        "data/106.risk_score_DEG/data_dds.tsv"
    )
    data_dds
}

run_risk_enrich <- function(data_dds) {
    data_dds <- data_dds %>%
        as_tibble() %>%
        dplyr::select(-.sample, -counts_scaled_adjusted, -risk, -risk_score) %>%
        distinct()
    entrez <- bitr(data_dds$.feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # inner join
    data_dds <- inner_join(data_dds, entrez, by = c(".feature" = "SYMBOL"))
    gene_list <- list(
        up = data_dds$ENTREZID[data_dds$log2FoldChange > 0],
        down = data_dds$ENTREZID[data_dds$log2FoldChange < 0]
    )
    gene_all <- data_dds$ENTREZID
    gene_gsea <- data_dds %>%
        pull(log2FoldChange) %>%
        setNames(data_dds$ENTREZID) %>%
        sort(decreasing = TRUE)
    # unique gene_gsea based on names
    gene_gsea <- gene_gsea[!duplicated(names(gene_gsea))]

    kk_up <- enrichKEGG(
        gene = gene_list$up,
        organism = "hsa",
        qvalueCutoff = 1
        # pvalueCutoff = 1
    )
    # if category == "Metabolism" is not NULL, filter it
    # kk_up <- kk_up %>% filter(category == "Metabolism")
    p <- dotplot(kk_up)
    ggsave("result/107.risk_score_enrich/KEGG_up.png", p)
    kk_down <- enrichKEGG(
        gene = gene_list$down,
        organism = "hsa",
        qvalueCutoff = 1
        # pvalueCutoff = 1
    )
    # kk_down <- kk_down %>% filter(category == "Metabolism")
    p <- dotplot(kk_down)
    ggsave("result/107.risk_score_enrich/KEGG_down.png", p)
    kk_gse <- gseKEGG(
        geneList = gene_gsea,
        nPermSimple = 10000,
        organism = "hsa",
        # by = "DOSE",
        # nPerm = 1000
        # scoreType = "pos",
    )
    # search Tryptophan  in kk_gse$Description
    # kk_gse <- kk_gse %>% filter(grepl(".*Tryptophan.*", Description))
    p <- dotplot(kk_gse)
    ggsave("result/107.risk_score_enrich/KEGG_gse.png", p)
    kk_up <- setReadable(kk_up, org.Hs.eg.db, keyType = "ENTREZID")
    kk_down <- setReadable(kk_down, org.Hs.eg.db, keyType = "ENTREZID")
    kk_gse <- setReadable(kk_gse, org.Hs.eg.db, keyType = "ENTREZID")
    write_tsv(kk_up %>% as.data.frame(), "data/107.risk_score_enrich/KEGG_up.tsv")
    write_tsv(kk_down %>% as.data.frame(), "data/107.risk_score_enrich/KEGG_down.tsv")
    write_tsv(kk_gse %>% as.data.frame(), "data/107.risk_score_enrich/KEGG_gse.tsv")

    p <- plotEnrichAdv(
        kk_up %>% as.data.frame() %>% filter(category == "Metabolism") %>% rename(FoldEnrich = FoldEnrichment),
        kk_down %>% as.data.frame() %>% filter(category == "Metabolism") %>% rename(FoldEnrich = FoldEnrichment),
        plot_type = "one",
        term_metric = "FoldEnrich",
        stats_metric = "p.adjust"
    )
    ggsave("result/107.risk_score_enrich/metabolism_enrich.png", p, width = 10)
    trp_pathway <- pathview(
        gene.data = gene_gsea, same.layer = FALSE, limit = list(gene = 4, cpd = 1),
        pathway.id = "hsa00380", species = "hsa", kegg.dir = "result/107.risk_score_enrich/"
    )
    # move hsa00380.pathview.png to result/107.risk_score_enrich/trp_pathway.png
    file.rename("./hsa00380.pathview.png", "result/107.risk_score_enrich/trp_pathway.png")
    metabolism_gsea <- gseaplot2(kk_gse, geneSetID = 19, title = kk_gse$Description[19])
    ggsave("result/107.risk_score_enrich/metabolism_gsea.png", metabolism_gsea)
    EA_trp <- kk_down %>% filter(grepl(".*ryptophan.*", Description))
    # trp_heatmap <- heatplot(EA_trp, foldChange = kk_gse@geneList[kk_gse@geneList < 0])
    # ggsave("result/107.risk_score_enrich/trp_heatmap.png", trp_heatmap, height = 3)
    # select try_pathway$plot.data.gene %>% as_tibble() %>% select(mol.data no NA)
    ansEA <- list(
        kk_up = kk_up,
        kk_down = kk_down,
        kk_gse = kk_gse,
        trp_pathway = trp_pathway
    )
}
