run_risk_DEG <- function(data_filt, data_risk_score) {
    data_dds <- data_filt %>%
        left_join(data_risk_score, by = c(".sample" = ".sample")) %>%
        # add risk column, if risk_score is above median, set as high, otherwise low
        mutate(risk = ifelse(risk_score > median(risk_score, na.rm = TRUE), "high", "low"))
    assay(data_dds, "counts_scaled_adjusted") <- round(assay(data_dds, "counts_scaled_adjusted"))
    data_dds <- data_dds %>%
        test_differential_abundance(
            .formula = ~ 0 + risk,
            .abundance = counts_scaled_adjusted,
            scaling_method = "none",
            contrasts = list(c("risk", "high", "low")),
            action = "add",
            method = "DESeq2"
        )
    data_dds <- data_dds %>%
        as_tibble() %>%
        dplyr::mutate(log2FoldChange = `log2FoldChange___risk high-low`) %>%
        dplyr::mutate(pvalue = `pvalue___risk high-low`) %>%
        dplyr::mutate(padj = `padj___risk high-low`) %>%
        dplyr::select(.sample, .feature, counts_scaled_adjusted, risk, risk_score, log2FoldChange, pvalue, padj) %>%
        filter(padj < 0.05 & abs(log2FoldChange) > 1)
    write_tsv(
        data_dds %>%
            dplyr::select(.feature, log2FoldChange, pvalue, padj) %>% distinct(),
        "data/106.risk_score_DEG/data_dds.tsv"
    )
    data_dds
}

run_risk_enrich <- function(data_dds, msigdb) {
    data_dds <- data_dds %>%
        as_tibble() %>%
        dplyr::select(-.sample, -counts_scaled_adjusted, -risk, -risk_score) %>%
        distinct()
    entrez <- bitr(data_dds$.feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # inner join
    data_dds <- inner_join(data_dds, entrez, by = c(".feature" = "SYMBOL"))
    gene_gsea <- data_dds %>%
        pull(log2FoldChange) %>%
        setNames(data_dds$ENTREZID) %>%
        sort(decreasing = TRUE)

    gs <- geneset::getMsigdb(org = "human")
    gs$geneset <- msigdb
    gsea_metabolism <- genGSEA(gene_gsea, gs)
    p <- plotGSEA(gsea_metabolism, plot_type = "ridge", stats_metric = "p.adjust", show_pathway = nrow(gsea_metabolism$gsea_df))
    ggsave("result/107.risk_score_enrich/gsea_metabolism.png", p)
    p <- plotGSEA(gsea_metabolism, plot_type = "bar", colour = c("navyblue", "orange"))
    ggsave("result/107.risk_score_enrich/gsea_metabolism_bar.png", p)
    write_tsv(gsea_metabolism$gsea_df, "data/107.risk_score_enrich/gsea_metabolism.tsv")

    gsea_metabolism
}
