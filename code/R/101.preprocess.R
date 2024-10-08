# library(tidyverse)
# library(TCGAbiolinks)
# library(SummarizedExperiment)
# library(tidySummarizedExperiment)
download_data <- function() {
    load("data/101.raw_data/KIRC.rda")
    data
}
download_gtex <- function() {
    # download https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_kidney_cortex.gct.gz to data/101.raw_data
    download.file(
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/counts-by-tissue/gene_reads_2017-06-05_v8_kidney_cortex.gct.gz",
        "data/101.raw_data/gene_tpm_2017-06-05_v8_kidney_cortex.gct.gz"
    )
    # download https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
    download.file(
        "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
        "data/101.raw_data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    )
    c("data/101.raw_data/gene_tpm_2017-06-05_v8_kidney_cortex.gct.gz", "data/101.raw_data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
}
download_msigdb <- function() {
    msigdb <- read.gmt("data/101.raw_data/msigdb.v2024.1.Hs.entrez.gmt")
    # filter "METABOLI*" in term
    msigdb <- msigdb %>%
        filter(grepl("METABOLI", term)) %>%
        rename(gs_name = term, entrez_gene = gene)
}

preprocess_data <- function(data, gtex_data_path) {
    gtex_data <- read_tsv(gtex_data_path[1], skip = 2) %>%
        dplyr::select(-id, -Name)
    # only keep the first two parts of colnames sep by -
    colnames(gtex_data)[2:ncol(gtex_data)] <- colnames(gtex_data)[2:ncol(gtex_data)] %>%
        str_split("-", simplify = TRUE) %>%
        `[`(, 1:2) %>%
        apply(., 1, paste, collapse = "-")
    # combine rows with the same Description, get the mean value
    gtex_data <- gtex_data %>%
        group_by(Description) %>%
        summarise(across(everything(), mean)) %>%
        ungroup()
    gtex_data <- gtex_data %>%
        column_to_rownames("Description")
    gtex_meta <- read_tsv(gtex_data_path[2])
    # filter only colnames in gtex_data
    gtex_meta <- gtex_meta %>%
        filter(SUBJID %in% colnames(gtex_data)) %>%
        column_to_rownames("SUBJID")
    # construct summarizedExperiment
    gtex_data <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(gtex_data)),
        colData = DataFrame(gtex_meta),
        rowData = data.frame(gene_name = rownames(gtex_data))
    )
    # get the assay data
    data_mRNA <- data %>%
        filter(gene_type == "protein_coding") %>%
        aggregate_duplicates(.transcript = gene_name)
    data_mRNA <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = assay(data_mRNA)),
        rowData = rowData(data_mRNA)
    )
    data_mRNA <- data_mRNA %>%
        select(.feature, .sample, counts)
    gtex_data <- gtex_data %>%
        select(.feature, .sample, counts)
    data_mRNA_rownames <- transId(rownames(assay(data_mRNA)), transTo = "symbol", unique = TRUE)
    # filter rows in data_mRNA with rownames in data_mRNA_rownames$input_id, then replace rownames with data_mRNA_rownames$symbol
    data_mRNA <- data_mRNA %>%
        dplyr::filter(.feature %in% data_mRNA_rownames$input_id)
    data_mRNA_assay <- assay(data_mRNA)
    rownames(data_mRNA_assay) <- data_mRNA_rownames$symbol
    gtex_rownames <- transId(rownames(assay(gtex_data)), transTo = "symbol", unique = TRUE)
    gtex_data <- gtex_data %>%
        dplyr::filter(.feature %in% gtex_rownames$input)
    gtex_data_assay <- assay(gtex_data)
    rownames(gtex_data_assay) <- gtex_rownames$symbol
        gene_intersect <- intersect(rownames(data_mRNA_assay), rownames(gtex_data_assay))
    data_mRNA_assay <- data_mRNA_assay[gene_intersect, ]
    gtex_data_assay <- gtex_data_assay[gene_intersect, ]
    data_combined <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = cbind(data_mRNA_assay, gtex_data_assay))
    )
    browser()
    data_combined <- data_combined %>%
        mutate(
            group = ifelse(startsWith(.sample, "TCGA") & as.numeric(str_sub(.sample, 14, 15)) < 10, "tumor", "normal"),
            batch = ifelse(startsWith(.sample, "TCGA"), "TCGA", "GTEX")
        )
    data_combined_filt <- data_combined %>%
        keep_abundant(factor_of_interest = group) %>%
        identify_abundant(factor_of_interest = group) %>%
        scale_abundance() %>%
        adjust_abundance(.factor_unwanted = batch, .factor_of_interest = group, .abundance = counts_scaled, method = "limma_remove_batch_effect")

    # normalize
    p <- ggplot(
        data_combined_filt, aes(x = .sample, y = counts_scaled_adjusted)
    ) +
        geom_boxplot() +
        theme_minimal() +
        theme(axis.text.x = element_blank())
    ggsave("result/101.preprocess/boxplot.png", p, width = 15, height = 15)
    data_old_pca <- data_combined_filt %>%
        reduce_dimensions(method = "PCA", .dims = 3, .abundance = counts_scaled)
    data_pca <- data_combined_filt %>% reduce_dimensions(method = "PCA", .dims = 3, .abundance = counts_scaled_adjusted)
    p1 <- data_old_pca %>%
        pivot_sample() %>%
        ggplot(aes(x = `PC1`, y = `PC2`, color = group, shape = batch)) +
        geom_point()
    p2 <- data_pca %>%
        pivot_sample() %>%
        ggplot(aes(x = `PC1`, y = `PC2`, color = group, shape = batch)) +
        geom_point()
    p <- p1 + p2
    ggsave("result/101.preprocess/pca.png", p, width = 14)
    data_combined_filt
}
