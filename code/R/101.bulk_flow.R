# library(tidyverse)
# library(TCGAbiolinks)
# library(SummarizedExperiment)
# library(tidySummarizedExperiment)
download_data <- function() {
    load("data/101.raw_data/KIRC.rda")
    data
}
preprocess_data <- function(data) {
    # get the assay data
    data_mRNA <- data[rowData(data)$gene_type == "protein_coding", ]
    outliers <- TCGAanalyze_Preprocessing(data_mRNA, cor.cut = 0.6, filename = "result/101.preprocess/preprocess.png")
    # turns out no outliers

    # normalize
    data_mRNA_norm <- TCGAanalyze_Normalization(tabDF = data_mRNA, geneInfo = geneInfoHT)

    # filter
    data_filt <- TCGAanalyze_Filtering(tabDF = data_mRNA_norm, method = "quantile", qnt.cut = 0.25)
    TCGAbatch_Correction <- function(tabDF, batch.factor = NULL, adjustment = NULL, ClinicalDF = data.frame(),
                                     UnpublishedData = FALSE, AnnotationDF = data.frame()) {
        if (UnpublishedData == TRUE) {
            if (!"Batch" %in% colnames(AnnotationDF)) {
                stop("AnnotationDF should have a Batch column")
            } else {
                batch.factor <- as.factor(AnnotationDF$Batch)
            }
            if (!"Condition" %in% colnames(AnnotationDF)) {
                mod <- model.matrix(~ as.factor(Condition), data = AnnotationDF)
            } else {
                mod <- NULL
            }
            batch_corr <- sva::ComBat(
                dat = tabDF, batch = batch.factor,
                mod = mod, par.prior = TRUE, prior.plots = TRUE
            )
        }
        if (UnpublishedData == FALSE) {
            if (length(batch.factor) == 0 & length(adjustment) ==
                0) {
                message("batch correction will be skipped")
            } else if (batch.factor %in% adjustment) {
                stop(paste0("Cannot adjust and correct for the same factor|"))
            }
            my_IDs <- get_IDs(tabDF)
            if (length(batch.factor) > 0 || length(adjustment) >
                0) {
                if ((nrow(ClinicalDF) > 0 & batch.factor == "Year") ||
                    ("Year" %in% adjustment == TRUE & nrow(ClinicalDF) >
                        0)) {
                    names(ClinicalDF)[names(ClinicalDF) == "bcr_patient_barcode"] <- "patient"
                    ClinicalDF$age_at_diag_year <- floor(ClinicalDF$age_at_diagnosis / 365)
                    ClinicalDF$diag_year <- ClinicalDF$age_at_diag_year +
                        ClinicalDF$year_of_birth
                    diag_yearDF <- ClinicalDF[, c("patient", "diag_year")]
                    Year <- merge(my_IDs, diag_yearDF, by = "patient")
                    Year <- Year$diag_year
                    Year <- as.factor(Year)
                } else if (nrow(ClinicalDF) == 0 & batch.factor ==
                    "Year") {
                    stop("Cannot extract Year data. Clinical data was not provided")
                }
            }
            Plate <- as.factor(my_IDs$plate)
            Condition <- as.factor(my_IDs$condition)
            TSS <- as.factor(my_IDs$tss)
            Portion <- as.factor(my_IDs$portion)
            Sequencing.Center <- as.factor(my_IDs$center)
            design.mod.combat <- model.matrix(~Condition)
            options <- c("Plate", "TSS", "Year", "Portion", "Sequencing Center")
            if (length(batch.factor) > 1) {
                stop("Combat can only correct for one batch variable. Provide one batch factor")
            }
            if (batch.factor %in% options == FALSE) {
                stop(paste0(o, " is not a valid batch correction factor"))
            }
            for (o in adjustment) {
                if (o %in% options == FALSE) {
                    stop(paste0(o, " is not a valid adjustment factor"))
                }
            }
            adjustment.data <- c()
            for (a in adjustment) {
                if (a == "Sequencing Center") {
                    a <- Sequencing.Center
                }
                adjustment.data <- cbind(eval(parse(text = a)), adjustment.data)
            }
            if (batch.factor == "Sequencing Center") {
                batch.factor <- Sequencing.Center
            }
            batchCombat <- eval(parse(text = batch.factor))
            if (length(adjustment) > 0) {
                adjustment.formula <- paste(adjustment, collapse = "+")
                adjustment.formula <- paste0("+", adjustment.formula)
                adjustment.formula <- paste0("~Condition", adjustment.formula)
                print(adjustment.formula)
                model <- data.frame(batchCombat, row.names = colnames(tabDF))
                design.mod.combat <- model.matrix(eval(parse(text = adjustment.formula)),
                    data = model
                )
            }
            print(unique(batchCombat))
            batch_corr <- sva::ComBat(
                dat = tabDF, batch = batchCombat,
                mod = design.mod.combat, par.prior = TRUE, prior.plots = FALSE
            )
        }
        return(batch_corr)
    }
    # data_filt <- TCGAbatch_Correction(data_filt, batch.factor = "Plate")
    data_filt
}
run_DEG <- function(data_filt) {
    data_mRNA_tumor_index <- as.numeric(substr(colnames(data_filt), 14, 15)) < 10
    data_mRNA_tumor <- data_filt[, data_mRNA_tumor_index]
    data_mRNA_normal <- data_filt[, !data_mRNA_tumor_index]

    DEGs <- TCGAanalyze_DEA(
        mat1 = data_mRNA_normal,
        mat2 = data_mRNA_tumor,
        Cond1 = "Normal",
        Cond2 = "Tumor",
        method = "glmLRT",
        pipeline = "edgeR",
        batch.factors = "Plate"
    )
    DEGs <- DEGs %>% filter(FDR < 0.05 & abs(logFC) > 1)
    # save the DEGs table to data/102.DEG
    write_tsv(DEGs, "data/102.DEG/DEGs.tsv")
    # DEGs table with expression values in normal and tumor samples
    DEGs_filt_level <- TCGAanalyze_LevelTab(
        FC_FDR_table_mRNA = DEGs,
        typeCond1 = "Normal",
        typeCond2 = "Tumor",
        TableCond1 = data_mRNA_normal,
        TableCond2 = data_mRNA_tumor
    )
    write_tsv(DEGs_filt_level, "data/102.DEG/DEGs_filt_level.tsv")
    DEG_res <- list(
        DEGs = DEGs,
        DEGs_filt_level = DEGs_filt_level
    )
}

run_enrich <- function(DEG_res) {
    DEGs_filt_level <- DEG_res$DEGs_filt_level
    entrez <- bitr(DEGs_filt_level$mRNA, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # inner join
    DEGs_filt_level <- inner_join(DEGs_filt_level, entrez, by = c("mRNA" = "ENSEMBL"))
    gene_list <- list(
        up = DEGs_filt_level$ENTREZID[DEGs_filt_level$logFC > 0],
        down = DEGs_filt_level$ENTREZID[DEGs_filt_level$logFC < 0]
    )
    gene_all <- DEGs_filt_level$ENTREZID
    gene_gsea <- DEGs_filt_level %>%
        pull(logFC) %>%
        setNames(DEGs_filt_level$ENTREZID) %>%
        sort(decreasing = TRUE)
    # unique gene_gsea based on names
    gene_gsea <- gene_gsea[!duplicated(names(gene_gsea))]

    kk_up <- enrichKEGG(
        gene = gene_list$up,
        organism = "hsa",
        qvalueCutoff = 1,
        # pvalueCutoff = 1,
        universe = gene_all
    )
    # if category == "Metabolism" is not NULL, filter it
    # kk_up <- kk_up %>% filter(category == "Metabolism")
    p <- dotplot(kk_up)
    ggsave("result/103.enrich/KEGG_up.png", p)
    kk_down <- enrichKEGG(
        gene = gene_list$down,
        organism = "hsa",
        qvalueCutoff = 1,
        # pvalueCutoff = 1,
        universe = gene_all
    )
    # kk_down <- kk_down %>% filter(category == "Metabolism")
    p <- dotplot(kk_down)
    ggsave("result/103.enrich/KEGG_down.png", p)
    kk_gse <- gseKEGG(
        geneList = gene_gsea,
        organism = "hsa",
        by = "DOSE",
        nPerm = 1000
        # scoreType = "pos",
        # nPermSimple = 10000
    )
    # search Tryptophan  in kk_gse$Description
    # kk_gse <- kk_gse %>% filter(grepl(".*Tryptophan.*", Description))
    p <- dotplot(kk_gse)
    ggsave("result/103.enrich/KEGG_gse.png", p)
    write_tsv(kk_up %>% as.data.frame(), "data/103.enrich/KEGG_up.tsv")
    write_tsv(kk_down %>% as.data.frame(), "data/103.enrich/KEGG_down.tsv")
    write_tsv(kk_gse %>% as.data.frame(), "data/103.enrich/KEGG_gse.tsv")

    gene_gsea <- DEGs_filt_level %>%
        pull(logFC) %>%
        setNames(DEGs_filt_level$ENTREZID) %>%
        sort(decreasing = TRUE)
    # unique gene_gsea based on names
    gene_gsea <- gene_gsea[!duplicated(names(gene_gsea))]
    gene_list <- list(
        up = gene_gsea[gene_gsea > 0],
        down = gene_gsea[gene_gsea < 0]
    )
    kk_up <- setReadable(kk_up, org.Hs.eg.db, keyType = "ENTREZID")
    kk_down <- setReadable(kk_down, org.Hs.eg.db, keyType = "ENTREZID")
    kk_gse <- setReadable(kk_gse, org.Hs.eg.db, keyType = "ENTREZID")
    trp_pathway <- pathview(
        gene.data = c(gene_list$up, gene_list$down), same.layer = FALSE, limit = list(gene = 4, cpd = 1),
        pathway.id = "hsa00380", species = "hsa", kegg.dir = "result/103.enrich/"
    )
    # move hsa00380.pathview.png to result/103.enrich/trp_pathway.png
    file.rename("./hsa00380.pathview.png", "result/103.enrich/trp_pathway.png")
    metabolism_gsea <- gseaplot2(kk_gse, geneSetID = 19, title = kk_gse$Description[19])
    ggsave("result/103.enrich/metabolism_gsea.png", metabolism_gsea)
    EA_trp <- kk_down %>% filter(grepl(".*ryptophan.*", Description))
    # trp_heatmap <- heatplot(EA_trp, foldChange = kk_gse@geneList[kk_gse@geneList < 0])
    # ggsave("result/103.enrich/trp_heatmap.png", trp_heatmap, height = 3)
    ansEA <- list(
        kk_up = kk_up,
        kk_down = kk_down,
        kk_gse = kk_gse
    )
}
