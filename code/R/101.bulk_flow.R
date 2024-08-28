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
    # select try_pathway$plot.data.gene %>% as_tibble() %>% select(mol.data no NA)
    ansEA <- list(
        kk_up = kk_up,
        kk_down = kk_down,
        kk_gse = kk_gse,
        trp_pathway = trp_pathway
    )
}

run_WGCNA <- function(data_filt) {
    # mdata <- data_filt %>%
    #     as.data.frame() %>%
    #     pivot_longer(cols = everything(), names_to = "sample", values_to = "value") %>%
    #     mutate(group = ifelse(as.numeric(substr(sample, 14, 15)) < 10, "tumor", "normal"))
    input_mat <- data_filt %>%
        as.data.frame() %>%
        t()
    allowWGCNAThreads()
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    sft <- pickSoftThreshold(
        input_mat,
        # blockSize = 30,
        powerVector = powers,
        verbose = 5
    )
    pdf("data/104.WGCNA/soft_threshold.pdf")
    par(mfrow = c(1, 2))
    cex1 <- 0.9

    plot(sft$fitIndices[, 1],
        -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        xlab = "Soft Threshold (power)",
        ylab = "Scale Free Topology Model Fit, signed R^2",
        main = paste("Scale independence")
    )
    text(sft$fitIndices[, 1],
        -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        labels = powers, cex = cex1, col = "red"
    )
    abline(h = 0.90, col = "red")
    plot(sft$fitIndices[, 1],
        sft$fitIndices[, 5],
        xlab = "Soft Threshold (power)",
        ylab = "Mean Connectivity",
        type = "n",
        main = paste("Mean connectivity")
    )
    text(sft$fitIndices[, 1],
        sft$fitIndices[, 5],
        labels = powers,
        cex = cex1, col = "red"
    )
    dev.off()

    picked_power <- 9
    temp_cor <- cor
    cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)
    netwk <- blockwiseModules(input_mat, # <= input here

        # == Adjacency Function ==
        power = picked_power, # <= power here
        networkType = "signed",

        # == Tree and Block Options ==
        deepSplit = 2,
        pamRespectsDendro = F,
        # detectCutHeight = 0.75,
        minModuleSize = 30,
        maxBlockSize = 4000,

        # == Module Adjustments ==
        reassignThreshold = 0,
        mergeCutHeight = 0.25,

        # == Output Options
        numericLabels = T,
        verbose = 3
    )
    # Convert labels to colors for plotting
    mergedColors <- labels2colors(netwk$colors)
    # Plot the dendrogram and the module colors underneath
    pdf("data/104.WGCNA/dendrogram.pdf")
    plotDendroAndColors(
        netwk$dendrograms[[1]],
        mergedColors[netwk$blockGenes[[1]]],
        "Module colors",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05
    )
    dev.off()
    module_df <- data.frame(
        gene_id = names(netwk$colors),
        colors = labels2colors(netwk$colors)
    )
    write_delim(module_df, "data/104.WGCNA/module_colors.tsv", delim = "\t")
    WGCNA_res <- list(
        netwk = netwk,
        module_df = module_df,
        input_mat = input_mat
    )
}

plot_WGCNA <- function(WGCNA_res, ansEA, data_filt, DEG_res) {
    EA_genes <- ansEA$trp_pathway$plot.data.gene %>%
        as_tibble() %>%
        dplyr::filter(!is.na(mol.data)) %>%
        pull(labels) %>%
        unique()
    # use bitr to convert SYMBOL to ensembl
    EA_genes <- bitr(EA_genes, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
    # filter rownames(data_filt) %in% EA_genes
    data_EA <- data_filt[rownames(data_filt) %in% EA_genes$ENSEMBL, ] %>%
        as.data.frame() %>%
        rownames_to_column(var = "gene_id") %>%
        as_tibble() %>%
        pivot_longer(cols = starts_with("TCGA"), names_to = "sample", values_to = "value") %>%
        mutate(group = ifelse(as.numeric(substr(sample, 14, 15)) < 10, "tumor", "normal")) %>%
        left_join(., EA_genes, by = c("gene_id" = "ENSEMBL")) %>%
        left_join(., DEG_res$DEGs, by = c("SYMBOL" = "gene_name"))
    p <- grouped_ggbetweenstats(
        data = data_EA,
        x = group,
        y = value,
        grouping.var = SYMBOL,
        # caption = "Trypthophan Pathway Genes Expression",
        xlab = "Group",
        ylab = "Expression",
        bf.message = FALSE
    )
    ggsave("result/103.enrich/trp_pathway_expression.png", p, width = 18, height = 15)

    netwk <- WGCNA_res$netwk
    module_df <- WGCNA_res$module_df
    input_mat <- WGCNA_res$input_mat
    # Get Module Eigengenes per cluster
    mergedColors <- labels2colors(netwk$colors)
    MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

    # Reorder modules so similar modules are next to each other
    MEs0 <- orderMEs(MEs0)
    module_order <- names(MEs0) %>% gsub("ME", "", .)

    # separate MEs0 by normal and tumor, into two dataframes
    MEs0_normal <- MEs0[as.numeric(substr(rownames(MEs0), 14, 15)) >= 10, ]
    MEs0_tumor <- MEs0[as.numeric(substr(rownames(MEs0), 14, 15)) < 10, ]
    # hclust cols of MEs0_tumor, get the order
    MEs0_tumor_order <- MEs0_tumor %>%
        t() %>%
        dist() %>%
        hclust()
    # use pheatmap to plot MEs0_tumor and MEs0_normal respectively
    p1 <- pheatmap(MEs0_tumor %>% t(),
        cluster_rows = MEs0_tumor_order, cluster_cols = TRUE, show_colnames = FALSE, silent = TRUE, legend = FALSE, border_color = NA,
        show_rownames = FALSE, breaks = seq(-1, 1, by = 0.01), color = colorRampPalette(c("blue", "white", "red"))(200)
    )
    p2 <- pheatmap(MEs0_normal %>% t(),
        cluster_rows = MEs0_tumor_order, cluster_cols = TRUE, show_colnames = FALSE, silent = TRUE, border_color = NA,
        treeheight_row = 0, breaks = seq(-1, 1, by = 0.01), color = colorRampPalette(c("blue", "white", "red"))(200)
    )
    # use patchwork to combine p1 and p2
    p <- wrap_plots(p1$gtable, p2$gtable)
    ggsave("result/104.WGCNA/module_eigengenes.png", p)

    write_tsv(
        data_EA %>% dplyr::select(gene_id, SYMBOL, logFC, PValue, FDR) %>% distinct() %>%
            left_join(., module_df, by = "gene_id"),
        "result/103.enrich/trp_pathway_expression.tsv"
    )
    c("result/103.enrich/trp_pathway_expression.png", "result/104.WGCNA/module_eigengenes.png", "result/103.enrich/trp_pathway_expression.tsv")
}
