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
    gene_intersect <- intersect(rowData(data_mRNA)$gene_name, rownames(gtex_data))
    gtex_data <- gtex_data %>%
        filter(gene_name %in% gene_intersect)
    data_mRNA <- data_mRNA %>% filter(gene_name %in% gene_intersect)
    data_mRNA <- data_mRNA %>%
        select(.feature, .sample, counts)
    gtex_data <- gtex_data %>%
        select(.feature, .sample, counts)
    data_combined <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = cbind(assay(data_mRNA), assay(gtex_data)))
    )
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
run_DEG <- function(data_filt) {
    data_dds <- data_filt %>%
        test_differential_abundance(
            .formula = ~ 0 + group,
            .abundance = counts_scaled_adjusted,
            scaling_method = "none",
            contrasts = c("grouptumor - groupnormal"),
            action = "add",
            method = "limma_voom"
        )
    data_dds <- data_dds %>%
        dplyr::mutate(log2FoldChange = `logFC___grouptumor - groupnormal`) %>%
        dplyr::mutate(pvalue = `P.Value___grouptumor - groupnormal`) %>%
        dplyr::mutate(padj = `adj.P.Val___grouptumor - groupnormal`) %>%
        select(.sample, .feature, counts_scaled_adjusted, group, log2FoldChange, pvalue, padj) %>%
        filter(padj < 0.05 & abs(log2FoldChange) > 1)
    write_tsv(
        data_dds %>% as_tibble() %>%
            select(-.sample, -counts_scaled_adjusted, -group) %>% distinct(),
        "data/102.DEG/data_dds.tsv"
    )
    data_dds
}

run_enrich <- function(data_dds) {
    data_dds <- data_dds %>%
        as_tibble() %>%
        select(-.sample, -counts_scaled_adjusted, -group)
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
    ggsave("result/103.enrich/KEGG_up.png", p)
    kk_down <- enrichKEGG(
        gene = gene_list$down,
        organism = "hsa",
        qvalueCutoff = 1
        # pvalueCutoff = 1
    )
    # kk_down <- kk_down %>% filter(category == "Metabolism")
    p <- dotplot(kk_down)
    ggsave("result/103.enrich/KEGG_down.png", p)
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
    ggsave("result/103.enrich/KEGG_gse.png", p)
    kk_up <- setReadable(kk_up, org.Hs.eg.db, keyType = "ENTREZID")
    kk_down <- setReadable(kk_down, org.Hs.eg.db, keyType = "ENTREZID")
    kk_gse <- setReadable(kk_gse, org.Hs.eg.db, keyType = "ENTREZID")
    write_tsv(kk_up %>% as.data.frame(), "data/103.enrich/KEGG_up.tsv")
    write_tsv(kk_down %>% as.data.frame(), "data/103.enrich/KEGG_down.tsv")
    write_tsv(kk_gse %>% as.data.frame(), "data/103.enrich/KEGG_gse.tsv")

    p <- plotEnrichAdv(
        kk_up %>% as.data.frame() %>% filter(category == "Metabolism") %>% rename(FoldEnrich = FoldEnrichment),
        kk_down %>% as.data.frame() %>% filter(category == "Metabolism") %>% rename(FoldEnrich = FoldEnrichment),
        plot_type = "one",
        term_metric = "FoldEnrich",
        stats_metric = "p.adjust"
    )
    ggsave("result/103.enrich/metabolism_enrich.png", p, width = 10)
    trp_pathway <- pathview(
        gene.data = gene_gsea, same.layer = FALSE, limit = list(gene = 4, cpd = 1),
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
        assay(., "counts_scaled_adjusted") %>%
        as.data.frame() %>%
        t()
    allowWGCNAThreads()
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    sft <- pickSoftThreshold(
        input_mat,
        RsquaredCut = 0.82,
        networkType = "signed",
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

    picked_power <- sft$fitIndices %>%
        dplyr::filter(SFT.R.sq > 0.8) %>%
        pull(Power) %>%
        min()
    temp_cor <- cor
    cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)
    netwk <- blockwiseModules(input_mat, # <= input here

        # == Adjacency Function ==
        power = picked_power, # <= power here
        networkType = "signed",
        TOMType = "signed",
        # == Tree and Block Options ==
        deepSplit = 2,
        pamRespectsDendro = F,
        # detectCutHeight = 0.75,
        minModuleSize = 30,
        maxBlockSize = 4000,

        # == Module Adjustments ==
        # reassignThreshold = 0,
        mergeCutHeight = 0.6,

        # == Output Options
        numericLabels = TRUE,
        verbose = 3,
        nThreads = 50
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
        input_mat = input_mat,
        picked_power = picked_power
    )
}

plot_WGCNA <- function(WGCNA_res, ansEA, data_filt, data_dds) {
    EA_genes <- ansEA$trp_pathway$plot.data.gene %>%
        as_tibble() %>%
        dplyr::filter(!is.na(mol.data)) %>%
        pull(labels) %>%
        unique()
    # filter rownames(data_filt) %in% EA_genes
    data_EA <- data_dds %>%
        filter(.feature %in% EA_genes)
    p <- grouped_ggbetweenstats(
        data = data_EA,
        x = group,
        y = counts_scaled_adjusted,
        grouping.var = .feature,
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
    tumor_samples <- data_dds %>%
        filter(group == "tumor") %>%
        pull(.sample)
    MEs0_normal <- MEs0 %>%
        as_tibble(rownames = "sample") %>%
        filter(!sample %in% tumor_samples)
    MEs0_tumor <- MEs0 %>%
        as_tibble(rownames = "sample") %>%
        filter(sample %in% tumor_samples)
    # hclust cols of MEs0_tumor, get the order
    MEs0_tumor_order <- MEs0_tumor %>%
        select(-sample) %>%
        t() %>%
        dist() %>%
        hclust()
    # use pheatmap to plot MEs0_tumor and MEs0_normal respectively
    p1 <- pheatmap(MEs0_tumor %>% column_to_rownames("sample") %>% t(),
        cluster_rows = MEs0_tumor_order, cluster_cols = TRUE, show_colnames = FALSE, silent = TRUE, legend = FALSE, border_color = NA,
        show_rownames = FALSE, breaks = seq(-1, 1, by = 0.01), color = colorRampPalette(c("blue", "white", "red"))(200)
    )
    p2 <- pheatmap(MEs0_normal %>% column_to_rownames("sample") %>% t(),
        cluster_rows = MEs0_tumor_order, cluster_cols = TRUE, show_colnames = FALSE, silent = TRUE, border_color = NA,
        treeheight_row = 0, breaks = seq(-1, 1, by = 0.01), color = colorRampPalette(c("blue", "white", "red"))(200)
    )
    # use patchwork to combine p1 and p2
    p <- wrap_plots(p1$gtable, p2$gtable)
    ggsave("result/104.WGCNA/module_eigengenes.png", p)

    data_EA_tidy <- data_EA %>%
        dplyr::select(.feature, log2FoldChange, pvalue, padj) %>%
        distinct() %>%
        left_join(., module_df, by = c(".feature" = "gene_id"))
    write_tsv(
        data_EA_tidy, "result/103.enrich/trp_pathway_expression.tsv"
    )
    data_EA_tidy
}

get_network <- function(WGCNA_res, data_dds, data_EA_tidy) {
    netwk <- WGCNA_res$netwk
    module_df <- WGCNA_res$module_df
    input_mat <- WGCNA_res$input_mat
    picked_power <- WGCNA_res$picked_power
    TOM <- TOMsimilarityFromExpr(t(data_dds %>% assay(., "counts_scaled_adjusted")),
        networkType = "signed",
        power = picked_power
    )
    rownames(TOM) <- rownames(data_dds %>% assay(., "counts_scaled_adjusted"))
    colnames(TOM) <- rownames(data_dds %>% assay(., "counts_scaled_adjusted"))
    # only keep the upper triangular part of the TOM:
    TOM[upper.tri(TOM)] <- NA

    # cast the network from wide to long format
    cur_network <- TOM %>%
        reshape2::melt() %>%
        dplyr::rename(gene1 = Var1, gene2 = Var2, weight = value) %>%
        subset(!is.na(weight)) %>%
        # remove gene1 == gene2
        subset(gene1 != gene2)

    # get the module & color info for gene1
    temp1 <- dplyr::inner_join(
        cur_network,
        module_df %>%
            dplyr::select(c(gene_id, colors)) %>%
            dplyr::rename(gene1 = gene_id, color1 = colors),
        by = "gene1"
    ) %>% dplyr::select(color1)

    # get the module & color info for gene2
    temp2 <- dplyr::inner_join(
        cur_network,
        module_df %>%
            dplyr::select(c(gene_id, colors)) %>%
            dplyr::rename(gene2 = gene_id, color2 = colors),
        by = "gene2"
    ) %>% dplyr::select(color2)

    # add the module & color info
    cur_network <- cbind(cur_network, temp1, temp2)

    # set the edge color to the module's color if they are the two genes are in the same module
    cur_network$edge_color <- ifelse(
        cur_network$color1 == cur_network$color2,
        as.character(cur_network$color1),
        "grey"
    )

    # keep this network before subsetting
    cur_network_full <- cur_network

    # keep the top 10% of edges
    edge_percent <- 0.1
    cur_network <- cur_network_full %>%
        dplyr::slice_max(
            order_by = weight,
            n = round(nrow(cur_network) * edge_percent)
        )

    # make the graph object with tidygraph
    graph <- cur_network %>%
        igraph::graph_from_data_frame() %>%
        tidygraph::as_tbl_graph(directed = FALSE) %>%
        tidygraph::activate(nodes)

    # add the module name to the graph:
    graph <- graph %>%
        left_join(module_df, by = c("name" = "gene_id")) %>%
        dplyr::rename(module = colors)

    # get the top 25 hub genes for each module
    # Get Module Eigengenes per cluster
    mergedColors <- labels2colors(netwk$colors)
    MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

    # Reorder modules so similar modules are next to each other
    MEs0 <- orderMEs(MEs0)
    datKME <- signedKME(data_dds %>% assay(., "counts_scaled_adjusted") %>% t(),
        MEs0,
        outputColumnName = "MM."
    )
    # get the top 25 hub genes for each module
    hub_genes <- datKME %>%
        rownames_to_column(var = "gene_name") %>%
        pivot_longer(cols = starts_with("MM"), names_to = "module", values_to = "kME") %>%
        group_by(module) %>%
        top_n(10, kME) %>%
        pull(gene_name)
    trp_genes <- data_EA_tidy$.feature
    graph <- graph %>%
        activate(nodes) %>%
        filter(name %in% c(hub_genes, trp_genes))

    # add a column for name in data_EA_tidy$.feature, yellow

    # if (!is_connected(graph)) {
    #     components <- clusters(graph)
    #     largest_component <- which.max(components$csize)
    #     graph <- induced_subgraph(graph, components$membership == largest_component) %>%
    #         tidygraph::as_tbl_graph(directed = FALSE) %>%
    #         tidygraph::activate(nodes)
    # }
    # make the plot with gggraph
    plot_network <- function(trp_gene) {
        graph <- graph %>%
            activate(nodes) %>%
            mutate(
                depth = bfs_dist(root = which(name == trp_gene)),
                label = ifelse(depth %in% c(0, 1), name, ""),
                label_color = ifelse(name %in% data_EA_tidy$.feature, "blue", "black"),
                label_font_face = ifelse(name %in% data_EA_tidy$.feature, "bold", "italic")
            )
        p <- ggraph(graph, layout = "fr") +
            geom_node_voronoi(aes(fill = module, color = module), alpha = 0.3) +
            geom_edge_link0(aes(alpha = weight, color = edge_color, edge_width = weight * 3)) +
            geom_node_point(aes(color = module)) +
            geom_node_label(aes(label = label, color = label_color, fontface = label_font_face),
                repel = TRUE, max.overlaps = Inf
            ) +
            scale_colour_manual(values = graph %>% activate(nodes) %>% pull(module) %>% unique() %>% setNames(., .)) +
            scale_edge_colour_manual(values = graph %>% activate(edges) %>% pull(edge_color) %>% unique() %>% setNames(., .)) +
            scale_fill_manual(values = graph %>% activate(nodes) %>% pull(module) %>% unique() %>% setNames(., .))
        ggsave(paste0("result/105.network/network_", trp_gene, ".png"), p, width = 15, height = 15)
    }
    walk(trp_genes, plot_network)
    hub_genes_25 <- datKME %>%
        rownames_to_column(var = "gene_name") %>%
        pivot_longer(cols = starts_with("MM"), names_to = "module", values_to = "kME") %>%
        group_by(module) %>%
        top_n(25, kME) %>%
        arrange(module, desc(kME))
    write_tsv(hub_genes_25, "result/105.network/hub_genes.tsv")
    "result/105.network/network.png"
}
get_string_network <- function(WGCNA_res, data_EA_tidy) {
    netwk <- WGCNA_res$netwk
    module_df <- WGCNA_res$module_df
    input_mat <- WGCNA_res$input_mat
    module <- data_EA_tidy %>%
        pull(colors) %>%
        unique()
    modProbes <- module_df %>%
        filter(colors %in% module) %>%
        pull(gene_id)
    string_db <- STRINGdb$new(version = "12", species = 9606)
    data_mapped <- string_db$map(modProbes %>% bitr(fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db), "ENTREZID")
    data_links <- data_mapped$STRING_id %>% string_db$get_interactions()
    # 转换stringID为Symbol，只取前两列和最后一列
    links <- data_links %>%
        mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>%
        mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%
        dplyr::select(from, to, last_col()) %>%
        dplyr::rename(weight = combined_score)
    write_tsv(links, "data/105.network/links.tsv")
    # 节点数据
    nodes <- links %>%
        {
            data.frame(gene = c(.$from, .$to))
        } %>%
        distinct()
    # 创建网络图
    # 根据links和nodes创建
    net <- igraph::graph_from_data_frame(d = links, vertices = nodes, directed = F)
    # 添加一些参数信息用于后续绘图
    # V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links)
    igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
    igraph::V(net)$size <- igraph::degree(net) / 5 #
    igraph::E(net)$width <- igraph::E(net)$weight / 10
    # 使用ggraph绘图
    # ggraph是基于ggplot2的包，语法和常规ggplot2类似
    p <- ggraph(net, layout = "kk") +
        geom_edge_fan(aes(edge_width = width), color = "lightblue", show.legend = F) +
        geom_node_point(aes(size = size), color = "orange", alpha = 0.7) +
        geom_node_text(aes(filter = deg > 5, label = name), size = 5, repel = T) +
        scale_edge_width(range = c(0.2, 1)) +
        scale_size_continuous(range = c(1, 10)) +
        guides(size = F) +
        theme_graph()
    ggsave("result/104.WGCNA/module_network.png", p)
}

# get the correlation between IDO1 and BCL11A in data_filt
get_correlation <- function(WGCNA_res, data_filt, data_EA_tidy) {
    netwk <- WGCNA_res$netwk
    module_df <- WGCNA_res$module_df
    input_mat <- WGCNA_res$input_mat
    module <- data_EA_tidy %>%
        filter(.feature %in% "IDO1") %>%
        pull(colors) %>%
        unique()
    modProbes <- module_df %>%
        filter(colors %in% module) %>%
        pull(gene_id)
    expr_of_interest <- data_filt %>% filter(.feature %in% modProbes)
    tidy_expr <- expr_of_interest %>%
        as.data.frame() %>%
        rownames_to_column(var = "gene_id") %>%
        as_tibble() %>%
        pivot_longer(cols = starts_with("TCGA"), names_to = "sample", values_to = "value")
}

get_survival <- function(data, data_filt, data_EA_tidy) {
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data),
        typesample = c("TP")
    )
    get_res <- function(gene) {
        data_clin <- data_filt %>%
            as_tibble() %>%
            filter(.feature %in% gene) %>%
            filter(.sample %in% samplesTP) %>%
            left_join(., colData(data) %>% as.data.frame() %>% rownames_to_column(".sample"), by = ".sample")
        gene_values <- data_clin %>%
            pull(counts_scaled_adjusted)
        data_clin_final <- data_clin %>%
            dplyr::mutate(
                survival_time = ifelse(vital_status == "Dead", days_to_death,
                    days_to_last_follow_up
                ) / 365.25,
                vital_status = ifelse(vital_status == "Dead", 1, 0)
            ) %>%
            mutate(
                gene_group = ifelse(counts_scaled_adjusted >= quantile(gene_values, 0.75), "High",
                    ifelse(counts_scaled_adjusted <= quantile(gene_values, 0.25), "Low", NA)
                )
            ) %>%
            dplyr::filter(!is.na(gene_group)) %>%
            dplyr::select(.sample, survival_time, vital_status, gene_group)
        fit <- survfit(Surv(survival_time, vital_status) ~ gene_group, data = data_clin_final)
        p <- ggsurvplot(fit,
            data = data_clin_final, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
            legend.title = paste0("Gene: ", gene),
            legend.labs = c("Low", "High")
        )$plot
        ggsave(paste0("result/106.survival/survival_", gene, ".png"), p)
    }
    walk(data_EA_tidy$.feature, get_res)
    "result/106.survival"
}

get_module_trait <- function(WGCNA_res, data_filt, data, data_EA_tidy) {
    netwk <- WGCNA_res$netwk
    module_df <- WGCNA_res$module_df
    input_mat <- WGCNA_res$input_mat
    mergedColors <- labels2colors(netwk$colors)
    MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

    # Reorder modules so similar modules are next to each other
    MEs0 <- orderMEs(MEs0)
    data_clin <- data_filt %>%
        pivot_sample() %>%
        left_join(., colData(data) %>% as.data.frame() %>% rownames_to_column(".sample"), by = ".sample") %>%
        dplyr::select(
            .sample, longest_dimension, ajcc_pathologic_stage, ajcc_pathologic_t,
            ajcc_pathologic_n, ajcc_pathologic_m, vital_status, days_to_death, days_to_last_follow_up
        ) %>%
        mutate(
            survival_time = ifelse(vital_status == "Dead", days_to_death,
                days_to_last_follow_up
            ) / 365.25,
            vital_status = ifelse(vital_status == "Dead", 1, 0),
            pathologic_stage = case_when(
                grepl("Stage I", ajcc_pathologic_stage, ignore.case = TRUE) ~ 1,
                grepl("Stage II", ajcc_pathologic_stage, ignore.case = TRUE) ~ 2,
                grepl("Stage III", ajcc_pathologic_stage, ignore.case = TRUE) ~ 3,
                grepl("Stage IV", ajcc_pathologic_stage, ignore.case = TRUE) ~ 4,
                TRUE ~ 0
            ),
            pathologic_t = factor(ajcc_pathologic_t %>% replace_na("N")),
            pathologic_n = factor(ajcc_pathologic_n %>% replace_na("N")),
            pathologic_m = factor(ajcc_pathologic_m %>% replace_na("N"))
        ) %>%
        dplyr::select(.sample, survival_time, vital_status, pathologic_stage, pathologic_t, pathologic_n, pathologic_m)
    # make rows of MEs0 order like data_clin$.sample order
    MEs0 <- MEs0[match(data_clin$.sample, rownames(MEs0)), ]
    # calculate correlation between MEs0 and pathologic_stage
    cor_stage <- corr.test(MEs0, data_clin$pathologic_stage %>% as.numeric(), method = "spearman")
    cor_t <- corr.test(MEs0, data_clin$pathologic_t %>% as.numeric(), method = "spearman")
    cor_n <- corr.test(MEs0, data_clin$pathologic_n %>% as.numeric(), method = "spearman")
    cor_m <- corr.test(MEs0, data_clin$pathologic_m %>% as.numeric(), method = "spearman")
    cor_survival <- corr.test(MEs0, data_clin$survival_time, method = "spearman")
    # combine the correlation results
    cor_data <- list(
        stage = cor_stage$r,
        t = cor_t$r,
        n = cor_n$r,
        m = cor_m$r,
        survival = cor_survival$r
    ) %>%
        purrr::reduce(cbind)
    colnames(cor_data) <- c("stage", "t", "n", "m", "survival")
    cor_data <- cor_data %>%
        as_tibble(rownames = "module") %>%
        pivot_longer(cols = -module, names_to = "trait", values_to = "correlation")
    p_data <- list(
        stage = cor_stage$p,
        t = cor_t$p,
        n = cor_n$p,
        m = cor_m$p,
        survival = cor_survival$p
    ) %>%
        purrr::reduce(cbind)
    colnames(p_data) <- c("stage", "t", "n", "m", "survival")
    p_data <- p_data %>%
        as_tibble(rownames = "module") %>%
        pivot_longer(cols = -module, names_to = "trait", values_to = "p_value")
    # combine cor_data and p_data
    cor_data <- cor_data %>%
        left_join(., p_data, by = c("module", "trait")) %>%
        # set correlation to NA if p_value > 0.05
        mutate(correlation = ifelse(p_value > 0.05, NA, correlation))
    p <- tidyheatmap(cor_data,
        rows = module,
        columns = trait,
        values = correlation,
        color_na = "grey",
        scale = "none",
        main = "Correlation between Module Eigengenes and Clinical Traits",
        angle_col = "0",
        silent = TRUE
    )
    ggsave("result/106.survival/module_trait_correlation.png", p)

    data_clin <- data_filt %>%
        dplyr::filter(.feature %in% data_EA_tidy$.feature) %>%
        left_join(., colData(data) %>% as.data.frame() %>% rownames_to_column(".sample"), by = ".sample") %>%
        dplyr::select(
            .sample, .feature, counts_scaled_adjusted, longest_dimension, ajcc_pathologic_stage, ajcc_pathologic_t,
            ajcc_pathologic_n, ajcc_pathologic_m, vital_status, days_to_death, days_to_last_follow_up
        ) %>%
        mutate(
            survival_time = ifelse(vital_status == "Dead", days_to_death,
                days_to_last_follow_up
            ) / 365.25,
            vital_status = ifelse(vital_status == "Dead", 1, 0),
            pathologic_stage = case_when(
                grepl("Stage I", ajcc_pathologic_stage, ignore.case = TRUE) ~ 1,
                grepl("Stage II", ajcc_pathologic_stage, ignore.case = TRUE) ~ 2,
                grepl("Stage III", ajcc_pathologic_stage, ignore.case = TRUE) ~ 3,
                grepl("Stage IV", ajcc_pathologic_stage, ignore.case = TRUE) ~ 4,
                TRUE ~ 0
            ),
            pathologic_t = factor(ajcc_pathologic_t %>% replace_na("N")),
            pathologic_n = factor(ajcc_pathologic_n %>% replace_na("N")),
            pathologic_m = factor(ajcc_pathologic_m %>% replace_na("N"))
        ) %>%
        dplyr::select(.feature, counts_scaled_adjusted, survival_time, vital_status, pathologic_stage, pathologic_t, pathologic_n, pathologic_m)
    get_cor_res <- function(gene) {
        data_clin_final <- data_clin %>%
            filter(.feature == gene)
        get_res <- function(trait) {
            res <- corr.test(data_clin_final$counts_scaled_adjusted, data_clin_final %>% pull(trait) %>% as.numeric(),
                method = "spearman"
            )
            res
        }
        traits <- c("survival_time", "pathologic_stage", "pathologic_t", "pathologic_n", "pathologic_m")
        cor_res <- map(traits, get_res)
        cor_res_tidy <- cor_res %>%
            map_dfr(~ tibble(correlation = .x$r, p_value = .x$p)) %>%
            mutate(gene = gene, trait = traits)
    }
    cor_res <- map_dfr(data_EA_tidy$.feature, get_cor_res)
    cor_res <- cor_res %>%
        mutate(correlation = ifelse(p_value > 0.05, NA, correlation))
    p <- tidyheatmap(cor_res,
        rows = gene,
        columns = trait,
        values = correlation,
        color_na = "grey",
        scale = "none",
        main = "Correlation between Gene Expression and Clinical Traits",
        angle_col = "0",
        silent = TRUE
    )
    ggsave("result/106.survival/gene_trait_correlation.png", p)

    # # survival analysis for module
    # data_clin <- data_clin %>%
    #     left_join(., MEs0 %>% as_tibble() %>% rownames_to_column(".sample"), by = ".sample")
    # fit <- survfit(Surv(survival_time, vital_status) ~ module, data = data_clin)
    # p <- ggsurvplot(fit,
    #     data = data_clin, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
    #     legend.title = "Module",
    #     legend.labs = data_clin$module %>% unique()
    # )$plot
    # ggsave("result/106.survival/module_survival.png", p)
    "result/106.survival"
}

deconv <- function(data_filt, data) {
    data_cell <- data_filt %>%
        deconvolve_cellularity(action = "get", cores = 30, prefix = "cibersort__") %>%
        pivot_sample()
    data_cell <- data_cell %>%
        dplyr::select(.sample, starts_with("cibersort__"))
    data_cell <- data_cell %>%
        pivot_longer(
            names_to = "Cell_type_inferred",
            values_to = "proportion",
            names_prefix = "cibersort__",
            cols = contains("cibersort__")
        )
    p <- data_cell %>%
        ggplot(aes(x = Cell_type_inferred, y = proportion)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio = 1 / 5)
    ggsave("result/107.deconv/cellularity.png", p)
    data_filt %>% test_differential_cellularity(. ~ group )
}
