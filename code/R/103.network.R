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