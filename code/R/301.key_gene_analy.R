get_data_full <- function(data, data_filt) {
    # add colData(data) to data_filt, based on .sample
    data_full <- data_filt %>%
        left_join(colData(data) %>% as.data.frame() %>%
            rownames_to_column(".sample"), by = c(".sample" = ".sample"))
    data_full
}


kg_DEG_single <- function(data_full, data_dds, kgene) {
    format_group_labels <- function(breaks, data, value_col) {
        sapply(breaks, function(b) {
            median_val <- data %>%
                filter(get(quo_name(enquo(value_col))) == b) %>%
                pull(counts_scaled_adjusted) %>%
                median() %>%
                round(2)
            paste0(b, "\nmedian: ", median_val)
        })
    }
    # kgene <- "GBP5"
    data_dds_features <- data_dds %>% pull(.feature)
    if (kgene %in% data_dds_features) {
        data_gene_dds <- data_dds %>%
            filter(.feature %in% kgene) %>%
            as_tibble() %>%
            distinct()
        padj <- format(data_gene_dds$padj %>% unique(), digits = 3)
        log2FoldChange <- format(data_gene_dds$log2FoldChange %>% unique(), digits = 3)
    } else {
        data_gene_dds <- data_full %>%
            filter(.feature %in% kgene) %>%
            as_tibble() %>%
            distinct()
        padj <- "none"
        log2FoldChange <- "none"
    }
    p1 <- data_gene_dds %>%
        tidyplot(x = group, y = counts_scaled_adjusted, color = group) %>%
        add_violin() %>%
        add_data_points_beeswarm() %>%
        adjust_x_axis(
            labels = function(x) format_group_labels(x, data_gene_dds, group),
            rotate_labels = TRUE
        ) %>%
        add_title(paste0(kgene, " expression in tumor and normal")) %>%
        add_caption(paste("padj: ", padj, "log2FoldChange: ", log2FoldChange))
    dir.create(paste0("result/301.key_gene_analy/", kgene), showWarnings = FALSE)
    save_plot(p1, filename = paste0("result/301.key_gene_analy/", kgene, "/group_violin.png"))

    # plot across ajcc_pathologic_stage
    data_full <- data_full %>%
        as_tibble()
    data_gene <- data_full %>%
        filter(group == "tumor") %>%
        filter(.feature %in% kgene) %>%
        distinct() %>%
        filter(!is.na(ajcc_pathologic_stage) & !is.na(ajcc_pathologic_t) & !is.na(ajcc_pathologic_n) & !is.na(ajcc_pathologic_m)) %>%
        mutate(
            ajcc_pathologic_stage = factor(ajcc_pathologic_stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
            ajcc_pathologic_t = factor(ajcc_pathologic_t, levels = c("T1", "T1a", "T1b", "T2", "T2a", "T2b", "T3", "T3a", "T3b", "T3c", "T4")),
            ajcc_pathologic_n = factor(ajcc_pathologic_n, levels = c("N0", "N1", "NX")),
            ajcc_pathologic_m = factor(ajcc_pathologic_m, levels = c("M0", "M1", "MX"))
        )
    p2 <- data_gene %>%
        tidyplot(x = ajcc_pathologic_stage, y = counts_scaled_adjusted, color = ajcc_pathologic_stage) %>%
        add_violin() %>%
        add_data_points_beeswarm() %>%
        adjust_x_axis(
            labels = function(x) format_group_labels(x, data_gene, ajcc_pathologic_stage),
            rotate_labels = TRUE
        ) %>%
        add_test_pvalue(method = "wilcoxon", p.adjust.method = "BH", ref.group = "Stage I", hide.ns = TRUE, hide_info = TRUE) %>%
        add_title(paste0(kgene, " expression in different ajcc_pathologic_stage"))
    save_plot(p2, filename = paste0("result/301.key_gene_analy/", kgene, "/ajcc_pathologic_stage_violin.png"))
    p3 <- data_gene %>%
        tidyplot(x = ajcc_pathologic_t, y = counts_scaled_adjusted, color = ajcc_pathologic_t) %>%
        add_violin() %>%
        add_data_points_beeswarm() %>%
        adjust_x_axis(
            labels = function(x) format_group_labels(x, data_gene, ajcc_pathologic_t),
            rotate_labels = TRUE
        ) %>%
        add_test_pvalue(method = "wilcoxon", p.adjust.method = "BH", ref.group = "T1", hide.ns = TRUE, hide_info = TRUE) %>%
        add_title(paste0(kgene, " expression in different ajcc_pathologic_t")) %>%
        adjust_size(width = NA, height = NA)
    save_plot(p3, filename = paste0("result/301.key_gene_analy/", kgene, "/ajcc_pathologic_t_violin.png"))
    data_gene_c <- data_gene %>%
        mutate(ajcc_pathologic_t_c = case_when(
            grepl("T1", ajcc_pathologic_t) ~ "T1",
            grepl("T2", ajcc_pathologic_t) ~ "T2",
            grepl("T3", ajcc_pathologic_t) ~ "T3",
            grepl("T4", ajcc_pathologic_t) ~ "T4",
            TRUE ~ ajcc_pathologic_t
        ))
    p3_c <- data_gene_c %>%
        tidyplot(x = ajcc_pathologic_t_c, y = counts_scaled_adjusted, color = ajcc_pathologic_t_c) %>%
        add_violin() %>%
        add_data_points_beeswarm() %>%
        adjust_x_axis(
            labels = function(x) format_group_labels(x, data_gene_c, ajcc_pathologic_t_c),
            rotate_labels = TRUE
        ) %>%
        add_test_pvalue(method = "wilcoxon", p.adjust.method = "BH", ref.group = "T1", hide.ns = TRUE, hide_info = TRUE) %>%
        add_title(paste0(kgene, " expression in different combined ajcc_pathologic_t")) %>%
        adjust_size(width = NA, height = NA)
    save_plot(p3, filename = paste0("result/301.key_gene_analy/", kgene, "/ajcc_pathologic_t_c_violin.png"))
    p4 <- data_gene %>%
        tidyplot(x = ajcc_pathologic_n, y = counts_scaled_adjusted, color = ajcc_pathologic_n) %>%
        add_violin() %>%
        add_data_points_beeswarm() %>%
        adjust_x_axis(
            labels = function(x) format_group_labels(x, data_gene, ajcc_pathologic_n),
            rotate_labels = TRUE
        ) %>%
        add_test_pvalue(method = "wilcoxon", p.adjust.method = "BH", ref.group = "N0", hide.ns = TRUE, hide_info = TRUE) %>%
        add_title(paste0(kgene, " expression in different ajcc_pathologic_n"))
    save_plot(p4, filename = paste0("result/301.key_gene_analy/", kgene, "/ajcc_pathologic_n_violin.png"))
    p5 <- data_gene %>%
        tidyplot(x = ajcc_pathologic_m, y = counts_scaled_adjusted, color = ajcc_pathologic_m) %>%
        add_violin() %>%
        add_data_points_beeswarm() %>%
        adjust_x_axis(
            labels = function(x) format_group_labels(x, data_gene, ajcc_pathologic_m),
            rotate_labels = TRUE
        ) %>%
        add_test_pvalue(method = "wilcoxon", p.adjust.method = "BH", ref.group = "M0", hide.ns = TRUE, hide_info = TRUE) %>%
        add_title(paste0(kgene, " expression in different ajcc_pathologic_m"))
    save_plot(p5, filename = paste0("result/301.key_gene_analy/", kgene, "/ajcc_pathologic_m_violin.png"))
    p_list <- list(p1, p2, p4, p5, p3_c, p3)
    # use patchwork to combine plots
    p <- wrap_plots(p_list[1:5], nrow = 2) / p_list[[6]] + plot_layout(heights = c(2, 1))
    ggsave(paste0("result/301.key_gene_analy/", kgene, "/combine_violin.png"), p, width = 13, height = 10)
    p_list
}

kg_sc <- function(sc_opt, kgene) {
    # kgene <- "GBP5"
    gene_expression <- FetchData(sc_opt, vars = kgene)

    # 定义表达量大于0的细胞
    cells.highlight <- rownames(gene_expression)[gene_expression[[kgene]] > 0]
    p1 <- DimPlot2(sc_opt,
        features = c(kgene, "cell_type_opt"), label = TRUE,
        # 高亮表达量大于0的细胞
        cells.highlight = cells.highlight,
        nrow = 1
    )
    p2 <- VlnPlot2(sc_opt, features = c(kgene), group.by = "cell_type_opt")
    toplot <- CalcStats(sc_opt, features = kgene, method = "zscore", order = "p", group.by = "cell_type_opt")
    p3 <- Heatmap(toplot, lab_fill = "zscore")
    toplot <- CalcStats(
        sc_opt %>%
            filter(cell_type_opt == "RCC"),
        features = kgene, method = "mean", order = "p", group.by = "batch", p.threshold = 1
    )
    p4 <- Heatmap(toplot, lab_fill = "mean")
    p <- (p1 / p2) / p3 / p4 + plot_layout(heights = c(1, 1, 0.2, 0.2))

    ggsave(paste0("result/301.key_gene_analy/", kgene, "/sc_", kgene, ".png"), p, width = 20, height = 14)

    # see percentage of AGER expression in each cell type
    high_cells <- cells.highlight
    sc_opt <- sc_opt %>%
        mutate(gene_exp = ifelse(.cell %in% high_cells, paste0(kgene, "_exp"), paste0(kgene, "_not_exp")))

    p <- ClusterDistrBar(origin = sc_opt$cell_type_opt, cluster = sc_opt$gene_exp)
    ggsave(paste0("result/301.key_gene_analy/", kgene, "/sc_", kgene, "_cluster.png"), p, width = 10, height = 6)
    "result/301.key_gene_analy"
}

# survival
kg_survival <- function(data_full, kgene) {
    # kgene = "GBP5"
    data_clin <- data_full %>%
        as_tibble() %>%
        filter(.feature %in% kgene) %>%
        filter(group == "tumor")
    gene_values <- data_clin %>%
        pull(counts_scaled_adjusted)

    # Cox regression analysis
    data_clin_final_cox <- data_clin %>%
        dplyr::mutate(
            survival_time = ifelse(vital_status == "Dead", days_to_death,
                days_to_last_follow_up
            ) / 365.25,
            vital_status = ifelse(vital_status == "Dead", 1, 0),
            gender = ifelse(gender == "male", 1, 0)
        ) %>%
        dplyr::select(.sample, survival_time, vital_status, counts_scaled_adjusted, gender, age_at_index)
    formula <- as.formula(paste("Surv(survival_time, vital_status) ~", "counts_scaled_adjusted + gender + age_at_index"))
    fit_cox <- coxph(formula, data = data_clin_final_cox)
    fit_summary <- summary(fit_cox)
    p.value <- signif(fit_summary$coefficients["counts_scaled_adjusted", "Pr(>|z|)"], digits = 2)
    beta <- signif(fit_summary$coefficients["counts_scaled_adjusted", "coef"], digits = 2)
    HR <- signif(fit_summary$coefficients["counts_scaled_adjusted", "exp(coef)"], digits = 2)
    cox_res <- c(beta, HR, p.value)
    names(cox_res) <- c("coef", "HR", "p.value")

    # KM analysis
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
    p_survival <- survdiff(Surv(survival_time, vital_status) ~ gene_group, data = data_clin_final)$pvalue

    # Add Cox results to plot
    cox_label <- sprintf("Cox analysis:\nHR = %.2f\np = %.2e", HR, p.value)
    p_list <- ggsurvplot(fit,
        data = data_clin_final, risk.table = TRUE, pval = FALSE, conf.int = TRUE,
        legend.title = paste0("Gene: ", kgene),
        # legend.labs = c("Low", "High"),
        cumevents = TRUE, cumcensor = TRUE, combine = TRUE
    )
    p1 <- p_list$plot + annotate("text",
        x = max(data_clin_final$survival_time) * 0.3,
        y = 0.25,
        label = sprintf("%s\nSurvival p = %.2e", cox_label, p_survival),
        hjust = 1
    )
    p <- p1 / p_list$table + plot_layout(heights = c(2, 1))

    ggsave(paste0("result/301.key_gene_analy/", kgene, "/survival.png"), p, width = 10, height = 10)
    "result/301.key_gene_analy"
}

kg_deconv_tumor <- function(data_full, ref_deconv) {
    # read data/201.load_sc/reference.h5ad
    # Extract the bulk RNA-seq counts
    data_full <- data_full %>%
        filter(group == "tumor")
    bulk_counts <- assay(data_full, "counts_scaled_adjusted")
    bulk_counts <- as.matrix(bulk_counts)

    # Extract single cell data and metadata from Seurat object
    ref_deconv <- ref_deconv %>%
        filter(summaryDescription == "Tumour")
    sc_counts <- as.matrix(ref_deconv@assays$RNA@data)
    cell_type_annotations <- ref_deconv$broad_type
    # batch_ids <- ref_deconv$dataset # Assuming this is your batch/sample identifier
    # Subsample cells (optional but recommended for computational efficiency)
    max_cells_per_cell_type_dtl <- 500

    sampled_metadata <- ref_deconv@meta.data %>%
        as_tibble() %>%
        group_by(broad_type) %>%
        nest() %>%
        mutate(n = map_dbl(data, nrow)) %>%
        filter(n >= 100) %>%
        mutate(n = pmin(n, max_cells_per_cell_type_dtl)) %>%
        ungroup() %>%
        mutate(samp = map2(data, n, ~ slice_sample(.x, n = .y))) %>%
        select(-data) %>%
        # filter out n < 100
        unnest(samp)

    ref_deconv_sampled <- ref_deconv %>%
        filter(`_index` %in% sampled_metadata$`_index`)

    # Extract counts matrix from subsampled data
    sc_counts_sampled <- as.matrix(ref_deconv_sampled@assays$RNA@counts)
    cell_type_annotations_sampled <- ref_deconv_sampled$broad_type
    # batch_ids_sampled <- ref_deconv_sampled$dataset
    # # Perform deconvolution using DWLS
    # # First build the signature matrix
    # signature_matrix_dwls <- omnideconv::build_model(
    #     single_cell_object = sc_counts_sampled,
    #     cell_type_annotations = cell_type_annotations_sampled,
    #     method = "dwls",
    #     dwls_method = "mast_optimized"
    # )

    # # Then perform deconvolution
    # deconv_results_dwls <- deconvolute_dwls(
    #     bulk_gene_expression = bulk_counts,
    #     signature = signature_matrix_dwls,
    #     dwls_submethod = "DampenedWLS"
    # )

    # Perform deconvolution using BayesPrism
    deconv_results_bayesprism <- omnideconv::deconvolute_bayesprism(
        bulk_gene_expression = bulk_counts,
        single_cell_object = sc_counts_sampled,
        cell_type_annotations = cell_type_annotations_sampled,
        # signature = NULL,
        n_cores = 50 # Adjust based on your system
    )

    deconv_results_bayesprism
}
kg_deconv_normal <- function(data_full, ref_deconv) {
    # read data/201.load_sc/reference.h5ad
    # Extract the bulk RNA-seq counts
    data_full <- data_full %>%
        filter(group == "normal")
    bulk_counts <- assay(data_full, "counts_scaled_adjusted")
    bulk_counts <- as.matrix(bulk_counts)

    # Extract single cell data and metadata from Seurat object
    ref_deconv <- ref_deconv %>%
        filter(summaryDescription == "Normal kidney")
    sc_counts <- as.matrix(ref_deconv@assays$RNA@data)
    cell_type_annotations <- ref_deconv$broad_type
    # batch_ids <- ref_deconv$dataset # Assuming this is your batch/sample identifier
    # Subsample cells (optional but recommended for computational efficiency)
    max_cells_per_cell_type_dtl <- 500

    sampled_metadata <- ref_deconv@meta.data %>%
        as_tibble() %>%
        group_by(broad_type) %>%
        nest() %>%
        mutate(n = map_dbl(data, nrow)) %>%
        filter(n >= 100) %>%
        mutate(n = pmin(n, max_cells_per_cell_type_dtl)) %>%
        ungroup() %>%
        mutate(samp = map2(data, n, ~ slice_sample(.x, n = .y))) %>%
        select(-data) %>%
        # filter out n < 100
        unnest(samp)

    ref_deconv_sampled <- ref_deconv %>%
        filter(`_index` %in% sampled_metadata$`_index`)

    # Extract counts matrix from subsampled data
    sc_counts_sampled <- as.matrix(ref_deconv_sampled@assays$RNA@counts)
    cell_type_annotations_sampled <- ref_deconv_sampled$broad_type
    # batch_ids_sampled <- ref_deconv_sampled$dataset
    # # Perform deconvolution using DWLS
    # # First build the signature matrix
    # signature_matrix_dwls <- omnideconv::build_model(
    #     single_cell_object = sc_counts_sampled,
    #     cell_type_annotations = cell_type_annotations_sampled,
    #     method = "dwls",
    #     dwls_method = "mast_optimized"
    # )

    # # Then perform deconvolution
    # deconv_results_dwls <- deconvolute_dwls(
    #     bulk_gene_expression = bulk_counts,
    #     signature = signature_matrix_dwls,
    #     dwls_submethod = "DampenedWLS"
    # )

    # Perform deconvolution using BayesPrism
    deconv_results_bayesprism <- omnideconv::deconvolute_bayesprism(
        bulk_gene_expression = bulk_counts,
        single_cell_object = sc_counts_sampled,
        cell_type_annotations = cell_type_annotations_sampled,
        # signature = NULL,
        n_cores = 50 # Adjust based on your system
    )

    deconv_results_bayesprism
}

plot_deconv <- function(kg_deconv_tumor_res, kg_deconv_normal_res, data_full) {
    deconv_res_tumor <- kg_deconv_tumor_res$theta %>%
        as.data.frame() %>%
        rownames_to_column("Sample") %>%
        pivot_longer(-Sample, names_to = "celltype", values_to = "cell_fraction")
    deconv_res_normal <- kg_deconv_normal_res$theta %>%
        as.data.frame() %>%
        rownames_to_column("Sample") %>%
        pivot_longer(-Sample, names_to = "celltype", values_to = "cell_fraction")
    # find cell type that exist only in tumor or normal, fill the cell type of the other with 0
    deconv_res_tumor$type <- "tumor"
    deconv_res_normal$type <- "normal"
    deconv_res <- rbind(deconv_res_tumor, deconv_res_normal)

    # Get all unique cell types
    all_celltypes <- unique(c(
        unique(deconv_res_tumor$celltype),
        unique(deconv_res_normal$celltype)
    ))

    # Create missing combinations for tumor
    missing_tumor <- expand.grid(
        Sample = unique(deconv_res_tumor$Sample),
        celltype = setdiff(all_celltypes, unique(deconv_res_tumor$celltype)),
        stringsAsFactors = FALSE
    ) %>%
        mutate(
            cell_fraction = 0,
            type = "tumor"
        ) %>%
        as_tibble()

    # Create missing combinations for normal
    missing_normal <- expand.grid(
        Sample = unique(deconv_res_normal$Sample),
        celltype = setdiff(all_celltypes, unique(deconv_res_normal$celltype)),
        stringsAsFactors = FALSE
    ) %>%
        mutate(
            cell_fraction = 0,
            type = "normal"
        ) %>%
        as_tibble()

    deconv_res_full <- rbind(deconv_res, missing_tumor, missing_normal)

    p <- tidyplot(deconv_res_full, x = celltype, y = cell_fraction, color = type) %>%
        add_boxplot() %>%
        # add_data_points_beeswarm() %>%
        add_title("Cell fraction in tumor and normal") %>%
        adjust_size(width = 300, height = 100) %>%
        adjust_x_axis(rotate_labels = 45) %>%
        add_test_asterisks(method = "wilcoxon", p.adjust.method = "BH", hide.ns = TRUE, hide_info = TRUE)
    save_plot(p, filename = "result/301.key_gene_analy/deconv.png", view_plot = FALSE)
    # join data_full and deconv_res, by Sample
    # Pivot deconv_res wider to get one column per cell type
    deconv_wide <- deconv_res_full %>%
        select(-type) %>%
        pivot_wider(
            id_cols = Sample,
            names_from = celltype,
            values_from = cell_fraction,
            names_prefix = "celltype_"
        )

    # Convert Sample names to match data_full colnames
    coldata <- colData(data_full) %>%
        as.data.frame() %>%
        rownames_to_column(".sample")

    # Join deconv results with colData
    new_coldata <- coldata %>%
        left_join(deconv_wide, by = c(".sample" = "Sample")) %>%
        column_to_rownames(".sample")

    # Update colData in data_full
    colData(data_full) <- DataFrame(new_coldata)

    list(deconv_res_full, data_full)
}

evaluate_gene_celltype_correlation <- function(deconv_res_list, kgene) {
    # kgene <- "GBP5"
    data_full <- deconv_res_list[[2]]
    # Extract gene expression data
    gene_expr <- assay(data_full)[kgene, ]

    # Extract cell type proportions and sample type
    metadata <- colData(data_full) %>%
        as.data.frame() %>%
        select(group, matches("^celltype_"))

    # Calculate correlations separately for tumor and normal
    cor_results <- lapply(c("tumor", "normal"), function(sample_type) {
        # Filter for specific sample type
        sample_idx <- metadata$group == sample_type
        sample_expr <- gene_expr[sample_idx]
        sample_props <- metadata[sample_idx, -1] # Exclude group column

        # Calculate correlations for this sample type
        sample_cors <- apply(sample_props, 2, function(x) {
            cor.test(sample_expr, x, method = "spearman", exact = FALSE)
        })

        # Format results
        data.frame(
            celltype = names(sample_cors),
            correlation = sapply(sample_cors, function(x) x$estimate),
            pvalue = sapply(sample_cors, function(x) x$p.value),
            sample_type = sample_type
        )
    })

    # Combine results
    cor_df <- do.call(rbind, cor_results) %>%
        mutate(
            celltype = gsub("celltype_", "", celltype),
            padj = p.adjust(pvalue, method = "BH"),
            significance = case_when(
                padj < 0.001 ~ "***",
                padj < 0.01 ~ "**",
                padj < 0.05 ~ "*",
                TRUE ~ ""
            )
        )

    # Create correlation plot
    p <- ggplot(cor_df, aes(x = celltype, y = correlation, fill = sample_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = significance), position = position_dodge(width = 0.9), vjust = -0.5) +
        scale_fill_manual(values = c("normal" = "turquoise", "tumor" = "salmon")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(plot.background = element_rect(fill = "white")) +
        labs(
            title = paste0("Correlation between ", kgene, " expression and cell type proportions"),
            x = "Cell Type",
            y = "Spearman Correlation"
        )

    save_plot(p, filename = paste0("result/301.key_gene_analy/", kgene, "/celltype_correlation.png"))

    # 添加分组分析
    # Only analyze tumor samples
    tumor_idx <- metadata$group == "tumor"
    tumor_expr <- gene_expr[tumor_idx]
    tumor_metadata <- metadata[tumor_idx, ]

    quartile_analysis <- lapply(unique(cor_df$celltype), function(ct) {
        # Get cell proportions for current cell type (tumor only)
        cell_prop <- tumor_metadata[[paste0("celltype_", ct)]]

        # Calculate quartiles
        q25 <- quantile(cell_prop, 0.25)
        q75 <- quantile(cell_prop, 0.75)

        # Group samples
        high_group <- cell_prop >= q75
        low_group <- cell_prop <= q25

        # Get gene expression for each group
        expr_high <- tumor_expr[high_group]
        expr_low <- tumor_expr[low_group]

        # Wilcoxon test
        test_result <- wilcox.test(expr_high, expr_low)

        # Create plot data
        plot_data <- data.frame(
            expression = c(expr_high, expr_low),
            group = factor(
                c(
                    rep("High", sum(high_group)),
                    rep("Low", sum(low_group))
                ),
                levels = c("Low", "High")
            ),
            celltype = ct
        )

        list(
            celltype = ct,
            wilcox_pvalue = test_result$p.value,
            fold_change = mean(expr_high) / mean(expr_low),
            plot_data = plot_data
        )
    })

    # 整合分析结果
    quartile_df <- do.call(rbind, lapply(quartile_analysis, function(x) {
        data.frame(
            celltype = x$celltype,
            wilcox_pvalue = x$wilcox_pvalue,
            fold_change = x$fold_change
        )
    })) %>%
        # drop celltype containing "Epi"
        filter(!grepl("Epi", celltype)) %>%
    mutate(
        wilcox_pvalue_adj = p.adjust(wilcox_pvalue, method = "BH")
    )

    # Create dot plot
    p_quartile <- ggplot(quartile_df, aes(x = log2(fold_change), y = -log10(wilcox_pvalue_adj))) +
        geom_point(aes(color = case_when(
            wilcox_pvalue_adj <= 0.05 & log2(fold_change) >= 1 ~ "Up",
            wilcox_pvalue_adj <= 0.05 & log2(fold_change) <= -1 ~ "Down",
            TRUE ~ "NS"
        )), size = 3, show.legend = FALSE) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
        ggrepel::geom_text_repel(aes(label = celltype)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "red") +
        theme_minimal() +
        labs(
            y = "-log10(adjusted p-value)", 
            x = "log2(Fold Change)",
            title = paste0(kgene, " expression in high vs low cell proportion groups")
        ) +
        theme(plot.background = element_rect(fill = "white"))

    save_plot(p_quartile, filename = paste0("result/301.key_gene_analy/", kgene, "/celltype_compare_gene.png"))
    return(cor_df)
}
