deconv <- function(data_filt, data, data_risk_score) {
        samplesTP <- TCGAquery_SampleTypes(
            barcode = colnames(data_filt),
            typesample = c("TP")
        )
    data_cell <- data_filt %>%
        filter(.sample %in% samplesTP) %>%
        deconvolve_cellularity(.abundance = counts_scaled_adjusted, action = "get", cores = 30, prefix = "cibersort__") %>%
        pivot_sample() %>%
        left_join(data_risk_score, by = c(".sample" = ".sample")) %>%
        mutate(risk = ifelse(risk_score > median(risk_score, na.rm = TRUE), "high", "low"))
    data_cell_only <- data_cell %>%
        dplyr::select(.sample, starts_with("cibersort__")) %>%
        pivot_longer(
            names_to = "Cell_type_inferred",
            values_to = "proportion",
            names_prefix = "cibersort__",
            cols = contains("cibersort__")
        )
    p <- data_cell_only %>%
        ggplot(aes(x = Cell_type_inferred, y = proportion)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("result/107.deconv/cellularity.png", p)
    # cell_dds <- test_differential_cellularity(data_cell, . ~ risk, .abundance = counts_scaled_adjusted)
    data_cell_prep <- data_cell %>%
        dplyr::select(starts_with("cibersort__"), risk, risk_score) %>%
        # convert risk to numeric, low is 1, high is 2
        mutate(risk = ifelse(risk == "low", 1, 2)) %>%
        # round columns starting with cibersort__ to 5 digits
        mutate_at(vars(starts_with("cibersort__")), ~ round(., 5)) %>%
        # remove cibersort__ prefix
        rename_with(~ gsub("cibersort__", "", .x))
    cell_cor <- data_cell_prep %>%
        binarize(n_bins = 6, thresh_infreq = 0.01)
    cell_cor_res <- cell_cor %>%
        correlationfunnel::correlate(target = risk__1, method = "kendall")
    p <- plot_correlation_funnel(cell_cor_res)
    ggsave("result/107.deconv/cellularity_correlation_funnel.png", p, width = 10, height = 10)
    cor <- corrr::correlate(data_cell_prep, method = "kendall")
    p <- corrr::rplot(cor) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("result/107.deconv/cellularity_correlation.png", p, width = 14)

    data_cell_long <- data_cell %>%
        dplyr::select(starts_with("cibersort__"), risk, risk_score) %>%
        rename_with(~ gsub("cibersort__", "", .x)) %>%
        pivot_longer(
            cols = -c(risk, risk_score),
            names_to = "Cell_type_inferred",
            values_to = "proportion"
        ) %>%
        # add a column as the broad cell type, which is the first part of the cell type (sep by ".")
        mutate(broad_cell_type = str_split(Cell_type_inferred, "\\.") %>% map_chr(~ .x[1]))
    # p <- grouped_ggbetweenstats(
    #     data = data_cell_long,
    #     x = risk,
    #     y = proportion,
    #     grouping.var = Cell_type_inferred,
    #     type = "p"
    # )
    broad_cell_types <- data_cell_long$broad_cell_type %>% unique()
    walk(broad_cell_types, ~ {
        cell_type <- .x
        # fig width is 5 * number of cell types
        fig_width <- data_cell_long %>%
            filter(broad_cell_type == cell_type) %>%
            pull(Cell_type_inferred) %>%
            unique() %>%
            length() * 5
        fig_height <- 7
        if (cell_type == "T") {
            fig_width <- 15
            fig_height <- 15
        }
        data_cell_long %>%
            filter(broad_cell_type == cell_type) %>%
            grouped_ggbetweenstats(
                x = risk,
                y = proportion,
                grouping.var = Cell_type_inferred,
                type = "p"
            ) %>%
            ggsave(paste0("result/107.deconv/cellularity_", cell_type, "_test.png"), ., width = fig_width, height = fig_height, device = grDevices::png)
    })
    "result/107.deconv"
    # ggsave("result/107.deconv/cellularity_test.png", p, width = 20, height = 20)
    # ggsave("result/107.deconv/cellularity_test.pdf", p, width = 20, height = 20)
}
