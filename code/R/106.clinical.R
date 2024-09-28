get_corr_clinical <- function(data, data_risk_score) {
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data),
        typesample = c("TP")
    )
    data_clinical <- data %>%
        as_tibble() %>%
        left_join(data_risk_score, by = c(".sample" = ".sample")) %>%
        filter(.sample %in% samplesTP) %>%
        dplyr::select(
            intermediate_dimension, ajcc_pathologic_stage, tumor_grade,
            age_at_diagnosis, ajcc_pathologic_t, ajcc_pathologic_n, ajcc_pathologic_m, gender, age_at_index, risk_score, risk
        ) %>%
        distinct() %>%
        mutate(across(contains("ajcc_pathologic") | contains("gender") | contains("risk"), ~ as.numeric(as.factor(.)))) %>%
        drop_na()

    data_clinical_bin <- data_clinical %>%
        binarize(n_bins = 4, thresh_infreq = 0.01)
    # get the colname in data_clinical_bin which contains both "risk_score" and "_Inf"
    max_risk_score <- colnames(data_clinical_bin) %>%
        str_subset("risk_score") %>%
        str_subset("_Inf") 
    data_cor <- data_clinical_bin %>%
        correlationfunnel::correlate(target = max_risk_score, method = "kendall")
    p <- plot_correlation_funnel(data_cor)
    ggsave("result/108.clinical/clinical_correlation_funnel.png", p)

    # use data_clinical, test the correlation between risk_score and ajcc_pathologic_stage (convert to factor)
    # cor <- cor.test(data_clinical$risk_score, data_clinical$ajcc_pathologic_stage %>%
    #     as.factor() %>% as.numeric(), method = "kendall")
    cor <- corrr::correlate(data_clinical, method = "kendall")
    p <- corrr::rplot(cor)
    ggsave("result/108.clinical/clinical_correlation.png", p, width = 14)

    # pivot_longer
    data_clinical_long <- data %>%
        as_tibble() %>%
        left_join(data_risk_score, by = c(".sample" = ".sample")) %>%
        filter(.sample %in% samplesTP) %>%
        dplyr::select(
            ajcc_pathologic_stage, 
            ajcc_pathologic_t, ajcc_pathologic_n, ajcc_pathologic_m, risk_score
        ) %>%
        distinct() %>%
        pivot_longer(
            cols = -c(risk_score),
            names_to = "feature",
            values_to = "value"
        )
    p <- grouped_ggbetweenstats(
        data = data_clinical_long,
        x = value,
        y = risk_score,
        grouping.var = feature,
        type = "np"
    )
    ggsave("result/108.clinical/clinical_test.png", p, width = 20, height = 20)
    data_clinical_risk <- data_clinical_long
}
