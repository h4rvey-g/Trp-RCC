get_corr_clinical <- function(data, data_risk_score){
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data),
        typesample = c("TP")
    )
    data_clinical <- data %>%
        as_tibble() %>%
        left_join(data_risk_score, by = c(".sample" = ".sample")) %>%
        filter(.sample %in% samplesTP) %>%
        dplyr::select(
            sample_type, intermediate_dimension, ajcc_pathologic_stage, last_known_disease_status,
            age_at_diagnosis, ajcc_pathologic_t, ajcc_pathologic_n, ajcc_pathologic_m, gender, age_at_index, risk_score, risk
        ) %>%
        distinct()
    # drop rows containing NA
    data_clinical <- data_clinical %>%
        drop_na()
    data_clinical_bin <- data_clinical %>%
        binarize(n_bins = 4, thresh_infreq = 0.01)
    data_cor <- data_clinical_bin %>%
        correlate(target = "risk_score__0.422483067225408_Inf")
    p <- plot_correlation_funnel(data_cor)
    ggsave("result/108.clinical/clinical_correlation_funnel.png", p)
}