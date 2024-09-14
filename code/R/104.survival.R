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

get_cox <- function(data, data_filt, data_dds) {
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data),
        typesample = c("TP")
    )
    get_cox_res <- function(gene) {
        library(tidySummarizedExperiment)
        data_clin <- data_filt %>%
            as_tibble() %>%
            filter(.feature %in% gene) %>%
            filter(.sample %in% samplesTP) %>%
            left_join(., colData(data) %>% as.data.frame() %>% rownames_to_column(".sample"), by = ".sample")
        data_clin_final <- data_clin %>%
            dplyr::mutate(
                survival_time = ifelse(vital_status == "Dead", days_to_death,
                    days_to_last_follow_up
                ) / 365.25,
                vital_status = ifelse(vital_status == "Dead", 1, 0)
            ) %>%
            dplyr::select(.sample, survival_time, vital_status, counts_scaled_adjusted)
        formula <- as.formula(paste("Surv(survival_time, vital_status) ~", "counts_scaled_adjusted"))
        fit <- coxph(formula, data = data_clin_final)
        fit_summary <- summary(fit)
        p.value <- signif(fit_summary$wald["pvalue"], digits = 2)
        beta <- signif(fit_summary$coef[1], digits = 2)
        HR <- signif(fit_summary$coef[2], digits = 2)
        HR.confint.lower <- signif(fit_summary$conf.int[, "lower .95"], 2)
        HR.confint.upper <- signif(fit_summary$conf.int[, "upper .95"], 2)
        HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
        res <- c(beta, HR, p.value)
        names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
        res
    }
    options(future.globals.maxSize = 9000 * 1024^2)
    plan(multisession)
    dds_feature <- data_dds %>%
        pull(.feature) %>%
        unique()
    cox_list <- future_map(dds_feature, get_cox_res)
    # convert to tibble
    cox_res <- cox_list %>%
        purrr::reduce(rbind) %>%
        as_tibble() %>%
        mutate(gene = dds_feature) %>%
        filter(p.value < 0.05)
   cox_res 
}

get_lasso <- function(data, data_filt, data_dds) {
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data),
        typesample = c("TP")
    )
    gene <- data_dds %>%
        pull(.feature) %>%
        .[1]
    data_clin <- data_filt %>%
        as_tibble() %>%
        filter(.feature %in% gene) %>%
        filter(.sample %in% samplesTP) %>%
        left_join(., colData(data) %>% as.data.frame() %>% rownames_to_column(".sample"), by = ".sample")
    data_clin_final <- data_clin %>%
        dplyr::mutate(
            survival_time = ifelse(vital_status == "Dead", days_to_death,
                days_to_last_follow_up
            ) / 365.25,
            vital_status = ifelse(vital_status == "Dead", 1, 0)
        ) %>%
        dplyr::select(.sample, survival_time, vital_status, counts_scaled_adjusted)
    x <- assay(data_dds, "counts_scaled_adjusted") %>% as.matrix()
    y <- Surv(data_clin_final$survival_time, data_clin_final$vital_status) %>% as.matrix()
    fit <- cv.glmnet(
        x = x,
        y = y, family = "cox",
        nlambda = 100
    )
    pdf("result/106.survival/lasso.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()

    fit.cv <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "cox")
    pdf("result/106.survival/lasso_cv.pdf")
    plot(fit.cv)
    dev.off()

    feature_all <- as.data.frame(as.matrix(coef(fit.cv, s = "lambda.min")))
    feature_s <- rownames(feature_all)[feature_all != 0]
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
    data_filt %>% test_differential_cellularity(. ~ group)
}
