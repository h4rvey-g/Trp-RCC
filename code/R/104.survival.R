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
}

get_cox <- function(data, data_filt, data_dds) {
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data),
        typesample = c("TP")
    )
    data_clin <- data_filt %>%
        as_tibble() %>%
        filter(.sample %in% samplesTP) %>%
        left_join(., colData(data) %>% as.data.frame() %>% rownames_to_column(".sample"), by = ".sample") %>%
        dplyr::select(.sample, .feature, counts_scaled_adjusted, days_to_death, days_to_last_follow_up, vital_status)
    get_cox_res <- function(gene, p) {
        p()
        data_clin <- data_clin %>%
            filter(.feature %in% gene)
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
    plan(multisession, workers = 30)
    dds_feature <- data_dds %>%
        pull(.feature) %>%
        unique()
    with_progress({
        p <- progressor(steps = length(dds_feature))
        cox_list <- future_map(dds_feature, get_cox_res, p = p)
    })
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
        pull(.feature)
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
        # remove survival time = 0
        dplyr::filter(survival_time > 0) %>%
        dplyr::select(.sample, survival_time, vital_status, counts_scaled_adjusted)
    x <- assay(
        data_dds %>%
            filter(.sample %in% data_clin_final$.sample),
        "counts_scaled_adjusted"
    ) %>%
        as.matrix() %>%
        t()
    y <- data_clin_final %>%
        dplyr::select(-counts_scaled_adjusted) %>%
        distinct() %>%
        {
            Surv(time = .$survival_time, event = .$vital_status)
        }
    fit.cv <- cv.glmnet(
        x = x,
        y = y, family = "cox",
        alpha = 1
    )
    pdf("result/106.survival/lasso_cv.pdf")
    plot(fit.cv, xvar = "lambda", label = TRUE)
    dev.off()

    fit <- glmnet(x, y, family = "cox", nlambda = 100)
    pdf("result/106.survival/lasso.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()

    feature_all <- as_tibble(as.matrix(coef(fit, s = fit.cv$lambda.min)), rownames = "feature") %>%
        filter(`1` != 0) %>%
        arrange(desc(abs(`1`))) %>%
        dplyr::rename(coef = `1`)

    p <- ggplot(
        feature_all, aes(x = reorder(feature, coef), y = coef, fill = factor(sign(coef)))
    ) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("blue", "red")) +
        labs(x = "Gene", y = "Coefficients") +
        theme(legend.position = "none")
    ggsave("result/106.survival/lasso.png", p)
    feature_all
}

get_risk_score <- function(lasso_res, data_filt, data) {
    gene <- lasso_res %>%
        pull(feature)
    data_gene <- data_filt %>%
        as_tibble() %>%
        filter(.feature %in% gene) %>%
        dplyr::select(.sample, .feature, counts_scaled_adjusted)
    # Merge data_risk_score and lasso_res by the gene feature
    merged_data <- left_join(data_gene, lasso_res, by = c(".feature" = "feature"))

    # Calculate the risk score for each row (counts_scaled_adjusted * coef)
    merged_data$risk_contrib <- merged_data$counts_scaled_adjusted * merged_data$coef
    # print the formula for risk score, coef is from lasso_res
    cat("Risk score = ", paste0(round(lasso_res$coef, 4), "*", lasso_res$feature, " + "), "\n")

    # Summarize risk score for each sample
    data_risk_score <- merged_data %>%
        group_by(.sample) %>%
        summarize(risk_score = sum(risk_contrib)) %>%
        # Add a new column to data_risk_score to indicate the risk group, if risk_score is above median, set as high, otherwise low
        mutate(risk = ifelse(risk_score > median(risk_score, na.rm = TRUE), "high", "low"))
    risk_score_vs_group <- data_filt %>%
        as_tibble() %>%
        dplyr::select(.sample, group) %>%
        distinct() %>%
        left_join(., data_risk_score, by = ".sample")
    # draw stacked barplot on group and risk, risk is the fill color
    p <- risk_score_vs_group %>%
        group_by(group, risk) %>%
        summarize(counts = n()) %>%
        ggplot(aes(x = group, y = counts, fill = risk)) +
        geom_col(position = "stack") +
        labs(x = "Group", y = "Counts", fill = "Risk") +
        theme(legend.position = "top")
    ggsave("result/106.survival/risk_score_vs_group.png", p)

    # do survival analysis on risk score
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data_filt),
        typesample = c("TP")
    )
    data_clin_final <- data_filt %>%
        as_tibble() %>%
        filter(.sample %in% samplesTP) %>%
        left_join(., colData(data) %>% as.data.frame() %>% rownames_to_column(".sample"), by = ".sample") %>%
        left_join(., data_risk_score, by = ".sample") %>%
        dplyr::mutate(
            survival_time = ifelse(vital_status == "Dead", days_to_death,
                days_to_last_follow_up
            ) / 365.25,
            vital_status = ifelse(vital_status == "Dead", 1, 0)
        ) %>%
        dplyr::select(.sample, survival_time, vital_status, risk)
    fit <- survfit(Surv(survival_time, vital_status) ~ risk, data = data_clin_final)
    p <- ggsurvplot(fit,
        data = data_clin_final, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        legend.title = "Risk score",
        legend.labs = c("Low", "High")
    )$plot
    ggsave("result/106.survival/risk_score_survival.png", p)
    data_risk_score
}

get_random_forest <- function(data, data_filt, data_dds) {
    samplesTP <- TCGAquery_SampleTypes(
        barcode = colnames(data),
        typesample = c("TP")
    )
    gene <- data_dds %>%
        pull(.feature)
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
        # remove survival time = 0
        dplyr::filter(survival_time > 0) %>%
        dplyr::select(.sample, survival_time, vital_status, counts_scaled_adjusted)
    x <- assay(
        data_dds %>%
            filter(.sample %in% data_clin_final$.sample),
        "counts_scaled_adjusted"
    ) %>%
        as.matrix() %>%
        t()
    y <- data_clin_final %>%
        dplyr::select(-counts_scaled_adjusted) %>%
        distinct() %>%
        {
            Surv(time = .$survival_time, event = .$vital_status)
        }
    rf <- randomForestSRC::rfsrc(
        x = x,
        y = y,
        ntree = 1000,
        nodedepth = 5,
        seed = 123
    )
    varimp <- rf$importance %>%
        as_tibble() %>%
        rownames_to_column("feature") %>%
        arrange(desc(importance))
    p <- ggplot(
        varimp, aes(x = reorder(feature, importance), y = importance)
    ) +
        geom_col() +
        coord_flip() +
        labs(x = "Gene", y = "Importance") +
        theme(legend.position = "none")
    ggsave("result/106.survival/random_forest.png", p)
    varimp
}
