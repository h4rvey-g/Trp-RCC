library(targets)
library(crew)
source("code/R/101.preprocess.R")
source("code/R/102.DEG_and_enrich.R")
source("code/R/103.network.R")
source("code/R/104.survival.R")
source("code/R/105.risk_score_DEG.R")
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "TCGAbiolinks", "SummarizedExperiment", "tidySummarizedExperiment", "clusterProfiler", "org.Hs.eg.db", "pathview",
        "enrichplot", "DOSE", "WGCNA", "ggstatsplot", "pheatmap", "patchwork", "igraph", "limma", "tidybulk", "DESeq2", "tidygraph",
        "ggraph", "genekitr", "survival", "survminer", "psych", "tidyheatmaps", "furrr", "progressr", "glmnet"
    ),
    controller = crew_controller_local(workers = 20, seconds_timeout = 36000),
    format = "qs",
    storage = "worker", retrieval = "worker"
)
tar_config_set(
    script = "code/_targets.R",
    store = "data/_targets"
)
list(
    tar_target(data, download_data()),
    tar_target(gtex_data_path, download_gtex()),
    tar_target(data_filt, preprocess_data(data, gtex_data_path)),
    tar_target(data_dds, run_DEG(data_filt)),
    tar_target(ansEA, run_enrich(data_dds), deployment = "main"),
    tar_target(WGCNA_res, run_WGCNA(data_filt)),
    tar_target(data_EA_tidy, plot_WGCNA(WGCNA_res, ansEA, data_filt, data_dds)),
    tar_target(network_res, get_network(WGCNA_res, data_dds, data_EA_tidy)),
    tar_target(survival_res, get_survival(data, data_filt, data_EA_tidy)),
    tar_target(cox_res, get_cox(data, data_filt, data_dds), deployment = "main"),
    tar_target(lasso_res, get_lasso(data, data_filt, data_dds)),
    tar_target(data_risk_score, get_risk_score(lasso_res, data_filt, data)),
    tar_target(data_risk_dds, run_risk_DEG(data_filt, data_risk_score)),
    tar_target(data_risk_EA, run_risk_enrich(data_risk_dds)),
    tar_target(trait_res, get_module_trait(WGCNA_res, data_filt, data, data_EA_tidy), format = "file")
)
