library(targets)
library(crew)
source("code/R/101.bulk_flow.R")
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "TCGAbiolinks", "SummarizedExperiment", "tidySummarizedExperiment", "clusterProfiler", "org.Hs.eg.db", "pathview",
        "enrichplot", "DOSE", "WGCNA", "ggstatsplot", "pheatmap", "patchwork", "igraph", "limma", "tidybulk", "DESeq2", "tidygraph",
        "ggraph", "genekitr", "survival", "survminer", "psych","tidyheatmaps"
    ),
    controller = crew_controller_local(workers = 2, seconds_timeout = 36000),
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
    tar_target(surv_res, get_survival(data, data_filt, data_EA_tidy), format = "file"),
    tar_target(trait_res, get_module_trait(WGCNA_res, data_filt, data, data_EA_tidy), format = "file")
)
