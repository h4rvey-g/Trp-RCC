library(targets)
library(crew)
source("code/R/101.bulk_flow.R")
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "TCGAbiolinks", "SummarizedExperiment", "tidySummarizedExperiment", "clusterProfiler", "org.Hs.eg.db", "pathview",
        "enrichplot", "DOSE", "WGCNA"
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
    tar_target(data_filt, preprocess_data(data)),
    tar_target(DEG_res, run_DEG(data_filt)),
    tar_target(ansEA, run_enrich(DEG_res), deployment = "main")
)
