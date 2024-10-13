library(targets)
library(crew)
source("code/R/101.preprocess.R")
source("code/R/102.DEG_and_enrich.R")
source("code/R/103.network.R")
source("code/R/104.survival.R")
source("code/R/105.risk_score_DEG.R")
source("code/R/106.clinical.R")
source("code/R/107.deconv.R")
source("code/R/108.direct_enrich.R")
source("code/R/201.load_sc.R")
source("code/R/202.annotation.R")
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "TCGAbiolinks", "SummarizedExperiment", "tidySummarizedExperiment", "clusterProfiler",
        "org.Hs.eg.db", "pathview", "enrichplot", "DOSE", "WGCNA", "ggstatsplot", "pheatmap", "patchwork", "igraph",
        "limma", "tidybulk", "DESeq2", "tidygraph", "ggraph", "genekitr", "survival", "survminer", "psych",
        "tidyheatmaps", "furrr", "progressr", "glmnet", "msigdb", "ggstatsplot", "correlationfunnel", "corrr",
        "EnhancedVolcano", "ggupset", "writexl", "tidyseurat", "SeuratDisk", "Seurat", "anndata"
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
    tar_target(msigdb, download_msigdb()),
    tar_target(gtex_data_path, download_gtex()),
    tar_target(data_filt, preprocess_data(data, gtex_data_path)),
    tar_target(data_dds, run_DEG(data_filt)),
    tar_target(ansEA, run_enrich(data_dds), deployment = "main"),
    tar_target(WGCNA_res, run_WGCNA(data_filt)),
    tar_target(data_EA_tidy, plot_WGCNA(WGCNA_res, ansEA, data_filt, data_dds)),
    tar_target(network_res, get_network(WGCNA_res, data_dds, data_EA_tidy)),
    tar_target(survival_res, get_survival(data, data_filt, data_EA_tidy)),
    tar_target(cox_res, get_cox(data, data_filt, data_dds), deployment = "main"),
    tar_target(uni_EA_res, uni_enrich(cox_res, data_dds, msigdb)),
    tar_target(uni_EA_plot, plot_uni_enrich(uni_EA_res, data_dds, cox_res)),
    tar_target(uni_EA_intersect, select_intersect_genes(uni_EA_res, data_dds, cox_res)),
    tar_target(lasso_res, get_lasso(data, data_filt, data_dds)),
    tar_target(data_risk_score, get_risk_score(lasso_res, data_filt, data)),
    tar_target(data_risk_dds, run_risk_DEG(data_filt, data_risk_score)),
    tar_target(data_risk_EA, run_risk_enrich(data_risk_dds, msigdb)),
    tar_target(data_risk_dds_tumor, run_risk_DEG_tumor(data_filt, data_risk_score)),
    tar_target(data_risk_EA_tumor, run_risk_enrich_tumor(data_risk_dds_tumor, msigdb)),
    tar_target(data_clinical_risk, get_corr_clinical(data, data_risk_score)),
    tar_target(deconv_res, deconv(data_filt, data, data_risk_score), format = "file"),
    tar_target(trait_res, get_module_trait(WGCNA_res, data_filt, data, data_EA_tidy), format = "file"),
    # Single cell
    tar_target(h5ad_path, "data/101.raw_data/after_integration.h5ad", format = "file"),
    tar_target(h5seurat_path, "data/101.raw_data/sce_pca.h5seurat", format = "file"),
    tar_target(sc_pre, load_sc_pre(h5seurat_path)), 
    tar_target(latent, run_integration(), format = "file"),
    tar_target(sc, import_integration(latent, sc_pre)),
    tar_target(sc_pro, preprocess_sc(sc)),
    tar_target(predicted_labels_path, run_annotation_train(), format = "file"),
    tar_target(sc_anno, combine_annotation(sc, predicted_labels_path))
)
