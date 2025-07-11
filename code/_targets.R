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
source("code/R/204.cluster.R")
source("code/R/205.sub_cluster_temp.R")
source("code/R/206.Tex.R")
source("code/R/301.key_gene_analy.R")
log_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
log_dir <- paste0("logs/", log_timestamp)
tar_option_set(
    error = "trim",
    packages = c(
        "tidyverse",
        # "TCGAbiolinks",
        # "SummarizedExperiment",
        # "tidySummarizedExperiment",
        # "clusterProfiler",
        # "org.Hs.eg.db",
        # "pathview",
        # "enrichplot",
        # "DOSE",
        # "WGCNA",
        # "ggstatsplot",
        # "pheatmap",
        "patchwork",
        # "igraph",
        # "limma",
        # "tidybulk",
        # "DESeq2",
        # "tidygraph",
        # "ggraph",
        # "genekitr",
        # "survival",
        # "survminer",
        # "psych",
        # "tidyheatmaps",
        "furrr",
        "progressr",
        # "glmnet",
        # "msigdb",
        # "ggstatsplot",
        # "corrr",
        # "EnhancedVolcano",
        # "ggupset",
        # "writexl",
        "tidyseurat",
        # "SeuratDisk",
        "Seurat",
        "SeuratExtend",
        # "tidyplots",
        "scCustomize",
        "scDblFinder",
        "SingleCellExperiment",
        # "omnideconv",
        "reticulate",
        "schard"
    ),
    controller = crew_controller_local(
        workers = 8,
        seconds_timeout = 10800,
        options_local = crew_options_local(log_directory = log_dir)
    ),
    format = "auto",
    resources = tar_resources(
        qs = tar_resources_qs(nthreads = 8L)
    )
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
    tar_target(
        cox_res,
        get_cox(data, data_filt, data_dds),
        deployment = "main"
    ),
    tar_target(uni_EA_res, uni_enrich(cox_res, data_dds, msigdb)),
    tar_target(uni_EA_plot, plot_uni_enrich(uni_EA_res, data_dds, cox_res)),
    tar_target(
        uni_EA_intersect,
        select_intersect_genes(uni_EA_res, data_dds, cox_res)
    ),
    tar_target(lasso_res, get_lasso(data, data_filt, data_dds)),
    tar_target(data_risk_score, get_risk_score(lasso_res, data_filt, data)),
    tar_target(data_risk_dds, run_risk_DEG(data_filt, data_risk_score)),
    tar_target(data_risk_EA, run_risk_enrich(data_risk_dds, msigdb)),
    tar_target(
        data_risk_dds_tumor,
        run_risk_DEG_tumor(data_filt, data_risk_score)
    ),
    tar_target(
        data_risk_EA_tumor,
        run_risk_enrich_tumor(data_risk_dds_tumor, msigdb)
    ),
    tar_target(data_clinical_risk, get_corr_clinical(data, data_risk_score)),
    tar_target(
        deconv_res,
        deconv(data_filt, data, data_risk_score),
        format = "file"
    ),
    tar_target(
        trait_res,
        get_module_trait(WGCNA_res, data_filt, data, data_EA_tidy),
        format = "file"
    ),
    # Single cell - Raw Data Loading and Merging
    tar_target(sc_kourtis, load_kourtis("data/201.load_sc/raw")),
    tar_target(sc_krishna, load_krishna("data/201.load_sc/raw")),
    tar_target(sc_li, load_li("data/201.load_sc/raw")),
    tar_target(sc_obradovic, load_obradovic("data/201.load_sc/raw")),
    tar_target(sc_saout, load_saout("data/201.load_sc/raw")),
    tar_target(sc_yu, load_yu("data/201.load_sc/raw")),
    tar_target(
        sc_merged_raw,
        merge_sc_objects(
            kourtis = sc_kourtis,
            krishna = sc_krishna,
            li = sc_li,
            obradovic = sc_obradovic,
            saout = sc_saout,
            yu = sc_yu
        )
    ),
    tar_target(
        sc_merged_qc,
        qc_merged_sc(sc_merged_raw)
    ),
    tar_target(
        sc_merged_no_doublets,
        remove_doublets(sc_merged_qc),
        deployment = "main",
        packages = c("scDblFinder")
    ),
    tar_target(
        sc_merged_h5ad,
        save_merged_h5ad(
            sc_merged_no_doublets,
            "data/201.load_sc/merged_sc_qc.h5ad"
        ),
        format = "file",
        packages = c("reticulate", tar_option_get("packages"))
    ),

    # Single cell - Integration
    tar_target(
        python_integration_script,
        "code/Python/scvi_integrate.py",
        format = "file"
    ),
    tar_target(
        sc_integrated_h5ad,
        run_scvi_integration(
            python_script_path = python_integration_script,
            input_h5ad = sc_merged_h5ad,
            output_h5ad = "data/202.annotation/integrated.h5ad",
            host_project_root = "/data0/work/guozhonghao/Trp-RCC",
            host_venv_name = ".venv_scvi_gpu",
            host_user = "guozhonghao", # Your username on the host machine
            host_address = "host.docker.internal" # Or your host's IP
        ),
        format = "file",
        deployment = "main" # Run on main process due to venv
    ),
    tar_target(
        sc_integrated,
        load_integrated_h5ad(sc_integrated_h5ad)
    )

    # (Deprecated) Single cell - Integration and Downstream Analysis
    # tar_target(
    #     h5ad_path,
    #     "data/101.raw_data/after_integration.h5ad",
    #     format = "file",
    #     cue = tar_cue("never")
    # ),
    # tar_target(
    #     h5seurat_path,
    #     "data/101.raw_data/sce_pca.h5seurat",
    #     format = "file",
    #     cue = tar_cue("never")
    # ),
    # tar_target(sc_pre, load_sc_pre(h5seurat_path), cue = tar_cue("never")),
    # tar_target(latent, run_integration(), format = "file"),
    # tar_target(sc, import_integration(latent, sc_pre)),
    # tar_target(sc_pro, preprocess_sc(sc)),
    # tar_target(predicted_labels_path, run_annotation_train(), format = "file"),
    # tar_target(sc_anno, combine_annotation(sc, predicted_labels_path)),
    # tar_target(sc_cluster, run_clusters(sc_anno)),
    # tar_target(cluster_tree_path, cluster_tree(sc_cluster)),
    # tar_target(
    #     final_annotation,
    #     optimize_clusters(sc_cluster),
    #     format = "file"
    # ),
    # tar_target(sc_opt, get_final_annotation(sc_cluster, final_annotation)),
    # # tar_target(sub_res, c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5)),
    # # tar_target(T_sub_cluster, get_sub_cluster(sc_opt, sub_res), pattern = map(sub_res)),
    # tar_target(
    #     predicted_labels_T,
    #     "data/205.sub_cluster/predicted_labels_T.csv",
    #     format = "file"
    # ),
    # tar_target(sc_T_sub_cluster, get_sub_cluster(sc_opt)),
    # tar_target(sc_T_anno, anno_Tex(sc_T_sub_cluster, predicted_labels_T)),
    # tar_target(Tex_score, get_module_score_T(sc_T_anno)),
    # tar_target(Tex_marker, DEG_Tex(sc_T_anno, sc_opt, Tex_score)),
    # # key gene analysis
    # tar_target(
    #     kgenes,
    #     c("GBP5", "CD27", "ANXA7", "HAPLN3", "CHI3L2", "NCKAP5", "FIBCD1")
    # ),
    # tar_target(data_full, get_data_full(data, data_filt)),
    # tar_target(
    #     kg_DEG_single_plot,
    #     kg_DEG_single(data_full, data_dds, kgenes),
    #     pattern = map(kgenes),
    #     iteration = "list",
    #     error = "null"
    # ),
    # tar_target(
    #     kg_survival_plot,
    #     kg_survival(data_full, kgenes),
    #     pattern = map(kgenes),
    #     iteration = "list",
    #     error = "null"
    # ),
    # tar_target(
    #     kg_sc_path,
    #     kg_sc(sc_opt, kgenes),
    #     pattern = map(kgenes),
    #     iteration = "list",
    #     error = "null"
    # ),
    # tar_target(
    #     ref_deconv,
    #     schard::h5ad2seurat("data/201.load_sc/reference.h5ad")
    # ),
    # tar_target(kg_deconv_tumor_res, kg_deconv_tumor(data_full, ref_deconv)),
    # tar_target(kg_deconv_normal_res, kg_deconv_normal(data_full, ref_deconv)),
    # tar_target(
    #     deconv_res_list,
    #     plot_deconv(kg_deconv_tumor_res, kg_deconv_normal_res, data_full)
    # ),
    # tar_target(
    #     kg_deconv_gene_res,
    #     evaluate_gene_celltype_correlation(deconv_res_list, kgenes),
    #     pattern = map(kgenes),
    #     iteration = "list",
    #     error = "null"
    # )
)
