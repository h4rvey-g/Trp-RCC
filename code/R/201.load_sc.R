# ------------------------------------------------------------------------------
# 201.load_sc.R - Functions for Single-Cell Data Loading
#
# This script provides a set of functions to load 6 different single-cell
# RNA-seq datasets. It is designed to be sourced by a `targets` pipeline.
#
# Contains:
# - Individual functions to load each of the 6 datasets.
# - A main function `load_and_merge_sc` that orchestrates the loading,
#   merging, and saving of the final Seurat object.
# ------------------------------------------------------------------------------

#' Load Kourtis et al. 2022 Dataset
#' @param data_path Base path to the raw data directory.
#' @return A Seurat object.
load_kourtis <- function(data_path) {
    message("Loading Kourtis dataset..")
    data_dir <- file.path(data_path, "kourtisSinglecellMapDynamic2022")

    matrix_file <- file.path(
        data_dir,
        "GSE181061_ccRCC_4pt_scRNAseq_CD45plus_matrix.mtx.gz"
    )
    features_file <- file.path(
        data_dir,
        "GSE181061_ccRCC_4pt_scRNAseq_CD45plus_genes.tsv.gz"
    )
    barcodes_file <- file.path(
        data_dir,
        "GSE181061_ccRCC_4pt_scRNAseq_CD45plus_barcodes.tsv.gz"
    )

    counts <- ReadMtx(
        mtx = matrix_file,
        features = features_file,
        cells = barcodes_file,
        feature.column = 2
    )

    meta_path <- file.path(
        data_dir,
        "GSE181061_ccRCC_4pt_scRNAseq_Tcells_final_Metadata.txt.gz"
    )
    meta <- read.delim(meta_path, row.names = 1)

    common_barcodes <- intersect(colnames(counts), rownames(meta))
    counts <- counts[, common_barcodes]
    meta <- meta[common_barcodes, , drop = FALSE]

    obj <- CreateSeuratObject(
        counts = counts,
        project = "Kourtis",
        meta.data = meta
    )
    obj$dataset <- "kourtis"
    obj$sample_id <- "kourtis_GSE181061"
    return(obj)
}

#' Load Krishna et al. 2021 Dataset
#' @param data_path Base path to the raw data directory.
#' @return A Seurat object.
load_krishna <- function(data_path) {
    message("Loading Krishna dataset..")
    data_dir <- file.path(data_path, "krishnaSinglecellSequencingLinks2021")
    obj <- readRDS(file.path(data_dir, "ccRCC_6pat_Seurat.rds"))

    # The loaded object is an old Seurat v2 object. Update it to the current version.
    if (packageVersion("Seurat") >= "3.0.0" && class(obj) == "seurat") {
        message("Updating Seurat object from a previous version...")
        obj <- UpdateSeuratObject(obj)
    }

    annot <- read.table(
        file.path(data_dir, "ccRCC_6pat_cell_annotations.txt"),
        header = TRUE,
        sep = "\t",
        row.names = 1
    )

    # Ensure the annotation rows match the cell names in the object
    common_cells <- intersect(colnames(obj), rownames(annot))
    annot_filtered <- annot[common_cells, , drop = FALSE]
    obj <- AddMetaData(obj, metadata = annot_filtered)

    obj$dataset <- "krishna"
    # Use the standard 'orig.ident' column, which should contain the patient info,
    # instead of the non-existent 'patient' column.
    if ("orig.ident" %in% colnames(obj@meta.data)) {
        obj$sample_id <- paste("krishna", obj$orig.ident, sep = "_")
    } else {
        message(
            "'orig.ident' metadata not found for Krishna dataset, using a default sample_id."
        )
        obj$sample_id <- "krishna_dataset"
    }
    return(obj)
}

#' Load Li et al. 2022 Dataset
#' @param data_path Base path to the raw data directory.
#' @return A Seurat object.
load_li <- function(data_path) {
    message("Loading Li dataset..")
    h5ad_file <- file.path(
        data_path,
        "liMappingSinglecellTranscriptomes2022",
        "RCC_upload_final_raw_counts.h5ad"
    )

    if (!requireNamespace("schard", quietly = TRUE)) {
        stop("Package 'schard' is required. Please install it.")
    }

    obj <- schard::h5ad2seurat(h5ad_file)

    obj$dataset <- "li"
    if ("batch" %in% colnames(obj@meta.data)) {
        obj$sample_id <- paste("li", obj$batch, sep = "_")
    } else {
        message(
            "'batch' metadata not found for Li dataset, using a default sample_id."
        )
        obj$sample_id <- "li_dataset"
    }
    return(obj)
}

#' Load Obradovic et al. 2021 Dataset
#' @param data_path Base path to the raw data directory.
#' @return A merged Seurat object.
load_obradovic <- function(data_path) {
    message("Loading Obradovic dataset..")
    base_dir <- file.path(data_path, "obradovicSinglecellProteinActivity2021")
    all_dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
    matrix_dirs <- grep(
        "filtered_(gene|feature)_bc_matrices(/GRCh38)?$|filtered_feature_bc_matrix$",
        all_dirs,
        value = TRUE
    )

    obj_list <- lapply(matrix_dirs, function(d) {
        path_parts <- strsplit(d, .Platform$file.sep)[]
        sample_part <- grep("^CN[0-9]+b?$", path_parts, value = TRUE)

        if (length(sample_part) == 1) {
            sample_name <- sample_part
        } else {
            # Fallback for paths that don't match the CN pattern
            sample_name <- basename(dirname(dirname(dirname(d))))
        }

        # Robustly find file paths, checking for both gzipped and unzipped versions
        barcodes_f <- file.path(d, "barcodes.tsv.gz")
        features_f <- file.path(d, "features.tsv.gz")
        matrix_f <- file.path(d, "matrix.mtx.gz")

        if (!file.exists(features_f)) {
            features_f <- file.path(d, "genes.tsv.gz")
        }

        if (!file.exists(barcodes_f)) {
            barcodes_f <- file.path(d, "barcodes.tsv")
            features_f <- file.path(d, "features.tsv")
            if (!file.exists(features_f)) {
                features_f <- file.path(d, "genes.tsv")
            }
            matrix_f <- file.path(d, "matrix.mtx")
        }

        if (!all(file.exists(barcodes_f, features_f, matrix_f))) {
            warning(paste(
                "Skipping directory:",
                d,
                "- required files not found."
            ))
            return(NULL)
        }

        counts <- ReadMtx(
            mtx = matrix_f,
            cells = barcodes_f,
            features = features_f,
            feature.column = 2
        )
        obj <- CreateSeuratObject(counts = counts, project = "Obradovic")
        obj$sample_id <- paste("obradovic", sample_name, sep = "_")
        # Rename cells before merging to avoid clashes
        obj <- RenameCells(obj, add.cell.id = sample_name)
        return(obj)
    })

    obj_list <- obj_list[!sapply(obj_list, is.null)]
    if (length(obj_list) == 0) {
        return(NULL)
    }

    # Correctly merge the list of Seurat objects
    merged_obj <- if (length(obj_list) > 1) {
        merge(x = obj_list[[1]], y = obj_list[-1])
    } else {
        obj_list[]
    }
    merged_obj$dataset <- "obradovic"
    return(merged_obj)
}

#' Load Saout et al. 2023 Dataset
#' @param data_path Base path to the raw data directory.
#' @return A merged Seurat object.
load_saout <- function(data_path) {
    message("Loading Saout dataset..")
    data_dir <- file.path(data_path, "saoutSinglecellDeconvolutionSpecific2023")
    barcodes_files <- list.files(
        data_dir,
        pattern = "_barcodes\\.tsv\\.gz$",
        full.names = TRUE
    )
    samples <- gsub("_barcodes\\.tsv\\.gz$", "", basename(barcodes_files))

    obj_list <- lapply(samples, function(s) {
        matrix_file <- file.path(data_dir, paste0(s, "_matrix.mtx.gz"))
        features_file <- file.path(data_dir, paste0(s, "_features.tsv.gz"))
        barcodes_file <- file.path(data_dir, paste0(s, "_barcodes.tsv.gz"))

        if (!file.exists(features_file)) {
            features_file <- file.path(data_dir, paste0(s, "_genes.tsv.gz"))
        }
        if (!all(file.exists(matrix_file, features_file, barcodes_file))) {
            return(NULL)
        }

        counts <- ReadMtx(
            mtx = matrix_file,
            features = features_file,
            cells = barcodes_file,
            feature.column = 2
        )
        obj <- CreateSeuratObject(counts = counts, project = "Saout")
        obj$sample_id <- paste("saout", s, sep = "_")
        # Rename cells before merging to avoid clashes
        obj <- RenameCells(obj, add.cell.id = s)
        return(obj)
    })

    obj_list <- obj_list[!sapply(obj_list, is.null)]
    if (length(obj_list) == 0) {
        return(NULL)
    }

    # Correctly merge the list of Seurat objects
    merged_obj <- if (length(obj_list) > 1) {
        merge(x = obj_list[[1]], y = obj_list[-1])
    } else {
        obj_list[]
    }
    merged_obj$dataset <- "saout"
    return(merged_obj)
}

#' Load Yu et al. 2023 Dataset
#' @param data_path Base path to the raw data directory.
#' @return A merged Seurat object.
load_yu <- function(data_path) {
    message("Loading Yu dataset..")
    data_dir <- file.path(data_path, "yuIntegrativeSingleCellAnalysis2023")
    sample_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
    sample_dirs <- sample_dirs[grepl("RCC", basename(sample_dirs))]

    obj_list <- lapply(sample_dirs, function(d) {
        counts <- Read10X(data.dir = d)
        obj <- CreateSeuratObject(counts = counts, project = "Yu")
        obj$sample_id <- paste("yu", basename(d), sep = "_")
        # Rename cells before merging to avoid clashes
        obj <- RenameCells(obj, add.cell.id = basename(d))
        return(obj)
    })

    obj_list <- obj_list[!sapply(obj_list, is.null)]
    if (length(obj_list) == 0) {
        return(NULL)
    }

    # Correctly merge the list of Seurat objects
    merged_obj <- if (length(obj_list) > 1) {
        merge(x = obj_list[[1]], y = obj_list[-1])
    } else {
        obj_list[]
    }
    merged_obj$dataset <- "yu"
    return(merged_obj)
}

#' Merge multiple Seurat objects
#'
#' This function takes multiple Seurat objects, merges them, and saves the result.
#' It's designed to be called from a `targets` pipeline, where each input object
#' is a dependency.
#'
#' @param output_file The path to save the final merged RDS file.
#' @param ... Named Seurat objects to be merged. The names will be used as prefixes
#'   for cell barcodes.
#' @return The path to the saved RDS file.
merge_sc_objects <- function(...) {
    all_data_list <- list(...)
    all_data_list <- all_data_list[!sapply(all_data_list, is.null)]
    if (length(all_data_list) < 1) {
        stop("No valid Seurat objects were provided to merge.")
    }

    message("Merging all datasets into a single Seurat object..")

    merged_seurat_obj <- if (length(all_data_list) > 1) {
        merge(
            x = all_data_list[[1]],
            y = all_data_list[-1],
            add.cell.ids = names(all_data_list)
        )
    } else {
        # If only one object, just add the prefix manually
        obj <- all_data_list[]
        obj <- RenameCells(
            obj,
            add.cell.id = names(all_data_list)
        )
        obj
    }
    return(merged_seurat_obj)
}

#' Perform QC on the merged Seurat object
#'
#' @param merged_obj The merged Seurat object.
#' @return A Seurat object after QC.
qc_merged_sc <- function(merged_obj) {
    message("Performing QC on merged Seurat object..")

    # Set up output directory
    plot_dir <- "result/201.load_sc"
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

    # Calculate mitochondrial gene percentage
    merged_obj[["percent.mt"]] <- PercentageFeatureSet(
        merged_obj,
        pattern = "^MT-"
    )

    # Create before QC plots
    message("Creating before QC plots...")
    before_plot <- VlnPlot(
        merged_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.1
    ) +
        plot_annotation(title = "Before QC")

    ggsave(
        filename = file.path(plot_dir, "before_qc_violin.png"),
        plot = before_plot,
        width = 12,
        height = 6,
        dpi = 300
    )

    # Apply QC filters
    merged_obj_qc <- subset(
        merged_obj,
        subset = nFeature_RNA > 200 & percent.mt < 20
    )

    # Create after QC plots
    message("Creating after QC plots...")
    after_plot <- VlnPlot(
        merged_obj_qc,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.1
    ) +
        plot_annotation(title = "After QC")

    ggsave(
        filename = file.path(plot_dir, "after_qc_violin.png"),
        plot = after_plot,
        width = 12,
        height = 6,
        dpi = 300
    )

    # Summary statistics
    message(paste(
        "Cells before QC:",
        ncol(merged_obj),
        "| Cells after QC:",
        ncol(merged_obj_qc)
    ))

    return(merged_obj_qc)
}

#' Remove doublets using scDblFinder
#'
#' @param qc_obj A Seurat object that has passed QC.
#' @return A Seurat object with doublets removed.
remove_doublets <- function(qc_obj) {
    # Set up logging
    log_dir <- "result/201.load_sc"
    dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
    log_file <- file.path(log_dir, "doublet_removal.log")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        log_entry <- paste("[", timestamp, "]", msg)
        message(msg)
        cat(log_entry, "\n", file = log_file, append = TRUE)
    }

    log_message("Removing doublets with scDblFinder...")
    if (!requireNamespace("scDblFinder", quietly = TRUE)) {
        stop("Package 'scDblFinder' is required. Please install it.")
    }

    log_message("Extracting counts from Seurat v5 object for scDblFinder..")
    qc_obj <- JoinLayers(qc_obj)
    counts <- Seurat::GetAssayData(qc_obj, layer = "counts")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = qc_obj[[]]
    )

    bpparam <- BiocParallel::MulticoreParam(
        workers = 3
    )
    sce <- scDblFinder::scDblFinder(
        sce,
        samples = "sample_id",
        clusters = TRUE,
        multiSampleMode = "split",
        BPPARAM = bpparam
    )

    # Add doublet info to Seurat object metadata and filter
    qc_obj$scDblFinder_class <- sce$scDblFinder.class
    qc_obj$scDblFinder_score <- sce$scDblFinder.score

    n_doublets <- sum(qc_obj$scDblFinder_class == "doublet", na.rm = TRUE)
    log_message(paste("Found and removing", n_doublets, "doublets."))

    obj_no_doublets <- qc_obj[, qc_obj$scDblFinder_class == "singlet"]

    log_message(paste(
        "Cells before doublet removal:",
        ncol(qc_obj),
        "| Cells after doublet removal:",
        ncol(obj_no_doublets)
    ))

    return(obj_no_doublets)
}

#' Save Seurat object to H5AD file
#'
#' @param qc_obj The Seurat object after QC.
#' @param output_file The path to save the H5AD file.
#' @return The path to the saved H5AD file.
save_merged_h5ad <- function(qc_obj, output_file) {
    reticulate::use_virtualenv("./.venv_anndata", required = TRUE)
    message("Saving QC'd Seurat object to H5AD file...")
    if (!requireNamespace("scCustomize", quietly = TRUE)) {
        stop("Package 'scCustomize' is required. Please install it.")
    }
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    # Force-regenerate the 'data' slot from 'counts' to fix inconsistencies
    message("Re-normalizing data to ensure consistency before saving...")
    qc_obj <- NormalizeData(qc_obj, verbose = FALSE)

    scCustomize::as.anndata(
        x = qc_obj,
        file_path = dirname(output_file),
        file_name = basename(output_file)
    )
    return(output_file)
}


#' Run scVI Integration Python Script
#'
#' This function serves as a wrapper to execute a Python integration script
#' from within R. It ensures the correct Python virtual environment is used
#' and passes the necessary file paths. The path to the python script is
#' passed as an argument to make it a trackable dependency in the `targets`
#' pipeline.
#'
#' @param python_script_path The file path to the python integration script.
#' @param input_h5ad The path to the H5AD file to be integrated.
#' @param output_h5ad The path where the integrated H5AD file will be saved.
#' @return The path to the integrated H5AD file, for use in `targets`.
run_scvi_integration <- function(python_script_path, input_h5ad, output_h5ad) {
    message("--- Running scVI integration script ---")

    # Define the path to the virtual environment
    venv_path <- "./.venv_scvi"

    # Ensure the virtual environment and python script exist
    if (!dir.exists(venv_path)) {
        stop(paste("Python virtual environment not found at:", venv_path))
    }
    if (!file.exists(python_script_path)) {
        stop(paste("Python script not found at:", python_script_path))
    }
    reticulate::use_virtualenv(venv_path, required = TRUE)

    # Construct the command and execute it
    message("Executing Python script for scVI integration...")
    result <- system2(
        "python",
        args = c(shQuote(python_script_path), shQuote(input_h5ad), shQuote(output_h5ad)),
        stdout = TRUE,
        stderr = TRUE
    )

    # Check for errors during execution
    status <- attr(result, "status")
    if (!is.null(status) && status != 0) {
        stop("Python script execution failed. Error:\n", paste(result, collapse = "\n"))
    } else {
        message("Python script executed successfully.")
        print(result)
    }

    # Return the output file path for the `targets` pipeline
    return(output_h5ad)
}
