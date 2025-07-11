# ------------------------------------------------------------------------------
# 202.annotation.R - Functions for Single-Cell Annotation
#
# This script provides functions for annotating single-cell datasets after
# integration. It is designed to be sourced by a `targets` pipeline.
#
# Contains:
# - A function to load an integrated H5AD file as a Seurat object.
# ------------------------------------------------------------------------------

#' Load Integrated H5AD as Seurat Object
#'
#' This function loads a Python-generated H5AD file containing integrated
#' single-cell data and converts it into a Seurat object for downstream
#' analysis in R.
#'
#' @param h5ad_path The file path to the integrated H5AD file.
#' @return A Seurat object with the integrated data.
load_integrated_h5ad <- function(h5ad_path) {
    message(paste("Loading integrated H5AD file from:", h5ad_path))

    if (!requireNamespace("schard", quietly = TRUE)) {
        stop("Package 'schard' is required. Please install it to load H5AD files.")
    }

    # Load the H5AD file and convert it to a Seurat object
    seurat_obj <- schard::h5ad2seurat(h5ad_path)

    message("Successfully converted H5AD to Seurat object.")
    return(seurat_obj)
}