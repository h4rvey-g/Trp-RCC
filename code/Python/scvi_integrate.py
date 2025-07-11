import argparse
import os
import scanpy as sc
import scvi
import torch
import scipy.sparse

def main(input_file, output_file, gpu):
    """
    Integrates single-cell data using scVI, generates UMAP plots before and
    after integration, and saves the integrated data.

    Args:
        input_file (str): Path to the input H5AD file (e.g., 'sc_merged_qc.h5ad').
        output_file (str): Path to save the integrated H5AD file.
        gpu (int): GPU device to use.
    """
    # --------------------------------------------------------------------------
    # 0. Setup
    # --------------------------------------------------------------------------
    print("--- Setting up environment ---")
    
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    plot_dir = "result/201.load_sc"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # Configure scvi-tools to use 3/4 of available CPU threads
    available_threads = torch.get_num_threads()
    threads_to_use = max(1, int(available_threads * 0.75))
    scvi.settings.num_threads = threads_to_use
    print(f"Using {threads_to_use} out of {available_threads} CPU threads.")


    # --------------------------------------------------------------------------
    # 1. Load Data
    # --------------------------------------------------------------------------
    print(f"--- Loading data from {input_file} ---")
    adata = sc.read_h5ad(input_file)
    
    # Store raw data right away
    adata.raw = adata
    
    # Ensure the 'counts' layer exists for scVI, using 'counts_RNA' as a source
    if "counts_RNA" in adata.layers:
        adata.layers["counts"] = adata.layers["counts_RNA"].copy()
    else:
        # Fallback if 'counts_RNA' is not present
        adata.layers["counts"] = adata.X.copy()
    
    # Use normalized/log-transformed data for PCA/UMAP. Start from raw counts.
    adata.X = adata.layers['counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # --------------------------------------------------------------------------
    # 2. Plot UMAP Before Integration
    # --------------------------------------------------------------------------
    print("--- Generating UMAP before integration ---")
    # Find highly variable genes before PCA
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=4000
    )
    # Run PCA on the highly variable genes
    sc.pp.pca(adata, mask_var="highly_variable")
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color="dataset",
        title="UMAP Before Integration (by Dataset)",
        save="_before_integration_dataset.png",
        show=False
    )
    # Move the plot to the correct directory
    os.rename(
        os.path.join("figures", "umap_before_integration_dataset.png"),
        os.path.join(plot_dir, "umap_before_integration_dataset.png")
    )


    # --------------------------------------------------------------------------
    # 3. Integrate with scVI
    # --------------------------------------------------------------------------
    print("--- Setting up and training scVI model ---")

    # Convert to CSR for faster training
    if scipy.sparse.issparse(adata.layers["counts"]):
        adata.layers["counts"] = adata.layers["counts"].tocsr()

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")
    
    model = scvi.model.SCVI(adata)
    if gpu == -1:
        print("--- Training on CPU ---")
        model.train(accelerator="cpu")
    else:
        print(f"--- Training on GPU device {gpu} ---")
        model.train(accelerator="gpu", devices=[gpu])

    # --------------------------------------------------------------------------
    # 4. Plot UMAP After Integration
    # --------------------------------------------------------------------------
    print("--- Generating UMAP after integration ---")
    adata.obsm["X_scVI"] = model.get_latent_representation()
    
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    
    sc.pl.umap(
        adata,
        color="dataset",
        title="UMAP After scVI Integration (by Dataset)",
        save="_after_scvi_integration_dataset.png",
        show=False
    )
    # Move the plot to the correct directory
    os.rename(
        os.path.join("figures", "umap_after_scvi_integration_dataset.png"),
        os.path.join(plot_dir, "umap_after_scvi_integration_dataset.png")
    )


    # --------------------------------------------------------------------------
    # 5. Save Integrated Data
    # --------------------------------------------------------------------------
    print(f"--- Saving integrated data to {output_file} ---")
    adata.write_h5ad(output_file, compression="gzip")
    
    print("--- Integration complete ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Integrate single-cell data using scVI.")
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input H5AD file."
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path where the integrated H5AD file will be saved."
    )
    parser.add_argument(
        "--gpu",
        type=int,
        default=0,
        help="GPU device to use (e.g., 0). Set to -1 to force CPU. Defaults to GPU 0."
    )
    
    args = parser.parse_args()
    main(args.input_file, args.output_file, args.gpu)