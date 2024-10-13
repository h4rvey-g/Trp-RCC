load_sc_pre <- function(h5seurat_path) {
  sc_pre <- LoadH5Seurat(h5seurat_path)
  sc_pre <- sc_pre %>%
    # create col batch based on dataset, if dataset starts with "PD", batch is "PD",
    # if dataset starts with "RCC", batch is "RCC"
    # if dataset starts with "GSM", batch is the first 7 characters of dataset
    mutate(batch = case_when(
      str_detect(dataset, "^PD") ~ "PD",
      str_detect(dataset, "^RCC") ~ "RCC",
      str_detect(dataset, "^GSM") ~ str_sub(dataset, 1, 7),
      TRUE ~ "unknown"
    ))
  sc_pre <- RunUMAP(sc_pre, dims = 1:20, reduction = "pca", reduction.name = "umap_on_pca")
  sc_anndata <- AnnData(
    X = sc_pre[["RNA"]]$counts %>% t(),
    obs = sc_pre[[]],
    obsm = list(
      pca = as.matrix(Embeddings(sc_pre, reduction = "pca"))
    )
  )
  write_h5ad(sc_anndata, "data/201.load_sc/sc_pre.h5ad")
  sc_pre
}

check_batch_effect <- function(sc_pre) {
  p1 <- DimPlot(sc_pre, reduction = "umap_on_pca", group.by = "type")
  p2 <- DimPlot(sc_pre, reduction = "umap_on_pca", group.by = "type", split.by = "batch")
  p <- patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1, 8))
  ggsave("result/201.load_sc/batch_effect.png", p, width = 30, height = 7)
  p_list <- lapply(unique(sc_pre$batch), function(one_batch) {
    sc_batch <- sc_pre %>% filter(batch == one_batch)
    DimPlot(sc_batch, reduction = "umap_on_pca", group.by = "dataset") + 
      ggtitle(paste("Batch:", one_batch))
  })
  p_combined <- patchwork::wrap_plots(p_list, nrow = 3)
  ggsave("result/201.load_sc/batch_effect_dataset.png", p_combined, width = 20, height = 20)
  "result/201.load_sc/batch_effect.png"
}

run_integration <- function() {
  # run code/Python/201.integration.ipynb
  # use_python("/home/vscode/.pyenv/shims/python", required = TRUE)
  # py_run_file("code/Python/201.integration.py")
  "data/201.load_sc/scvi_latent.tsv"
}

import_integration <- function(latent, sc_pre) {
  scvi <- read_tsv(latent) %>% as.matrix()
  rownames(scvi) <- colnames(sc_pre[["RNA"]])
  colnames(scvi) <- paste0("scVI_", 1:ncol(scvi))
  sc <- sc_pre
  sc[["scvi"]] <- CreateDimReducObject(embeddings = scvi, key = "scVI_", assay = "RNA")
  sc
}

preprocess_sc <- function(sc) {
  # mkdir if not exist
  dir.create("result/201.load_sc", showWarnings = FALSE)
  sc_before <- FindNeighbors(sc, dims = 1:10, reduction = "pca")
  sc_before <- FindClusters(sc_before, resolution = 0.5)
  sc_before <- RunUMAP(sc_before, dims = 1:10, reduction = "pca", reduction.name = "umap_on_pca")
  p1 <- DimPlot(sc_before, reduction = "umap_on_pca", group.by = "type")
  p2 <- DimPlot(sc_before, reduction = "umap_on_pca", group.by = "RNA_snn_res.0.5")
  p3 <- DimPlot(sc_before, reduction = "umap_on_pca", split.by = "type")
  p <- (p1 + p2) / p3
  ggsave("result/201.load_sc/umap_after_pca.png", p, width = 10, height = 10)

  sc <- FindNeighbors(sc, dims = 1:10, reduction = "scvi")
  sc <- FindClusters(sc, resolution = 0.5)
  sc <- RunUMAP(sc, dims = 1:10, reduction = "scvi")
  p1 <- DimPlot(sc, reduction = "umap", group.by = "type")
  p2 <- DimPlot(sc, reduction = "umap", group.by = "RNA_snn_res.0.5")
  p3 <- DimPlot(sc, reduction = "umap", split.by = "type")
  p <- (p1 + p2) / p3
  ggsave("result/201.load_sc/umap_after_scvi.png", p, width = 10, height = 10)

  sc_pro <- sc
  sc_pro
}
