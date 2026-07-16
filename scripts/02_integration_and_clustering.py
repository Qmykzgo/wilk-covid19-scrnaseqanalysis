import scanpy as sc
import scanpy.external as sce
import harmonypy as hm
import celltypist
from celltypist import models
import pandas as pd
import os

sc.settings.verbosity = 1

def main():
    print("Loading normalized data...")
    adata = sc.read_h5ad("results/wilk_covid_qc_normalized.h5ad")
    
    # Extract highly variable genes
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    
    print("Scaling and running PCA...")
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver="arpack")
    
    # Run Harmony Batch Integration directly to avoid wrapper shape errors
    print("Integrating batches with Harmony...")
    harmony_out = hm.run_harmony(adata_hvg.obsm['X_pca'], adata_hvg.obs, 'sample')
    
    # Handle version compatibility of harmonypy output shape
    if harmony_out.Z_corr.shape[1] == adata_hvg.shape[0]:
        adata_hvg.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
    else:
        adata_hvg.obsm['X_pca_harmony'] = harmony_out.Z_corr
    
    N_PCS = 30
    print(f"Building neighborhood graph with {N_PCS} PCs...")
    sc.pp.neighbors(adata_hvg, use_rep="X_pca_harmony", n_pcs=N_PCS)
    
    print("Running UMAP and Leiden clustering...")
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")
    
    # Transfer representations back to full object
    print("Transferring embeddings back to full object...")
    adata.obs["leiden"] = adata_hvg.obs["leiden"].values
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]
    
    print("Annotating cell types using CellTypist...")
    models.download_models(model="Immune_All_Low.pkl")
    model = models.Model.load("Immune_All_Low.pkl")
    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    adata.obs["celltypist_label"] = predictions.predicted_labels.majority_voting
    
    # In a full run, we would dynamically map the leiden clusters based on markers,
    # but here we use majority voting direct mapping as a strong baseline
    adata.obs["cell_type"] = adata.obs["celltypist_label"].copy()
    
    out_path = "results/wilk_covid_clustered.h5ad"
    adata.write(out_path)
    print(f"Saved clustered and annotated object to {out_path}")

if __name__ == "__main__":
    main()
