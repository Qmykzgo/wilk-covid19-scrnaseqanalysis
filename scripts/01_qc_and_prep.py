import numpy as np
import pandas as pd
import scanpy as sc
import os

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor="white")

def main():
    print("Loading raw processed data...")
    DATA_PATH = "wilk_covid19_pbmc.h5ad" # This assumes data has been downloaded from CZ CELLxGENE
    
    if not os.path.exists(DATA_PATH):
        print(f"File {DATA_PATH} not found. Please download from CZ CELLxGENE first.")
        return

    adata = sc.read_h5ad(DATA_PATH)
    
    # Map Ensembl IDs to Gene Symbols
    adata.var_names = adata.var["feature_name"].astype(str)
    adata.var_names_make_unique()
    adata.var.index.name = None
    
    # 1. Standardize Metadata
    print("Standardizing metadata...")
    # Map from clinical_metadata.csv
    meta_df = pd.read_csv("clinical_metadata.csv").set_index("sample_id")
    
    sample_col = "donor_id" 
    condition_col = "disease"
    
    adata.obs["sample"] = adata.obs[sample_col].astype(str)
    
    # CZ CELLxGENE labels for healthy are often 'normal'
    adata.obs["condition"] = adata.obs[condition_col].astype(str).str.lower().map(
        lambda x: "healthy" if ("normal" in x or "healthy" in x) else "covid19"
    )
    
    # Apply ventilation map
    adata.obs["severity"] = adata.obs["sample"].map(meta_df["ventilation_status"]).fillna("unknown")
    
    # 2. QC Metrics Calculation
    print("Calculating QC metrics...")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P|E)]")
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True)
    
    # 3. Filtering
    print("Filtering cells and genes...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Seq-Well has lower depth; 300 is a safer bound than 10x's 500
    adata = adata[adata.obs["total_counts"] > 300, :]
    adata = adata[adata.obs["total_counts"] < 25000, :]
    adata = adata[adata.obs["n_genes_by_counts"] > 250, :]
    adata = adata[adata.obs["n_genes_by_counts"] < 4000, :]
    adata = adata[adata.obs["pct_counts_mt"] < 20, :]
    
    # 4. Doublet Detection
    print("Running Scrublet for doublet detection per sample...")
    import scanpy.external as sce
    sce.pp.scrublet(adata, batch_key="sample", threshold=0.25)
    
    # Filter doublets
    adata = adata[adata.obs["predicted_doublet"] == False].copy()
    
    # 5. Normalization & Checkpoint
    print("Normalizing and finding HVGs...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample")
    
    out_path = "results/wilk_covid_qc_normalized.h5ad"
    os.makedirs("results", exist_ok=True)
    adata.write(out_path)
    print(f"Saved QC'd and normalized object to {out_path}")

if __name__ == "__main__":
    main()
