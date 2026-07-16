import scanpy as sc
import pandas as pd
try:
    import milopy
    import milopy.core as milo
    HAS_MILOPY = True
except ImportError:
    HAS_MILOPY = False
    print("Warning: milopy not found. Milo differential abundance analysis will be skipped.")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import liana as li

sc.settings.verbosity = 1

def run_milopy(adata, n_pcs=30):
    print("--- Running Differential Abundance (Milo) ---")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=n_pcs)
    milo.make_nhoods(adata)
    milo.count_nhoods(adata, sample_col="sample")
    
    # Filter to known severities for testing
    milo.DA_nhoods(adata, design="~severity")
    milo_results = adata.uns["nhood_adata"].obs
    print(milo_results.sort_values("SpatialFDR").head(10))
    return milo_results

def run_pydeseq2(adata):
    print("\n--- Running Pseudobulk Differential Expression (PyDESeq2) ---")
    # Subsetting to monocytes to test the paper's main claim
    mono = adata[adata.obs["cell_type"].str.contains("monocyte|Monocyte", case=False)]
    
    if len(mono) == 0:
        print("Warning: No monocyte cluster found for DESeq2 analysis.")
        return
    
    counts = pd.DataFrame(
        mono.layers["counts"].toarray() if hasattr(mono.layers["counts"], "toarray") else mono.layers["counts"],
        index=mono.obs_names, columns=mono.var_names,
    )
    counts["sample"] = mono.obs["sample"].values
    pb_counts = counts.groupby("sample").sum()
    
    pb_meta = mono.obs[["sample", "condition"]].drop_duplicates().set_index("sample").loc[pb_counts.index]
    
    try:
        # Compatible across newer and older pydeseq2 versions
        dds = DeseqDataSet(counts=pb_counts.astype(int), metadata=pb_meta, design_factors="condition")
    except TypeError:
        # Fallback to older clinical / design arg setup
        dds = DeseqDataSet(counts=pb_counts.astype(int), clinical=pb_meta, design="~condition")
    
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", "covid19", "healthy"])
    stat_res.summary()
    
    hla_genes = [g for g in ["HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DPA1"] if g in stat_res.results_df.index]
    cytokine_genes = [g for g in ["IL1B", "IL6", "TNF", "CXCL8"] if g in stat_res.results_df.index]
    
    print("\nHLA Class II and Cytokine Gene Log2FC (COVID-19 vs Healthy Monocytes):")
    print(stat_res.results_df.loc[hla_genes + cytokine_genes])

def run_liana(adata):
    print("\n--- Running Cell-Cell Communication (LIANA) ---")
    results = {}
    for cond in adata.obs["condition"].unique():
        sub = adata[adata.obs["condition"] == cond].copy()
        try:
            li.mt.rank_aggregate(sub, groupby="cell_type",
                                 use_raw=False,
                                 expr_prop=0.1,
                                 resource_name="consensus", verbose=False)
            results[cond] = sub.uns["liana_res"]
        except Exception as e:
            print(f"Skipping LIANA for {cond} due to error: {e}")
            
    if "covid19" in results and "healthy" in results:
        covid_pairs = set(zip(results["covid19"]["ligand_complex"], results["covid19"]["receptor_complex"],
                               results["covid19"]["source"], results["covid19"]["target"]))
        healthy_pairs = set(zip(results["healthy"]["ligand_complex"], results["healthy"]["receptor_complex"],
                                 results["healthy"]["source"], results["healthy"]["target"]))
        
        covid_specific = covid_pairs - healthy_pairs
        print(f"\n{len(covid_specific)} signaling pairs detected ONLY in COVID-19 samples.")

def main():
    print("Loading clustered data...")
    adata = sc.read_h5ad("results/wilk_covid_clustered.h5ad")
    
    if HAS_MILOPY:
        run_milopy(adata)
    else:
        print("\n--- Skipping Milo (milopy not installed) ---")
    run_pydeseq2(adata)
    run_liana(adata)
    
if __name__ == "__main__":
    main()
