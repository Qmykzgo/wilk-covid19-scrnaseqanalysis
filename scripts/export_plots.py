import os
import scanpy as sc
import matplotlib.pyplot as plt

def main():
    print("Loading clustered data...")
    adata_path = "/home/yer_kanat/Downloads/wilk covid 19/results/wilk_covid_clustered.h5ad"
    if not os.path.exists(adata_path):
        print(f"Error: {adata_path} does not exist.")
        return
        
    adata = sc.read_h5ad(adata_path)
    
    # Configure high quality plotting settings
    sc.set_figure_params(dpi=150, facecolor='white', frameon=False)
    
    print("Generating UMAP subplots...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Left subplot: CellTypist Annotated Cell Types
    sc.pl.umap(
        adata, 
        color='celltypist_label', 
        ax=axes[0], 
        show=False, 
        legend_fontsize=8, 
        legend_loc='right margin', 
        title='Cell Types (CellTypist Predictions)'
    )
    
    # Right subplot: Disease Condition
    sc.pl.umap(
        adata, 
        color='condition', 
        ax=axes[1], 
        show=False, 
        legend_fontsize=10, 
        legend_loc='right margin', 
        title='Condition (COVID-19 vs. Healthy)'
    )
    
    plt.suptitle('PBMC Single-Cell UMAP (Harmony Integrated)', y=0.98, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = "/home/yer_kanat/Downloads/wilk covid 19/results/umap_plots.png"
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Successfully saved UMAP plot to: {output_path}")

if __name__ == "__main__":
    main()
