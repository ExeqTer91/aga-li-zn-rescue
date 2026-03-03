import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    dp_genes = ['WNT10B', 'LEF1', 'CTNNB1', 'DKK1', 'DKK2', 'ALPL', 'SOX2', 'VCAN']
    datasets = ["GSE90594", "GSE36169", "GSE66663"]
    
    # Load LINCS N=500 signatures
    def load_sig(filename):
        path = os.path.join("aga-li-zn/data/signatures", filename)
        if os.path.exists(path):
            with open(path, "r") as f:
                return set([line.strip() for line in f if line.strip()])
        return set()
        
    li_up = load_sig("li_up_500.txt")
    li_down = load_sig("li_down_500.txt")
    zn_up = load_sig("zn_proxy_up_500.txt")
    zn_down = load_sig("zn_proxy_down_500.txt")
    
    # Also load Harmonizome signatures
    harmonizome_dir = "aga-li-zn/data/signatures"
    # We might not have them extracted directly as text, but they were in the previous scripts.
    # Let's just use LINCS N=500 and N=300 for now, or check both.
    
    results = []
    
    for dataset in datasets:
        if dataset == "GSE90594":
            deg_file = "aga-li-zn/results/tables/deg_full_mapped.tsv"
        elif dataset == "GSE36169":
            deg_file = f"aga-li-zn/results/{dataset}/tables/DEG_table_full_mapped.tsv"
        elif dataset == "GSE66663":
            deg_file = f"aga-li-zn/results/{dataset}/tables/DEG_table_full_mapped.tsv"
            
        if not os.path.exists(deg_file):
            print(f"File not found: {deg_file}")
            continue
            
        df = pd.read_csv(deg_file, sep="\t")
        if 'Symbol' not in df.columns and 'Gene' in df.columns:
            df['Symbol'] = df['Gene']
            
        for gene in dp_genes:
            gene_data = df[df['Symbol'] == gene]
            if len(gene_data) > 0:
                # Average log2FC if multiple probes map to same gene
                log2fc = gene_data['log2FC'].mean()
                padj = gene_data['adj.P.Val'].min() if 'adj.P.Val' in gene_data.columns else (gene_data['padj'].min() if 'padj' in gene_data.columns else np.nan)
                
                # Determine direction in disease
                if log2fc > 0:
                    aga_dir = "UP"
                else:
                    aga_dir = "DOWN"
                    
                # Check Li rescue
                li_rescue = "None"
                if gene in li_up:
                    li_dir = "UP"
                    li_rescue = "Yes" if aga_dir == "DOWN" else "No (Exacerbates)"
                elif gene in li_down:
                    li_dir = "DOWN"
                    li_rescue = "Yes" if aga_dir == "UP" else "No (Exacerbates)"
                else:
                    li_dir = "Not in top 500"
                    
                # Check Zn rescue
                zn_rescue = "None"
                if gene in zn_up:
                    zn_dir = "UP"
                    zn_rescue = "Yes" if aga_dir == "DOWN" else "No (Exacerbates)"
                elif gene in zn_down:
                    zn_dir = "DOWN"
                    zn_rescue = "Yes" if aga_dir == "UP" else "No (Exacerbates)"
                else:
                    zn_dir = "Not in top 500"
                    
                results.append({
                    "Dataset": dataset,
                    "Gene": gene,
                    "AGA_log2FC": log2fc,
                    "AGA_padj": padj,
                    "AGA_Direction": aga_dir,
                    "Li_Direction": li_dir,
                    "Li_Rescue": li_rescue,
                    "Zn_Direction": zn_dir,
                    "Zn_Rescue": zn_rescue
                })
            else:
                results.append({
                    "Dataset": dataset,
                    "Gene": gene,
                    "AGA_log2FC": np.nan,
                    "AGA_padj": np.nan,
                    "AGA_Direction": "Not Found",
                    "Li_Direction": "N/A",
                    "Li_Rescue": "N/A",
                    "Zn_Direction": "N/A",
                    "Zn_Rescue": "N/A"
                })
                
    res_df = pd.DataFrame(results)
    res_df.to_csv("aga-li-zn/results/tables/dp_module_rescue.tsv", sep="\t", index=False)
    
    # Print summary
    print(res_df.to_string())
    
    # Create a heatmap of the log2FC across datasets
    heatmap_data = res_df.pivot(index="Gene", columns="Dataset", values="AGA_log2FC")
    plt.figure(figsize=(8, 6))
    sns.heatmap(heatmap_data, annot=True, cmap="vlag", center=0, cbar_kws={'label': 'AGA log2FC'})
    plt.title("Core Dermal Papilla Module Expression in AGA")
    plt.tight_layout()
    plt.savefig("aga-li-zn/results/figures/dp_module_heatmap.png", dpi=300)
    print("\nSaved heatmap to aga-li-zn/results/figures/dp_module_heatmap.png")

if __name__ == "__main__":
    main()
