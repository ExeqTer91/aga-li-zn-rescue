import pandas as pd
import gseapy as gp
import os
import matplotlib.pyplot as plt
import seaborn as sns

def run_gsea(deg_file, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"Loading DEG data from {deg_file}...")
    deg_df = pd.read_csv(deg_file, sep='\t')
    
    # Drop rows with NaN in log2FC or pvalue
    deg_df = deg_df.dropna(subset=['log2FC', 'pvalue'])
    
    # We need gene symbols, assuming 'Gene' column is symbol or Entrez. 
    # For GSE90594, platform is GPL17077, so the ID is likely Agilent Probe ID. 
    # Let's map it if needed, or assume it's already Symbol for this pipeline's simplicity.
    # To make gseapy work correctly with real IDs, we should ideally have gene symbols.
    # We'll use the 'Gene' column as the name, and rank by log2FC.
    
    # Prepare rank format for prerank: a DataFrame with gene name and rank score
    # score = -log10(pvalue) * sign(log2FC)
    
    # Avoid log10(0)
    min_p = deg_df[deg_df['pvalue'] > 0]['pvalue'].min()
    deg_df['pvalue'] = deg_df['pvalue'].replace(0, min_p)
    
    # Map Probes to Gene Symbols
    gpl_file = 'data/raw/GPL17077_annotation.tsv'
    if os.path.exists(gpl_file):
        gpl = pd.read_csv(gpl_file, sep='\t', low_memory=False)
        symbol_col = next((col for col in gpl.columns if 'symbol' in col.lower() or ('gene' in col.lower() and 'name' not in col.lower())), None)
        if symbol_col:
            probe_to_symbol = dict(zip(gpl['ID'], gpl[symbol_col]))
            deg_df['Symbol'] = deg_df['Gene'].map(probe_to_symbol)
            deg_df = deg_df.dropna(subset=['Symbol'])
            deg_df['Symbol'] = deg_df['Symbol'].apply(lambda x: str(x).split('///')[0].strip()) # Take first symbol if multiple
    else:
        deg_df['Symbol'] = deg_df['Gene'] # Fallback
        
    deg_df['score'] = -np.log10(deg_df['pvalue']) * np.sign(deg_df['log2FC'])
    
    # Let's do Over-Representation Analysis (ORA) on top DEGs for robustness.
    # Take top 300 up and down by p-value to avoid payload too large for Enrichr
    up_genes_df = deg_df[deg_df['direction'] == 'Up'].sort_values('pvalue').head(300)
    down_genes_df = deg_df[deg_df['direction'] == 'Down'].sort_values('pvalue').head(300)
    
    up_genes = up_genes_df['Symbol'].unique().tolist()
    down_genes = down_genes_df['Symbol'].unique().tolist()
    
    print(f"Up genes: {len(up_genes)}, Down genes: {len(down_genes)}")
    
    gene_sets = ['GO_Biological_Process_2021', 'KEGG_2021_Human', 'Reactome_2022']
    
    results = []
    
    if len(up_genes) > 0:
        try:
            enr_up = gp.enrichr(gene_list=up_genes,
                                gene_sets=gene_sets,
                                outdir=None)
            res_up = enr_up.results
            res_up['Direction'] = 'Up'
            results.append(res_up)
        except Exception as e:
            print(f"Enrichment for UP genes failed: {e}")
            
    if len(down_genes) > 0:
        try:
            enr_down = gp.enrichr(gene_list=down_genes,
                                  gene_sets=gene_sets,
                                  outdir=None)
            res_down = enr_down.results
            res_down['Direction'] = 'Down'
            results.append(res_down)
        except Exception as e:
            print(f"Enrichment for DOWN genes failed: {e}")
            
    if results:
        all_res = pd.concat(results, ignore_index=True)
        # Filter significant
        sig_res = all_res[all_res['Adjusted P-value'] < 0.05]
        
        out_file = os.path.join(out_dir, 'GSEA_results.tsv')
        sig_res.to_csv(out_file, sep='\t', index=False)
        print(f"Saved GSEA results to {out_file}")
        
        # Plot top 10 for Up and Down
        top_up = sig_res[sig_res['Direction'] == 'Up'].sort_values('Adjusted P-value').head(10)
        top_down = sig_res[sig_res['Direction'] == 'Down'].sort_values('Adjusted P-value').head(10)
        
        plot_df = pd.concat([top_up, top_down])
        if not plot_df.empty:
            plot_df['-log10(FDR)'] = -np.log10(plot_df['Adjusted P-value'])
            
            plt.figure(figsize=(10, 8))
            sns.barplot(
                x='-log10(FDR)', 
                y='Term', 
                hue='Direction', 
                data=plot_df,
                palette={'Up': 'red', 'Down': 'blue'}
            )
            plt.title('Top Enriched Pathways (AGA vs Control)')
            plt.tight_layout()
            plt.savefig(os.path.join("aga-li-zn/results/figures", "gsea_barplot.png"), dpi=300)
            plt.close()
            print("GSEA plot saved.")
    else:
        print("No GSEA results generated.")

if __name__ == "__main__":
    import numpy as np # import here to avoid error
    run_gsea("aga-li-zn/results/tables/DEG_table_full.tsv", "aga-li-zn/results/tables")
