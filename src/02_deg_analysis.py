import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def run_deg(expr_file, pheno_file, out_dir):
    os.makedirs(os.path.join(out_dir, "tables"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "figures"), exist_ok=True)
    
    print(f"Loading preprocessed data from {expr_file}...")
    expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    
    group_col = 'Group'
    if group_col not in pheno_df.columns:
        print("Error: 'Group' column not found in phenotype data.")
        return
        
    aga_samples = pheno_df[pheno_df[group_col] == 'AGA']['sample_id'].values
    ctrl_samples = pheno_df[pheno_df[group_col] == 'Control']['sample_id'].values
    
    aga_samples = [s for s in aga_samples if s in expr_df.columns]
    ctrl_samples = [s for s in ctrl_samples if s in expr_df.columns]
    
    if len(aga_samples) == 0 or len(ctrl_samples) == 0:
        print("Not enough samples for comparison.")
        return
        
    print(f"Comparing {len(aga_samples)} AGA vs {len(ctrl_samples)} Control...")
    
    aga_data = expr_df[aga_samples].values
    ctrl_data = expr_df[ctrl_samples].values
    
    log2fc = np.mean(aga_data, axis=1) - np.mean(ctrl_data, axis=1)
    
    t_stat, p_vals = stats.ttest_ind(aga_data, ctrl_data, axis=1, nan_policy='omit')
    
    mask = ~np.isnan(p_vals)
    fdr = np.full(p_vals.shape, np.nan)
    if np.sum(mask) > 0:
        _, fdr[mask], _, _ = multipletests(p_vals[mask], method='fdr_bh')
    
    deg_df = pd.DataFrame({
        'Gene': expr_df.index,
        'log2FC': log2fc,
        'pvalue': p_vals,
        'FDR': fdr
    })
    
    deg_df['direction'] = 'Not Sig'
    deg_df.loc[(deg_df['pvalue'] < 0.05) & (deg_df['log2FC'] > 0.3), 'direction'] = 'Up'
    deg_df.loc[(deg_df['pvalue'] < 0.05) & (deg_df['log2FC'] < -0.3), 'direction'] = 'Down'
    
    out_table = os.path.join(out_dir, "tables", "DEG_table_full.tsv")
    deg_df.to_csv(out_table, sep='\t', index=False)
    print(f"DEG table saved to {out_table}")
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x='log2FC', y=-np.log10(deg_df['pvalue']),
        hue='direction',
        palette={'Up': 'red', 'Down': 'blue', 'Not Sig': 'grey'},
        data=deg_df,
        s=10, alpha=0.5
    )
    plt.axvline(x=0.3, color='k', linestyle='--')
    plt.axvline(x=-0.3, color='k', linestyle='--')
    plt.axhline(y=-np.log10(0.05), color='k', linestyle='--')
    plt.title('Volcano Plot (AGA vs Control)')
    plt.savefig(os.path.join(out_dir, "figures", "volcano.png"), dpi=300)
    plt.close()
    
    top50 = deg_df.sort_values('pvalue').head(50)
    top50_genes = top50['Gene'].values
    
    heatmap_data = expr_df.loc[top50_genes, list(aga_samples) + list(ctrl_samples)]
    
    z_scores = stats.zscore(heatmap_data.values, axis=1, nan_policy='omit')
    heatmap_data = pd.DataFrame(z_scores, index=heatmap_data.index, columns=heatmap_data.columns)
    heatmap_data = heatmap_data.dropna()
    
    plt.figure(figsize=(10, 12))
    if not heatmap_data.empty:
        sns.heatmap(heatmap_data, cmap='viridis', center=0, yticklabels=False, xticklabels=False)
    plt.title('Top 50 DEGs Heatmap')
    plt.savefig(os.path.join(out_dir, "figures", "heatmap_top50.png"), dpi=300)
    plt.close()
    
    print("Plots generated successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--expr", required=True)
    parser.add_argument("--pheno", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    run_deg(args.expr, args.pheno, args.outdir)
