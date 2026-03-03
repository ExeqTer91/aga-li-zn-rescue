import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

datasets = ['GSE90594', 'GSE36169', 'GSE66663']
all_results = []

for ds in datasets:
    if ds == 'GSE90594':
        path = f'aga-li-zn/results/tables/modal_convergence_stats_{ds}.tsv'
    else:
        path = f'aga-li-zn/results/{ds}/tables/modal_convergence_stats_{ds}.tsv'
        
    if os.path.exists(path):
        df = pd.read_csv(path, sep='\t')
        df['Dataset'] = ds
        all_results.append(df)

if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    out_dir = 'aga-li-zn/results/figures'
    os.makedirs(out_dir, exist_ok=True)
    
    plt.figure(figsize=(10, 6))
    sns.barplot(data=combined, x='Dataset', y='Z_Score', hue='Treatment', palette='Set2')
    plt.axhline(1.96, ls='--', color='red', alpha=0.7, label='p=0.05 (Z=1.96)')
    plt.title('Convergence Z-Scores Across AGA Datasets (Modal Permutations)')
    plt.ylabel('Z-Score (higher = more rescue)')
    plt.legend()
    plt.tight_layout()
    
    out_path = os.path.join(out_dir, 'cross_dataset_zscores.png')
    plt.savefig(out_path, dpi=300)
    print(f"Saved plot to {out_path}")
    
    combined.to_csv('aga-li-zn/results/tables/cross_dataset_modal_summary.tsv', sep='\t', index=False)
    print("Saved combined summary to aga-li-zn/results/tables/cross_dataset_modal_summary.tsv")
