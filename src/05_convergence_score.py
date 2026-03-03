import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

def compute_convergence(deg_file, signature_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    
    deg_df = pd.read_csv(deg_file.replace('.tsv', '_mapped.tsv'), sep="\t")

    if 'Symbol' not in deg_df.columns:
        deg_df['Symbol'] = deg_df['Gene']

    aga_up = set(deg_df[deg_df['direction'] == 'Up']['Symbol'])
    aga_down = set(deg_df[deg_df['direction'] == 'Down']['Symbol'])
    
    def load_sig(filename):
        path = os.path.join(signature_dir, filename)
        if os.path.exists(path):
            with open(path, 'r') as f:
                return set([line.strip() for line in f if line.strip()])
        return set()
        
    li_up = load_sig('li_up.txt')
    li_down = load_sig('li_down.txt')
    zn_up = load_sig('zn_up.txt')
    zn_down = load_sig('zn_down.txt')
    
    # Rescue Score Directional:
    # How much drug reverses the disease signature
    # Rescue = (AGA_Up intersected with Drug_Down) + (AGA_Down intersected with Drug_Up)
    
    def calc_rescue(drug_up, drug_down):
        rescue_up = aga_up.intersection(drug_down)
        rescue_down = aga_down.intersection(drug_up)
        return len(rescue_up) + len(rescue_down)
        
    li_score = calc_rescue(li_up, li_down)
    zn_score = calc_rescue(zn_up, zn_down)
    
    combo_up = li_up.union(zn_up)
    combo_down = li_down.union(zn_down)
    combo_score = calc_rescue(combo_up, combo_down)
    
    results = pd.DataFrame({
        'Treatment': ['Li', 'Zn', 'Li+Zn Combo'],
        'Directional_Rescue_Score': [li_score, zn_score, combo_score]
    })
    
    out_file = os.path.join(out_dir, 'Connectivity_scores.tsv')
    results.to_csv(out_file, sep='\t', index=False)
    
    plt.figure(figsize=(6, 5))
    sns.barplot(x='Treatment', y='Directional_Rescue_Score', data=results, palette='Set2')
    plt.title('Convergence / Rescue Score')
    plt.ylabel('Number of Rescued Genes')
    plt.tight_layout()
    plt.savefig(os.path.join("aga-li-zn/results/figures", 'convergence_barplot.png'), dpi=300)
    plt.close()
    
    print(f"Convergence scores saved to {out_file}")

if __name__ == "__main__":
    compute_convergence(
        "aga-li-zn/results/tables/DEG_table_full.tsv",
        "aga-li-zn/data/raw",
        "aga-li-zn/results/tables"
    )
