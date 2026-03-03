import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def run_null_model(deg_file, signature_dir, out_dir, n_permutations=10000):
    os.makedirs(out_dir, exist_ok=True)
    
    print("Loading data for null model...")
    deg_df = pd.read_csv(deg_file, sep='\t')
    
    # Map Probes to Gene Symbols if needed
    gpl_file = 'data/raw/GPL17077_annotation.tsv'
    if os.path.exists(gpl_file):
        gpl = pd.read_csv(gpl_file, sep='\t', low_memory=False)
        symbol_col = next((col for col in gpl.columns if 'symbol' in col.lower() or ('gene' in col.lower() and 'name' not in col.lower())), None)
        if symbol_col:
            probe_to_symbol = dict(zip(gpl['ID'], gpl[symbol_col]))
            deg_df['Symbol'] = deg_df['Gene'].map(probe_to_symbol)
            deg_df = deg_df.dropna(subset=['Symbol'])
            deg_df['Symbol'] = deg_df['Symbol'].apply(lambda x: str(x).split('///')[0].strip())
    else:
        deg_df['Symbol'] = deg_df['Gene']

    # Focus only on significant DEGs
    aga_up = set(deg_df[deg_df['direction'] == 'Up']['Symbol'])
    aga_down = set(deg_df[deg_df['direction'] == 'Down']['Symbol'])
    background_genes = list(deg_df['Symbol'].unique())
    
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
    
    # Restrict to background
    li_up = li_up.intersection(background_genes)
    li_down = li_down.intersection(background_genes)
    zn_up = zn_up.intersection(background_genes)
    zn_down = zn_down.intersection(background_genes)
    
    def calc_rescue(drug_up, drug_down):
        rescue_up = aga_up.intersection(drug_down)
        rescue_down = aga_down.intersection(drug_up)
        return len(rescue_up) + len(rescue_down)
        
    actual_li = calc_rescue(li_up, li_down)
    actual_zn = calc_rescue(zn_up, zn_down)
    actual_combo = calc_rescue(li_up.union(zn_up), li_down.union(zn_down))
    
    # Permutations
    print(f"Running {n_permutations} permutations...")
    li_nulls = []
    zn_nulls = []
    combo_nulls = []
    
    len_li_up = len(li_up)
    len_li_down = len(li_down)
    len_zn_up = len(zn_up)
    len_zn_down = len(zn_down)
    
    bg_array = np.array(background_genes)
    
    for i in range(n_permutations):
        if i % 1000 == 0:
            print(f"Permutation {i}/{n_permutations}")
        # Shuffle background to pick random sets
        np.random.shuffle(bg_array)
        
        rand_li_up = set(bg_array[:len_li_up])
        rand_li_down = set(bg_array[len_li_up:len_li_up+len_li_down])
        
        rand_zn_up = set(bg_array[:len_zn_up])
        rand_zn_down = set(bg_array[len_zn_up:len_zn_up+len_zn_down])
        
        li_nulls.append(calc_rescue(rand_li_up, rand_li_down))
        zn_nulls.append(calc_rescue(rand_zn_up, rand_zn_down))
        combo_nulls.append(calc_rescue(rand_li_up.union(rand_zn_up), rand_li_down.union(rand_zn_down)))
        
    # Calculate p-values
    li_pval = (np.sum(np.array(li_nulls) >= actual_li) + 1) / (n_permutations + 1)
    zn_pval = (np.sum(np.array(zn_nulls) >= actual_zn) + 1) / (n_permutations + 1)
    combo_pval = (np.sum(np.array(combo_nulls) >= actual_combo) + 1) / (n_permutations + 1)
    
    li_z = (actual_li - np.mean(li_nulls)) / (np.std(li_nulls) + 1e-9)
    zn_z = (actual_zn - np.mean(zn_nulls)) / (np.std(zn_nulls) + 1e-9)
    combo_z = (actual_combo - np.mean(combo_nulls)) / (np.std(combo_nulls) + 1e-9)
    
    res_df = pd.DataFrame({
        'Treatment': ['Li', 'Zn', 'Li+Zn Combo'],
        'Actual_Score': [actual_li, actual_zn, actual_combo],
        'Z_Score': [li_z, zn_z, combo_z],
        'P_Value': [li_pval, zn_pval, combo_pval]
    })
    
    out_table = os.path.join("aga-li-zn/results/tables", 'null_models_results.tsv')
    res_df.to_csv(out_table, sep='\t', index=False)
    print(f"Null model results saved to {out_table}")
    print(res_df)
    
    # Plot Null distribution for Combo
    plt.figure(figsize=(8, 6))
    sns.histplot(combo_nulls, bins=range(max(combo_nulls)+5), kde=True, color='gray')
    plt.axvline(actual_combo, color='red', linestyle='--', label=f'Actual Score ({actual_combo})')
    plt.title(f'Null Distribution (Combo) - p={combo_pval:.4f}')
    plt.xlabel('Rescue Score')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(os.path.join("aga-li-zn/results/figures", 'null_distribution_combo.png'), dpi=300)
    plt.close()

if __name__ == "__main__":
    run_null_model(
        "aga-li-zn/results/tables/DEG_table_full.tsv",
        "aga-li-zn/data/raw",
        "aga-li-zn/results/figures"
    )
