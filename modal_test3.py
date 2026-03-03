import modal
import os
import pandas as pd
import numpy as np

app = modal.App("aga-li-zn-robustness")

image = modal.Image.debian_slim().pip_install("pandas", "numpy", "scipy")

@app.function(image=image, timeout=3600)
def run_leave_10pct_out_permutations(n_iterations: int, n_permutations: int, deg_symbols: list, li_up: set, li_down: set, zn_up: set, zn_down: set, aga_up: set, aga_down: set, symbol_to_fc: dict):
    import numpy as np
    import random
    
    bg_array = np.array(deg_symbols)
    
    def calc_weighted_rescue(drug_up, drug_down):
        rescue_up = aga_up.intersection(drug_down)
        rescue_down = aga_down.intersection(drug_up)
        
        score = 0
        for gene in rescue_up:
            score += abs(symbol_to_fc.get(gene, 0))
        for gene in rescue_down:
            score += abs(symbol_to_fc.get(gene, 0))
        return score
    
    results = []
    
    for i in range(n_iterations):
        def drop_10pct(sig_set):
            sig_list = list(sig_set)
            drop_count = max(1, int(len(sig_list) * 0.1))
            random.shuffle(sig_list)
            return set(sig_list[drop_count:])
            
        cur_li_up = drop_10pct(li_up)
        cur_li_down = drop_10pct(li_down)
        cur_zn_up = drop_10pct(zn_up)
        cur_zn_down = drop_10pct(zn_down)
        
        cur_combo_up = cur_li_up.union(cur_zn_up)
        cur_combo_down = cur_li_down.union(cur_zn_down)
        
        actual_combo = calc_weighted_rescue(cur_combo_up, cur_combo_down)
        
        combo_nulls = []
        c_up_len = len(cur_combo_up)
        c_down_len = len(cur_combo_down)
        
        for _ in range(n_permutations):
            np.random.shuffle(bg_array)
            rand_up = set(bg_array[:c_up_len])
            rand_down = set(bg_array[c_up_len:c_up_len+c_down_len])
            combo_nulls.append(calc_weighted_rescue(rand_up, rand_down))
            
        combo_pval = (np.sum(np.array(combo_nulls) >= actual_combo) + 1) / (n_permutations + 1)
        combo_z = (actual_combo - np.mean(combo_nulls)) / (np.std(combo_nulls) + 1e-9)
        
        results.append({
            "iteration": i,
            "actual_combo": actual_combo,
            "combo_z": combo_z,
            "combo_pval": combo_pval
        })
        
    return results

@app.local_entrypoint()
def main(dataset: str = "GSE90594"):
    if dataset == "GSE90594":
        deg_file = "aga-li-zn/results/tables/deg_full.tsv"
    else:
        deg_file = f"aga-li-zn/results/{dataset}/tables/DEG_table_full.tsv"
        
    deg_df = pd.read_csv(deg_file.replace('.tsv', '_mapped.tsv'), sep='\t')
    if 'Symbol' not in deg_df.columns:
        deg_df['Symbol'] = deg_df['Gene']

    aga_up = set(deg_df[deg_df["direction"] == "Up"]["Symbol"])
    aga_down = set(deg_df[deg_df["direction"] == "Down"]["Symbol"])
    symbol_to_fc = deg_df.groupby("Symbol")["log2FC"].mean().to_dict()
    background_genes = list(deg_df["Symbol"].unique())
    
    def load_sig(filename):
        path = os.path.join("aga-li-zn/data/raw", filename)
        if os.path.exists(path):
            with open(path, "r") as f:
                return set([line.strip() for line in f if line.strip()])
        return set()
        
    li_up = load_sig("li_up.txt").intersection(background_genes)
    li_down = load_sig("li_down.txt").intersection(background_genes)
    zn_up = load_sig("zn_up.txt").intersection(background_genes)
    zn_down = load_sig("zn_down.txt").intersection(background_genes)
    
    print(f"Running robustness Test 3 (Leave-10%-out) for {dataset}...")
    results = run_leave_10pct_out_permutations.remote(
        n_iterations=50, 
        n_permutations=1000, 
        deg_symbols=background_genes,
        li_up=li_up, li_down=li_down, zn_up=zn_up, zn_down=zn_down,
        aga_up=aga_up, aga_down=aga_down, symbol_to_fc=symbol_to_fc
    )
    
    res_df = pd.DataFrame(results)
    out_dir = f"aga-li-zn/results/{dataset}/tables" if dataset != "GSE90594" else "aga-li-zn/results/tables"
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"robustness_test3_leave10pct_{dataset}.tsv")
    res_df.to_csv(out_file, sep='\t', index=False)
    
    print("====================================")
    print(f"RESULTS FOR {dataset} (Leave-10%-out)")
    print(f"Mean Z-score across 100 iterations: {res_df['combo_z'].mean():.3f}")
    print(f"Min Z-score observed: {res_df['combo_z'].min():.3f}")
    print(f"Max Z-score observed: {res_df['combo_z'].max():.3f}")
    print(f"Mean P-value: {res_df['combo_pval'].mean():.4f}")
    print(f"Saved details to: {out_file}")
    print("====================================")
