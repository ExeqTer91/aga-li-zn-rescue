import modal
import os
import pandas as pd
import numpy as np

app = modal.App("aga-li-zn-null-models")

image = modal.Image.debian_slim().pip_install(
    "pandas",
    "numpy",
    "scipy",
    "scikit-learn"
)

@app.function(image=image, timeout=3600)
def run_null_models_on_modal(n_permutations: int, deg_symbols: list, log2fc_vals: list, li_up_len: int, li_down_len: int, zn_up_len: int, zn_down_len: int, actual_li_weighted: float, actual_zn_weighted: float, actual_combo_weighted: float, aga_up: set, aga_down: set, symbol_to_fc: dict):
    import numpy as np
    
    print(f"Running {n_permutations} permutations on Modal (Weighted)...")
    
    li_nulls = []
    zn_nulls = []
    combo_nulls = []
    
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
        
    for i in range(n_permutations):
        if i % 1000 == 0:
            print(f"Permutation {i}/{n_permutations}")
        
        np.random.shuffle(bg_array)
        
        rand_li_up = set(bg_array[:li_up_len])
        rand_li_down = set(bg_array[li_up_len:li_up_len+li_down_len])
        
        rand_zn_up = set(bg_array[:zn_up_len])
        rand_zn_down = set(bg_array[zn_up_len:zn_up_len+zn_down_len])
        
        li_nulls.append(calc_weighted_rescue(rand_li_up, rand_li_down))
        zn_nulls.append(calc_weighted_rescue(rand_zn_up, rand_zn_down))
        combo_nulls.append(calc_weighted_rescue(rand_li_up.union(rand_zn_up), rand_li_down.union(rand_zn_down)))
        
    li_pval = (np.sum(np.array(li_nulls) >= actual_li_weighted) + 1) / (n_permutations + 1)
    zn_pval = (np.sum(np.array(zn_nulls) >= actual_zn_weighted) + 1) / (n_permutations + 1)
    combo_pval = (np.sum(np.array(combo_nulls) >= actual_combo_weighted) + 1) / (n_permutations + 1)
    
    li_z = (actual_li_weighted - np.mean(li_nulls)) / (np.std(li_nulls) + 1e-9)
    zn_z = (actual_zn_weighted - np.mean(zn_nulls)) / (np.std(zn_nulls) + 1e-9)
    combo_z = (actual_combo_weighted - np.mean(combo_nulls)) / (np.std(combo_nulls) + 1e-9)
    
    return {
        "li_pval": li_pval, "zn_pval": zn_pval, "combo_pval": combo_pval,
        "li_z": li_z, "zn_z": zn_z, "combo_z": combo_z,
        "actual_li": actual_li_weighted, "actual_zn": actual_zn_weighted, "actual_combo": actual_combo_weighted
    }

@app.local_entrypoint()
def main(dataset: str = "GSE90594", n_permutations: int = 10000):
    import pandas as pd
    
    if dataset == "GSE90594":
        deg_file = "aga-li-zn/results/tables/deg_full.tsv"
        out_table = "aga-li-zn/results/tables/modal_convergence_stats_GSE90594.tsv"
        gpl_file = "data/raw/GPL17077_annotation.tsv"
    elif dataset == "GSE36169":
        deg_file = f"aga-li-zn/results/{dataset}/tables/DEG_table_full.tsv"
        out_table = f"aga-li-zn/results/{dataset}/tables/modal_convergence_stats_{dataset}.tsv"
        gpl_file = "data/GSE36169/GPL96_annotation.tsv"
    elif dataset == "GSE66663":
        deg_file = f"aga-li-zn/results/{dataset}/tables/DEG_table_full.tsv"
        out_table = f"aga-li-zn/results/{dataset}/tables/modal_convergence_stats_{dataset}.tsv"
        gpl_file = "data/GSE66663/GPL10558_annotation.tsv"
        
    deg_df = pd.read_csv(deg_file.replace('.tsv', '_mapped.tsv'), sep="\t")

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
    
    def calc_weighted_rescue(drug_up, drug_down):
        rescue_up = aga_up.intersection(drug_down)
        rescue_down = aga_down.intersection(drug_up)
        
        score = 0
        for gene in rescue_up:
            score += abs(symbol_to_fc.get(gene, 0))
        for gene in rescue_down:
            score += abs(symbol_to_fc.get(gene, 0))
        return score
        
    actual_li = calc_weighted_rescue(li_up, li_down)
    actual_zn = calc_weighted_rescue(zn_up, zn_down)
    actual_combo = calc_weighted_rescue(li_up.union(zn_up), li_down.union(zn_down))
    
    print(f"Calling Modal Function for {dataset}...")
    result = run_null_models_on_modal.remote(
        n_permutations=n_permutations, 
        deg_symbols=background_genes, 
        log2fc_vals=list(symbol_to_fc.values()),
        li_up_len=len(li_up), li_down_len=len(li_down), 
        zn_up_len=len(zn_up), zn_down_len=len(zn_down), 
        actual_li_weighted=actual_li, actual_zn_weighted=actual_zn, actual_combo_weighted=actual_combo,
        aga_up=aga_up, aga_down=aga_down,
        symbol_to_fc=symbol_to_fc
    )
    
    print(f"Modal computation finished for {dataset}.")
    
    res_df = pd.DataFrame({
        "Treatment": ["Li", "Zn", "Li+Zn Combo"],
        "Actual_Weighted_Score": [result["actual_li"], result["actual_zn"], result["actual_combo"]],
        "Z_Score": [result["li_z"], result["zn_z"], result["combo_z"]],
        "P_Value": [result["li_pval"], result["zn_pval"], result["combo_pval"]]
    })
    
    res_df.to_csv(out_table, sep="\t", index=False)
    print(f"Saved stats to {out_table}")
    print(res_df)
