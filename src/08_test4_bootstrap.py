import pandas as pd
import numpy as np
import os
import random
from scipy import stats
from statsmodels.stats.multitest import multipletests

def run_bootstrap_test(expr_file, pheno_file, gpl_file, out_dir, n_bootstraps=1000):
    print("Running Dataset Imbalance Control (Bootstrap) for GSE66663...")
    os.makedirs(out_dir, exist_ok=True)
    
    expr = pd.read_csv(expr_file, sep='\t', index_col=0)
    pheno = pd.read_csv(pheno_file, sep='\t')
    
    aga_samples = pheno[pheno['Group'] == 'AGA']['sample_id'].tolist()
    ctrl_samples = pheno[pheno['Group'] == 'Control']['sample_id'].tolist()
    
    gpl = pd.read_csv(gpl_file, sep='\t', low_memory=False, skiprows=0, on_bad_lines='skip')
    symbol_col = next((col for col in gpl.columns if 'symbol' in str(col).lower() or ('gene' in str(col).lower() and 'name' not in str(col).lower())), None)
    probe_to_symbol = dict(zip(gpl['ID'], gpl[symbol_col])) if symbol_col else {}
    
    def load_sig(filename):
        path = os.path.join("aga-li-zn/data/raw", filename)
        if os.path.exists(path):
            with open(path, "r") as f:
                return set([line.strip() for line in f if line.strip()])
        return set()
        
    li_up = load_sig("li_up.txt")
    li_down = load_sig("li_down.txt")
    zn_up = load_sig("zn_up.txt")
    zn_down = load_sig("zn_down.txt")
    combo_up = li_up.union(zn_up)
    combo_down = li_down.union(zn_down)
    
    bootstrap_results = []
    
    for i in range(n_bootstraps):
        if i % 100 == 0:
            print(f"Bootstrap {i}/{n_bootstraps}")
            
        boot_aga = random.choices(aga_samples, k=len(aga_samples))
        boot_ctrl = random.choices(ctrl_samples, k=len(ctrl_samples))
        
        aga_expr = expr[boot_aga]
        ctrl_expr = expr[boot_ctrl]
        
        log2fc = aga_expr.mean(axis=1) - ctrl_expr.mean(axis=1)
        
        t_stat, p_val = stats.ttest_ind(aga_expr, ctrl_expr, axis=1, equal_var=False)
        
        boot_deg = pd.DataFrame({'Gene': expr.index, 'log2FC': log2fc, 'pvalue': p_val})
        boot_deg = boot_deg.dropna()
        _, fdr, _, _ = multipletests(boot_deg['pvalue'], alpha=0.05, method='fdr_bh')
        boot_deg['FDR'] = fdr
        
        boot_deg['Symbol'] = boot_deg['Gene'].map(probe_to_symbol)
        boot_deg = boot_deg.dropna(subset=['Symbol'])
        boot_deg['Symbol'] = boot_deg['Symbol'].apply(lambda x: str(x).split('///')[0].strip())
        
        boot_deg['direction'] = 'Not Sig'
        boot_deg.loc[(boot_deg['pvalue'] < 0.05) & (boot_deg['log2FC'] > 0), 'direction'] = 'Up'
        boot_deg.loc[(boot_deg['pvalue'] < 0.05) & (boot_deg['log2FC'] < 0), 'direction'] = 'Down'
        
        aga_up_set = set(boot_deg[boot_deg['direction'] == 'Up']['Symbol'])
        aga_down_set = set(boot_deg[boot_deg['direction'] == 'Down']['Symbol'])
        
        symbol_to_fc = boot_deg.groupby("Symbol")["log2FC"].mean().to_dict()
        
        def calc_weighted_rescue(drug_up, drug_down):
            rescue_up = aga_up_set.intersection(drug_down)
            rescue_down = aga_down_set.intersection(drug_up)
            
            score = 0
            for gene in rescue_up:
                score += abs(symbol_to_fc.get(gene, 0))
            for gene in rescue_down:
                score += abs(symbol_to_fc.get(gene, 0))
            return score
            
        actual_combo = calc_weighted_rescue(combo_up, combo_down)
        bootstrap_results.append(actual_combo)
        
    res_df = pd.DataFrame({'Bootstrap': range(n_bootstraps), 'Weighted_Rescue_Score': bootstrap_results})
    out_file = os.path.join(out_dir, "robustness_test4_bootstrap_GSE66663.tsv")
    res_df.to_csv(out_file, sep='\t', index=False)
    
    mean_score = np.mean(bootstrap_results)
    std_score = np.std(bootstrap_results)
    ci_lower = np.percentile(bootstrap_results, 2.5)
    ci_upper = np.percentile(bootstrap_results, 97.5)
    
    print("====================================")
    print(f"RESULTS FOR GSE66663 BOOTSTRAP (1000 iter)")
    print(f"Mean Bootstrapped Rescue Score: {mean_score:.3f}")
    print(f"Std Dev: {std_score:.3f}")
    print(f"95% Confidence Interval: [{ci_lower:.3f}, {ci_upper:.3f}]")
    print("If CI > 0, the rescue is stable despite small sample size.")
    print("====================================")

if __name__ == "__main__":
    run_bootstrap_test(
        "data/GSE66663/expr_raw.tsv",
        "data/GSE66663/pheno.tsv",
        "data/GSE66663/GPL10558_annotation.tsv",
        "aga-li-zn/results/GSE66663/tables"
    )
