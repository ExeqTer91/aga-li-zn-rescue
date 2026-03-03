import modal
import os

app = modal.App("aga-random-controls")

image = modal.Image.debian_slim().apt_install("wget", "gzip").pip_install(
    "pandas", "numpy", "scipy", "cmapPy", "requests", "h5py", "matplotlib", "seaborn"
)

@app.function(image=image, timeout=14400, memory=16384)
def run_random_controls_all_datasets(N_list, random_pert_dict, dataset_data):
    import pandas as pd
    import numpy as np
    import urllib.request
    from cmapPy.pandasGEXpress.parse import parse
    import subprocess
    
    random_pert_ids = list(random_pert_dict.keys())
    
    targets = {
        "Li": ["BRD-M74254599"],
        "Zn_proxy": ["BRD-K09255212"],
        "GSK3": ["BRD-K16189898", "BRD-K37312348"]
    }
    
    all_target_ids = []
    for v in targets.values():
        all_target_ids.extend(v)
        
    all_pert_ids = random_pert_ids + all_target_ids
    
    print("Downloading sig_info and gene_info...")
    sig_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz"
    
    import time
    def download_with_retry(url, filename, max_retries=5):
        for attempt in range(max_retries):
            try:
                urllib.request.urlretrieve(url, filename)
                return True
            except Exception as e:
                print(f"Download failed (attempt {attempt+1}/{max_retries}): {e}")
                time.sleep(5)
        return False
        
    download_with_retry(sig_url, "sig_info.txt.gz")
    sig_info = pd.read_csv("sig_info.txt.gz", sep="\t", compression="gzip", low_memory=False)
    
    gene_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz"
    download_with_retry(gene_url, "gene_info.txt.gz")
    gene_info = pd.read_csv("gene_info.txt.gz", sep="\t", compression="gzip")
    gene_id_to_symbol = dict(zip(gene_info['pr_gene_id'].astype(str), gene_info['pr_gene_symbol']))
    
    relevant_sigs = sig_info[sig_info['pert_id'].isin(all_pert_ids) & (sig_info['pert_type'] == 'trt_cp')]
    sig_ids_to_extract = relevant_sigs['sig_id'].tolist()
    
    print(f"Found {len(sig_ids_to_extract)} signatures to extract from GCTX.")
    
    gctx_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"
    gctx_file = "level5.gctx"
    
    print("Downloading GCTX (21GB)...")
    subprocess.run(["wget", "-q", "-O", gctx_file + ".gz", gctx_url])
    print("Unzipping GCTX...")
    subprocess.run(["gunzip", gctx_file + ".gz"])
    
    print("Parsing relevant columns from GCTX using cmapPy...")
    gct = parse(gctx_file, cid=sig_ids_to_extract)
    data = gct.data_df
    
    print("Mapping gene symbols...")
    data.index = data.index.astype(str).map(gene_id_to_symbol)
    data = data[data.index.notna()]
    
    print("Computing consensus signatures...")
    pert_signatures = {}
    for pert in all_pert_ids:
        group_sigs = relevant_sigs[relevant_sigs['pert_id'] == pert]['sig_id'].tolist()
        group_sigs = [s for s in group_sigs if s in data.columns]
        if group_sigs:
            group_matrix = data[group_sigs]
            consensus = group_matrix.median(axis=1)
            pert_signatures[pert] = consensus
            
    datasets = ["GSE90594", "GSE36169", "GSE66663"]
    all_results_long = []
    all_results_summary = []
    
    for dataset_name in datasets:
        print(f"--- Processing {dataset_name} ---")
        
        aga_up = set(dataset_data[dataset_name]["aga_up"])
        aga_down = set(dataset_data[dataset_name]["aga_down"])
        symbol_to_fc = dataset_data[dataset_name]["symbol_to_fc"]
        background_genes = dataset_data[dataset_name]["background_genes"]
        
        def calc_weighted_rescue(drug_up, drug_down):
            rescue_up = aga_up.intersection(drug_down)
            rescue_down = aga_down.intersection(drug_up)
            
            score = 0
            for gene in rescue_up:
                score += abs(symbol_to_fc.get(gene, 0))
            for gene in rescue_down:
                score += abs(symbol_to_fc.get(gene, 0))
            return score

        for N in N_list:
            print(f"  N = {N}")
            # Target scores
            target_up_down = {}
            for label, pert_list in targets.items():
                target_cons = []
                for p in pert_list:
                    if p in pert_signatures:
                        target_cons.append(pert_signatures[p])
                if target_cons:
                    combo = pd.concat(target_cons, axis=1).median(axis=1)
                    target_up_down[label] = {
                        "up": set(combo.sort_values(ascending=False).head(N).index).intersection(background_genes),
                        "down": set(combo.sort_values(ascending=True).head(N).index).intersection(background_genes)
                    }
            
            if "Li" in target_up_down and "Zn_proxy" in target_up_down:
                target_up_down["Li_Zn"] = {
                    "up": target_up_down["Li"]["up"].union(target_up_down["Zn_proxy"]["up"]),
                    "down": target_up_down["Li"]["down"].union(target_up_down["Zn_proxy"]["down"])
                }
                
            if "GSK3" in target_up_down and "Zn_proxy" in target_up_down:
                target_up_down["GSK3_Zn"] = {
                    "up": target_up_down["GSK3"]["up"].union(target_up_down["Zn_proxy"]["up"]),
                    "down": target_up_down["GSK3"]["down"].union(target_up_down["Zn_proxy"]["down"])
                }
                
            target_scores = {}
            for label, ud in target_up_down.items():
                target_scores[label] = calc_weighted_rescue(ud["up"], ud["down"])
                all_results_long.append({
                    "dataset": dataset_name,
                    "N": N,
                    "pert_id": label,
                    "pert_iname": label,
                    "score": target_scores[label],
                    "is_target": 1
                })
                
            # Random scores
            random_scores = []
            for p in random_pert_ids:
                if p in pert_signatures:
                    cons = pert_signatures[p]
                    up = set(cons.sort_values(ascending=False).head(N).index).intersection(background_genes)
                    down = set(cons.sort_values(ascending=True).head(N).index).intersection(background_genes)
                    
                    score = calc_weighted_rescue(up, down)
                    random_scores.append(score)
                    
                    iname = random_pert_dict.get(p, "Unknown")
                    
                    all_results_long.append({
                        "dataset": dataset_name,
                        "N": N,
                        "pert_id": p,
                        "pert_iname": iname,
                        "score": score,
                        "is_target": 0
                    })
                    
            random_scores = np.array(random_scores)
            rand_median = float(np.median(random_scores))
            rand_mean = float(np.mean(random_scores))
            rand_sd = float(np.std(random_scores))
            
            for label, score in target_scores.items():
                percentile = float((np.sum(random_scores <= score) / len(random_scores)) * 100)
                emp_p = float((1 + np.sum(random_scores >= score)) / (len(random_scores) + 1))
                rank = int(np.sum(random_scores >= score) + 1)
                
                all_results_summary.append({
                    "dataset": dataset_name,
                    "N": N,
                    "condition": label,
                    "K_random": len(random_scores),
                    "score_target": score,
                    "random_median": rand_median,
                    "random_mean": rand_mean,
                    "random_sd": rand_sd,
                    "percentile": percentile,
                    "empirical_p": emp_p,
                    "rank": rank
                })
                
    return all_results_long, all_results_summary

@app.local_entrypoint()
def main():
    import pandas as pd
    import os
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    print("Loading local data to pass to Modal...")
    random_perts_df = pd.read_csv("aga-li-zn/data/processed/random_pert_ids.tsv", sep="\t")
    random_pert_ids = random_perts_df['pert_id'].tolist()
    random_pert_inames = random_perts_df['pert_iname'].tolist()
    random_pert_dict = dict(zip(random_pert_ids, random_pert_inames))
    
    datasets = ["GSE90594", "GSE36169", "GSE66663"]
    dataset_data = {}
    
    for dataset_name in datasets:
        if dataset_name == "GSE90594":
            deg_file = "aga-li-zn/results/tables/deg_full_mapped.tsv"
        else:
            deg_file = f"aga-li-zn/results/{dataset_name}/tables/DEG_table_full_mapped.tsv"
            
        df = pd.read_csv(deg_file, sep="\t")
        if 'Symbol' not in df.columns:
            df['Symbol'] = df['Gene']
            
        aga_up = list(df[df["direction"] == "Up"]["Symbol"])
        aga_down = list(df[df["direction"] == "Down"]["Symbol"])
        symbol_to_fc = df.groupby("Symbol")["log2FC"].mean().to_dict()
        background_genes = list(df["Symbol"].unique())
        
        dataset_data[dataset_name] = {
            "aga_up": aga_up,
            "aga_down": aga_down,
            "symbol_to_fc": symbol_to_fc,
            "background_genes": background_genes
        }
    
    print("Submitting Modal job for K=200 random compounds comparison (N=100, 300)...")
    res_long, res_summary = run_random_controls_all_datasets.remote([100, 300], random_pert_dict, dataset_data)
    
    df_long = pd.DataFrame(res_long)
    df_summary = pd.DataFrame(res_summary)
    
    os.makedirs("aga-li-zn/results/tables", exist_ok=True)
    os.makedirs("aga-li-zn/results/figures", exist_ok=True)
    
    df_long.to_csv("aga-li-zn/results/tables/random_drug_controls_scores_long.tsv", sep="\t", index=False)
    df_summary.to_csv("aga-li-zn/results/tables/random_drug_controls_summary.tsv", sep="\t", index=False)
    
    print("Saved tables. Generating plots...")
    
    # 1. Violin plots
    for dataset in datasets:
        for N in [100, 300]:
            sub_long = df_long[(df_long['dataset'] == dataset) & (df_long['N'] == N)]
            random_only = sub_long[sub_long['is_target'] == 0]
            targets_only = sub_long[sub_long['is_target'] == 1]
            
            plt.figure(figsize=(8, 6))
            sns.violinplot(data=random_only, y='score', color='lightgray', inner='quartile')
            sns.stripplot(data=random_only, y='score', color='gray', alpha=0.5, jitter=True, size=3)
            
            colors = {"Li": "blue", "Zn_proxy": "orange", "Li_Zn": "red", "GSK3_Zn": "purple", "GSK3": "green"}
            
            for _, row in targets_only.iterrows():
                plt.axhline(row['score'], color=colors.get(row['pert_id'], 'black'), linestyle='--', 
                            label=f"{row['pert_id']} (p={df_summary[(df_summary['dataset']==dataset) & (df_summary['N']==N) & (df_summary['condition']==row['pert_id'])]['empirical_p'].values[0]:.3f})")
                
            plt.title(f"Random Controls Distribution - {dataset} (N={N})")
            plt.ylabel("Weighted Rescue Score")
            plt.legend(loc="upper right")
            plt.tight_layout()
            plt.savefig(f"aga-li-zn/results/figures/random_controls_violin_{dataset}_N{N}.png", dpi=300)
            plt.close()
            
    # 2. Forest plot style for percentiles
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df_summary, x='percentile', y='dataset', hue='condition', style='N', s=100, palette='Set1')
    plt.axvline(95, color='red', linestyle=':', label='Top 5%')
    plt.axvline(99, color='darkred', linestyle='--', label='Top 1%')
    plt.xlim(0, 105)
    plt.title("Percentile Rank of Target Treatments vs 200 Random LINCS Compounds")
    plt.xlabel("Percentile (higher is better rescue)")
    plt.ylabel("Dataset")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("aga-li-zn/results/figures/random_controls_percentile.png", dpi=300)
    plt.close()
    
    print("Done! All results and plots saved.")

if __name__ == "__main__":
    main()
