import pandas as pd
import numpy as np
import os
import random

def main():
    print("Loading sig_info.txt.gz...")
    # Read the signature info file
    sig_info_path = "sig_info.txt.gz"
    if not os.path.exists(sig_info_path):
        import urllib.request
        print("Downloading sig_info.txt.gz...")
        sig_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz"
        urllib.request.urlretrieve(sig_url, sig_info_path)
        
    sig_info = pd.read_csv(sig_info_path, sep="\t", compression="gzip", low_memory=False)
    
    # Target compounds to exclude from random pool
    targets = {
        "li": "BRD-M74254599",
        "zn_proxy": "BRD-K09255212",
        "gsk3_1": "BRD-K16189898",
        "gsk3_2": "BRD-K37312348"
    }
    target_ids = list(targets.values())
    
    print("Filtering for trt_cp compounds...")
    trt_cp_info = sig_info[sig_info['pert_type'] == 'trt_cp']
    
    # Count sig_ids per pert_id
    pert_counts = trt_cp_info['pert_id'].value_counts()
    
    # Filter: at least 10 sig_ids
    min_sig_ids = 10
    valid_perts = pert_counts[pert_counts >= min_sig_ids].index.tolist()
    
    # Exclude targets
    valid_perts = [p for p in valid_perts if p not in target_ids]
    
    print(f"Found {len(valid_perts)} compounds with >= {min_sig_ids} signatures.")
    
    # Sample K=200 compounds with fixed seed
    K = 200
    random.seed(42)
    
    if len(valid_perts) < K:
        print(f"Warning: only found {len(valid_perts)} valid compounds, using all.")
        sampled_perts = valid_perts
    else:
        sampled_perts = random.sample(valid_perts, K)
        
    # Get pert_inames for the sampled compounds
    pert_inames = {}
    for p in sampled_perts:
        iname = trt_cp_info[trt_cp_info['pert_id'] == p]['pert_iname'].iloc[0]
        pert_inames[p] = iname
        
    # Create dataframe
    df = pd.DataFrame({
        'pert_id': sampled_perts,
        'pert_iname': [pert_inames[p] for p in sampled_perts]
    })
    
    os.makedirs("aga-li-zn/data/processed", exist_ok=True)
    out_file = "aga-li-zn/data/processed/random_pert_ids.tsv"
    df.to_csv(out_file, sep="\t", index=False)
    print(f"Saved {len(df)} random compounds to {out_file}")

if __name__ == "__main__":
    main()
