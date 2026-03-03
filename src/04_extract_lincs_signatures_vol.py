import modal
import os

app = modal.App("extract-lincs-signatures-vol")
vol = modal.Volume.from_name("lincs-signatures-vol", create_if_missing=True)

image = modal.Image.debian_slim().apt_install("wget", "gzip").pip_install("pandas", "cmapPy", "requests", "h5py")

@app.function(image=image, timeout=7200, memory=16384, volumes={"/data": vol})
def extract_signatures():
    import pandas as pd
    import urllib.request
    from cmapPy.pandasGEXpress.parse import parse
    import subprocess
    
    print("Downloading sig_info.txt.gz...")
    sig_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz"
    urllib.request.urlretrieve(sig_url, "sig_info.txt.gz")
    
    sig_info = pd.read_csv("sig_info.txt.gz", sep="\t", compression="gzip", low_memory=False)
    
    targets = {
        "li": ["BRD-M74254599"],
        "zn_proxy": ["BRD-K09255212"],
        "gsk3_proxy": ["BRD-K16189898", "BRD-K37312348"]
    }
    
    print("Downloading gene_info.txt.gz...")
    gene_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz"
    urllib.request.urlretrieve(gene_url, "gene_info.txt.gz")
    gene_info = pd.read_csv("gene_info.txt.gz", sep="\t", compression="gzip")
    gene_id_to_symbol = dict(zip(gene_info['pr_gene_id'].astype(str), gene_info['pr_gene_symbol']))
    
    all_target_ids = []
    for v in targets.values():
        all_target_ids.extend(v)
        
    relevant_sigs = sig_info[(sig_info['pert_id'].isin(all_target_ids)) & (sig_info['pert_type'] == 'trt_cp')]
    sig_ids = relevant_sigs['sig_id'].tolist()
    print(f"Found {len(sig_ids)} total signature profiles.")
    
    gctx_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"
    gctx_file = "level5.gctx"
    
    print("Downloading GCTX (21GB)...")
    subprocess.run(["wget", "-q", "-O", gctx_file + ".gz", gctx_url])
    print("Unzipping GCTX...")
    subprocess.run(["gunzip", gctx_file + ".gz"])
    
    print("Parsing relevant columns from GCTX using cmapPy...")
    gct = parse(gctx_file, cid=sig_ids)
    data = gct.data_df
    
    print("Mapping gene symbols...")
    data.index = data.index.astype(str).map(gene_id_to_symbol)
    data = data[data.index.notna()]
    
    for label, brd_list in targets.items():
        group_sigs = relevant_sigs[relevant_sigs['pert_id'].isin(brd_list)]['sig_id'].tolist()
        group_sigs = [s for s in group_sigs if s in data.columns]
        
        print(f"{label}: found {len(group_sigs)} matching profiles in the GCTX")
        if not group_sigs:
            continue
            
        group_matrix = data[group_sigs]
        consensus = group_matrix.median(axis=1)
        
        for N in [100, 300, 500]:
            up = consensus.sort_values(ascending=False).head(N).index.tolist()
            down = consensus.sort_values(ascending=True).head(N).index.tolist()
            
            with open(f"/data/{label}_up_{N}.txt", "w") as f:
                f.write("\n".join(up) + "\n")
            with open(f"/data/{label}_down_{N}.txt", "w") as f:
                f.write("\n".join(down) + "\n")
    
    vol.commit()
    print("All signatures extracted and saved to volume.")

@app.local_entrypoint()
def main():
    print("Extracting signatures via Modal Volume...")
    extract_signatures.remote()
    print("Done. Use 'modal volume get lincs-signatures-vol /* .' to download.")
