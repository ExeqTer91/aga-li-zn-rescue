import GEOparse
import pandas as pd
import os
import argparse

def download_geo(gse_id, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    print(f"Downloading {gse_id}...")
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=out_dir)
    except Exception as e:
        print(f"Error downloading {gse_id}: {e}")
        return
    
    pheno_data = []
    for gsm_name, gsm in gse.gsms.items():
        metadata = gsm.metadata
        row = {'sample_id': gsm_name}
        for key, values in metadata.items():
            row[key] = " | ".join(values)
        pheno_data.append(row)
    
    pheno_df = pd.DataFrame(pheno_data)
    pheno_file = os.path.join(out_dir, 'pheno_raw.tsv')
    pheno_df.to_csv(pheno_file, sep='\t', index=False)
    print(f"Phenotype data saved to {pheno_file}")
    
    expr_df = gse.pivot_samples('VALUE')
    expr_file = os.path.join(out_dir, 'expr_raw.tsv')
    expr_df.to_csv(expr_file, sep='\t')
    print(f"Expression data saved to {expr_file}")
    
    for gpl_name, gpl in gse.gpls.items():
        print(f"Platform identified: {gpl_name}")
        gpl.table.to_csv(os.path.join(out_dir, f"{gpl_name}_annotation.tsv"), sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", type=str, required=True)
    parser.add_argument("--out", type=str, required=True)
    args = parser.parse_args()
    
    download_geo(args.gse, args.out)
