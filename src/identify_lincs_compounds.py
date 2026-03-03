import pandas as pd
import gzip
import re

file_path = "aga-li-zn/data/raw/GSE92742_Broad_LINCS_sig_info.txt.gz"

try:
    sig_info = pd.read_csv(file_path, sep="\t", compression="gzip")
    print("Columns:", sig_info.columns.tolist())

    # lithium
    li_hits = sig_info[sig_info["pert_iname"].astype(str).str.contains("lithium", case=False, na=False)]
    print("\nLithium matches:")
    if not li_hits.empty:
        print(li_hits[["pert_iname", "pert_id"]].drop_duplicates())
    else:
        print("None found.")

    # zinc
    zn_hits = sig_info[sig_info["pert_iname"].astype(str).str.contains("zinc", case=False, na=False)]
    print("\nZinc matches:")
    if not zn_hits.empty:
        print(zn_hits[["pert_iname", "pert_id"]].drop_duplicates())
    else:
        print("None found.")
        
except Exception as e:
    print("Error:", e)
