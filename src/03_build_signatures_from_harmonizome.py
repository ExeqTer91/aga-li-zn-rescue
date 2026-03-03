import pandas as pd
import gzip
import os
import sys

file_path = "aga-li-zn/data/raw/Harmonizome_LINCS_L1000_CMAP_Chemical_Perturbation_Consensus_Signatures.tsv.gz"
out_dir = "aga-li-zn/data/signatures"
os.makedirs(out_dir, exist_ok=True)

print("Loading Harmonizome matrix...")
# The matrix has 2 rows of metadata (with `#`), and the 3rd row contains the perturbation names.
df = pd.read_csv(file_path, sep='\t', skiprows=2, index_col=0, low_memory=False)

print(f"Matrix shape: {df.shape}")

# Drop columns that are completely null or unnamed just in case
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

# Find columns containing lithium
li_cols = [c for c in df.columns if "lithium chloride" in str(c).lower() or "lithium carbonate" in str(c).lower() or "lithium acetate" in str(c).lower()]
print(f"Found {len(li_cols)} Lithium columns: {li_cols}")

# Find columns containing zinc
zn_cols = [c for c in df.columns if "zinc chloride" in str(c).lower() or "zinc sulfate" in str(c).lower() or "zinc acetate" in str(c).lower()]
print(f"Found {len(zn_cols)} Zinc columns: {zn_cols}")

if li_cols:
    li_scores = df[li_cols].mean(axis=1).sort_values(ascending=False)
    li_up = li_scores.head(500).index.dropna().tolist()
    li_down = li_scores.tail(500).index.dropna().tolist()
    
    with open(os.path.join(out_dir, "li_up.txt"), "w") as f:
        f.write("\n".join(str(g).split('///')[0].strip() for g in li_up) + "\n")
    with open(os.path.join(out_dir, "li_down.txt"), "w") as f:
        f.write("\n".join(str(g).split('///')[0].strip() for g in li_down) + "\n")
    print(f"Saved Li signatures: {len(li_up)} up, {len(li_down)} down")
else:
    print("WARNING: No Lithium columns found!")

if zn_cols:
    zn_scores = df[zn_cols].mean(axis=1).sort_values(ascending=False)
    zn_up = zn_scores.head(500).index.dropna().tolist()
    zn_down = zn_scores.tail(500).index.dropna().tolist()
    
    with open(os.path.join(out_dir, "zn_up.txt"), "w") as f:
        f.write("\n".join(str(g).split('///')[0].strip() for g in zn_up) + "\n")
    with open(os.path.join(out_dir, "zn_down.txt"), "w") as f:
        f.write("\n".join(str(g).split('///')[0].strip() for g in zn_down) + "\n")
    print(f"Saved Zn signatures: {len(zn_up)} up, {len(zn_down)} down")
else:
    print("WARNING: No Zinc columns found!")
