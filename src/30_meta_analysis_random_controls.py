import pandas as pd
import numpy as np
from scipy.stats import chi2, norm
import matplotlib.pyplot as plt
import os

df = pd.read_csv("aga-li-zn/results/tables/random_drug_controls_summary.tsv", sep="\t")

conditions_to_meta = ["Li_Zn", "GSK3_Zn", "Li", "Zn_proxy", "GSK3"]

meta_rows = []

for cond in conditions_to_meta:
    sub = df[df["condition"] == cond]
    if sub.empty:
        continue
    for N_val in sub["N"].unique():
        grp = sub[sub["N"] == N_val].copy()
        k = len(grp)
        if k < 2:
            continue

        percentiles = grp["percentile"].values
        p_values = grp["empirical_p"].values
        p_values = np.clip(p_values, 1e-10, 1 - 1e-10)

        median_pct = float(np.median(percentiles))
        mean_pct = float(np.mean(percentiles))

        chi2_stat = -2 * np.sum(np.log(p_values))
        dof = 2 * k
        fisher_p = float(1 - chi2.cdf(chi2_stat, dof))

        z_values = norm.ppf(1 - p_values)
        z_meta = float(np.sum(z_values) / np.sqrt(k))
        stouffer_p = float(2 * (1 - norm.cdf(abs(z_meta))))

        z_mean = np.mean(z_values)
        Q = float(np.sum((z_values - z_mean) ** 2))
        I2 = max(0, (Q - (k - 1)) / Q) if Q > 0 else 0.0

        meta_rows.append({
            "condition": cond,
            "N": N_val,
            "k_datasets": k,
            "median_percentile": round(median_pct, 2),
            "mean_percentile": round(mean_pct, 2),
            "fisher_chi2": round(chi2_stat, 4),
            "fisher_p": round(fisher_p, 6),
            "stouffer_z": round(z_meta, 4),
            "stouffer_p": round(stouffer_p, 6),
            "I2": round(I2, 4),
        })

meta_df = pd.DataFrame(meta_rows)
os.makedirs("aga-li-zn/results/tables", exist_ok=True)
os.makedirs("aga-li-zn/results/figures", exist_ok=True)
meta_df.to_csv("aga-li-zn/results/tables/meta_analysis_random_controls.tsv", sep="\t", index=False)
print("Meta-analysis results:")
print(meta_df.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

dataset_labels = ["GSE90594", "GSE36169", "GSE66663"]
colors_cond = {
    "Li_Zn": "#e41a1c",
    "GSK3_Zn": "#984ea3",
    "Li": "#377eb8",
    "Zn_proxy": "#ff7f00",
    "GSK3": "#4daf4a",
}

for ax_idx, N_val in enumerate([100, 300]):
    ax = axes[ax_idx]
    y_pos = 0
    yticks = []
    ytick_labels = []

    for cond in conditions_to_meta:
        sub = df[(df["condition"] == cond) & (df["N"] == N_val)]
        if sub.empty:
            continue

        for _, row in sub.iterrows():
            ax.barh(y_pos, row["percentile"], color=colors_cond.get(cond, "gray"), alpha=0.7, height=0.6)
            ax.text(row["percentile"] + 1, y_pos, f'{row["percentile"]:.1f}%', va='center', fontsize=8)
            yticks.append(y_pos)
            ytick_labels.append(f'{row["dataset"]} — {cond}')
            y_pos += 1

        meta_row = meta_df[(meta_df["condition"] == cond) & (meta_df["N"] == N_val)]
        if not meta_row.empty:
            mean_pct = meta_row["mean_percentile"].values[0]
            ax.axvline(mean_pct, color=colors_cond.get(cond, "gray"), linestyle="--", alpha=0.5, linewidth=1)

        y_pos += 0.5

    ax.axvline(95, color="red", linestyle=":", alpha=0.6, linewidth=1)
    ax.axvline(99, color="darkred", linestyle="--", alpha=0.6, linewidth=1)
    ax.text(95.5, y_pos - 1, "95th", fontsize=7, color="red")
    ax.text(99.2, y_pos - 1, "99th", fontsize=7, color="darkred")

    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels, fontsize=8)
    ax.set_xlabel("Percentile among 200 random LINCS compounds")
    ax.set_title(f"N = {N_val}")
    ax.set_xlim(0, 108)
    ax.invert_yaxis()

plt.suptitle("Meta-Analysis: Target Treatment Percentiles vs Random Drug Controls", fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig("aga-li-zn/results/figures/meta_random_controls_forest.png", dpi=300, bbox_inches="tight")
plt.close()
print("\nForest plot saved to results/figures/meta_random_controls_forest.png")
