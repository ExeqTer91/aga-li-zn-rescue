# AGA Li+Zn Transcriptomic Rescue

**Systems pharmacology identifies dermal papilla-specific transcriptomic reversal of androgenetic alopecia by lithium-zinc pathway convergence**

## Overview

This repository contains the complete computational pipeline and results for an *in silico* study testing whether Lithium (Li) and Zinc (Zn) treatment signatures can "rescue" the androgenetic alopecia (AGA) disease transcriptome. The analysis integrates public GEO microarray datasets with LINCS L1000 perturbational signatures to compute weighted directional convergence scores, validated by extensive robustness testing.

## Key Findings

| Experiment | Result | Status |
|---|---|---|
| Cross-dataset permutation (10,000 iterations) | Li+Zn: GSE90594 z=3.42, p=0.009; GSE36169 z=3.99, p=0.017; GSE66663 z=9.79, p<0.001 | Passed |
| Leave-10%-out stress test | Signal distributed, not single-gene driven | Passed |
| GSE66663 bootstrap (1,000 iterations) | Mean=14.09, 95% CI [13.34, 15.19] | Passed |
| Random drug benchmarking (K=200) | GSE66663 Li+Zn percentile 99.5%, p=0.010 | Passed |
| Meta-analysis | Li+Zn N=300: Fisher p=0.0088, Stouffer Z=2.58, I²=0 | Passed |
| Mechanistic positive control (GSK3i+Zn) | N=300: Fisher p=0.0017, Stouffer Z=3.18, I²=0 | Passed |

## Datasets

| Dataset | Tissue | Comparison | Role |
|---|---|---|---|
| GSE90594 | Vertex scalp biopsies | 14 AGA vs 14 control | Primary bulk discovery |
| GSE36169 | Paired bald/haired scalp | 5 bald vs 5 haired | Bulk replication |
| GSE66663 | Immortalized DP cells | 3 balding vs 3 non-balding | DP-specific amplification |

## LINCS Compounds

- **Lithium chloride**: BRD-M74254599 (direct)
- **Zinc proxy**: Clioquinol (BRD-K09255212, zinc ionophore)
- **GSK3 inhibitors** (positive control): BRD-K16189898, BRD-K37312348

## Repository Structure

```
aga-li-zn/
├── src/                          # Analysis scripts (numbered pipeline)
│   ├── 00_download_geo.py        # Download GEO datasets
│   ├── 01_preprocess_qc.py       # QC and PCA
│   ├── 02_deg_analysis.py        # Differential expression
│   ├── 03_gsea.py                # Gene set enrichment
│   ├── 04_extract_lincs_*.py     # LINCS signature extraction (Modal)
│   ├── 05_convergence_score.py   # Weighted directional rescue
│   ├── 06_null_models.py         # 10,000-permutation null models
│   ├── 07_summary.py             # Cross-dataset summary
│   ├── 08_test4_bootstrap.py     # Bootstrap resampling (GSE66663)
│   ├── 09_test2_signature_length.py  # Signature length sensitivity
│   ├── 10_dp_module_rescue.py    # DP module gene-level rescue
│   ├── 20_sample_random_compounds.py # Sample K=200 random LINCS drugs
│   ├── 21_submit_modal_random_controls.py # Random drug controls (Modal)
│   └── 30_meta_analysis_random_controls.py # Cross-dataset meta-analysis
├── results/
│   ├── tables/                   # TSV result tables
│   │   ├── cross_dataset_modal_summary.tsv
│   │   ├── meta_analysis_random_controls.tsv
│   │   ├── random_drug_controls_summary.tsv
│   │   ├── test2_signature_length_results.tsv
│   │   ├── robustness_test3_leave10pct_GSE90594.tsv
│   │   └── ...
│   ├── figures/                  # Publication-quality figures
│   │   ├── meta_random_controls_forest.png
│   │   ├── random_controls_percentile.png
│   │   ├── random_controls_violin_*.png
│   │   ├── volcano.png
│   │   ├── heatmap_top50.png
│   │   └── ...
│   ├── GSE36169/                 # Dataset-specific results
│   ├── GSE66663/                 # Dataset-specific results
│   └── limitations_section.md    # Submission-ready limitations text
└── data/
    └── processed/
        └── random_pert_ids.tsv   # K=200 sampled compound IDs
```

## Methods Summary

1. **DEG Analysis**: Gene-wise group comparison with Benjamini-Hochberg FDR correction
2. **LINCS Extraction**: Level 5 consensus profiles parsed from GSE92742 GCTX (21GB)
3. **Rescue Scoring**: Weighted directional convergence (fold-change magnitude weighting)
4. **Null Models**: 10,000 gene-set permutations per dataset/condition (Modal cloud)
5. **Robustness**: Leave-10%-out, bootstrap resampling, signature length sensitivity
6. **Random Controls**: Benchmarking against 200 randomly sampled LINCS compounds
7. **Meta-Analysis**: Fisher combined p-values, Stouffer Z-scores, I² heterogeneity

## Infrastructure

- **Replit**: Primary development and analysis environment
- **Modal**: Cloud compute for permutation testing and LINCS GCTX processing
- **Python**: pandas, numpy, scipy, cmapPy, matplotlib, seaborn

## Citation

Manuscript in preparation for *Frontiers in Pharmacology* (Systems Pharmacology section).

## License

MIT
