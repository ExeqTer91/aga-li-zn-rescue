import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from scipy import stats
import os
import argparse

def preprocess_and_qc(expr_file, pheno_file, out_dir):
    os.makedirs(os.path.join(out_dir, "figures"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "tables"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "processed"), exist_ok=True)
    
    print(f"Loading data from {expr_file}...")
    expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    
    expr_df = expr_df.dropna()
    
    if expr_df.max().max() > 100:
        print("Applying log2 transformation...")
        expr_df = np.log2(expr_df + 1)
        
    print("Performing PCA...")
    scaler = StandardScaler()
    scaled_expr = scaler.fit_transform(expr_df.T)
    
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_expr)
    
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df['sample_id'] = expr_df.columns
    
    group_col = 'Group'
    if group_col not in pheno_df.columns:
        if 'title' in pheno_df.columns:
            pheno_df[group_col] = pheno_df['title'].apply(lambda x: 'AGA' if 'alopecia' in str(x).lower() or 'bab' in str(x).lower() else ('Control' if 'healthy' in str(x).lower() or 'ban' in str(x).lower() else 'Unknown'))
        else:
            pheno_df[group_col] = 'Unknown'
            
    pca_df = pca_df.merge(pheno_df[['sample_id', group_col]], on='sample_id', how='left')
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', hue=group_col, data=pca_df, palette='Set1')
    plt.title('PCA of Expression Data')
    plt.savefig(os.path.join(out_dir, 'figures', 'qc_pca.png'), dpi=300)
    plt.close()
    
    print("Calculating outlier scores...")
    distances = squareform(pdist(scaled_expr, metric='euclidean'))
    mean_dist = np.mean(distances, axis=1)
    z_scores = stats.zscore(mean_dist)
    
    outliers_df = pd.DataFrame({
        'sample_id': expr_df.columns,
        'mean_distance': mean_dist,
        'z_score': z_scores
    })
    
    outliers_df['is_outlier'] = outliers_df['z_score'].abs() > 2.5
    outliers_df.to_csv(os.path.join(out_dir, 'tables', 'outliers_report.tsv'), sep='\t', index=False)
    
    processed_expr_file = os.path.join(out_dir, 'processed', 'expr_processed.tsv')
    expr_df.to_csv(processed_expr_file, sep='\t')
    
    processed_pheno_file = os.path.join(out_dir, 'processed', 'pheno_processed.tsv')
    pheno_df[['sample_id', group_col]].to_csv(processed_pheno_file, sep='\t', index=False)
    
    print(f"Preprocessing complete. Files saved in {out_dir}/processed/")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--expr", required=True)
    parser.add_argument("--pheno", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    preprocess_and_qc(args.expr, args.pheno, args.outdir)
