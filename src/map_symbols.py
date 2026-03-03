import pandas as pd
import numpy as np
import os
import argparse

def map_symbols(deg_file, gpl_file, out_file):
    print(f"Mapping symbols for {deg_file} using {gpl_file}...")
    deg_df = pd.read_csv(deg_file, sep='\t')
    
    if gpl_file and os.path.exists(gpl_file):
        header_idx = 0
        with open(gpl_file, 'r', encoding='utf8', errors='ignore') as f:
            for i, line in enumerate(f):
                if line.startswith('ID') or line.startswith('Probe'):
                    header_idx = i
                    break
                    
        gpl = pd.read_csv(gpl_file, sep='\t', low_memory=False, on_bad_lines='skip', skiprows=header_idx)
        
        symbol_col = None
        for col in gpl.columns:
            c = str(col).lower()
            if c in ['gene symbol', 'symbol', 'gene_symbol']:
                symbol_col = col
                break
        if not symbol_col:
            symbol_col = next((col for col in gpl.columns if 'symbol' in str(col).lower() or ('gene' in str(col).lower() and 'name' not in str(col).lower())), None)
            
        print(f"Found symbol column: {symbol_col}")
        
        if symbol_col:
            probe_to_symbol = dict(zip(gpl.iloc[:, 0], gpl[symbol_col]))
            deg_df['Symbol'] = deg_df['Gene'].map(probe_to_symbol)
            deg_df = deg_df.dropna(subset=['Symbol'])
            deg_df['Symbol'] = deg_df['Symbol'].apply(lambda x: str(x).split('///')[0].strip())
        else:
            deg_df['Symbol'] = deg_df['Gene']
    else:
        deg_df['Symbol'] = deg_df['Gene']

    deg_df.to_csv(out_file, sep='\t', index=False)
    print(f"Saved mapped DEG to {out_file}")
    print(f"Mapped genes count: {len(deg_df['Symbol'].unique())}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--deg', required=True)
    parser.add_argument('--gpl', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    map_symbols(args.deg, args.gpl, args.out)
