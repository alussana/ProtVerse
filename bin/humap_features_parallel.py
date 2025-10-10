#!/usr/bin/env python3

import sys
import pandas as pd
from multiprocessing import Pool


def create_lookup_vectorized(humap3_table):
    # Create both directions
    forward = humap3_table[['gene_1', 'gene_2', 'humap_score']].copy()
    backward = humap3_table[['gene_2', 'gene_1', 'humap_score']].copy()
    backward.columns = ['gene_1', 'gene_2', 'humap_score']

    # Combine and get max scores
    combined = pd.concat([forward, backward], ignore_index=True)
    max_scores = combined.groupby(['gene_1', 'gene_2'])['humap_score'].max()
    
    return max_scores.to_dict()


def compute_features(humap3_dict, gene_pairs): 
    gene_1 = gene_pairs.loc['gene_1']
    gene_2 = gene_pairs.loc['gene_2']
    score = round(humap3_dict.get((gene_1, gene_2), -1), 3)

    feat_vec = [
        gene_1,
        gene_2,
        score,
    ]
    return(feat_vec)


def main():
    """
    gene_pairs_file = 'gene_pairs.tsv'
    humap_table_file = 'input/humap_table.tsv'
    """
    gene_pairs_file = sys.argv[1]
    humap3_table_file = sys.argv[2]
    n_proc = int(sys.argv[3])

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['gene_1', 'gene_2']
    humap3_table = pd.read_csv(humap3_table_file, sep='\t')

    humap3_dict = create_lookup_vectorized(humap3_table)
        
    df_cols = [
        'gene_1',
        'gene_2',
        'Hu.MAP 3.0 Score',
    ]

    gene_pairs_dfs = [gene_pairs.iloc[i,] for i in range(len(gene_pairs))]

    arg1 = [humap3_dict for i in range(len(gene_pairs_dfs))]
    arg2 = gene_pairs_dfs

    with Pool(processes=n_proc) as pool:
        df_rows = pool.starmap(compute_features, zip(arg1, arg2) )

    features_df = pd.DataFrame(df_rows)
    features_df.columns = df_cols

    print(features_df.to_csv(sep='\t', index=False))

if __name__ == '__main__':
    main()