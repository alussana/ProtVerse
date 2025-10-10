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


def main():
    """
    gene_pairs_file = 'input/gene_pairs.tsv'
    humap3_table_file = 'input/humap3_table.tsv'
    n_proc = 8
    """
    gene_pairs_file = sys.argv[1]
    humap3_table_file = sys.argv[2]

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['label', 'reactome', 'gene_1', 'gene_2']
    humap3_table = pd.read_csv(humap3_table_file, sep='\t')

    humap3_dict = create_lookup_vectorized(humap3_table)

    humap3_scores = []
    for i in range(len(gene_pairs)):
        gene_1 = gene_pairs.loc[i, 'gene_1']
        gene_2 = gene_pairs.loc[i, 'gene_2']
        humap3_scores.append(round(humap3_dict.get((gene_1, gene_2), -1), 3))
    
    gene_pairs['Hu.MAP 3.0 Score'] = humap3_scores

    print(gene_pairs.to_csv(sep='\t', index=False))

if __name__ == '__main__':
    main()