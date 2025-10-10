#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from binary_similarity import compute_binary_similarity
from multiprocessing import Pool


def return_missing_variables(n:int, default=-1):
    return [default for i in range(n)]


def p_to_activation(x, a=3):
    if x == -1:
        raise ValueError("f(x) is undefined at x = -1 (division by zero).")
    return (1 / ( x + 1 ) ** a) - x * (1 / ( x + 1 ) ** a)


def map_sign(x):
    return (x + 1) // 2


def compute_features(data_df, gene_pair_df):
    data_available = 1
    overlap_available = 1
    gene_1 = gene_pair_df['gene_1']
    gene_2 = gene_pair_df['gene_2']
    gene_1_expr = data_df.loc[data_df.index == gene_1, ].transpose()
    gene_2_expr = data_df.loc[data_df.index == gene_2, ].transpose()
    # if either gene is not found, return missing features
    if len(gene_1_expr.columns) == 0 or len(gene_2_expr.columns) == 0:
        nmi, cs, cs_sign, j, h = return_missing_variables(5)
        data_available = 0
    else:
        a = gene_1_expr[gene_1].values
        b = gene_2_expr[gene_2].values
        # if either gene is found but has only missing values, return missing features
        if all(np.isnan(a)) or all(np.isnan(b)) or (a == 0).all() or (b == 0).all():
            nmi, cs, cs_sign, j, h = return_missing_variables(5)
            data_available = 0
        else:
            # compute similarity
            similarity = compute_binary_similarity(a, b)
            cs = round(similarity['cosine'], 3)
            cs_sign = map_sign(np.sign(cs))
            j = round(similarity['jaccard'], 3)
            h = round(similarity['hamming_similarity'], 3)
            nmi = round(similarity['NMI_avg'], 3)
    feat_vec = [
        gene_1,
        gene_2,
        nmi,
        cs**2,
        cs_sign,
        j,
        h,
        #data_available,
    ]
    return(feat_vec)


def main():
    """
    gene_pairs_file = 'gene_pairs.tsv'
    ubiquitination_table_file = 'input/gene_tpm.tsv'
    n_proc = 10
    """
    gene_pairs_file = sys.argv[1]
    ubiquitination_table_file = sys.argv[2]
    n_proc = int(sys.argv[3])

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['gene_1', 'gene_2']

    ubiquitination = pd.read_csv(ubiquitination_table_file, sep='\t', index_col=0, header=None)

    df_cols = [
        'gene_1',
        'gene_2',
        'Ubiquitination NMI',
        r'Ubiquitination $\cos^2\theta$',
        r'Ubiquitination $\mathrm{sign}\,(\cos\theta)$',
        'Ubiquitination Jaccard',
        'Ubiquitination Hamming',
    ]

    gene_pairs_dfs = [gene_pairs.iloc[i,] for i in range(len(gene_pairs))]

    arg1 = [ubiquitination for i in range(len(gene_pairs_dfs))]
    arg2 = gene_pairs_dfs

    with Pool(processes=n_proc) as pool:
        df_rows = pool.starmap(compute_features, zip(arg1, arg2) )

    features_df = pd.DataFrame(df_rows)
    features_df.columns = df_cols

    print(features_df.to_csv(sep='\t', index=False))


if __name__ == '__main__':
    main()