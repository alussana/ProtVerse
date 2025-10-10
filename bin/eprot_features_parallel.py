#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from ccc.coef import ccc
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity
from information_theory import information_correlation_knn
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
        p, p_sign, p_pval, s, s_sign, s_pval, cs, cs_sign, nmi, c, max_cr, min_cr, tot_cr = return_missing_variables(13)
        data_available, overlap_available = return_missing_variables(2, 0)
    else:
        a = gene_1_expr[gene_1].values
        b = gene_2_expr[gene_2].values
        # if either gene is found but has only missing values, return missing features
        if all(np.isnan(a)) or all(np.isnan(b)) or (a == 0).all() or (b == 0).all():
            p, p_sign, p_pval, s, s_sign, s_pval, cs, cs_sign, nmi, c, max_cr, min_cr, tot_cr = return_missing_variables(13)
            data_available, overlap_available = return_missing_variables(2, 0)
        else:
            # create arrays where both genes have no missing values
            mask = np.logical_and(~np.isnan(a), ~np.isnan(b))
            a_filtered = a[mask]
            b_filtered = b[mask]
            # if no data overlap is found for the two genes, return missing features, but indicate that the data is available
            if len(a_filtered) == 0 or len(b_filtered) == 0:
                p, p_sign, p_pval, s, s_sign, s_pval, cs, cs_sign, nmi, c, max_cr, min_cr = return_missing_variables(12)
                overlap_available, tot_cr = return_missing_variables(2, 0)
            else:
                # codetection ratios
                den_min = min(sum(~np.isnan(a)), sum(~np.isnan(b)))
                den_max = max(sum(~np.isnan(a)), sum(~np.isnan(b)))
                codetection = len(a_filtered)
                max_cr = round(codetection / den_min, 3)
                min_cr = round(codetection / den_max, 3)
                tot_cr = round(codetection / len(a), 3)
                # if the data overlap for the two genes is greater than 3 points, compute features
                if codetection > 3 and np.std(a_filtered) != 0 and np.std(b_filtered) != 0:
                    # pearson correlation
                    p, p_pval = pearsonr(a_filtered, b_filtered)
                    p_sign = map_sign(np.sign(p))
                    p = round(p**2, 3)
                    p_pval = round(p_to_activation(p_pval), 3)
                    # spearman correlation
                    s, s_pval = spearmanr(a_filtered, b_filtered)
                    s_sign = map_sign(np.sign(s))
                    s = round(s**2, 3)
                    s_pval = round(p_to_activation(s_pval), 3)
                    # cosine similarity
                    #cos_dist = cosine(a_filtered, b_filtered)
                    #cs = 1 - cos_dist
                    cs = cosine_similarity(a_filtered.reshape(1, -1), b_filtered.reshape(1, -1))[0][0]
                    cs_sign = map_sign(np.sign(cs))
                    cs = round((cs) ** 2, 3)
                    # normalized mutual information
                    X = a_filtered.reshape([a_filtered.shape[0],1])
                    Y = b_filtered.reshape([b_filtered.shape[0],1])
                    nmi = max(0, round(information_correlation_knn(X, Y, add_noise=1e-10), 3))
                    # clutermatch correlation coefficient
                    c = round(ccc(gene_1_expr[gene_1], gene_2_expr[gene_2]), 3)
                else:
                   p, p_sign, p_pval, s, s_sign, s_pval, cs, cs_sign, nmi, c = return_missing_variables(10)
                   overlap_available = 0
    feat_vec = [
        gene_1,
        gene_2,
        nmi,
        c,
        cs,
        cs_sign,
        p,
        p_sign,
        p_pval,
        s,
        s_sign,
        s_pval,
        max_cr,
        min_cr,
        tot_cr,
        #data_available,
        #overlap_available,
    ]
    return(feat_vec)


def main():
    """
    gene_pairs_file = 'gene_pairs.tsv'
    eprot_table_file = 'input/eprot_table.tsv'
    n_proc = 1
    """
    gene_pairs_file = sys.argv[1]
    eprot_table_file = sys.argv[2]
    n_proc = int(sys.argv[3])

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['gene_1', 'gene_2']

    eprot = pd.read_csv(eprot_table_file, sep='\t', index_col=0, header=None)

    df_cols = [
        'gene_1',
        'gene_2',
        'PRIDE (tissues) Linfoot MI',
        'PRIDE (tissues) CCC',
        r'PRIDE (tissues) $\cos^2\theta$',
        r'PRIDE (tissues) $\mathrm{sign}\,(\cos\theta)$',
        r'PRIDE (tissues) $r^2$',
        r'PRIDE (tissues) $\mathrm{sign}\,(r)$',
        r'PRIDE (tissues) $r^2$ Sig.',
        r'PRIDE (tissues) $r_s^2$',
        r'PRIDE (tissues) $\mathrm{sign}\,(r_s)$',
        r'PRIDE (tissues) $r_s^2$ Sig.',
        'PRIDE (tissues) Max CR',
        'PRIDE (tissues) Min CR',
        'PRIDE (tissues) Tot CR',
        #'PRIDE (tissues) Data',
        #'PRIDE (tissues) Overlap'
    ]

    gene_pairs_dfs = [gene_pairs.iloc[i,] for i in range(len(gene_pairs))]

    arg1 = [eprot for i in range(len(gene_pairs_dfs))]
    arg2 = gene_pairs_dfs

    with Pool(processes=n_proc) as pool:
        df_rows = pool.starmap(compute_features, zip(arg1, arg2) )

    features_df = pd.DataFrame(df_rows)
    features_df.columns = df_cols

    print(features_df.to_csv(sep='\t', index=False))

if __name__ == '__main__':
    main()