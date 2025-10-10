#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from ccc.coef import ccc
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity
from information_theory import normalized_mutual_information
from multiprocessing import Pool


def return_missing_variables(n:int, default=-1):
    return [default for i in range(n)]


def p_to_activation(x, a=3):
    if x == -1:
        raise ValueError("f(x) is undefined at x = -1 (division by zero).")
    return (1 / ( x + 1 ) ** a) - x * (1 / ( x + 1 ) ** a)


def map_sign(x):
    return (x + 1) // 2


def compute_features(ptmdb, gene_pairs, index_gene):
    r2_opt_cover, r2_sign_opt_cover, r2_pval_opt_cover, s_opt_cover, s_sign_opt_cover, s_pval_opt_cover, cs_opt_cover, cs_sign_opt_cover, c_opt_cover = return_missing_variables(9)
    r2_opt_score, r2_sign_opt_score, r2_pval_opt_score, s_opt_score, s_sign_opt_score, s_pval_opt_score, cs_opt_score, cs_sign_opt_score, c_opt_score = return_missing_variables(9)
    max_cr, min_cr, tot_cr, max_cr_opt_cover, min_cr_opt_cover, tot_cr_opt_cover = return_missing_variables(6)
    nmi_detection, nmi_direction = return_missing_variables(2)
    max_cr_r2, min_cr_r2, tot_cr_r2, max_cr_s, min_cr_s, tot_cr_s, max_cr_cs, min_cr_cs, tot_cr_cs, max_cr_c, min_cr_c, tot_cr_c = return_missing_variables(12)
    data_available, overlap_available = return_missing_variables(2, 0)
    gene_1 = gene_pairs.loc['gene_1']
    gene_2 = gene_pairs.loc['gene_2']
    gene_1_indices = [i for i, gene in enumerate(index_gene) if gene_1==gene]
    gene_2_indices = [i for i, gene in enumerate(index_gene) if gene_2==gene]
    if len(gene_1_indices) > 0 and len(gene_2_indices) > 0:
        # normalized mutual information (direction of regulation)
        X = np.sign(ptmdb.iloc[gene_1_indices]).fillna(0).astype(int).transpose().values
        Y = np.sign(ptmdb.iloc[gene_2_indices]).fillna(0).astype(int).transpose().values
        nmi_direction = max(0, round(normalized_mutual_information(X, Y, normalization="avg"), 3))
        # normalized mutual information (detection)
        X = ptmdb.iloc[gene_1_indices].notna().astype('int8').transpose().values
        Y = ptmdb.iloc[gene_2_indices].notna().astype('int8').transpose().values
        nmi_detection = max(0, round(normalized_mutual_information(X, Y, normalization="avg"), 3))
        # total codetection ratios
        X_bitwise = np.bitwise_or.reduce(X, axis=1)
        Y_bitwise = np.bitwise_or.reduce(Y, axis=1)
        codetect_arr = np.logical_and(X_bitwise, Y_bitwise)
        max_cr = round(codetect_arr.sum() / X_bitwise.sum(), 3)
        min_cr = round(codetect_arr.sum() / Y_bitwise.sum(), 3)
        tot_cr = round(codetect_arr.sum() / len(codetect_arr), 3)
        for pos_1 in gene_1_indices:
            a = ptmdb.iloc[[pos_1]].transpose()
            phos_1 = a.columns[0]
            a = a.values.reshape(len(a))
            # if phosphosite a is all nan or all zeroes, skip it
            if all(np.isnan(a)) or (a == 0).all():
                continue
            for pos_2 in gene_2_indices:
                b = ptmdb.iloc[[pos_2]].transpose()
                phos_2 = b.columns[0]
                b = b.values.reshape(len(b))
                # if phosphosite b is all nan or all zeroes, skip it
                if all(np.isnan(b)) or (b == 0).all():
                    continue
                # create arrays where both phosphosites have no missing or infinite values
                #mask = np.logical_and(~np.isnan(a), ~np.isnan(b))
                mask = np.isfinite(a) & np.isfinite(b)
                a_filtered = a[mask]
                b_filtered = b[mask]
                # if no data overlap is found for the two phosphosites, skip the pair but indicate that the data is available
                if len(a_filtered) == 0 or len(b_filtered) == 0:
                    data_available = 1
                    max_cr_opt_cover = 0
                    min_cr_opt_cover = 0
                    tot_cr_opt_cover = 0
                    continue
                else:
                    # codetection ratios
                    den_min = min(sum(~np.isnan(a)), sum(~np.isnan(b)))
                    den_max = max(sum(~np.isnan(a)), sum(~np.isnan(b)))
                    codetection = len(a_filtered)
                    max_cr_opt_score_new = round(codetection / den_min, 3)
                    min_cr_opt_score_new = round(codetection / den_max, 3)
                    tot_cr_opt_score_new = round(codetection / len(a), 3)
                    # if the data overlap for the two phosphosites is greater than 3 points, compute features and indicate that overlap is available
                    if codetection > 3 and np.std(a_filtered) != 0 and np.std(b_filtered) != 0:
                        overlap_available = 1
                        # pearson correlation
                        r2_opt_score_new, r2_pval_opt_score_new = pearsonr(a_filtered, b_filtered)
                        r2_sign_opt_score_new = map_sign(np.sign(r2_opt_score_new))
                        r2_opt_score_new = round(r2_opt_score_new**2, 3)
                        r2_pval_opt_score_new = round(p_to_activation(r2_pval_opt_score_new), 3)
                        if r2_opt_score_new > r2_opt_score:
                            r2_opt_score = r2_opt_score_new
                            r2_sign_opt_score = r2_sign_opt_score_new
                            r2_pval_opt_score = r2_pval_opt_score_new
                            max_cr_r2 = max_cr_opt_score_new
                            min_cr_r2 = min_cr_opt_score_new
                            tot_cr_r2 = tot_cr_opt_score_new
                        # spearman correlation
                        s_opt_score_new, s_pval_opt_score_new = spearmanr(a_filtered, b_filtered)
                        s_sign_opt_score_new = map_sign(np.sign(s_opt_score_new))
                        s_opt_score_new = round(s_opt_score_new**2, 3)
                        s_pval_opt_score_new = round(p_to_activation(s_pval_opt_score_new), 3)
                        if s_opt_score_new > s_opt_score:
                            s_opt_score = s_opt_score_new
                            s_sign_opt_score = s_sign_opt_score_new
                            s_pval_opt_score = s_pval_opt_score_new
                            max_cr_s = max_cr_opt_score_new
                            min_cr_s = min_cr_opt_score_new
                            tot_cr_s = tot_cr_opt_score_new
                        # cosine similarity
                        #cos_dist = cosine(a_filtered, b_filtered)
                        #cs = 1 - cos_dist
                        cs = cosine_similarity(a_filtered.reshape(1, -1), b_filtered.reshape(1, -1))[0][0]
                        cs_sign_opt_score_new = map_sign(np.sign(cs))
                        cs_opt_score_new = round((cs) ** 2, 3)
                        if cs_opt_score_new > cs_opt_score:
                            cs_opt_score = cs_opt_score_new
                            cs_sign_opt_score = cs_sign_opt_score_new
                            max_cr_cs = max_cr_opt_score_new
                            min_cr_cs = min_cr_opt_score_new
                            tot_cr_cs = tot_cr_opt_score_new
                        # clutermatch correlation coefficient
                        c_opt_score_new = round(ccc(a, b), 3)
                        if c_opt_score_new > c_opt_score:
                            c_opt_score = c_opt_score_new
                            max_cr_c = max_cr_opt_score_new
                            min_cr_c = min_cr_opt_score_new
                            tot_cr_c = tot_cr_opt_score_new
                    else:
                        r2_opt_score_new, r2_sign_opt_score_new, r2_pval_opt_score_new, s_opt_score_new, s_sign_opt_score_new,s_pval_opt_score_new, cs_opt_score_new, cs_sign_opt_score_new, c_opt_score_new = return_missing_variables(9)
                    # overlap ratio
                    tot_cr_opt_cover_new = round(len(a_filtered) / len(a), 3)
                    # if better overlap has been discovered then update metrics
                    if tot_cr_opt_cover_new > tot_cr_opt_cover:
                        max_cr_opt_cover = max_cr_opt_score_new
                        min_cr_opt_cover = min_cr_opt_score_new
                        tot_cr_opt_cover = tot_cr_opt_cover_new
                        r2_opt_cover = r2_opt_score_new
                        r2_sign_opt_cover = r2_sign_opt_score_new
                        r2_pval_opt_cover = r2_pval_opt_score_new
                        s_opt_cover = s_opt_score_new
                        s_sign_opt_cover = s_sign_opt_score_new
                        s_pval_opt_cover = s_pval_opt_score_new
                        cs_opt_cover = cs_opt_score_new
                        cs_sign_opt_cover = cs_sign_opt_score_new
                        c_opt_cover = c_opt_score_new


    feat_vec = [
        gene_1,
        gene_2,
        nmi_direction,
        nmi_detection,
        max_cr,
        min_cr,
        tot_cr,
        c_opt_score,
        cs_opt_score,
        cs_sign_opt_score,
        r2_opt_score,
        r2_sign_opt_score,
        r2_pval_opt_score,
        s_opt_score,
        s_sign_opt_score,
        s_pval_opt_score,
        c_opt_cover,
        cs_opt_cover,
        cs_sign_opt_cover,
        r2_opt_cover,
        r2_sign_opt_cover,
        r2_pval_opt_cover,
        s_opt_cover,
        s_sign_opt_cover,
        s_pval_opt_cover,
        max_cr_r2,
        min_cr_r2,
        tot_cr_r2,
        max_cr_s,
        min_cr_s,
        tot_cr_s,
        max_cr_cs,
        min_cr_cs,
        tot_cr_cs,
        max_cr_c,
        min_cr_c,
        tot_cr_c,
        max_cr_opt_cover,
        min_cr_opt_cover,
        tot_cr_opt_cover,
        #data_available,
        #overlap_available,
    ]
    return(feat_vec)

def main():
    """
    gene_pairs_file = 'gene_pairs.tsv'
    ptmdb_table_file = 'all_logFC.tsv'
    """
    gene_pairs_file = sys.argv[1]
    ptmdb_table_file = sys.argv[2]
    n_proc = int(sys.argv[3])

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['gene_1', 'gene_2']

    ptmdb = pd.read_csv(ptmdb_table_file, sep='\t', index_col=0, header=None)

    index_gene = [item.split('_')[0] for item in list(ptmdb.index)]

    df_cols = [
        'gene_1',
        'gene_2',
        'PTMDB NMI (direction)',
        'PTMDB NMI (detection)',
        'PTMDB Max CR (global)',
        'PTMDB Min CR (global)',
        'PTMDB Tot CR (global)',
        'PTMDB CCC (Max)',
        r'PTMDB $\cos^2\theta$ (Max)',
        r'PTMDB $\mathrm{sign}\,(\cos\theta)$ (Max)',
        r'PTMDB $r^2$ (Max)',
        r'PTMDB $\mathrm{sign}\,(r)$ (Max)',
        r'PTMDB $r^2$ Sig. (Max)',
        r'PTMDB $r_s^2$ (Max)',
        r'PTMDB $\mathrm{sign}\,(r_s)$ (Max)',
        r'PTMDB $r_s^2$ Sig. (Max)',
        'PTMDB CCC (Cover)',
        r'PTMDB $\cos^2\theta$ (Cover)',
        r'PTMDB $\mathrm{sign}\,(\cos\theta)$ (Cover)',
        r'PTMDB $r^2$ (Cover)',
        r'PTMDB $\mathrm{sign}\,(r)$ (Cover)',
        r'PTMDB $r^2$ Sig. (Cover)',
        r'PTMDB $r_s^2$ (Cover)',
        r'PTMDB $\mathrm{sign}\,(r_s)$ (Cover)',
        r'PTMDB $r_s^2$ Sig. (Cover)',
        r'PTMDB Max CR ($r^2$)',
        r'PTMDB Min CR ($r^2$)',
        r'PTMDB Tot CR ($r^2$)',
        r'PTMDB Max CR ($r_s^2$)',
        r'PTMDB Min CR ($r_s^2$)',
        r'PTMDB Tot CR ($r_s^2$)',
        r'PTMDB Max CR ($\cos\theta$)',
        r'PTMDB Min CR ($\cos\theta$)',
        r'PTMDB Tot CR ($\cos\theta$)',
        'PTMDB Max CR (CCC)',
        'PTMDB Min CR (CCC)',
        'PTMDB Tot CR (CCC)',
        'PTMDB Max CR (site)',
        'PTMDB Min CR (site)',
        'PTMDB Tot CR (site)',
        #'PTMDB Data',
        #'PTMDB Overlap',
    ]

    gene_pairs_dfs = [gene_pairs.iloc[i,] for i in range(len(gene_pairs))]

    arg1 = [ptmdb for i in range(len(gene_pairs_dfs))]
    arg2 = gene_pairs_dfs
    arg3 = [index_gene for i in range(len(gene_pairs_dfs))]

    with Pool(processes=n_proc) as pool:
        df_rows = pool.starmap(compute_features, zip(arg1, arg2, arg3) )

    features_df = pd.DataFrame(df_rows)
    features_df.columns = df_cols

    print(features_df.to_csv(sep='\t', index=False))

if __name__ == '__main__':
    main()
