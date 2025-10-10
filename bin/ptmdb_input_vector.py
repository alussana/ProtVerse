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


def return_missing_variables(n:int, default=-1):
    return [default for i in range(n)]


def p_to_activation(x, a=3):
    if x == -1:
        raise ValueError("f(x) is undefined at x = -1 (division by zero).")
    return (1 / ( x + 1 ) ** a) - x * (1 / ( x + 1 ) ** a)


def map_sign(x):
    return (x + 1) // 2


def main():
    """
    gene_pairs_file = 'input/gene_pairs.tsv'
    ptmdb_table_file = 'all_logFC.tsv'
    """
    gene_pairs_file = sys.argv[1]
    ptmdb_table_file = sys.argv[2]

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['label', 'reactome', 'gene_1', 'gene_2']

    ptmdb = pd.read_csv(ptmdb_table_file, sep='\t', index_col=0, header=None)

    pearson_opt_score = []
    pearson_sign_opt_score = []
    pearson_pval_opt_score = []
    spearman_opt_score = []
    spearman_sign_opt_score = []
    spearman_pval_opt_score = []
    cos_sim_opt_score = []
    cos_sim_sign_opt_score = []
    ccc_coef_opt_score = []
    
    norm_mutual_info_categorical = []
    norm_mutual_info_boolean = []

    max_codetection_ratio_pearson = []
    min_codetection_ratio_pearson = []
    tot_codetection_ratio_pearson = []
    max_codetection_ratio_spearman = []
    min_codetection_ratio_spearman = []
    tot_codetection_ratio_spearman = []
    max_codetection_ratio_cos_sim = []
    min_codetection_ratio_cos_sim = []
    tot_codetection_ratio_cos_sim = []
    max_codetection_ratio_ccc = []
    min_codetection_ratio_ccc = []
    tot_codetection_ratio_ccc = []

    pearson_opt_cover = []
    pearson_sign_opt_cover = []
    pearson_pval_opt_cover = []
    spearman_opt_cover = []
    spearman_sign_opt_cover = []
    spearman_pval_opt_cover = []
    cos_sim_opt_cover = []
    cos_sim_sign_opt_cover = []
    ccc_coef_opt_cover = []

    max_cr_global = []
    min_cr_global = []
    tot_cr_global = []

    max_site_cr = []
    min_site_cr = []
    tot_site_cr = []

    data_flag = []
    overlap_flag = []

    index_gene = [item.split('_')[0] for item in list(ptmdb.index)]

    for i in range(len(gene_pairs)):
        r2_opt_cover, r2_sign_opt_cover, r2_pval_opt_cover, s_opt_cover, s_sign_opt_cover, s_pval_opt_cover, cs_opt_cover, cs_sign_opt_cover, c_opt_cover = return_missing_variables(9)
        r2_opt_score, r2_sign_opt_score, r2_pval_opt_score, s_opt_score, s_sign_opt_score, s_pval_opt_score, cs_opt_score, cs_sign_opt_score, c_opt_score = return_missing_variables(9)
        max_cr, min_cr, tot_cr, max_cr_opt_cover, min_cr_opt_cover, tot_cr_opt_cover = return_missing_variables(6)
        nmi_detection, nmi_direction = return_missing_variables(2)
        max_cr_r2, min_cr_r2, tot_cr_r2, max_cr_s, min_cr_s, tot_cr_s, max_cr_cs, min_cr_cs, tot_cr_cs, max_cr_c, min_cr_c, tot_cr_c = return_missing_variables(12)
        data_available, overlap_available = return_missing_variables(2, 0)
        gene_1 = gene_pairs.loc[i,]['gene_1']
        gene_2 = gene_pairs.loc[i,]['gene_2']
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

        pearson_opt_score.append(r2_opt_score)
        pearson_sign_opt_score.append(r2_sign_opt_score)
        pearson_pval_opt_score.append(r2_pval_opt_score)
        spearman_opt_score.append(s_opt_score)
        spearman_sign_opt_score.append(s_sign_opt_score)
        spearman_pval_opt_score.append(s_pval_opt_score)
        cos_sim_opt_score.append(cs_opt_score)
        cos_sim_sign_opt_score.append(cs_sign_opt_score)
        ccc_coef_opt_score.append(c_opt_score)

        norm_mutual_info_categorical.append(nmi_direction)
        norm_mutual_info_boolean.append(nmi_detection)

        max_codetection_ratio_pearson.append(max_cr_r2)
        min_codetection_ratio_pearson.append(min_cr_r2)
        tot_codetection_ratio_pearson.append(tot_cr_r2)
        max_codetection_ratio_spearman.append(max_cr_s)
        min_codetection_ratio_spearman.append(min_cr_s)
        tot_codetection_ratio_spearman.append(tot_cr_s)
        max_codetection_ratio_cos_sim.append(max_cr_cs)
        min_codetection_ratio_cos_sim.append(min_cr_cs)
        tot_codetection_ratio_cos_sim.append(tot_cr_cs)
        max_codetection_ratio_ccc.append(max_cr_c)
        min_codetection_ratio_ccc.append(min_cr_c)
        tot_codetection_ratio_ccc.append(tot_cr_c)

        pearson_opt_cover.append(r2_opt_cover)
        pearson_sign_opt_cover.append(r2_sign_opt_cover)
        pearson_pval_opt_cover.append(r2_pval_opt_cover)
        spearman_opt_cover.append(s_opt_cover)
        spearman_sign_opt_cover.append(s_sign_opt_cover)
        spearman_pval_opt_cover.append(s_pval_opt_cover)
        cos_sim_opt_cover.append(cs_opt_cover)
        cos_sim_sign_opt_cover.append(cs_sign_opt_cover)
        ccc_coef_opt_cover.append(c_opt_cover)

        max_cr_global.append(max_cr)
        min_cr_global.append(min_cr)
        tot_cr_global.append(tot_cr)

        max_site_cr.append(max_cr_opt_cover)
        min_site_cr.append(min_cr_opt_cover)
        tot_site_cr.append(tot_cr_opt_cover)

        data_flag.append(data_available)
        overlap_flag.append(overlap_available)


    gene_pairs['PTMDB NMI (direction)'] = norm_mutual_info_categorical
    gene_pairs['PTMDB NMI (detection)'] = norm_mutual_info_boolean

    gene_pairs['PTMDB Max CR (global)'] = max_cr_global
    gene_pairs['PTMDB Min CR (global)'] = min_cr_global
    gene_pairs['PTMDB Tot CR (global)'] = tot_cr_global

    gene_pairs['PTMDB CCC (Max)'] = ccc_coef_opt_score
    gene_pairs[r'PTMDB $\cos^2\theta$ (Max)'] = cos_sim_opt_score
    gene_pairs[r'PTMDB $\mathrm{sign}\,(\cos\theta)$ (Max)'] = cos_sim_sign_opt_score
    gene_pairs[r'PTMDB $r^2$ (Max)'] = pearson_opt_score
    gene_pairs[r'PTMDB $\mathrm{sign}\,(r)$ (Max)'] = pearson_sign_opt_score
    gene_pairs[r'PTMDB $r^2$ Sig. (Max)'] = pearson_pval_opt_score
    gene_pairs[r'PTMDB $r_s^2$ (Max)'] = spearman_opt_score
    gene_pairs[r'PTMDB $\mathrm{sign}\,(r_s)$ (Max)'] = spearman_sign_opt_score
    gene_pairs[r'PTMDB $r_s^2$ Sig. (Max)'] = spearman_pval_opt_score

    gene_pairs[r'PTMDB CCC (Cover)'] = ccc_coef_opt_cover
    gene_pairs[r'PTMDB $\cos^2\theta$ (Cover)'] = cos_sim_opt_cover
    gene_pairs[r'PTMDB $\mathrm{sign}\,(\cos\theta)$ (Cover)'] = cos_sim_sign_opt_cover
    gene_pairs[r'PTMDB $r^2$ (Cover)'] = pearson_opt_cover
    gene_pairs[r'PTMDB $\mathrm{sign}\,(r)$ (Cover)'] = pearson_sign_opt_cover
    gene_pairs[r'PTMDB $r^2$ Sig. (Cover)'] = pearson_pval_opt_cover
    gene_pairs[r'PTMDB $r_s^2$ (Cover)'] = spearman_opt_cover
    gene_pairs[r'PTMDB $\mathrm{sign}\,(r_s)$ (Cover)'] = spearman_sign_opt_cover
    gene_pairs[r'PTMDB $r_s^2$ Sig. (Cover)'] = spearman_pval_opt_cover

    gene_pairs[r'PTMDB Max CR ($r^2$)'] = max_codetection_ratio_pearson
    gene_pairs[r'PTMDB Min CR ($r^2$)'] = min_codetection_ratio_pearson
    gene_pairs[r'PTMDB Tot CR ($r^2$)'] = tot_codetection_ratio_pearson
    gene_pairs[r'PTMDB Max CR ($r_s^2$)'] = max_codetection_ratio_spearman
    gene_pairs[r'PTMDB Min CR ($r_s^2$)'] = min_codetection_ratio_spearman
    gene_pairs[r'PTMDB Tot CR ($r_s^2$)'] = tot_codetection_ratio_spearman
    gene_pairs[r'PTMDB Max CR ($\cos\theta$)'] = max_codetection_ratio_cos_sim
    gene_pairs[r'PTMDB Min CR ($\cos\theta$)'] = min_codetection_ratio_cos_sim
    gene_pairs[r'PTMDB Tot CR ($\cos\theta$)'] = tot_codetection_ratio_cos_sim
    gene_pairs['PTMDB Max CR (CCC)'] = max_codetection_ratio_ccc
    gene_pairs['PTMDB Min CR (CCC)'] = min_codetection_ratio_ccc
    gene_pairs['PTMDB Tot CR (CCC)'] = tot_codetection_ratio_ccc

    gene_pairs['PTMDB Max CR (site)'] = max_site_cr
    gene_pairs['PTMDB Min CR (site)'] = min_site_cr
    gene_pairs['PTMDB Tot CR (site)'] = tot_site_cr

    gene_pairs['PTMDB Data'] = data_flag
    gene_pairs['PTMDB Overlap'] = overlap_flag


    print(gene_pairs.to_csv(sep='\t', index=False))


if __name__ == '__main__':
    main()