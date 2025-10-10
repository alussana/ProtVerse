#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from binary_similarity import compute_binary_similarity


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
    gene_pairs_file = 'gene_pairs.tsv'
    ubiquitination_table_file = 'input/ubiquitination_table.tsv'
    """
    gene_pairs_file = sys.argv[1]
    ubiquitination_table_file = sys.argv[2]

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['label', 'reactome', 'gene_1', 'gene_2']

    ubiquitination = pd.read_csv(ubiquitination_table_file, sep='\t', index_col=0, header=None)

    cos_sim = []
    cos_sign = []
    jaccard = []
    hamming = []
    norm_mutual_info = []
    data_flag = []

    for i in range(len(gene_pairs)):
        data_available = 1
        gene_1 = gene_pairs.loc[i,]['gene_1']
        gene_2 = gene_pairs.loc[i,]['gene_2']
        gene_1_expr = ubiquitination.loc[ubiquitination.index == gene_1, ].transpose()
        gene_2_expr = ubiquitination.loc[ubiquitination.index == gene_2, ].transpose()
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
        cos_sim.append(cs**2)
        cos_sign.append(cs_sign)
        jaccard.append(j)
        hamming.append(h)
        norm_mutual_info.append(nmi)
        data_flag.append(data_available)


    gene_pairs['Ubiquitination NMI'] = norm_mutual_info
    gene_pairs[r'Ubiquitination $\cos^2\theta$'] = cos_sim
    gene_pairs[r'Ubiquitination $\mathrm{sign}\,(\cos\theta)$'] = cos_sign
    gene_pairs['Ubiquitination Jaccard'] = jaccard
    gene_pairs['Ubiquitination Hamming'] = hamming
    gene_pairs['Ubiquitination Data'] = data_flag
            
    

    print(gene_pairs.to_csv(sep='\t', index=False))


if __name__ == '__main__':
    main()