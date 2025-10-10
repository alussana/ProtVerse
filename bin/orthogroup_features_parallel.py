#!/usr/bin/env python3

import sys
import pandas as pd
from multiprocessing import Pool

def compute_features(orthogroup_pcs, gene_pairs, orthogroup_table):
    try: 
        gene_1 = gene_pairs.loc['gene_1']
        gene_2 = gene_pairs.loc['gene_2']
        orthogroups_gene_1 = orthogroup_table.loc[
            orthogroup_table.index == gene_1,
            'Associated_Orthogroups'].iloc[0].split(',')
        orthogroups_gene_1 = [int(o) for o in orthogroups_gene_1]
        orthogroups_gene_2 = orthogroup_table.loc[
            orthogroup_table.index == gene_2,
            'Associated_Orthogroups'].iloc[0].split(',')
        orthogroups_gene_2 = [int(o) for o in orthogroups_gene_2]
        max_pcs = 0
        is_same_orthogroup = 0
        for x in range(len(orthogroups_gene_1)):
            for y in range(len(orthogroups_gene_2)):
                if orthogroups_gene_1[x] == orthogroups_gene_2[y]:
                    is_same_orthogroup = 1
                pcs_table = orthogroup_pcs.loc[
                    (orthogroup_pcs.Orthogroup1 == orthogroups_gene_1[x]) |
                    (orthogroup_pcs.Orthogroup1 == orthogroups_gene_2[y])
                ]
                pcs_table = pcs_table.loc[
                    (pcs_table.Orthogroup2 == orthogroups_gene_1[x]) |
                    (pcs_table.Orthogroup2 == orthogroups_gene_2[y])
                ]
                if len(pcs_table) > 0:
                    if float(pcs_table.PCS.values[0]) > max_pcs:
                        max_pcs = pcs_table.PCS.values[0]
        pcs = max_pcs
        same_orthogroup = is_same_orthogroup
    except:
        pcs = 0
        same_orthogroup = 0

    feat_vec = [
        gene_1,
        gene_2,
        pcs,
        same_orthogroup
    ]
    return(feat_vec)

def main():
    gene_pairs_file = sys.argv[1]
    orthogroup_table_file = sys.argv[2]
    orthogroup_pcs_file = sys.argv[3]
    n_proc = int(sys.argv[4])

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['gene_1', 'gene_2']

    orthogroup_table = pd.read_csv(
        orthogroup_table_file, sep='\t', index_col=0, header=None
    )
    orthogroup_table.columns = [
        'NCBI_Gene_ID',
        'Gene_Index',
        'Associated_Orthogroups'
    ]
    orthogroup_pcs = pd.read_csv(orthogroup_pcs_file, sep='\t', index_col=None)

    df_cols = [
        'gene_1',
        'gene_2',
        'Orthogroup PCS',
        'Orthogroup Identity'
    ]

    gene_pairs_dfs = [gene_pairs.iloc[i,] for i in range(len(gene_pairs))]

    arg1 = [orthogroup_pcs for i in range(len(gene_pairs_dfs))]
    arg2 = gene_pairs_dfs
    arg3 = [orthogroup_table for i in range(len(gene_pairs_dfs))]

    with Pool(processes=n_proc) as pool:
        df_rows = pool.starmap(compute_features, zip(arg1, arg2, arg3) )

    features_df = pd.DataFrame(df_rows)
    features_df.columns = df_cols

    print(features_df.to_csv(sep='\t', index=False))

if __name__ == '__main__':
    main()