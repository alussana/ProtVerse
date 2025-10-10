#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def main():
    gene_pairs_file = sys.argv[1]
    orthogroup_table_file = sys.argv[2]
    orthogroup_pcs_file = sys.argv[3]

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['label', 'reactome', 'gene_1', 'gene_2']

    orthogroup_table = pd.read_csv(
        orthogroup_table_file, sep='\t', index_col=0, header=None
    )
    orthogroup_table.columns = [
        'NCBI_Gene_ID',
        'Gene_Index',
        'Associated_Orthogroups'
    ]
    orthogroup_pcs = pd.read_csv(orthogroup_pcs_file, sep='\t', index_col=None)

    pcs = []
    same_orthogroup = []
    orthogroup_data = []

    for i in range(len(gene_pairs)):
        try:
            gene_1 = gene_pairs.loc[i,]['gene_1']
            gene_2 = gene_pairs.loc[i,]['gene_2']
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
            max_pcs = max_pcs / 100
            max_pcs = round(max_pcs, 3)
            if max_pcs < 0:
                max_pcs = 0
            data_avail = 1
        except:
            max_pcs = -1
            is_same_orthogroup = -1
            data_avail = 0
        pcs.append(max_pcs)
        same_orthogroup.append(is_same_orthogroup)
        orthogroup_data.append(data_avail)

    gene_pairs['Orthogroup PCS'] = pcs
    gene_pairs['Orthogroup Identity'] = same_orthogroup
    gene_pairs['Orthogroup Data'] = orthogroup_data

    print(gene_pairs.to_csv(sep='\t', index=False))

if __name__ == '__main__':
    main()