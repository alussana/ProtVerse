#!/usr/bin/env python3

import sys
import pandas as pd

"""
gene_dict_file = 'input/gene_id_dict.tsv'
dependency_genes_file = 'input/dependency_genes.txt'
eprot_genes_file = 'input/eprot_genes.txt'
gtex_genes_file = 'input/gtex_genes.txt'
orthogroup_genes_file = 'input/orthogroup_genes.txt'
ptmdb_genes_file = 'input/ptmdb_genes.txt'
reactome_genes_file = 'input/reactome_sig_genes.txt'
"""

def main():
    gene_dict_file = sys.argv[1]
    dependency_genes_file = sys.argv[5]
    eprot_genes_file = sys.argv[4]
    gtex_genes_file = sys.argv[3]
    orthogroup_genes_file = sys.argv[6]
    ptmdb_genes_file = sys.argv[7]
    reactome_genes_file = sys.argv[2]
    ivkaphe_genes_file = sys.argv[8]
    humap_genes_file = sys.argv[9]

    with open(dependency_genes_file) as dependency_genes:
        dependency_ids = []
        for line in dependency_genes:
            dependency_ids.append(line.strip())
    dependency = [set(dependency_ids), 'genesymbol']

    with open(eprot_genes_file) as eprot_genes:
        eprot_ids = []
        for line in eprot_genes:
            eprot_ids.append(line.strip())
    eprot = [set(eprot_ids), 'ensembl_g']

    with open(gtex_genes_file) as gtex_genes:
        gtex_ids = []
        for line in gtex_genes:
            gtex_ids.append(line.strip())
    gtex = [set(gtex_ids), 'ensembl_g']

    with open(orthogroup_genes_file) as orthogroup_genes:
        orthogroup_ids = []
        for line in orthogroup_genes:
            orthogroup_ids.append(line.strip())
    #orthogroup = [set(orthogroup_ids), 'ncbi']
    orthogroup = [set(orthogroup_ids), 'genesymbol']

    with open(ptmdb_genes_file) as ptmdb_genes:
        ptmdb_ids = []
        for line in ptmdb_genes:
            ptmdb_ids.append(line.strip())
    ptmdb = [set(ptmdb_ids), 'uniprot']

    with open(reactome_genes_file) as reactome_genes:
        reactome_ids = []
        for line in reactome_genes:
            reactome_ids.append(line.strip())
    reactome = [set(reactome_ids), 'ensembl_g']

    with open(ivkaphe_genes_file) as ivkaphe_genes:
        ivkaphe_ids = []
        for line in ivkaphe_genes:
            ivkaphe_ids.append(line.strip())
    ivkaphe = [set(ivkaphe_ids), 'uniprot']

    with open(humap_genes_file) as humap_genes:
        humap_ids = []
        for line in humap_genes:
            humap_ids.append(line.strip())
    humap = [set(humap_ids), 'genesymbol']

    gene_dict = pd.read_csv(gene_dict_file, sep='\t', header=None, dtype=str)
    gene_dict.columns = ['ensembl_g','ncbi','uniprot','genesymbol', 'ensembl_p']

    genes = list(set(gene_dict.genesymbol.dropna()))
    data = {'reactome': reactome,
            'gtex': gtex,
            'eprot': eprot,
            'dependency': dependency,
            'orthogroup': orthogroup,
            'ptmdb': ptmdb,
            'ivkaphe': ivkaphe,
            'humap': humap}

    df = pd.DataFrame()
    for gene in genes:
        mapping = gene_dict.loc[gene_dict['genesymbol']==gene]
        df_row = []
        for dataset_name in list(data.keys()):
            dataset_genes = data[dataset_name][0]
            dataset_nomenclature = data[dataset_name][1]
            translations = set(mapping.loc[:, dataset_nomenclature].unique())
            intersection = list(translations.intersection(dataset_genes))
            if len(intersection) == 0:
                value = 'NA'
            else:
                value = intersection[0]
            df_row.append(value)
        new_row = pd.DataFrame([df_row], index=[gene], columns=list(data.keys()))
        df = pd.concat([df, new_row])

    print(df.to_csv(sep='\t', index=True, header=True, na_rep='NA'))

if __name__ == '__main__':
    main()