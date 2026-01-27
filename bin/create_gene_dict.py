#!/usr/bin/env python3

import sys
import pandas as pd
from multiprocessing import Pool
from collections import ChainMap

def create_dict_col(dataset_name, ref_genes, idmapping, nomenclatures, gene_lists):
    dataset_col = []
    for gene in ref_genes:
        id_type = nomenclatures[dataset_name]
        # look for identical correspondence of gene name
        if id_type == 'Gene_Name' and gene in gene_lists[dataset_name]:
            match = gene
        else:
            # get uniprot ids for the reference gene name
            uniprot_ids = list(
                idmapping.loc[idmapping['value']==gene]['uniprot']
            )
            if id_type == 0:
                # look for a match with gene names in the considered dataset
                matches = [
                    name for name in gene_lists[dataset_name] if name in uniprot_ids
                ]
            else:
                # get ids in the appropriate nomenclature for the retrieved
                # uniprot ids
                values = idmapping.loc[
                    (pd.Series([a in uniprot_ids for a in idmapping['uniprot']])) &
                    (idmapping['nomenclature'] == id_type),
                    'value'
                ].values
                # look for a match with gene names in the considered dataset
                matches = [
                    name for name in gene_lists[dataset_name] if name in values
                ]
                if len(matches) == 0 and id_type == "Gene_Name":
                    # look for "Gene_Synonym"s
                    values = idmapping.loc[
                        (pd.Series([a in uniprot_ids for a in idmapping['uniprot']])) &
                        (idmapping['nomenclature'] == "Gene_synonym"),
                        'value'
                    ].values
                    # look for a match with gene names in the considered dataset
                    matches = [
                        name for name in gene_lists[dataset_name] if name in values
                    ]
            # define translations as 'NA' if no match
            # otherwise, take the first match of the list
            # (only 1 match expected in the vast majority of cases)
            if len(matches) == 0:
                match = 'NA'
            else:
                match = matches[0]
        dataset_col.append(match)
    dataset_col = {dataset_name: dataset_col}
    return(dataset_col)

def main():
    """
    ref_genes_file = 'input/ref_genes.txt'
    id_mapping_file = 'input/idmapping.tsv'
    signalling_genes_file = 'input/signalling_genes.txt'
    gtex_genes_file = 'input/gtex_genes.txt'
    eprot_genes_file = 'input/eprot_genes.txt'
    orthogroup_genes_file = 'input/orthogroups_genes.txt'
    depmap_genes_file = 'input/depmap_genes.txt'
    ubiquitination_genes_file = 'input/ubiquitination_genes.txt'
    humap3_genes_file = 'input/humap3_genes.txt'
    ptmdb_genes_file = 'input/ptmdb_genes.txt'
    proteomehd_genes_file = 'input/proteomehd_genes.txt'
    mitchell2023_genes_file = 'input/mitchell2023_genes.txt'
    lopit2025_genes_file = 'input/lopit2025_genes.txt'
    metabolism_gene_file = 'input/metabolism_genes.txt'
    reactome_genes_file = 'input/reactome_genes.txt'
    n_proc = 2
    """
    ref_genes_file = sys.argv[1]
    id_mapping_file = sys.argv[2]
    signalling_genes_file = sys.argv[3]
    gtex_genes_file = sys.argv[4]
    eprot_genes_file = sys.argv[5]
    orthogroup_genes_file = sys.argv[6]
    depmap_genes_file = sys.argv[7]
    ubiquitination_genes_file = sys.argv[8]
    humap3_genes_file = sys.argv[9]
    ptmdb_genes_file = sys.argv[10]
    proteomehd_genes_file = sys.argv[11]
    mitchell2023_genes_file = sys.argv[12]
    lopit2025_genes_file = sys.argv[13]
    metabolism_genes_file = sys.argv[14]
    reactome_genes_file = sys.argv[15]
    n_proc = int(sys.argv[16])

    ref_genes = []
    with open(ref_genes_file) as genes_fh:
        for gene in genes_fh:
            ref_genes.append(gene.strip())

    signalling_genes = []
    with open(signalling_genes_file) as genes_fh:
        for gene in genes_fh:
            signalling_genes.append(gene.strip())

    gtex_genes = []
    with open(gtex_genes_file) as genes_fh:
        for gene in genes_fh:
            gtex_genes.append(gene.strip())

    eprot_genes = []
    with open(eprot_genes_file) as genes_fh:
        for gene in genes_fh:
            eprot_genes.append(gene.strip())

    orthogroup_genes = []
    with open(orthogroup_genes_file) as genes_fh:
        for gene in genes_fh:
            orthogroup_genes.append(gene.strip())

    depmap_genes = []
    with open(depmap_genes_file) as genes_fh:
        for gene in genes_fh:
            depmap_genes.append(gene.strip())

    ubiquitination_genes = []
    with open(ubiquitination_genes_file) as genes_fh:
        for gene in genes_fh:
            ubiquitination_genes.append(gene.strip())

    humap3_genes = []
    with open(humap3_genes_file) as genes_fh:
        for gene in genes_fh:
            humap3_genes.append(gene.strip())

    ptmdb_genes = []
    with open(ptmdb_genes_file) as genes_fh:
        for gene in genes_fh:
            ptmdb_genes.append(gene.strip())
            
    proteomehd_genes = []
    with open(proteomehd_genes_file) as genes_fh:
        for gene in genes_fh:
            proteomehd_genes.append(gene.strip())
            
    mitchell2023_genes = []
    with open(mitchell2023_genes_file) as genes_fh:
        for gene in genes_fh:
            mitchell2023_genes.append(gene.strip())

    lopit2025_genes = []
    with open(lopit2025_genes_file) as genes_fh:
        for gene in genes_fh:
            lopit2025_genes.append(gene.strip())    

    metabolism_genes = []
    with open(metabolism_genes_file) as genes_fh:
        for gene in genes_fh:
            metabolism_genes.append(gene.strip())

    reactome_genes = []
    with open(reactome_genes_file) as genes_fh:
        for gene in genes_fh:
            reactome_genes.append(gene.strip())

    idmapping = pd.read_csv(id_mapping_file, sep='\t', header=None, index_col=None)
    idmapping.columns = ['uniprot', 'nomenclature', 'value']
    idmapping.loc[idmapping['nomenclature'] == 'Gene_Synonym', 'nomenclature'] = 'Gene_Name'

    # specify nomenclatures used in each dataset to be translated
    # flag 0 means that the ids are UniProt IDs
    nomenclatures = {
        'signalling': 'Ensembl',
        'metabolism': 'Ensembl',
        'gtex': 'Ensembl',
        'eprot': 'Ensembl',
        'orthogroup': 'Gene_Name',
        'depmap': 'Gene_Name',
        'ubiquitination': 'Gene_Name',
        'humap3': 0,
        'ptmdb': 0,
        'proteomehd': 'Gene_Name',
        'mitchell2023': 'Gene_Name',
        'lopit2025': 0,
        'reactome': 'Ensembl',
    }

    gene_lists = {
        'signalling': signalling_genes,
        'metabolism': metabolism_genes,
        'gtex': gtex_genes,
        'eprot': eprot_genes,
        'orthogroup': orthogroup_genes,
        'depmap': depmap_genes,
        'ubiquitination': ubiquitination_genes,
        'humap3': humap3_genes,
        'ptmdb': ptmdb_genes,
        'proteomehd': proteomehd_genes,
        'mitchell2023': mitchell2023_genes,
        'lopit2025': lopit2025_genes,
        'reactome': reactome_genes
    }

    arg1 = [name for name in gene_lists.keys()]
    arg2 = [ref_genes for i in range(len(arg1))]
    arg3 = [idmapping for i in range(len(arg1))]
    arg4 = [nomenclatures for i in range(len(arg1))]
    arg5 = [gene_lists for i in range(len(arg1))]
    with Pool(processes=n_proc) as pool:
        dataset_cols = pool.starmap(create_dict_col, zip(arg1,arg2,arg3,arg4,arg5) )

    gene_dict = pd.DataFrame(dict(ChainMap(*dataset_cols)))
    gene_dict.index = ref_genes

    print(gene_dict.to_csv(sep='\t', header=True, index=True, index_label='gene_name'))

if __name__ == '__main__':
    main()
