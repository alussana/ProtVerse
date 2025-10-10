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
    metabolism_genes_file = 'input/metabolism_genes.txt'
    n_proc = 2
    """
    ref_genes_file = sys.argv[1]
    id_mapping_file = sys.argv[2]
    metabolism_genes_file = sys.argv[3]
    n_proc = int(sys.argv[4])

    ref_genes = []
    with open(ref_genes_file) as genes_fh:
        for gene in genes_fh:
            ref_genes.append(gene.strip())

    metabolism_genes = []
    with open(metabolism_genes_file) as genes_fh:
        for gene in genes_fh:
            metabolism_genes.append(gene.strip())

    idmapping = pd.read_csv(id_mapping_file, sep='\t', header=None, index_col=None)
    idmapping.columns = ['uniprot', 'nomenclature', 'value']
    idmapping.loc[idmapping['nomenclature'] == 'Gene_Synonym', 'nomenclature'] = 'Gene_Name'

    # specify nomenclatures used in each dataset to be translated
    # flag 0 means that the ids are UniProt IDs
    nomenclatures = {
        'metabolism': 'Ensembl'
    }

    gene_lists = {
        'metabolism': metabolism_genes
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
