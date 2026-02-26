#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    """
    map_file = 'input/mapping.tsv.gz'
    map_file = 'mapping.tsv.gz'
    IDa = 'Ensembl_PRO'
    IDa = 'Gene_Synonym'
    IDb = 'HGNC'
    IDb = 'Gene_Name'
    """
    map_file = sys.argv[1]
    IDa = sys.argv[2]
    IDb = sys.argv[3]
    
    all_map = pd.read_csv(map_file, sep='\t', header=None, index_col=None)
    all_map.columns = ['UniProt_AC', 'ID_type', 'ID']

    uniprot2IDa = all_map.loc[all_map['ID_type']==f'{IDa}', ['UniProt_AC', 'ID']]
    uniprot2IDb = all_map.loc[all_map['ID_type']==f'{IDb}', ['UniProt_AC', 'ID']]

    genes = list(set(uniprot2IDa['ID']))

    IDa_column = []
    uniprot_column = []
    IDb_column = []

    for gene_name in genes:
        # get list of uniprot AC corresponding to gene
        uniprot_tr = list(uniprot2IDa.loc[uniprot2IDa['ID'] == gene_name, 'UniProt_AC'])
        for uniprot_id in uniprot_tr:
            # get list of IDb_column IDs corresponding to uniprot_id
            IDb_tr = list(uniprot2IDb.loc[uniprot2IDb['UniProt_AC'] == uniprot_id, 'ID'])
            for IDb_id in IDb_tr:
                IDa_column.append(gene_name)
                uniprot_column.append(uniprot_id)
                IDb_column.append(IDb_id)
            if len(IDb_tr) == 0:
                IDa_column.append(gene_name)
                uniprot_column.append(uniprot_id)
                IDb_column.append('NA')
        if len(uniprot_tr) == 0:
            IDa_column.append(gene_name)
            uniprot_column.append('NA')
            IDb_column.append('NA')
    
    id_map = pd.DataFrame({
        'IDa': IDa_column,
        'uniprot': uniprot_column,
        'IDb': IDb_column
    })
    
    # in rows where the IDa word appears among the IDb words,
    # map the IDa word to itself
    # this help avoiding annotation mistakes from UniProt
    IDb_words = set(list(id_map['IDb']))
    #id_map = id_map[~id_map['IDa'].isin(IDb_words)]
    id_map['IDb'] = id_map.apply(lambda x: x['IDa'] if x['IDa'] in IDb_words else x['IDb'], axis=1)

    print(id_map.to_csv(sep='\t', index=False, header=False))

if __name__ == '__main__':
    main()
