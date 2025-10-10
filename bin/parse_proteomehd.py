#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    """
    table_S1_csv = 'table_S1.csv'
    """
    table_S1_csv = sys.argv[1]

    table_S1 = pd.read_csv(table_S1_csv, sep=',')
    table_S1 = table_S1.loc[pd.isna(table_S1['Gene_names'])==False,]
    
    SILAC_ratios = table_S1.iloc[:,3:]
    
    """
    SILAC_ratios_multiple_gene_names = SILAC_ratios.loc[
        [';' in SILAC_ratios.iloc[i,0] for i in range(len(SILAC_ratios))],
    ]
    SILAC_ratios_multiple_gene_names_exploded = SILAC_ratios_multiple_gene_names.copy()
    SILAC_ratios_multiple_gene_names_exploded['Gene_names'] = SILAC_ratios_multiple_gene_names.apply(lambda x: x['Gene_names'].split(';'), 1)
    SILAC_ratios_multiple_gene_names_exploded = SILAC_ratios_multiple_gene_names_exploded.explode('Gene_names')
    """
        
    SILAC_ratios_exploded = SILAC_ratios.copy()
    SILAC_ratios_exploded['Gene_names'] = SILAC_ratios.apply(lambda x: x['Gene_names'].split(';'), 1)
    SILAC_ratios_exploded = SILAC_ratios_exploded.explode('Gene_names')

    print(SILAC_ratios_exploded.to_csv(sep='\t', index=False, header=True))

if __name__ == '__main__':
    main()