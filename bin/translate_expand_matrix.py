#!/usr/bin/env python

import sys
import pandas as pd

def main():
    """
    target_file = 'input/target.txt'
    dict_file = 'dictionary.txt'
    from_col = 1
    to_col = 2
    """   
    
    dict_file = sys.argv[1]
    target_file = sys.argv[2]
    from_col = int(sys.argv[3])
    to_col = int(sys.argv[4])
    
    word_dict = pd.read_csv(dict_file, sep='\t', header=None)
    target = pd.read_csv(target_file, sep='\t', header=None)
    
    rows = {}
    for gene_i in word_dict.index:
        gene_name = word_dict.iloc[gene_i, from_col - 1]
        gene_name_tr = word_dict.iloc[gene_i, to_col - 1]
        row = target.loc[target[0] == gene_name, ].iloc[0, 1:].values
        if len(row) != 0:
            rows[gene_name_tr] = row

    df = pd.DataFrame.from_dict(rows, orient='index')
    print(df.to_csv(sep='\t', header=False, index=True))
            
if __name__ == '__main__':
    main()