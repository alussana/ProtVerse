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

    target = []
    with open(target_file, 'r') as target_file_fh:
        for line in target_file_fh:
            target.append(line.strip())
        
    word_dict = pd.read_csv(dict_file, sep='\t', header=None)

    for gene_i in word_dict.index:
        gene_name = word_dict.iloc[gene_i, from_col - 1]
        if gene_name in target:
            gene_name_tr = word_dict.iloc[gene_i, to_col - 1]
            print(gene_name_tr)

if __name__ == '__main__':
    main()