#!/usr/bin/env python

import sys
import pandas as pd

def main():
    """
    target_file = 'input/target.txt'
    dict_file = 'dictionary.txt'
    from_col = 1
    to_col = 2
    keep_untranslated = False
    """   
    
    dict_file = sys.argv[1]
    target_file = sys.argv[2]
    from_col = int(sys.argv[3])
    to_col = int(sys.argv[4])
    keep_untranslated = bool(int(sys.argv[5]))
    
    word_dict = pd.read_csv(dict_file, sep='\t', header=None)

    ref_genes = set(word_dict.iloc[:, to_col - 1].values)

    with open(target_file, 'r') as target_fh:
        for line in target_fh:
            fields = line.strip().split('\t')
            if len(fields) > 1:
                translated = True
                # get gene names
                gene_0 = fields[0]
                gene_1 = fields[1]
                # translate gene names (multiple translations allowed)
                gene_0_tr = list(
                    word_dict.loc[
                        word_dict[from_col - 1] == gene_0,
                        ][
                            to_col - 1
                        ].values
                )
                gene_1_tr = list(
                    word_dict.loc[
                        word_dict[from_col - 1] == gene_1,
                    ][
                        to_col - 1
                    ].values
                )
                # remove translation that do not belong to
                # reference proteomes gene names
                gene_0_tr = [gene for gene in gene_0_tr if gene in ref_genes]
                gene_1_tr = [gene for gene in gene_1_tr if gene in ref_genes]
                # flag if translation was found for both genes
                # if not gene name is reversed to untranslated version
                if len(gene_0_tr) == 0:
                    gene_0_tr = [gene_0]
                    translated = False
                if len(gene_1_tr) == 0:
                    gene_1_tr = [gene_1]
                    translated = False
                # print line if appropriate
                # if multiple reference gene names are found per translation
                # print every combination (TODO)
                if translated or keep_untranslated:
                    for i in range(len(gene_0_tr)):
                        for j in range(len(gene_1_tr)):
                            fields[0] = gene_0_tr[i]
                            fields[1] = gene_1_tr[j]
                            print('\t'.join(fields))
            
if __name__ == '__main__':
    main()