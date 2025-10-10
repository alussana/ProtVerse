#!/usr/bin/env python

import sys
import pandas as pd

def translateWord(word: str, word_dict: pd.DataFrame, from_col:str, to_col:str):
    trs = word_dict.loc[word_dict[from_col] == word, ][to_col]
    if len(trs) == 0:
        return(None)
    else:
        return(trs.values[0])

def main():
    """
    target_file = 'input/tab.tsv'
    dict_file = 'input/dictionary.txt'
    from_col = "reactome"
    to_col = "gene_name"
    translate_on = 'all'
    keep_untranslated_fields = False
    """
    dict_file = sys.argv[1]
    target_file = sys.argv[2]
    from_col = sys.argv[3]
    to_col = sys.argv[4]
    translate_on = sys.argv[5]
    fields_to_translate = str(translate_on).split(',')
    if translate_on != 'all':
        fields_to_translate = [int(i) for i in fields_to_translate]
    keep_untranslated_fields = bool(int(sys.argv[6]))

    word_dict = pd.read_csv(dict_file, sep='\t')

    with open(target_file, 'r') as target_fh:
        for line in target_fh:
            fields = line.strip().split('\t')
            translated = True
            untr_idx = []
            if translate_on == 'all':
                fields_to_translate = range(1, len(fields) + 1)
            for i in fields_to_translate:
                tr_word = translateWord(fields[i-1], word_dict, from_col=from_col, to_col=to_col)
                if tr_word == None:
                    translated = False
                    untr_idx.append(i)
                else:
                    fields[i-1] = tr_word
            if translated or keep_untranslated_fields:
                print('\t'.join(fields))
            else:
                c = 1
                for idx in untr_idx:
                    fields.pop(idx - c)
                    c = c + 1
                print('\t'.join(fields))
            
if __name__ == '__main__':
    main()