#!/usr/bin/env python

import sys
import pandas as pd

def translateWord(word: str, word_dict: pd.DataFrame, from_col=1, to_col=2):
    trs = word_dict.loc[word_dict[from_col - 1] == word, ][to_col - 1]
    if len(trs) == 0:
        return(None)
    else:
        return(trs.values[0])

def main():
    """
    target_file = 'input/target.txt'
    dict_file = 'dictionary.txt'
    from_col = 1
    to_col = 2
    translate_on = [1,2]
    keep_untranslated = False
    """
    dict_file = sys.argv[1]
    target_file = sys.argv[2]
    from_col = int(sys.argv[3])
    to_col = int(sys.argv[4])
    translate_on = str(sys.argv[5]).split(',')
    translate_on = [int(i) for i in translate_on]
    keep_untranslated = bool(int(sys.argv[6]))

    word_dict = pd.read_csv(dict_file, sep='\t', header=None)

    with open(target_file, 'r') as target_fh:
        for line in target_fh:
            fields = line.removesuffix('\n').split('\t')
            translated = True
            for i in translate_on:
                tr_word = translateWord(fields[i-1], word_dict, from_col=from_col, to_col=to_col)
                if tr_word == None:
                    translated = False
                    break
                else:
                    fields[i-1] = tr_word
            if translated or keep_untranslated:
                print('\t'.join(fields))
            
if __name__ == '__main__':
    main()