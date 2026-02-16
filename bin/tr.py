#!/usr/bin/env python3

# Usage: in the following example the words found in the second field
# (sys.argv[5]) of file.tsv (sys.argv[2]) that are specified in the first
# column (sys.argv[3]) of dict.tsv (sys.argv[1]) are translated with the
# corresponding word found in the third column (sys.argv[4]) of dict.tsv.
# Rows of file.tsv where no translation could occurr won't be discarded 
# (sys.argv[6]).
#
#   translator.py \
#       dict.tsv \
#       file.tsv \
#       1 \
#       3 \
#       2 \
#       1 \
#       > translated_file.tsv
#
# Specifying multiple fields on which to perform the translation is also
# possible, e.g. to translate on the first and the second field of file.tsv:
#
#   translator.py \
#       dict.tsv \
#       file.tsv \
#       1 \
#       3 \
#       1,2 \
#       1 \
#       > translated_file.tsv
#
# or you can also use 'all' to translate everywhere:
#
#   translator.py \
#       dict.tsv \
#       file.tsv \
#       1 \
#       3 \
#       all \
#       1 \
#       > translated_file.tsv
#

import sys
import pandas as pd


def translateWord(word: str, word_dict: pd.DataFrame, from_col=1, to_col=2):
    trs = word_dict.loc[word_dict[from_col - 1] == word,][to_col - 1]
    if len(trs) == 0:
        return None
    else:
        return trs.values[0]


def main():
    dict_file = sys.argv[1]
    target_file = sys.argv[2]
    from_col = int(sys.argv[3])
    to_col = int(sys.argv[4])
    translate_on = sys.argv[5]
    keep_untranslated = bool(int(sys.argv[6]))

    if translate_on != "all":
        translate_on = str(translate_on).split(",")
        try:
            translate_on = [int(i) for i in translate_on]
        except:
            translate_on = [int(translate_on)]

    word_dict = pd.read_csv(dict_file, sep="\t", header=None).astype('string')
    word_dict = word_dict.dropna()

    if translate_on != "all":
        with open(target_file, "r") as target_fh:
            for line in target_fh:
                fields = line.removesuffix("\n").split("\t")
                translated = True
                for i in translate_on:
                    tr_word = translateWord(
                        fields[i - 1], word_dict, from_col=from_col, to_col=to_col
                    )
                    if tr_word == None:
                        translated = False
                    else:
                        fields[i - 1] = tr_word
                if translated or keep_untranslated:
                    print("\t".join([str(field) for field in fields]))
    else:
        with open(target_file, "r") as target_fh:
            for line in target_fh:
                fields = line.removesuffix("\n").split("\t")
                translated = True
                for i in range(len(fields)):
                    tr_word = translateWord(
                        fields[i - 1], word_dict, from_col=from_col, to_col=to_col
                    )
                    if tr_word == None:
                        translated = False
                    else:
                        fields[i - 1] = tr_word
                if translated or keep_untranslated:
                    print("\t".join([str(field) for field in fields]))


if __name__ == "__main__":
    main()
