#!/usr/bin/env python3

import sys
import pandas as pd


def idmapping_to_list(idmapping_tsv, max_nan):
    df = pd.read_csv(idmapping_tsv, sep='\t', index_col=0)
    df.drop(columns=['reactome', 'metabolism', 'signalling'], inplace=True)
    df = df[df.isna().sum(axis=1) <= max_nan]
    l = list(df.index)
    return l


def main():
    idmapping_tsv = sys.argv[1]
    max_na_int = int(sys.argv[2])

    l = idmapping_to_list(idmapping_tsv, max_na_int)
    
    for i in range(len(l) - 1):
        item_a = l[i]
        for j in range(i + 1, len(l)):
            item_b = l[j]
            print(f'{item_a}\t{item_b}')


if __name__ == '__main__':
    main()