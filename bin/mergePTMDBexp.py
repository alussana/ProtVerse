#!/usr/bin/env python3

import pandas as pd
import sys

def main():
    ids_file = sys.argv[1]
    exp_dir = sys.argv[2]
    
    exp_ids = []
    with open(ids_file) as ids:
        for line in ids:
            exp_ids.append(line.strip())

    datasets = {}
    for exp_id in exp_ids:
        table = pd.read_csv(f'{exp_dir}/{exp_id}.tsv', sep='\t', index_col=0)
        table = table[~table.index.duplicated(keep='first')] # see modules/ptmdb.nf for duplicated phosphosites issues
        datasets[exp_id] = table

    merged_exp = pd.concat(list(datasets.values()), axis=1) # outer join

    print(merged_exp.to_csv(sep='\t', na_rep='NA'))

if __name__ == '__main__':
    main()