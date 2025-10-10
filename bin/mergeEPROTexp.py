#!/usr/bin/env python3

from email import header
import pandas as pd
import sys

def main():
    ids_file = sys.argv[1]
    exp_dir = sys.argv[2]
    counts_out = sys.argv[3]
    
    eprot_ids = []
    with open(ids_file) as ids:
        for line in ids:
            eprot_ids.append(line.strip())

    datasets = {}
    for eprot_id in eprot_ids:
        table = pd.read_csv(f'{exp_dir}/{eprot_id}.tsv', sep='\t', index_col=0)
        #table = table.loc[:, ~table.columns.str.contains('^Unnamed')]
        datasets[eprot_id] = table

    merged_exp = pd.concat(list(datasets.values()), axis=1) # outer join

    non_NA_counts = merged_exp.count(axis=1)
    non_NA_counts.to_csv(counts_out, sep='\t', header=False)

    print(merged_exp.to_csv(sep='\t', na_rep='NA'))


if __name__ == '__main__':
    main()