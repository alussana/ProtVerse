#!/usr/bin/env python3

import sys
import pandas as pd
import collections

def main():
    """
    eprot_file = 'training/eprot.tsv.gz'
    proteomehd_file = 'training/proteomehd.tsv.gz'
    mitchell2023_file = 'training/mitchell2023.tsv.gz'
    gtex_file = 'training/gtex.tsv.gz'
    ptmdb_file = 'training/ptmdb.tsv.gz'
    ubiquitination_file = 'training/ubiquitination.tsv.gz'
    dependency_file = 'training/dependency.tsv.gz'
    orthogroup_file = 'training/orthogroup.tsv.gz'
    humap3_file = 'training/humap3.tsv.gz'
    """
    eprot_file = sys.argv[1]
    proteomehd_file = sys.argv[2]
    mitchell2023_file = sys.argv[3]
    gtex_file = sys.argv[4]
    ptmdb_file = sys.argv[5]
    ubiquitination_file = sys.argv[6]
    dependency_file = sys.argv[7]
    orthogroup_file = sys.argv[8]
    human3_file = sys.argv[9]

    eprot = pd.read_csv(eprot_file, sep='\t', index_col=0)
    proteomehd = pd.read_csv(proteomehd_file, sep='\t', index_col=0)
    mitchell2023 = pd.read_csv(mitchell2023_file, sep='\t', index_col=0)
    gtex = pd.read_csv(gtex_file, sep='\t', index_col=0)
    ptmdb = pd.read_csv(ptmdb_file, sep='\t', index_col=0)
    ubiquitination = pd.read_csv(ubiquitination_file, sep='\t', index_col=0)
    dependency = pd.read_csv(dependency_file, sep='\t', index_col=0)
    orthogroup = pd.read_csv(orthogroup_file, sep='\t', index_col=0)
    humap3 = pd.read_csv(human3_file, sep='\t', index_col=0)

    dfs = [eprot, proteomehd, mitchell2023, gtex, ptmdb, ubiquitination, dependency, orthogroup, humap3]

    for i in range(len(dfs)):
        duplicated_idx = [item for item, count in collections.Counter(dfs[i].index).items() if count > 1]
        dfs[i] = dfs[i].loc[([idx not in duplicated_idx for idx in list(dfs[i].index)]) or (dfs[i].label == 1), ]

    labels = dfs[0]['label']
    for i in range(len(dfs)):
        dfs[i] = dfs[i].drop(columns=['label'])

    dfs.append(pd.DataFrame(labels))
    examples = pd.concat(dfs, axis=1, join='inner')

    print(examples.to_csv(sep='\t'))

if __name__ == '__main__':
    main()