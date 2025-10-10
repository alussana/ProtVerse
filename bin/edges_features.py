#!/usr/bin/env python3

import sys
import pandas as pd
#import dask.dataframe as dd

def main():
    """
    mkdir -p edges_features/
    for table in $(ls *gz); do zcat training/${table} > edges_features/${table::-3}; done

    eprot_file = 'edges_features/eprot.tsv'
    proteomehd_file = 'edges_features/proteomehd.tsv'
    gtex_file = 'edges_features/gtex.tsv'
    ptmdb_file = 'edges_features/ptmdb.tsv'
    dependency_file = 'edges_features/dependency.tsv'
    orthogroup_file = 'edges_features/orthogroup.tsv'
    ubiquitination_file = 'edges_features/ubiquitination.tsv'
    humap_file = 'edges_features/humap.tsv'
    """
    eprot_file = sys.argv[1]
    proteomehd_file = sys.argv[2]
    mitchell2023_file = sys.argv[3]
    gtex_file = sys.argv[4]
    ptmdb_file = sys.argv[5]
    dependency_file = sys.argv[6]
    orthogroup_file = sys.argv[7]
    ubiquitination_file = sys.argv[8]
    humap_file = sys.argv[9]

    # use pandas
    eprot = pd.read_csv(eprot_file, sep='\t').set_index('index')
    proteomehd = pd.read_csv(proteomehd_file, sep='\t').set_index('index')
    mitchell2023 = pd.read_csv(mitchell2023_file, sep='\t').set_index('index')
    gtex = pd.read_csv(gtex_file, sep='\t').set_index('index')
    ptmdb = pd.read_csv(ptmdb_file, sep='\t').set_index('index')
    dependency = pd.read_csv(dependency_file, sep='\t').set_index('index')
    orthogroup = pd.read_csv(orthogroup_file, sep='\t').set_index('index')
    ubiquitination = pd.read_csv(ubiquitination_file, sep='\t').set_index('index')
    humap = pd.read_csv(humap_file, sep='\t').set_index('index')
    """ use dask
    eprot = dd.read_csv(eprot_file, sep='\t').set_index('index')
    gtex = dd.read_csv(gtex_file, sep='\t').set_index('index')
    ptmdb = dd.read_csv(ptmdb_file, sep='\t').set_index('index')
    dependency = dd.read_csv(dependency_file, sep='\t').set_index('index')
    orthogroup = dd.read_csv(orthogroup_file, sep='\t').set_index('index')
    ubiquitination = dd.read_csv(ubiquitination_file, sep='\t').set_index('index')
    """

    """
    dfs = [eprot, gtex, ptmdb, dependency, orthogroup, ubiquitination]
    for i in range(len(dfs)):
        duplicated_idx = [item for item, count in collections.Counter(dfs[i].index).items() if count > 1]
        dfs[i] = dfs[i].loc[([idx not in duplicated_idx for idx in list(dfs[i].index)]) or (dfs[i].label == 1), ]
    labels = dfs[0]['label']
    for i in range(len(dfs)):
        dfs[i] = dfs[i].drop(columns=['label'])
    """

    # merge all partial input vectors
    # use pandas
    dfs = [eprot, proteomehd, mitchell2023, gtex, ptmdb, dependency, orthogroup, ubiquitination, humap]
    edges_features = pd.concat(dfs, axis=1)
    """ use dask
    edges_features = dd.merge(eprot, gtex)
    edges_features = dd.merge(edges_features, ptmdb)
    edges_features = dd.merge(edges_features, dependency)
    edges_features = dd.merge(edges_features, orthogroup)
    edges_features = dd.merge(edges_features, ubiquitination)
    edges_features.to_csv('test.tsv', sep='\t') # tmp
    """

    print(edges_features.to_csv(sep='\t'))

if __name__ == '__main__':
    main()