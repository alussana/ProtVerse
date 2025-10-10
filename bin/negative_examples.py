#!/usr/bin/env python3

import sys
from random import shuffle

def main():
    N = float(sys.argv[1])
    uniprotID_file = sys.argv[2]
    notUniprotID_file = sys.argv[3]

    ids = []
    with open(uniprotID_file) as uniprotIDs:
        for id in uniprotIDs:
            ids.append(id.strip())
    
    not_ids = []
    with open(notUniprotID_file) as notUniprotID:
        for id in notUniprotID:
            not_ids.append(id.strip())

    N_neg = (len(ids) * (len(ids) - 1)) * N

    ids_idx = list(range(len(ids)))
    not_ids_idx = list(range(len(not_ids)))
    shuffle(ids_idx)
    shuffle(not_ids_idx)
    ids_l = len(ids_idx)
    not_ids_l = len(not_ids_idx)
    ids_i = 0
    not_ids_i = 0
    c = 0
    ids_i = 0
    while c < N_neg:
        if ids_i >= ids_l:
            ids_i = 0
        not_ids_i = c % not_ids_l
        if not_ids_i == 0 and c > 0:
            x = not_ids_idx.pop(0)
            not_ids_idx.append(x)
            ids_i = 0
        gene_1 = ids[ids_idx[ids_i]]
        gene_2 = not_ids[not_ids_idx[not_ids_i]]
        print(f'0\tNA\t{gene_1}\t{gene_2}')
        c += 1
        ids_i += 1

if __name__ == '__main__':
    main()