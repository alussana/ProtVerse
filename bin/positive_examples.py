#!/usr/bin/env python3

import sys

def main():
    gene_file = sys.argv[1]
    react_id = sys.argv[2]
    ids = []
    with open(gene_file) as genes:
        for id in genes:
            ids.append(id.strip())
    if len(ids) < 2:
        print()
    else:
        idx = list(range(len(ids)))
        gene_1 = []
        gene_2 = []
        for i in range(len(ids)):
            idx.pop(0)
            for y in range(len(idx)):
                gene_1.append(ids[i])
                gene_2.append(ids[idx[y]])
                print(f'1\t{react_id}\t{ids[i]}\t{ids[idx[y]]}')

if __name__ == '__main__':
    main()