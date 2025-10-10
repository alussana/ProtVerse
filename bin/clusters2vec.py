#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    gene_list_file = sys.argv[1]
    clusters_file = sys.argv[2]

    genes = []
    with open(gene_list_file, 'r')  as gene_list:
        for line in gene_list:
            genes.append(line.strip())
    
    cl = []
    with open(clusters_file, 'r') as clusters:
        for line in clusters:
            items = line.strip().split('\t')
            cl.append(items)

    vecs = pd.DataFrame(index=genes, columns=[i for i in range(len(cl))])
    for i in range(len(genes)):
        gene_name = genes[i]
        vec = [gene_name in cl[c] for c in range(len(cl))]
        vecs.loc[gene_name, ] = vec

    vecs = vecs.astype(int)

    print(vecs.to_csv(sep='\t', index=True))            

if __name__ == '__main__':
    main()