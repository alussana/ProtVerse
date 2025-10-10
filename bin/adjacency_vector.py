#!/usr/bin/env python3

import dask.dataframe as dd
import sys


def read_gene_pairs(file_path:str) -> list:
    gene_pairs = []
    with open(file_path) as gene_pairs_fh:
        for line in gene_pairs_fh:
            gene_pair = line.strip().split("\t")
            gene_pairs.append((gene_pair[0], gene_pair[1]))
    return gene_pairs


def get_adjacency_for_gene(gene: str, df: dd, out_dir:str):
    gene_adj = df.loc[(df["source"] == gene) | (df["target"] == gene)].compute()
    gene_adj = gene_adj.apply(
        lambda x: (
            x.drop("target").rename({"source": "target"})
            if x["target"] == gene
            else x.drop("source")
        ),
        axis=1,
    )
    gene_name_path_friendly = gene.replace('/','-') # some gene names contain "/", which causes errors in the following line
    gene_adj.to_csv(f"{out_dir}/{gene_name_path_friendly}.tsv", sep="\t", header=True, index=False)


def main():
    
    edges_parquet = sys.argv[1]
    genes_txt = sys.argv[2]
    out_dir = sys.argv[3]
    
    genes = []
    with open(genes_txt) as gene_fh:
        for line in gene_fh:
            genes.append(line.strip())
  
    df = dd.read_parquet(
        edges_parquet,
        names=["source", "target", "score"],
        engine='pyarrow'
    )

    for gene in genes:
        get_adjacency_for_gene(gene, df, out_dir)
    

if __name__ == "__main__":
    main()
