#!/usr/bin/env python3

import sys
import pandas as pd


def read_gene_pairs(file_path:str) -> list:
    gene_pairs = []
    with open(file_path) as gene_pairs_fh:
        for line in gene_pairs_fh:
            gene_pair = line.strip().split("\t")
            gene_pairs.append((gene_pair[0], gene_pair[1]))
    return gene_pairs


def get_adj(adj_dir:str, gene:str):
    try:
        gene_name_path_friendly = gene.replace('/','-')
        adj_df = pd.read_csv(f"{adj_dir}/{gene_name_path_friendly}.tsv", sep="\t")
        return adj_df
    except Exception as e:
        print(e, file=sys.stderr)
        return None


def corr_adj(adj_dir:str, gene_pair:tuple):
    gene_a = gene_pair[0]
    gene_b = gene_pair[1]
    adj_a_df = get_adj(adj_dir, gene_a)
    adj_b_df = get_adj(adj_dir, gene_b)
    if adj_a_df is not None and adj_b_df is not None:
        scores_df = adj_a_df.merge(adj_b_df, on="target")
        corr = round(scores_df["score_x"].corr(scores_df["score_y"]), 3)
        print(f"{gene_a}\t{gene_b}\t{corr}")


def main():
    
    gene_pairs_txt = sys.argv[1]
    adj_dir = sys.argv[2]
    #n_cores = int(sys.argv[3])
    
    gene_pairs = read_gene_pairs(gene_pairs_txt)
    
    for gene_pair in gene_pairs:
        corr_adj(adj_dir, gene_pair)
    

if __name__ == "__main__":
    main()
