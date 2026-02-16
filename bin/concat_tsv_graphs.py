#!/usr/bin/env python3

import pandas as pd
import sys


def main():
    graph1_tsv = sys.argv[1]
    graph2_tsv = sys.argv[2]

    df1 = pd.read_csv(graph1_tsv, sep="\t")
    df2 = pd.read_csv(graph2_tsv, sep="\t")

    merged_df = pd.concat([df1, df2]).groupby(["source", "target"]).agg("max").reset_index()

    print(merged_df.to_csv(sep="\t", header=True, index=False))


if __name__ == '__main__':
    main()
