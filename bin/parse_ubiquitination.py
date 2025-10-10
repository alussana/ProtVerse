#!/usr/bin/env python3

import pandas as pd
import re
import sys


def main():
    
    source_tsv = sys.argv[1]
    
    df = pd.read_csv(source_tsv, sep="\t", header=0, low_memory=False)

    # keep only sites with identification FDR < 0.01 observed in at least one dataset
    df = df.loc[df["PTM_FLR_category"]!="Bronze"]

    # drop columns that are not needed
    columns_to_drop = ["UniProtAC", "PTM_residue", "Decoy_mod"]
    df = df.drop(columns=columns_to_drop)
    patterns = "|".join(["PSM_counts", "FLR", "BinomialScore", "PSMcount"])
    df = df.drop(columns=df.filter(regex=patterns).columns)

    # booleanize datasets
    columns = df.columns
    columns = [re.sub(r"_peptide.*pos$", "", s) for s in columns]

    df.columns = columns

    datasets = columns[2:]

    for dataset in datasets:
        df[dataset] = df[dataset].notna().astype(int)

    print(df.to_csv(sep="\t", index=False, header=True))

if __name__ == "__main__":
    main()
