#!/usr/bin/env python3

import dask.dataframe as dd
import sys


def main():

    df = dd.read_csv(
        "input/file.tsv",
        sep="\\t",
        names=["source", "target", "score"],
        dtype={"source": str, "target": str, "score": float},
    )

    df.to_parquet('edges.parquet', engine='pyarrow')
    

if __name__ == "__main__":
    main()
