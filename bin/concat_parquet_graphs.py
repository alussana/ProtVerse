#!/usr/bin/env python3

import dask
import dask.dataframe as dd
import sys


def cat_undirected_pq_graphs(
    edges_parquet_1: str,
    edges_parquet_2: str,
    out_path: str,
    *,
    engine: str = "pyarrow",
    columns=("source", "target", "score"),
    score_col: str = "score",
    source_col: str = "source",
    target_col: str = "target",
    split_out: int = 64,
    shuffle: str = "tasks",
):
    dask.config.set({
        "dataframe.shuffle.method": shuffle
    })

    e1 = dd.read_parquet(edges_parquet_1, columns=columns, engine=engine)
    e2 = dd.read_parquet(edges_parquet_2, columns=columns, engine=engine)

    # union edges (stack)
    all_edges = dd.concat([e1, e2], interleave_partitions=True)

    # drop exact duplicates before the expensive groupby
    all_edges = all_edges.drop_duplicates(subset=[source_col, target_col, score_col])

    # group to keep max score per undirected edge
    merged = (
        all_edges.groupby([source_col, target_col])["score"]
        .max(split_out=split_out)
        .reset_index()
    )

    merged.to_parquet(out_path, engine=engine, write_index=False)


def main():
    graph1_pq = sys.argv[1]
    graph2_pq = sys.argv[2]
    out_prefix = sys.argv[3]

    cat_undirected_pq_graphs(graph1_pq, graph2_pq, out_prefix)


if __name__ == '__main__':
    main()
