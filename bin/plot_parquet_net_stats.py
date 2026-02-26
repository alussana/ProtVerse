#!/usr/bin/env python3

import sys

import dask
import dask.dataframe as dd
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns


sns.set_theme(style="ticks")

plt.rcParams['axes.titlesize']     = 7
plt.rcParams['axes.labelsize']     = 6
plt.rcParams['xtick.labelsize']    = 6
plt.rcParams['ytick.labelsize']    = 6
plt.rcParams['legend.fontsize']    = 6
plt.rcParams['figure.titlesize']   = 7
plt.rcParams['font.size']          = 6
plt.rcParams['axes.linewidth']     = 0.6
plt.rcParams['xtick.major.width']  = 0.6
plt.rcParams['ytick.major.width']  = 0.6
plt.rcParams['xtick.minor.width']  = 0.4
plt.rcParams['ytick.minor.width']  = 0.4
plt.rcParams['xtick.major.size']   = 3
plt.rcParams['ytick.major.size']   = 3
plt.rcParams['xtick.minor.size']   = 2
plt.rcParams['ytick.minor.size']   = 2

plt.rcParams["figure.figsize"]     = (3, 2)

red_shade = "#d62728"
grey_shade = "#7f7f7f"

custom_palette = sns.color_palette([red_shade, grey_shade], n_colors=2)
sns.set_palette(custom_palette)


def plot_parquet_net_stats(
    edges_parquet: str,
    out_prefix: str,
    *,
    engine: str = "pyarrow",
    columns=("source", "target", "score"),
    score_col: str = "score",
    source_col: str = "source",
    target_col: str = "target",
    undirected: bool = True,
    canonicalize_undirected: bool = True,
    split_out: int = 64,
    shuffle: str = "tasks",
) -> None:
    """
    Generate plots for a weighted graph stored as parquet edges.
    

    Parameters
    ----------
    edges_parquet : str
        Input parquet path (directory or file glob).
    out_prefix : str
        Output path prefix to be used for pdf plots.
    undirected : bool
        If True, compute significance from both endpoints.
    canonicalize_undirected : bool
        If True, enforce a canonical ordering per edge (min(source,target), max(...))
        to avoid treating swapped pairs as distinct rows.
    split_out : int
        Controls parallelism for groupby reductions (tune for cluster size).
    shuffle : str
        Dask shuffle strategy ("tasks" is often safe; "p2p" if available can be faster).

    Returns
    -------
    None
    """
    dask.config.set({
        "dataframe.shuffle.method": shuffle
    })

    # Read edges out-of-core
    df = dd.read_parquet(edges_parquet, engine=engine, columns=list(columns))

    # Optionally canonicalize undirected edges so (u,v) and (v,u) are identical
    if undirected and canonicalize_undirected:
        # Use map_partitions to avoid creating a giant intermediate in memory
        def _canon(pdf: pd.DataFrame) -> pd.DataFrame:
            s = pdf[source_col].values
            t = pdf[target_col].values
            lo = np.minimum(s, t)
            hi = np.maximum(s, t)
            pdf[source_col] = lo
            pdf[target_col] = hi
            return pdf

        df = df.map_partitions(_canon, meta=df._meta)

    # Compute node degree (k) out-of-core
    # Each undirected edge contributes to both endpoints.
    a = df[[source_col, score_col]].rename(columns={source_col: "node", score_col: "w"})
    b = df[[target_col, score_col]].rename(columns={target_col: "node", score_col: "w"})
    nodes = dd.concat([a, b], interleave_partitions=True)

    node_stats = nodes.groupby("node").agg(
        k=("w", "count"),
        avg_w=("w", "mean"),
        split_out=split_out,
    ).reset_index()

    # Execute computation graph
    node_stats_df = node_stats.compute()
    k = node_stats["k"]
    avg_w = node_stats["avg_w"]
    
    # Number of nodes
    n = len(k)

    # Plot degree distribution on histogram with logarithmic axis
    bins = np.logspace(np.log10(k.min()), np.log10(k.max()), 50)
    plt.clf()
    sns.histplot(
        k,
        bins=bins
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Degree (k)")
    plt.ylabel("Number of nodes")
    plt.title(f"Degree distribution (log–log) (n={n})")
    plt.tight_layout()
    plt.savefig(f'{out_prefix}_degree_loglog_hist.pdf')

    bins = np.logspace(np.log10(k.min()), np.log10(k.max()), 50)
    plt.clf()
    sns.histplot(
        k,
        bins=bins
    )
    plt.xscale("log")
    plt.xlabel("Degree (k)")
    plt.ylabel("Number of nodes")
    plt.title(f"Degree distribution (log–log) (n={n})")
    plt.tight_layout()
    plt.savefig(f'{out_prefix}_degree_loglin_hist.pdf')

    # Plot degree complementary cumulative density function
    k_sorted = np.sort(node_stats["k"].values)
    ccdf = 1.0 - np.arange(1, len(k_sorted) + 1) / len(k_sorted)
    plt.clf()
    plt.figure()
    plt.loglog(k_sorted, ccdf, marker=".", linestyle="none")
    plt.xlabel("Degree (k)")
    plt.ylabel("P(K ≥ k)")
    plt.title(f"Degree CCDF (n={n})")
    plt.tight_layout()
    plt.savefig(f'{out_prefix}_degree_loglog_ccdf.pdf')

    # Plot per-node mean edge weight histogram
    plt.clf()
    sns.histplot(
        avg_w,
    )
    plt.xlabel("Mean weight)")
    plt.ylabel("Number of nodes")
    plt.title(f"Distribution of mean edge weight (n={n})")
    plt.tight_layout()
    plt.savefig(f'{out_prefix}_mean_weight_hist.pdf')

    # Plot correlation between mean edge weight and node degree
    plt.clf()
    sns.scatterplot(
        x=k,
        y=avg_w,
        alpha=0.5
    )
    plt.xlabel("Degree (k)")
    plt.ylabel("Mean edge weight")
    plt.title(f"Mean edge weight vs degree (n={n})")
    plt.tight_layout()
    plt.savefig(f'{out_prefix}_mean_weight_vs_degree.pdf')


def main():
    edges_pq = sys.argv[1]
    out_prefix = sys.argv[2]

    plot_parquet_net_stats(edges_pq, out_prefix)


if __name__ == "__main__":
    main()
