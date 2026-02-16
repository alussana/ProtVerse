#!/usr/bin/env python3

import sys

import dask
import dask.dataframe as dd
import pandas as pd
import numpy as np


def prune_disparity_filter_parquet(
    edges_parquet: str,
    out_tsv: str,
    *,
    alpha: float = 0.05,
    engine: str = "pyarrow",
    columns=("source", "target", "score"),
    score_col: str = "score",
    source_col: str = "source",
    target_col: str = "target",
    undirected: bool = True,
    canonicalize_undirected: bool = True,
    keep_rule: str = "either",   # "either" (min) or "both" (max)
    shuffle: str = "tasks",
):
    """
    Disparity filter pruning for a large weighted graph stored as parquet edges.
    Out-of-core via Dask.

    Parameters
    ----------
    edges_parquet : str
        Input parquet path (directory or file glob).
    out_tsv : str
        Output tsv path.
    alpha : float
        Significance threshold; keep edges with alpha_ij < alpha by rule below.
    undirected : bool
        If True, compute significance from both endpoints.
    canonicalize_undirected : bool
        If True, enforce a canonical ordering per edge (min(source,target), max(...))
        to avoid treating swapped pairs as distinct rows.
    keep_rule : {"either","both"}
        For undirected networks:
          - "either": keep if significant from at least one endpoint (min(alpha_s, alpha_t) < alpha)
          - "both":   keep only if significant from both endpoints (max(alpha_s, alpha_t) < alpha)
    shuffle : str
        Dask shuffle strategy ("tasks" is often safe; "p2p" if available can be faster).

    Returns
    -------
    dask.dataframe.DataFrame
        A lazy Dask DataFrame representing the pruned edges (also written to out_parquet).
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
            s = pdf[source_col].to_numpy()
            t = pdf[target_col].to_numpy()
            lo = np.minimum(s, t)
            hi = np.maximum(s, t)
            pdf[source_col] = lo
            pdf[target_col] = hi
            return pdf

        df = df.map_partitions(_canon, meta=df._meta)

    # ---- 1) Compute node degree (k) and strength (s) out-of-core ----
    # Each undirected edge contributes to both endpoints.
    a = df[[source_col, score_col]].rename(columns={source_col: "node", score_col: "w"})
    b = df[[target_col, score_col]].rename(columns={target_col: "node", score_col: "w"})
    nodes = dd.concat([a, b], interleave_partitions=True)

    node_stats = nodes.groupby("node").agg(
        k=("w", "count"),
        s=("w", "sum"),
    ).reset_index()

    # Persist is safe even if large (spills to disk if needed).
    # Comment out for fully lazy behavior.
    node_stats = node_stats.persist()

    # ---- 2) Join node stats onto edges (for both endpoints) ----
    # Join for source endpoint
    src_stats = node_stats.rename(columns={"node": source_col, "k": "k_s", "s": "s_s"})
    df2 = df.merge(src_stats, on=source_col, how="left")

    # Join for target endpoint
    tgt_stats = node_stats.rename(columns={"node": target_col, "k": "k_t", "s": "s_t"})
    df2 = df2.merge(tgt_stats, on=target_col, how="left")

    # ---- 3) Compute disparity alpha values and filter ----
    def _compute_alpha_and_filter(pdf: pd.DataFrame) -> pd.DataFrame:
        w = pdf[score_col].to_numpy(dtype="float64")

        ks = pdf["k_s"].to_numpy(dtype="int64")
        ss = pdf["s_s"].to_numpy(dtype="float64")
        kt = pdf["k_t"].to_numpy(dtype="int64")
        st = pdf["s_t"].to_numpy(dtype="float64")

        # p_ij = w / s_i  (guard against s=0)
        ps = np.divide(w, ss, out=np.zeros_like(w), where=ss > 0)
        pt = np.divide(w, st, out=np.zeros_like(w), where=st > 0)

        # alpha_ij = (1 - p_ij)^(k_i - 1)
        # For k<=1, disparity filter typically treats edges as always significant (alpha=0).
        one_minus_ps = np.clip(1.0 - ps, 0.0, 1.0)
        one_minus_pt = np.clip(1.0 - pt, 0.0, 1.0)

        alpha_s = np.where(ks <= 1, 0.0, np.power(one_minus_ps, ks - 1))
        alpha_t = np.where(kt <= 1, 0.0, np.power(one_minus_pt, kt - 1))

        if keep_rule == "both":
            alpha_edge = np.maximum(alpha_s, alpha_t)
        else:  # "either"
            alpha_edge = np.minimum(alpha_s, alpha_t)

        keep = alpha_edge < alpha

        out = pdf.loc[keep, [source_col, target_col, score_col]].copy()
        # If you want diagnostics, uncomment:
        # out["alpha_s"] = alpha_s[keep]
        # out["alpha_t"] = alpha_t[keep]
        # out["alpha_edge"] = alpha_edge[keep]
        return out

    meta = pd.DataFrame(
        {
            source_col: df2._meta[source_col],
            target_col: df2._meta[target_col],
            score_col: df2._meta[score_col],
        }
    )

    pruned = df2.map_partitions(_compute_alpha_and_filter, meta=meta)

    # ---- 4) Write out-of-core ----
    pruned_df = pruned.compute()

    pruned_df.to_csv(out_tsv, header=True, index=False, sep="\t")

    return pruned


def main():
    edges_pq = sys.argv[1]
    out_pq = sys.argv[2]
    alpha = float(sys.argv[3])

    prune_disparity_filter_parquet(edges_pq, out_pq, alpha=alpha)


if __name__ == "__main__":
    main()
