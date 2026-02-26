#!/usr/bin/env python3
"""
Network propagation (Random Walk with Restart) with degree-matched seed permutations.

Inputs
------
1) First argument  : two tab-separated columns describing an undirected graph
2) Second argument : one node id per line

Outputs
-------
A TSV with per-node:
- RWR score using real seeds
- empirical p-value from degree-matched seed permutations
- degree, and whether the node is a real seed

Notes
-----
- RWR is computed by power iteration:
    p_{t+1} = (1-r) * P^T p_t + r * p0
  where P is row-stochastic (random-walk transition) from the adjacency matrix.
- Permutations: sample the same number of seeds per degree as the real seed set.
  Exact degree matching is used. If a degree bucket is too small, sampling for that
  degree falls back to sampling *with replacement* (preserving the degree multiset).
- Parallelism: uses multiprocessing with fork where available.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from multiprocessing import get_context
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm


# --------------------------
# Logging
# --------------------------
def setup_logging(verbosity: int) -> None:
    level = logging.INFO
    if verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# --------------------------
# I/O helpers
# --------------------------
def read_seeds(path: str) -> List[str]:
    seeds: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                seeds.append(s)
    return seeds


def read_edges(path: str, has_header: bool, col_v: str, col_u: str) -> pd.DataFrame:
    """
    Read edges as a DataFrame with columns [v, u] (strings).
    """
    if has_header:
        df = pd.read_csv(
            path,
            sep="\t",
            compression="gzip" if path.endswith(".gz") else None,
            dtype=str,
        )
        if col_v not in df.columns or col_u not in df.columns:
            raise ValueError(
                f"Header read, but columns {col_v!r}/{col_u!r} not found. "
                f"Available columns: {list(df.columns)}"
            )
        df = df[[col_v, col_u]].rename(columns={col_v: "v", col_u: "u"})
    else:
        df = pd.read_csv(
            path,
            sep="\t",
            compression="gzip" if path.endswith(".gz") else None,
            header=None,
            names=["v", "u"],
            usecols=[0, 1],
            dtype=str,
        )
    df = df.dropna()
    df["v"] = df["v"].astype(str)
    df["u"] = df["u"].astype(str)
    return df


# --------------------------
# Graph + RWR
# --------------------------
@dataclass
class GraphData:
    nodes: np.ndarray           # index -> node_id (dtype=object/str)
    node_to_idx: Dict[str, int] # node_id -> index
    adj: sparse.csr_matrix      # adjacency (undirected, unweighted)
    P_T: sparse.csr_matrix      # transpose of row-stochastic transition matrix
    degree: np.ndarray          # degree per node (int64)


def build_graph(edges: pd.DataFrame) -> GraphData:
    """
    Build sparse adjacency and row-stochastic transition (then transpose).
    """
    logging.info("Factorizing node ids...")
    all_nodes = pd.concat([edges["v"], edges["u"]], ignore_index=True)
    codes, uniques = pd.factorize(all_nodes, sort=False)
    n_edges = len(edges)
    v_codes = codes[:n_edges].astype(np.int64, copy=False)
    u_codes = codes[n_edges:].astype(np.int64, copy=False)
    nodes = uniques.astype(object)

    n = len(nodes)
    logging.info("Building sparse adjacency (%d nodes, %d edges)...", n, n_edges)

    # Undirected: add both directions
    row = np.concatenate([v_codes, u_codes])
    col = np.concatenate([u_codes, v_codes])
    data = np.ones(row.shape[0], dtype=np.float32)

    adj = sparse.coo_matrix((data, (row, col)), shape=(n, n), dtype=np.float32).tocsr()
    adj.sum_duplicates()

    # Degree (for undirected simple graph this is row-sum)
    deg = np.asarray(adj.sum(axis=1)).reshape(-1).astype(np.int64)

    # Build P = D^{-1} A (row-stochastic). For isolated nodes (deg=0), keep row all zeros.
    inv_deg = np.zeros_like(deg, dtype=np.float32)
    nz = deg > 0
    inv_deg[nz] = 1.0 / deg[nz].astype(np.float32)
    D_inv = sparse.diags(inv_deg, offsets=0, format="csr", dtype=np.float32)
    P = D_inv @ adj
    P_T = P.transpose().tocsr()

    node_to_idx = {str(nodes[i]): i for i in range(n)}
    return GraphData(nodes=nodes, node_to_idx=node_to_idx, adj=adj, P_T=P_T, degree=deg)


def rwr_power_iteration(
    P_T: sparse.csr_matrix,
    seeds_idx: np.ndarray,
    n_nodes: int,
    restart: float,
    tol: float,
    max_iter: int,
) -> np.ndarray:
    """
    Random Walk with Restart via power iteration.

    Returns:
        p (float64) stationary-ish distribution (scores), length n_nodes
    """
    if len(seeds_idx) == 0:
        raise ValueError("Seed set is empty after filtering to nodes present in the graph.")

    p0 = np.zeros(n_nodes, dtype=np.float64)
    # If duplicates exist in seeds_idx, they naturally increase weight; normalize anyway:
    # (this matters only if degree buckets were sampled with replacement).
    for i in seeds_idx:
        p0[i] += 1.0
    p0 /= p0.sum()

    p = p0.copy()
    one_minus = 1.0 - restart

    for it in range(max_iter):
        p_next = one_minus * (P_T @ p) + restart * p0
        # L1 norm convergence is common for RWR
        diff = np.abs(p_next - p).sum()
        p = p_next
        if diff < tol:
            logging.debug("RWR converged at iter=%d (L1 diff=%.3e)", it + 1, diff)
            break
    else:
        logging.warning("RWR hit max_iter=%d without meeting tol=%.2e", max_iter, tol)

    return p


# --------------------------
# Degree-matched permutations
# --------------------------
def build_degree_buckets(degree: np.ndarray) -> Dict[int, np.ndarray]:
    buckets: Dict[int, List[int]] = defaultdict(list)
    for idx, d in enumerate(degree.tolist()):
        buckets[int(d)].append(idx)
    return {d: np.asarray(idxs, dtype=np.int64) for d, idxs in buckets.items()}


def degree_matched_sample(
    rng: np.random.Generator,
    seed_deg_counts: Counter,
    buckets: Dict[int, np.ndarray],
) -> np.ndarray:
    """
    Sample node indices to match the exact degree multiset of the real seeds.
    If a bucket is too small for sampling without replacement, sample with replacement.
    """
    sampled: List[np.ndarray] = []
    for d, k in seed_deg_counts.items():
        pool = buckets.get(int(d))
        if pool is None or len(pool) == 0:
            # This should not happen if degrees were computed on the same graph,
            # but keep a safe fallback.
            raise RuntimeError(f"No nodes available with degree={d} to match seeds.")
        replace = k > len(pool)
        sampled.append(rng.choice(pool, size=k, replace=replace))
    return np.concatenate(sampled)


# --------------------------
# Multiprocessing worker
# --------------------------
# Globals for forked workers (kept read-only)
_G_P_T = None
_G_NODES = None
_G_REAL = None
_G_BUCKETS = None
_G_SEED_DEG_COUNTS = None
_G_RWR_ARGS = None


def _init_worker(P_T, n_nodes, real_scores, buckets, seed_deg_counts, rwr_args):
    global _G_P_T, _G_NODES, _G_REAL, _G_BUCKETS, _G_SEED_DEG_COUNTS, _G_RWR_ARGS
    _G_P_T = P_T
    _G_NODES = n_nodes
    _G_REAL = real_scores
    _G_BUCKETS = buckets
    _G_SEED_DEG_COUNTS = seed_deg_counts
    _G_RWR_ARGS = rwr_args


def _worker_run(task: Tuple[int, int, int]) -> np.ndarray:
    """
    task = (n_perms_for_this_worker, base_seed, worker_id)
    Returns count_ge: int64 array length n_nodes
    """
    n_perms, base_seed, worker_id = task
    rng = np.random.default_rng(base_seed + 1000003 * worker_id)

    count_ge = np.zeros(_G_NODES, dtype=np.int64)

    restart, tol, max_iter = _G_RWR_ARGS

    for _ in range(n_perms):
        perm_seeds = degree_matched_sample(rng, _G_SEED_DEG_COUNTS, _G_BUCKETS)
        perm_scores = rwr_power_iteration(
            _G_P_T, perm_seeds, _G_NODES, restart=restart, tol=tol, max_iter=max_iter
        )
        count_ge += (perm_scores >= _G_REAL).astype(np.int64)

    return count_ge


# --------------------------
# Main
# --------------------------
def main():
    ap = argparse.ArgumentParser(
        description="RWR network propagation + degree-matched seed permutation null + empirical p-values"
    )
    ap.add_argument("--edges", default="input/edges.tsv.gz", help="Path to edges TSV(.gz)")
    ap.add_argument("--seeds", default="input/seeds.txt", help="Path to seeds.txt")
    ap.add_argument("--out", default="output/rwr_pvalues.tsv", help="Output TSV path")
    ap.add_argument("--n-perm", type=int, default=100, help="Number of degree-matched permutations")
    ap.add_argument("--restart", type=float, default=0.5, help="Restart probability r in [0,1]")
    ap.add_argument("--tol", type=float, default=1e-10, help="Convergence tolerance (L1 diff)")
    ap.add_argument("--max-iter", type=int, default=200, help="Maximum iterations for power method")
    ap.add_argument("--n-jobs", type=int, default=0, help="Parallel workers (0 => all cores)")
    ap.add_argument("--seed", type=int, default=1, help="Base RNG seed for permutations")

    ap.add_argument("--edges-has-header", action="store_true", help="Edges file has header row")
    ap.add_argument("--col-v", default="v", help="Column name for v (if header present)")
    ap.add_argument("--col-u", default="u", help="Column name for u (if header present)")

    ap.add_argument("-v", "--verbose", action="count", default=1, help="Increase logging verbosity")
    args = ap.parse_args()

    setup_logging(args.verbose)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    logging.info("Reading edges: %s", args.edges)
    edges = read_edges(args.edges, args.edges_has_header, args.col_v, args.col_u)
    logging.info("Edges loaded: %d", len(edges))

    graph = build_graph(edges)

    logging.info("Reading seeds: %s", args.seeds)
    seeds = read_seeds(args.seeds)
    logging.info("Seeds loaded: %d", len(seeds))

    # Keep only seeds present in graph
    seeds_in_graph = [s for s in seeds if s in graph.node_to_idx]
    missing = len(seeds) - len(seeds_in_graph)
    if missing:
        logging.warning("Dropped %d seeds not found in graph nodes.", missing)

    seeds_idx = np.asarray([graph.node_to_idx[s] for s in seeds_in_graph], dtype=np.int64)
    is_seed = np.zeros(len(graph.nodes), dtype=bool)
    is_seed[seeds_idx] = True

    # 1) Real RWR
    logging.info("Computing real RWR scores...")
    real_scores = rwr_power_iteration(
        graph.P_T,
        seeds_idx,
        n_nodes=len(graph.nodes),
        restart=args.restart,
        tol=args.tol,
        max_iter=args.max_iter,
    ).astype(np.float64)

    # 2) Null via degree-matched permutations
    seed_degrees = graph.degree[seeds_idx].astype(np.int64)
    seed_deg_counts = Counter(seed_degrees.tolist())
    buckets = build_degree_buckets(graph.degree)

    n_perm = int(args.n_perm)
    if n_perm < 1:
        raise ValueError("--n-perm must be >= 1")

    n_jobs = args.n_jobs if args.n_jobs and args.n_jobs > 0 else os.cpu_count() or 1
    n_jobs = max(1, int(n_jobs))

    # Split work across workers
    per_worker = [n_perm // n_jobs] * n_jobs
    for i in range(n_perm % n_jobs):
        per_worker[i] += 1

    logging.info("Running %d permutations on %d worker(s)...", n_perm, n_jobs)

    rwr_args = (float(args.restart), float(args.tol), int(args.max_iter))

    # Prefer fork for large shared sparse matrices (Linux). On macOS/Windows, fork may be unavailable.
    # We'll use get_context("fork") when possible; otherwise default context.
    try:
        ctx = get_context("fork")
    except ValueError:
        ctx = get_context()

    with ctx.Pool(
        processes=n_jobs,
        initializer=_init_worker,
        initargs=(graph.P_T, len(graph.nodes), real_scores, buckets, seed_deg_counts, rwr_args),
    ) as pool:
        tasks = [(per_worker[i], int(args.seed), i) for i in range(n_jobs)]
        partial_counts = list(tqdm(pool.imap(_worker_run, tasks), total=n_jobs, desc="Workers"))

    count_ge = np.sum(partial_counts, axis=0).astype(np.int64)

    # 3) Empirical p-values (one-sided: enrichment over null)
    # Add +1 smoothing to avoid p=0.
    pvals = (1.0 + count_ge.astype(np.float64)) / (n_perm + 1.0)

    # Output
    out_df = pd.DataFrame(
        {
            "node": graph.nodes.astype(str),
            "score": real_scores,
            "p_empirical": pvals,
            "degree": graph.degree,
            "is_seed": is_seed.astype(int),
        }
    )
    out_df = out_df.sort_values(["p_empirical", "score"], ascending=[True, False]).reset_index(drop=True)

    logging.info("Writing output: %s", args.out)
    out_df.to_csv(args.out, sep="\t", index=False)

    logging.info("Done.")
    logging.info("Top 10 nodes by smallest p-value:\n%s", out_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
