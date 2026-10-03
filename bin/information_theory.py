
#!/usr/bin/env python3
"""
Utilities for entropy-like dependence measures on (multi-dimensional) boolean
or categorical matrices X and Y (shape: n_samples x n_dims).
"""

from __future__ import annotations

import numpy as np
from collections import Counter
from math import log2, sqrt
from typing import Callable, Tuple, Optional

# ---------------------------------------------------------------------------
# Basic entropy utilities
# ---------------------------------------------------------------------------

def _entropy_from_probs(prob_dist) -> float:
    """Compute Shannon entropy (base-2) from an iterable of probabilities."""
    return -sum(p * log2(p) for p in prob_dist if p > 0)


def _column_entropy(col: np.ndarray) -> float:
    """Entropy of a single discrete column (1-D array)."""
    values, counts = np.unique(col, return_counts=True)
    ps = counts / counts.sum()
    return _entropy_from_probs(ps)


def marginal_entropy(X: np.ndarray) -> float:
    """
    Sum of entropies of each column in matrix X (n_samples x n_dims).

    This is *not* the joint entropy of X; it's \sum_j H(X_j).
    """
    if X.ndim != 2:
        raise ValueError("X must be a 2D array.")
    return float(sum(_column_entropy(X[:, j]) for j in range(X.shape[1])))


# ---------------------------------------------------------------------------
# Joint / vector entropy
# ---------------------------------------------------------------------------

def _rows_to_tuples(M: np.ndarray) -> list[tuple]:
    """Convert each row of M into a tuple for counting joint states."""
    if M.ndim != 2:
        raise ValueError("Input must be 2D.")
    # Ensure 2D even when M has a single column
    return [tuple(row.tolist()) for row in M]


def vector_entropy(M: np.ndarray) -> float:
    """
    Joint entropy of all columns in M, treating each row as one joint symbol.
    """
    rows = _rows_to_tuples(M)
    counts = Counter(rows)
    total = len(rows)
    ps = [c / total for c in counts.values()]
    return float(_entropy_from_probs(ps))


def joint_entropy(X: np.ndarray, Y: np.ndarray) -> float:
    """Joint entropy H([X,Y]) where rows are concatenated across columns."""
    if X.shape[0] != Y.shape[0]:
        raise ValueError("X and Y must have the same number of rows (samples).")
    XY = np.concatenate([X, Y], axis=1)
    return vector_entropy(XY)


# ---------------------------------------------------------------------------
# Total correlation (a.k.a. multivariate mutual information / multi-information)
# ---------------------------------------------------------------------------

def total_correlation(X: np.ndarray, Y: np.ndarray) -> float:
    """
    Total correlation across all columns of X and Y:
        TC(X,Y) = sum_j H(X_j) + sum_k H(Y_k) - H([X,Y]).

    Note: If you define Z = [X,Y] (column-wise concatenation), then this equals
    \sum_r H(Z_r) - H(Z), which is the standard Watanabe total correlation.
    """
    HX_sum = marginal_entropy(X)
    HY_sum = marginal_entropy(Y)
    HXY = joint_entropy(X, Y)
    return float(HX_sum + HY_sum - HXY)


def normalized_total_correlation(
    X: np.ndarray,
    Y: np.ndarray,
    normalization: str = "avg_marginal_sums",
) -> float:
    """
    Normalized total correlation. We put TC(X,Y) in the numerator and choose
    a denominator from the options below.

    Args
    ----
    normalization:
      - "avg_marginal_sums" (default): 0.5*(sum_j H(X_j) + sum_k H(Y_k))
      - "geometric_marginal_sums": sqrt( sum_j H(X_j) * sum_k H(Y_k) )
      - "max_marginal_sums": max( sum_j H(X_j), sum_k H(Y_k) )
      - "min_marginal_sums": min( sum_j H(X_j), sum_k H(Y_k) )
      - "joint": H([X,Y])  (less common, but sometimes useful)

    Returns
    -------
    A float in [0, 1] for the first four normalizations. The "joint"
    option can exceed 1 on highly dependent data.
    """
    HX_sum = marginal_entropy(X)
    HY_sum = marginal_entropy(Y)
    HXY = joint_entropy(X, Y)
    TC = HX_sum + HY_sum - HXY

    if normalization == "avg_marginal_sums":
        denom = 0.5 * (HX_sum + HY_sum)
    elif normalization == "geometric_marginal_sums":
        denom = sqrt(HX_sum * HY_sum) if HX_sum > 0 and HY_sum > 0 else 0.0
    elif normalization == "max_marginal_sums":
        denom = max(HX_sum, HY_sum)
    elif normalization == "min_marginal_sums":
        denom = min(HX_sum, HY_sum)
    elif normalization == "joint":
        denom = HXY
    else:
        raise ValueError(f"Unknown normalization: {normalization!r}")

    if denom == 0.0:
        return 0.0
    return float(TC / denom)


# ---------------------------------------------------------------------------
# Mutual information between multivariate vectors X and Y
# ---------------------------------------------------------------------------

def mutual_information(X: np.ndarray, Y: np.ndarray) -> float:
    """
    Mutual information between *vectors* X and Y:
        I(X;Y) = H_vec(X) + H_vec(Y) - H([X,Y])
    where H_vec(.) is the *joint* entropy of the corresponding multivariate
    vector (all columns together).
    """
    HX_vec = vector_entropy(X)
    HY_vec = vector_entropy(Y)
    HXY = joint_entropy(X, Y)
    return float(HX_vec + HY_vec - HXY)


def normalized_mutual_information(
    X: np.ndarray,
    Y: np.ndarray,
    normalization: str = "avg",
) -> float:
    """
    Normalized mutual information between *vectors* X and Y.

    normalization in {"avg", "geometric", "min", "max"}:
        - "avg":       I / (0.5*(H_vec(X) + H_vec(Y)))
        - "geometric": I / sqrt(H_vec(X) * H_vec(Y))
        - "min":       I / min(H_vec(X), H_vec(Y))
        - "max":       I / max(H_vec(X), H_vec(Y))

    Returns 0.0 if the denominator is 0.
    """
    HX = vector_entropy(X)
    HY = vector_entropy(Y)
    I = mutual_information(X, Y)

    if normalization == "geometric":
        denom = sqrt(HX * HY) if HX > 0 and HY > 0 else 0.0
    elif normalization == "avg":
        denom = 0.5 * (HX + HY)
    elif normalization == "min":
        denom = min(HX, HY)
    elif normalization == "max":
        denom = max(HX, HY)
    else:
        raise ValueError(f"Unknown normalization: {normalization!r}")

    if denom == 0.0:
        return 0.0
    return float(I / denom)


# ---------------------------------------------------------------------------
# Permutation test
# ---------------------------------------------------------------------------

def permutation_test_total_correlation(
    X: np.ndarray,
    Y: np.ndarray,
    stat_fn: Callable[[np.ndarray, np.ndarray], float] = total_correlation,
    n_permutations: int = 1000,
    seed: Optional[int] = None,
) -> Tuple[float, float, np.ndarray]:
    """
    Permutation test for dependence using a TC-like statistic.

    We permute the rows of Y relative to X, recompute the statistic, and form
    an empirical null. p-value is computed as
        (count(stat_perm >= observed) + 1) / (n_permutations + 1)
    which is a conservative upper-bound.
    """
    rng = np.random.default_rng(seed)
    observed = float(stat_fn(X, Y))
    null_stats = np.empty(n_permutations, dtype=float)
    for b in range(n_permutations):
        perm = rng.permutation(X.shape[0])
        null_stats[b] = float(stat_fn(X, Y[perm]))
    pval = (np.sum(null_stats >= observed) + 1.0) / (n_permutations + 1.0)
    return observed, float(pval), null_stats


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n_samples, d_X, d_Y = 5000, 3, 4

    # Sparse-ish X
    X = (rng.random((n_samples, d_X)) < 0.15).astype(int)

    # Y depends on some of X with a bit of noise
    Y = np.zeros((n_samples, d_Y), dtype=int)
    noise = (rng.random((n_samples, d_Y)) < 0.05).astype(int)
    for i in range(d_Y):
        Y[:, i] = X[:, i % d_X] | noise[:, i]

    tc_val = total_correlation(X, Y)
    n_tc_val = normalized_total_correlation(X, Y)
    mi_val = mutual_information(X, Y)
    nmi_avg = normalized_mutual_information(X, Y, "avg")
    nmi_geo = normalized_mutual_information(X, Y, "geometric")

    print(f"TC(X,Y): {tc_val:.4f}")
    print(f"Normalized TC (avg_marginal_sums): {n_tc_val:.4f}")
    print(f"MI(X;Y): {mi_val:.4f}")
    print(f"NMI avg: {nmi_avg:.4f} | NMI geom: {nmi_geo:.4f}")

    obs, p, null = permutation_test_total_correlation(
        X, Y, stat_fn=total_correlation, n_permutations=200, seed=42
    )
    print(f"Observed TC: {obs:.4f}; empirical p-value (>=): {p:.4f}")


# =============================================================================
# Continuous (k-NN) estimators (Kraskov / Kozachenko–Leonenko)
# =============================================================================

import math

_EULER_GAMMA = 0.5772156649015328606

def _try_import_scipy_ckdtree():
    try:
        from scipy.spatial import cKDTree  # type: ignore
        return cKDTree
    except Exception:
        return None

def _harmonic_numbers_upto(n: int):
    """Return array H where H[m] = sum_{r=1}^m 1/r, for m=0..n."""
    import numpy as np
    H = np.zeros(n + 1, dtype=float)
    if n >= 1:
        H[1:] = np.cumsum(1.0 / np.arange(1, n + 1, dtype=float))
    return H

def _psi_integers(vals, H):
    """Digamma for positive integers using harmonic numbers: psi(m) = H_{m-1} - gamma"""
    import numpy as np
    vals = np.asarray(vals, dtype=int)
    if (vals < 1).any():
        raise ValueError("psi_integers requires vals >= 1")
    return H[vals - 1] - _EULER_GAMMA

def _pairwise_chebyshev(X, Y=None):
    """Compute pairwise Chebyshev (infinity norm) distances."""
    import numpy as np
    X = np.asarray(X, float)
    Y = X if Y is None else np.asarray(Y, float)
    diffs = X[:, None, :] - Y[None, :, :]
    return (np.abs(diffs)).max(axis=2)

def _kth_neighbor_radius_supnorm(Z, k: int):
    """Return Chebyshev distance to k-th nearest neighbor for each point in Z."""
    import numpy as np
    n = Z.shape[0]
    cKDTree = _try_import_scipy_ckdtree()
    if cKDTree is not None and n >= 2:
        tree = cKDTree(Z)
        dists, _ = tree.query(Z, k=k+1, p=float('inf'))
        if k == 1:
            eps = dists
        else:
            eps = dists[:, -1]
        return eps, "ckdtree"
    # Fallback
    D = _pairwise_chebyshev(Z)
    np.fill_diagonal(D, np.inf)
    eps = np.partition(D, kth=k-1, axis=1)[:, k-1]
    return eps, "bruteforce"

def _counts_within_radius_supnorm(X, eps):
    """Counts of neighbors j != i with Chebyshev distance <= eps[i] - tiny."""
    import numpy as np
    n = X.shape[0]
    tiny = 1e-15
    cKDTree = _try_import_scipy_ckdtree()
    if cKDTree is not None and n >= 2:
        tree = cKDTree(X)
        idx_lists = tree.query_ball_point(X, r=(eps - tiny), p=float('inf'))
        return np.array([len(idx)-1 for idx in idx_lists], dtype=int)
    D = _pairwise_chebyshev(X)
    np.fill_diagonal(D, np.inf)
    return (D <= (eps[:, None] - tiny)).sum(axis=1).astype(int)

def mutual_information_knn(X, Y, k: int = 5, add_noise: float | None = None, base: float = 2.0) -> float:
    """
    KSG (algorithm 1) estimator of I(X;Y) with Chebyshev norm in joint space.

    Parameters
    ----------
    X, Y : arrays (n, d_x), (n, d_y)
    k : neighbors
    add_noise : if provided, add Gaussian jitter with std = add_noise * std(col)
    base : log base (2 bits, e nats)
    """
    import numpy as np
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    if X.ndim != 2 or Y.ndim != 2 or X.shape[0] != Y.shape[0]:
        raise ValueError("X and Y must be 2D with the same number of rows.")
    n = X.shape[0]
    if not (1 <= k < n):
        raise ValueError("k must satisfy 1 <= k < n_samples.")
    if add_noise:
        rng = np.random.default_rng(0)
        X = X + rng.normal(scale=np.std(X, axis=0, ddof=1, keepdims=True) * add_noise, size=X.shape)
        Y = Y + rng.normal(scale=np.std(Y, axis=0, ddof=1, keepdims=True) * add_noise, size=Y.shape)

    Z = np.concatenate([X, Y], axis=1)
    eps, _ = _kth_neighbor_radius_supnorm(Z, k=k)
    nx = _counts_within_radius_supnorm(X, eps)
    ny = _counts_within_radius_supnorm(Y, eps)

    H = _harmonic_numbers_upto(n)
    psi_k = float(_psi_integers([k], H)[0])
    psi_n = float(_psi_integers([n], H)[0])
    psi_nx1 = _psi_integers(nx + 1, H)
    psi_ny1 = _psi_integers(ny + 1, H)

    I = psi_k + psi_n - float((psi_nx1 + psi_ny1).mean())
    return float(I / math.log(base))

def differential_entropy_knn(X, k: int = 5, metric: str = "euclidean", add_noise: float | None = None, base: float = 2.0) -> float:
    """
    Kozachenko–Leonenko k-NN estimator of differential entropy H(X).

    h = psi(n) - psi(k) + ln(V_d) + d * E[ ln eps_i ],
    where V_d is the unit-ball volume in chosen norm:
      - euclidean:  V_d = pi^(d/2) / Gamma(d/2 + 1)
      - chebyshev:  V_d = 2^d
    """
    import numpy as np, math
    X = np.asarray(X, float)
    if X.ndim != 2:
        raise ValueError("X must be 2D (n_samples, n_dims).")
    n, d = X.shape
    if not (1 <= k < n):
        raise ValueError("k must satisfy 1 <= k < n_samples.")
    if add_noise:
        rng = np.random.default_rng(0)
        X = X + rng.normal(scale=np.std(X, axis=0, ddof=1, keepdims=True) * add_noise, size=X.shape)

    if metric in ("chebyshev", "sup", "infinity"):
        eps, _ = _kth_neighbor_radius_supnorm(X, k=k)
        ln_Vd = d * math.log(2.0)
    elif metric == "euclidean":
        try:
            from scipy.spatial import cKDTree  # type: ignore
            tree = cKDTree(X)
            dists, _ = tree.query(X, k=k+1, p=2)
            eps = dists if k == 1 else dists[:, -1]
        except Exception:
            diffs = X[:, None, :] - X[None, :, :]
            D = (diffs * diffs).sum(axis=2)**0.5
            np.fill_diagonal(D, np.inf)
            eps = np.partition(D, kth=k-1, axis=1)[:, k-1]
        ln_Vd = (d/2.0) * math.log(math.pi) - math.lgamma(d/2.0 + 1.0)
    else:
        raise ValueError("metric must be 'euclidean' or 'chebyshev'")

    if (eps <= 0).any():
        if add_noise is None:
            return differential_entropy_knn(X, k=k, metric=metric, add_noise=1e-10, base=base)
        eps = np.maximum(eps, 1e-300)

    H = _harmonic_numbers_upto(n)
    psi_n = float(_psi_integers([n], H)[0])
    psi_k = float(_psi_integers([k], H)[0])
    h = psi_n - psi_k + ln_Vd + d * float(np.log(eps).mean())
    return float(h / math.log(base))

def normalized_mutual_information_knn(
    X, Y, k: int = 5, normalization: str = "avg", base: float = 2.0, add_noise: float | None = None, entropy_metric: str = "euclidean"
) -> float:
    """
    Normalized MI for continuous data using KSG (for I) and KL (for entropies).

    normalization:
      - "avg": 2*I / (H(X)+H(Y))
      - "min": I / min(H(X), H(Y))
      - "max": I / max(H(X), H(Y))
      - "geometric": I / sqrt(H(X) * H(Y))
      - "joint": I / H([X,Y])

    Note: differential entropies can be negative; these ratios may exceed 1.
    """
    import numpy as np, math
    I = mutual_information_knn(X, Y, k=k, add_noise=add_noise, base=base)
    HX = differential_entropy_knn(X, k=k, metric=entropy_metric, add_noise=add_noise, base=base)
    HY = differential_entropy_knn(Y, k=k, metric=entropy_metric, add_noise=add_noise, base=base)

    if normalization == "avg":
        denom = 0.5 * (HX + HY)
        return 0.0 if denom == 0 else float(2.0 * I / (HX + HY))
    elif normalization == "min":
        denom = min(HX, HY)
    elif normalization == "max":
        denom = max(HX, HY)
    elif normalization == "geometric":
        denom = math.sqrt(HX * HY) if HX > 0 and HY > 0 else 0.0
    elif normalization == "joint":
        HXY = differential_entropy_knn(np.concatenate([X, Y], axis=1), k=k, metric=entropy_metric, add_noise=add_noise, base=base)
        denom = HXY
    else:
        raise ValueError(f"Unknown normalization: {normalization!r}")
    return 0.0 if denom == 0 else float(I / denom)

def total_correlation_knn(X, Y=None, k: int = 5, base: float = 2.0, add_noise: float | None = None, entropy_metric: str = "euclidean") -> float:
    """
    Total correlation for continuous data via KL entropy estimates.

    If Y is None:
        TC = sum_j H(X_j) - H(X)
    Else:
        TC = sum_j H(X_j) + sum_k H(Y_k) - H([X,Y])
    """
    import numpy as np
    X = np.asarray(X, float)
    if Y is not None:
        Y = np.asarray(Y, float)
        if X.shape[0] != Y.shape[0]:
            raise ValueError("X and Y must have the same number of rows.")
        M = np.concatenate([X, Y], axis=1)
    else:
        M = X

    HX_sum = 0.0
    for j in range(X.shape[1]):
        HX_sum += differential_entropy_knn(X[:, [j]], k=k, metric=entropy_metric, add_noise=add_noise, base=base)

    HY_sum = 0.0
    if Y is not None:
        for j in range(Y.shape[1]):
            HY_sum += differential_entropy_knn(Y[:, [j]], k=k, metric=entropy_metric, add_noise=add_noise, base=base)

    H_joint = differential_entropy_knn(M, k=k, metric=entropy_metric, add_noise=add_noise, base=base)
    return float(HX_sum + HY_sum - H_joint)

def normalized_total_correlation_knn(
    X, Y, k: int = 5, normalization: str = "avg_marginal_sums", base: float = 2.0, add_noise: float | None = None, entropy_metric: str = "euclidean"
) -> float:
    """
    Normalized total correlation for continuous data using KL entropy estimates.

      - "avg_marginal_sums" (default): TC / (0.5*(sum H(X_j) + sum H(Y_k)))
      - "geometric_marginal_sums":    TC / sqrt( sum H(X_j) * sum H(Y_k) )
      - "max_marginal_sums":          TC / max( sum H(X_j), sum H(Y_k) )
      - "min_marginal_sums":          TC / min( sum H(X_j), sum H(Y_k) )
      - "joint":                      TC / H([X,Y])
    """
    import numpy as np, math
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    if X.shape[0] != Y.shape[0]:
        raise ValueError("X and Y must have the same number of rows.")
    HX_sum = sum(differential_entropy_knn(X[:, [j]], k=k, metric=entropy_metric, add_noise=add_noise, base=base) for j in range(X.shape[1]))
    HY_sum = sum(differential_entropy_knn(Y[:, [j]], k=k, metric=entropy_metric, add_noise=add_noise, base=base) for j in range(Y.shape[1]))
    HXY = differential_entropy_knn(np.concatenate([X, Y], axis=1), k=k, metric=entropy_metric, add_noise=add_noise, base=base)
    TC = HX_sum + HY_sum - HXY

    if normalization == "avg_marginal_sums":
        denom = 0.5 * (HX_sum + HY_sum)
    elif normalization == "geometric_marginal_sums":
        denom = math.sqrt(HX_sum * HY_sum) if HX_sum > 0 and HY_sum > 0 else 0.0
    elif normalization == "max_marginal_sums":
        denom = max(HX_sum, HY_sum)
    elif normalization == "min_marginal_sums":
        denom = min(HX_sum, HY_sum)
    elif normalization == "joint":
        denom = HXY
    else:
        raise ValueError(f"Unknown normalization: {normalization!r}")

    return 0.0 if denom == 0 else float(TC / denom)


# =============================================================================
# Bounded dependence measures for comparability across pairs
# =============================================================================

def information_correlation_from_I(I: float, base: float = 2.0) -> float:
    """
    Linfoot's information correlation: rho_I = sqrt(1 - base**(-2 I)).
    - If base=2, I is in bits; if base=e, I is in nats.
    - Monotone in I, bounded in [0,1), invariant to invertible reparameterizations of X or Y.
    """
    import math
    if I < 0:
        I = 0.0  # numerical safeguard; true MI is nonnegative
    return math.sqrt(1.0 - (base ** (-2.0 * I)))

def information_correlation_knn(X, Y, k: int = 3, base: float = 2.0, add_noise: float | None = None) -> float:
    """
    Convenience: compute k-NN MI then convert to bounded [0,1) via Linfoot transform.
    """
    I = mutual_information_knn(X, Y, k=k, add_noise=add_noise, base=base)
    return information_correlation_from_I(I, base=base)

def distance_correlation(X, Y) -> float:
    """
    Distance correlation R in [0,1] using the (biased) Székely–Rizzo estimator.

    Works for 1D or multi-D arrays; for 1D, X and Y can be shape (n,) or (n,1).
    This equals 0 iff independence at the population level (R=0 <=> independent).
    """
    import numpy as np, math
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)

    if X.ndim == 1:
        X = X[:, None]
    if Y.ndim == 1:
        Y = Y[:, None]
    if X.shape[0] != Y.shape[0]:
        raise ValueError("X and Y must have the same number of rows.")

    # pairwise Euclidean distances
    DX = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(axis=2))
    DY = np.sqrt(((Y[:, None, :] - Y[None, :, :]) ** 2).sum(axis=2))

    # double-centering
    DX_row = DX.mean(axis=1, keepdims=True)
    DX_col = DX.mean(axis=0, keepdims=True)
    DX_all = DX.mean()
    AX = DX - DX_row - DX_col + DX_all

    DY_row = DY.mean(axis=1, keepdims=True)
    DY_col = DY.mean(axis=0, keepdims=True)
    DY_all = DY.mean()
    AY = DY - DY_row - DY_col + DY_all

    V2XY = (AX * AY).mean()
    V2X = (AX * AX).mean()
    V2Y = (AY * AY).mean()

    # numerical guards
    V2XY = max(0.0, float(V2XY))
    if V2X <= 0 or V2Y <= 0:
        return 0.0
    R2 = V2XY / math.sqrt(V2X * V2Y)
    R2 = min(1.0, max(0.0, R2))
    return math.sqrt(R2)
