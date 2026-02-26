#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
from typing import Literal, Tuple


def mann_whitney_u_test(
    x: np.ndarray,
    y: np.ndarray,
    alternative: Literal["two-sided", "less", "greater"] = "greater",
    method: Literal["auto", "asymptotic", "exact"] = "auto",
    continuity: bool = True,
) -> Tuple[float, float]:
    """
    Mann–Whitney U test (a.k.a. Wilcoxon rank-sum test) for two independent samples.

    Parameters
    ----------
    x, y : np.ndarray
        1D arrays of sample values. Lengths can differ.
    alternative : {"two-sided", "less", "greater"}
        - "two-sided": distributions differ
        - "less":      x tends to be smaller than y
        - "greater":   x tends to be larger than y
    method : {"auto", "asymptotic", "exact"}
        - "exact": exact p-value (only valid when there are no ties; practical for small samples)
        - "asymptotic": normal approximation with tie correction
        - "auto": uses exact when feasible (no ties and small), else asymptotic
    continuity : bool
        Apply continuity correction for asymptotic method.

    Returns
    -------
    U : float
        Mann–Whitney U statistic (for sample x).
    p_value : float
        P-value corresponding to the chosen alternative/method.

    Notes
    -----
    - Requires SciPy. If you truly need a SciPy-free implementation, tell me and I’ll provide one
      (asymptotic with tie correction).
    """
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()

    # Drop NaNs (common in real data); change this if you'd rather error on NaNs.
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    if x.size == 0 or y.size == 0:
        raise ValueError("Both samples must contain at least one non-NaN observation.")

    try:
        from scipy.stats import mannwhitneyu
    except ImportError as e:
        raise ImportError("SciPy is required for this implementation. Install with `pip install scipy`.") from e

    # SciPy supports method and alternative; 'use_continuity' applies for asymptotic.
    U, p = mannwhitneyu(
        x, y,
        alternative=alternative,
        method=method,
        use_continuity=continuity
    )
    return float(U), float(p)


def main():
    dep_tsv = sys.argv[1]
    rwr_tsv_gz = sys.argv[2]

    dep_df = pd.read_csv(dep_tsv, index_col=0).transpose()
    rwr_df = pd.read_csv(rwr_tsv_gz, sep='\t', index_col=0)

    # remove seeds
    rwr_df = rwr_df.loc[rwr_df["is_seed"]==0,]

    # -log10(p)
    neg_log_p = -np.log10(rwr_df['p_empirical'])

    # take nodes intersection
    df = pd.merge(neg_log_p, dep_df, right_index=True, left_index=True).dropna()
    df.columns=["Propagation Strength", "Dependency Probability"]

    # perform Wilcoxon rank-sum (Mann–Whitney U) test 
    propagated = df.loc[df['Propagation Strength']>2, "Dependency Probability"].values
    background = df.loc[df['Propagation Strength']<=2, "Dependency Probability"].values

    if len(propagated) != 0:
        U, p = mann_whitney_u_test(propagated, background, alternative="greater", method="auto")       
        print(f'{U}\t{p}')

    # perform Beta regression
    # [...]


if __name__ == "__main__":
    main()
