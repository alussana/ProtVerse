#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt


def qq_plot_pvalues(pvals, ax=None, log10=False, title="Q-Q plot of p-values",
                    point_kwargs=None, line_kwargs=None, out_path="out.pdf"):
    """
    Make a Q-Q plot for p-values against the Uniform(0,1) null.

    Parameters
    ----------
    pvals : array-like
        Iterable of p-values.
    ax : matplotlib.axes.Axes or None
        Axes to draw on. If None, creates a new figure/axes.
    log10 : bool
        If True, plot -log10(expected) vs -log10(observed) (common for p-values).
        If False, plot expected p vs observed p.
    title : str
        Plot title.
    point_kwargs : dict or None
        Keyword args passed to ax.scatter for points.
    line_kwargs : dict or None
        Keyword args passed to ax.plot for the y=x reference line.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes containing the plot.
    """
    p = np.asarray(pvals, dtype=float)

    # Keep only valid p-values in (0, 1]
    p = p[np.isfinite(p)]
    p = p[(p > 0) & (p <= 1)]
    n = p.size
    if n == 0:
        raise ValueError("No valid p-values found. Need values in (0, 1].")

    # Observed: sort ascending
    observed = np.sort(p)

    # Expected under Uniform(0,1): use (i - 0.5)/n to avoid 0 and 1 exactly
    i = np.arange(1, n + 1)
    expected = (i - 0.5) / n

    if log10:
        x = -np.log10(expected)
        y = -np.log10(observed)
        xlabel = r"Expected $-\log_{10}(p)$"
        ylabel = r"Observed $-\log_{10}(p)$"
    else:
        x = expected
        y = observed
        xlabel = "Expected p"
        ylabel = "Observed p"

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    if point_kwargs is None:
        point_kwargs = {}
    if line_kwargs is None:
        line_kwargs = {}

    ax.scatter(x, y, s=12, alpha=0.8, **point_kwargs)

    # Reference line y=x over plotting range
    lo = min(np.min(x), np.min(y))
    hi = max(np.max(x), np.max(y))
    ax.plot([lo, hi], [lo, hi], linewidth=1, **line_kwargs)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.set_aspect("equal", adjustable="box")

    plt.savefig(out_path)

    return ax


# Example:
# pvals = np.random.rand(10000)
# qq_plot_pvalues(pvals, log10=True, title="Q-Q plot (-log10)")
# plt.show()

def main():
    U_p_protverse_tsv = sys.argv[1]
    U_p_string_tsv = sys.argv[2]
    out_prefix = sys.argv[3]

    U_p_protverse = np.loadtxt(U_p_protverse_tsv, delimiter="\t") # (U, p)
    U_p_string = np.loadtxt(U_p_string_tsv, delimiter="\t") # (U, p)
    
    p_protverse = U_p_protverse[:, 1]
    p_string = U_p_string[:, 1]

    qq_plot_pvalues(p_protverse, out_path=f'{out_prefix}protverse_p_QQ.pdf')
    qq_plot_pvalues(p_string, out_path=f'{out_prefix}string_p_QQ.pdf')


if __name__ == "__main__":
    main()
