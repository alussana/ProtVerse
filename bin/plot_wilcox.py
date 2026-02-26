#!/usr/bin/env python3

import sys
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

red_shade = "#d62728"
grey_shade = "#7f7f7f"

custom_palette = sns.blend_palette([red_shade, grey_shade], n_colors=3)

sns.set_palette(custom_palette)


def qq_plot_pvalues(pvals, pvals_1, pvals_2, ax=None, log10=False, title="Q-Q plot of p-values",
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
    p_1 = np.asarray(pvals_1, dtype=float)
    p_2 = np.asarray(pvals_2, dtype=float)

    n = p.size
    n_1 = p_1.size
    n_2 = p_2.size

    # Observed: sort ascending
    observed = np.sort(p)
    observed_1 = np.sort(p_1)
    observed_2 = np.sort(p_2)

    # Expected under Uniform(0,1): use (i - 0.5)/n to avoid 0 and 1 exactly
    i = np.arange(1, n + 1)
    i_1 = np.arange(1, n_1 + 1)
    i_2 = np.arange(1, n_2 + 1)
    expected = (i - 0.5) / n
    expected_1 = (i_1 - 0.5) / n_1
    expected_2 = (i_2 - 0.5) / n_2

    if log10:
        x = -np.log10(expected)
        x_1 = -np.log10(expected_1)
        x_2 = -np.log10(expected_2)
        y = -np.log10(observed)
        y_1 = -np.log10(observed_1)
        y_2 = -np.log10(observed_2)
        xlabel = r"Expected $-\log_{10}(p)$"
        ylabel = r"Observed $-\log_{10}(p)$"
    else:
        x = expected
        x_1 = expected_1
        x_2 = expected_2
        y = observed
        y_1 = observed_1
        y_2 = observed_2
        xlabel = "Expected p"
        ylabel = "Observed p"

    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 3))

    if point_kwargs is None:
        point_kwargs = {}
    if line_kwargs is None:
        line_kwargs = {}

    # Reference line y=x over plotting range
    lo = min(np.min(x), np.min(y))
    hi = max(np.max(x), np.max(y))
    ax.plot([lo, hi], [lo, hi], linewidth=0.8, color='black', **line_kwargs)

    # QQ points
    ax.scatter(x, y, alpha=0.5, s=1, label="ProtVerse", **point_kwargs)
    ax.scatter(x_1, y_1, alpha=0.5, s=1, label="STRING", **point_kwargs)
    ax.scatter(x_2, y_2, alpha=0.5, s=1, label="Reactome", **point_kwargs)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.legend(frameon=False)
    #ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.set_aspect("equal", adjustable="box")
    
    sns.despine()
    
    fig.tight_layout()

    plt.savefig(out_path)

    return ax


# Example:
# pvals = np.random.rand(10000)
# qq_plot_pvalues(pvals, log10=True, title="Q-Q plot (-log10)")
# plt.show()

def main():
    U_p_protverse_tsv = sys.argv[1]
    U_p_string_tsv = sys.argv[2]
    U_p_reactome_tsv = sys.argv[3]
    out_prefix = sys.argv[4]

    U_p_protverse = np.loadtxt(U_p_protverse_tsv, delimiter="\t") # (U, p)
    U_p_string = np.loadtxt(U_p_string_tsv, delimiter="\t") # (U, p)
    U_p_reactome = np.loadtxt(U_p_reactome_tsv, delimiter="\t") # (U, p)
    
    p_protverse = U_p_protverse[:, 1]
    p_string = U_p_string[:, 1]
    p_reactome = U_p_reactome[:, 1]
    
    qq_plot_pvalues(p_protverse, p_string, p_reactome, out_path=f'{out_prefix}p_QQ.pdf')


if __name__ == "__main__":
    main()
