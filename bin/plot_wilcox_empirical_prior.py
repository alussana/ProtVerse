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


def qq_plot_pvalues(pvals, pvals_prior, pvals_1, pvals_1_prior, pvals_2, pvals_2_prior, ax=None, log10=False, title=r"Q-Q plot of Mann-Whitney U $p$ values",
                    point_kwargs=None, line_kwargs=None, out_path="out.pdf"):
    
    p = np.asarray(pvals, dtype=float)
    p_1 = np.asarray(pvals_1, dtype=float)
    p_2 = np.asarray(pvals_2, dtype=float)

    p_prior = np.asarray(pvals_prior, dtype=float)
    p_1_prior = np.asarray(pvals_1_prior, dtype=float)
    p_2_prior = np.asarray(pvals_2_prior, dtype=float)

    n = p.size
    n_1 = p_1.size
    n_2 = p_2.size

    # Observed: sort ascending
    observed = np.sort(p)
    observed_1 = np.sort(p_1)
    observed_2 = np.sort(p_2)

    # Expected under 1lof
    expected = np.quantile(p_prior, np.linspace(0, 1, n))
    expected_1 = np.quantile(p_1_prior, np.linspace(0, 1, n_1))
    expected_2 = np.quantile(p_2_prior, np.linspace(0, 1, n_2))

    if log10:
        x = -np.log10(expected)
        x_1 = -np.log10(expected_1)
        x_2 = -np.log10(expected_2)
        y = -np.log10(observed)
        y_1 = -np.log10(observed_1)
        y_2 = -np.log10(observed_2)
        xlabel = r"Observed $-\log_{10}(p)$, seeds with LOF allele dosage $= 1$"
        ylabel = r"Observed $-\log_{10}(p)$, seeds with LOF allele dosage $\geq 2$"
    else:
        x = expected
        x_1 = expected_1
        x_2 = expected_2
        y = observed
        y_1 = observed_1
        y_2 = observed_2
        xlabel = r"Observed $p$, seeds with LOF allele dosage $= 1$"
        ylabel = r"Observed $p$, seeds with LOF allele dosage $\geq 2$"

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
    U_p_protverse_1lof_tsv = sys.argv[1]
    U_p_protverse_2lof_tsv = sys.argv[2]
    U_p_string_1lof_tsv = sys.argv[3]
    U_p_string_2lof_tsv = sys.argv[4]
    U_p_reactome_1lof_tsv = sys.argv[5]
    U_p_reactome_2lof_tsv = sys.argv[6]
    out_prefix = sys.argv[7]

    U_p_protverse_1lof = np.loadtxt(U_p_protverse_1lof_tsv, delimiter="\t") # (U, p)
    U_p_string_1lof = np.loadtxt(U_p_string_1lof_tsv, delimiter="\t") # (U, p)
    U_p_reactome_1lof = np.loadtxt(U_p_reactome_1lof_tsv, delimiter="\t") # (U, p)

    U_p_protverse_2lof = np.loadtxt(U_p_protverse_2lof_tsv, delimiter="\t") # (U, p)
    U_p_string_2lof = np.loadtxt(U_p_string_2lof_tsv, delimiter="\t") # (U, p)
    U_p_reactome_2lof = np.loadtxt(U_p_reactome_2lof_tsv, delimiter="\t") # (U, p)
    
    p_protverse_1lof = U_p_protverse_1lof[:, 1]
    p_string_1lof = U_p_string_1lof[:, 1]
    p_reactome_1lof = U_p_reactome_1lof[:, 1]

    p_protverse_2lof = U_p_protverse_2lof[:, 1]
    p_string_2lof = U_p_string_2lof[:, 1]
    p_reactome_2lof = U_p_reactome_2lof[:, 1]
    
    qq_plot_pvalues(p_protverse_2lof, p_protverse_1lof, p_string_2lof, p_string_1lof, p_reactome_2lof, p_reactome_1lof, out_path=f'{out_prefix}p_QQ_empirical_prior.pdf')


if __name__ == "__main__":
    main()
