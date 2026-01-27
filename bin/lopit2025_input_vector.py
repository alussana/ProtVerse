#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

def fraction_covariation_residual_corr(
    df: pd.DataFrame,
    gene_a: str,
    gene_b: str,
    *,
    gene_col: str = "Gene",
    replicate_col: str = "Replicate",
    fraction_col: str = "Fraction",
    value_col: str = "Fold Change",
) -> float:
    """
    Compute correlation between residual Fold Change for two genes after regressing out:
      - gene mean (intercept)
      - replicate fixed effects (categorical)

    Residuals are computed per gene on its own rows, then aligned on (Replicate, Fraction),
    and correlated.
    """

    def _residuals_for_gene(gene_name: str) -> pd.DataFrame:
        # get data for this gene
        g = df.loc[df[gene_col] == gene_name, [replicate_col, fraction_col, value_col]].copy()
        g[replicate_col] = g[replicate_col].astype("category")

        # design matrix: intercept + replicate dummies (drop first to avoid collinearity with intercept)
        X = pd.get_dummies(g[replicate_col], prefix="rep", drop_first=True)
        X.insert(0, "Intercept", 1.0)

        # dependent variable
        y = g[value_col].astype(float).to_numpy()

        # OLS via least squares
        beta, *_ = np.linalg.lstsq(X.to_numpy(dtype=float), y, rcond=None)
        y_hat = X.to_numpy(dtype=float) @ beta
        resid = y - y_hat

        out = g[[replicate_col, fraction_col]].copy()
        out["residual"] = resid
        return out

    ra = _residuals_for_gene(gene_a)
    rb = _residuals_for_gene(gene_b)

    # align residuals on the same (Replicate, Fraction) observations
    merged = ra.merge(
        rb,
        on=[replicate_col, fraction_col],
        how="inner",
        suffixes=(f"_{gene_a}", f"_{gene_b}"),
    )

    x = merged[f"residual_{gene_a}"].to_numpy(dtype=float)
    y = merged[f"residual_{gene_b}"].to_numpy(dtype=float)

    r, p = pearsonr(x, y)

    return (float(r), float(p))


def return_missing_variables(n:int, default=-1):
    return [default for i in range(n)]


def p_to_activation(x, a=3):
    if x == -1:
        raise ValueError("f(x) is undefined at x = -1 (division by zero).")
    return (1 / ( x + 1 ) ** a) - x * (1 / ( x + 1 ) ** a)


def map_sign(x):
    return (x + 1) // 2


def main():
    """
    gene_pairs_file = 'input/gene_pairs.tsv'
    lopit2025_table_file = 'input/lopit2025_table.tsv'
    """
    gene_pairs_file = sys.argv[1]
    lopit2025_table_file = sys.argv[2]

    gene_pairs = pd.read_csv(gene_pairs_file, sep='\t', header=None)
    gene_pairs.columns = ['label', 'reactome', 'gene_1', 'gene_2']

    lopit2025 = pd.read_csv(lopit2025_table_file, sep='\t', header=None)
    lopit2025.columns = ["Gene", "Replicate", "Fraction", "Fold Change"]

    pearson = []
    pearson_sign = []
    pearson_pval = []
    data_flag = []

    for i in range(len(gene_pairs)):
        data_available = 1
        gene_1 = gene_pairs.loc[i,]['gene_1']
        gene_2 = gene_pairs.loc[i,]['gene_2']
        # if either gene is not found, return missing features
        if gene_1 not in list(lopit2025["Gene"]) or gene_2 not in list(lopit2025["Gene"]):
            p, p_sign, p_pval = return_missing_variables(3)
            data_available = return_missing_variables(1, 0).pop()
        else:
            p, p_pval = fraction_covariation_residual_corr(
                lopit2025,
                gene_1,
                gene_2,
                gene_col="Gene",
                replicate_col="Replicate",
                fraction_col="Fraction",
                value_col="Fold Change",
            )
            if np.isnan(p):
                p, p_sign, p_pval = return_missing_variables(3)
                data_available = return_missing_variables(1, 0).pop()
            else:
                p_sign = round(map_sign(p), 3)
                p_pval = round(p_to_activation(p_pval), 3)
                p = round(p**2, 3)
        pearson.append(p)
        pearson_sign.append(p_sign)
        pearson_pval.append(p_pval)
        data_flag.append(data_available)

    gene_pairs[r'Cotranslocation residuals $r^2$'] = pearson
    gene_pairs[r'Cotranslocation residuals $\mathrm{sign}\,(r)$'] = pearson_sign
    gene_pairs[r'Cotranslocation residuals $r^2$ Sig.'] = pearson_pval
    gene_pairs['Cotranslocation Data'] = data_flag
    
    print(gene_pairs.to_csv(sep='\t', index=False))


if __name__ == '__main__':
    main()