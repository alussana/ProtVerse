#!/usr/bin/env python

import pandas as pd
import sys


def main():
    df = pd.read_excel(sys.argv[1], sheet_name="Supp.Data3", header=2)
    meta_df = pd.read_excel(sys.argv[1], sheet_name="Sample_info", header=11)

    # remove unneeded columns
    to_drop_list = [
        "Gene",
        "Protein name",
        "Protein sequence coverage (%)",
        "Number of unique assigned peptides",
        "Subcellular markers",
        "Classification",
        "Classification probability",
        "Differential localisation score",
        "Outlier score",
        "Subcellular markers.1",
        "Classification.1",
        "Classification probability.1",
        "Differential localisation score.1",
        "Outlier score.1",
    ]
    df.drop(columns=to_drop_list, inplace=True)

    # unroll datapoints
    df = df.melt("Accession")
    df.columns = ["Accession", "TMT channel", "value"]

    # add replicate column
    df["Replicate"] = df["TMT channel"].astype(str).str.extract(r"\.([^\.]+)$")

    # add metadata: condition column
    df = df.merge(
        meta_df[["TMT channel", "Condition", "Spin speed"]],
        on="TMT channel",
        how="left",
    )

    # drop TMT channel column
    df.drop(columns=["TMT channel"], inplace=True)

    # compute protein fold changes between IR-treated and Control for each replicate and spin speed
    expr_df = df.loc[df["Condition"] == "IR-treated"].copy()
    ctrl_df = df.loc[df["Condition"] == "Control"].copy()
    expr_df.drop(columns=["Condition"], inplace=True)
    ctrl_df.drop(columns=["Condition"], inplace=True)
    expr_df.rename(columns={"value": "expr_value"}, inplace=True)
    ctrl_df.rename(columns={"value": "ctrl_value"}, inplace=True)
    merged_df = expr_df.merge(
        ctrl_df, on=["Accession", "Replicate", "Spin speed"], how="inner"
    )
    # check correlation between expr and ctrl values (same replicate and spin speed)
    # > from scipy.stats import pearsonr
    # > pearsonr(merged_df["expr_value"].values, merged_df["ctrl_value"].values)
    # PearsonRResult(statistic=0.9354219430500028, pvalue=0.0)
    merged_df["Fold Change"] = merged_df["expr_value"] / merged_df["ctrl_value"]
    merged_df.drop(columns=["expr_value", "ctrl_value"], inplace=True)

    merged_df.to_csv(sys.argv[2], index=False, sep="\t")


if __name__ == "__main__":
    main()
