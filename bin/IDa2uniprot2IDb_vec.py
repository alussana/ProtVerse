#!/usr/bin/env python3

import sys
import pandas as pd


def main():
    map_file = sys.argv[1]
    IDa = sys.argv[2]
    IDb = sys.argv[3]

    # Read as strings; name columns up front
    all_map = pd.read_csv(
        map_file,
        sep="\t",
        header=None,
        names=["UniProt_AC", "ID_type", "ID"],
        dtype="string",
        usecols=[0, 1, 2],
    )

    # Filter once, rename for clarity
    a = (
        all_map.loc[all_map["ID_type"] == IDa, ["UniProt_AC", "ID"]]
        .rename(columns={"ID": "IDa"})
        .drop_duplicates()
    )
    b = (
        all_map.loc[all_map["ID_type"] == IDb, ["UniProt_AC", "ID"]]
        .rename(columns={"ID": "IDb"})
        .drop_duplicates()
    )

    # Vectorized join: every (IDa, UniProt) gets all matching IDb rows
    id_map = a.merge(b, on="UniProt_AC", how="left")
    id_map["IDb"] = id_map["IDb"].fillna("NA")

    # Self-map if IDa appears among IDb words
    mask = id_map["IDa"].isin(id_map["IDb"])
    id_map.loc[mask, "IDb"] = id_map.loc[mask, "IDa"]

    # Match column names/order
    id_map = id_map.rename(columns={"UniProt_AC": "uniprot"})
    id_map.to_csv(sys.stdout, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()

