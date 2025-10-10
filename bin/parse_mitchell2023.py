#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    """
    table_S1_csv = 'input/table_S1.csv'
    """
    table_S1_csv = sys.argv[1]

    table_S1 = pd.read_csv(table_S1_csv, sep=',', encoding='unicode_escape')
    # remove 10 out of 9960 entries that do not have a reported gene name (they only have a UniProt AC)
    table_S1 = table_S1.loc[pd.isna(table_S1['Gene Name'])==False,]
    
    table_S1 = table_S1.drop('UniprotID', axis=1)

    print(table_S1.to_csv(sep='\t', index=False, header=True))

if __name__ == '__main__':
    main()