#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    rnk_file = sys.argv[1]
    df = pd.read_csv(rnk_file, sep='\t', header=None, index_col=0)
    df_groupedIndex_mean = df.groupby(df.index).mean()
    print(df_groupedIndex_mean.to_csv(sep='\t', index=True, header=False))

if __name__ == '__main__':
    main()