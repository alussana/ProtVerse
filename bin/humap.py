#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    """
    humap_file = 'input/net.tsv.gz'
    """
    humap_file = sys.argv[1]

    humap = pd.read_csv(humap_file, sep='\t', index_col=None, header=None)

    humap.columns = ['gene_a','gene_b','prob']

    humap.round(decimals={'prob': 5})

    print(humap.to_csv(sep='\t', index=False, header=True))

if __name__ == '__main__':
    main()