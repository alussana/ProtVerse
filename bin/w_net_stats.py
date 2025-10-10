#!/usr/bin/env python3

import sys
import seaborn as sns
import matplotlib.pyplot as plt


sns.set_theme(style="ticks")
sns.set_style("white")

plt.rcParams['axes.titlesize'] = 7        
plt.rcParams['axes.labelsize'] = 6        
plt.rcParams['xtick.labelsize'] = 6       
plt.rcParams['ytick.labelsize'] = 6       
plt.rcParams['legend.fontsize'] = 6       
plt.rcParams['figure.titlesize'] = 8   
plt.rcParams['font.size'] = 6

red_shade = "#d62728"
grey_shade = "#7f7f7f"

custom_palette = sns.color_palette([red_shade, grey_shade], n_colors=2)
sns.set_palette(custom_palette)


def main():
    """
    edges_tsv = 'input/net.tsv'
    net_id = 'wcsn'
    outFile = 'test.png'
    """
    edges_tsv = sys.argv[1]
    outFile = sys.argv[2]

    edge_weights = []
    with open(edges_tsv, 'r') as edges_fh:
        for line in edges_fh:
            w = float(line.strip().split('\t')[-1])
            edge_weights.append(w)

    sns.set_palette("blend:#009F4D,#287AE2")
    plt.clf()
    f, ax = plt.subplots()
    sns.despine(f)
    sns.histplot(edge_weights)
    ax.set_xlabel('edge weight')
    ax.set_ylabel('counts (edges)')
    plt.savefig(outFile)

if __name__ == '__main__':
    main()