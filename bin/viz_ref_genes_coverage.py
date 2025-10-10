#!/usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import missingno as msno
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


plt.rcParams['axes.titlesize'] = 7        
plt.rcParams['axes.labelsize'] = 6        
plt.rcParams['xtick.labelsize'] = 6       
plt.rcParams['ytick.labelsize'] = 6       
plt.rcParams['legend.fontsize'] = 6       
plt.rcParams['figure.titlesize'] = 8   
plt.rcParams['font.size'] = 6

red_shade = "#d62728"
grey_shade = "#7f7f7f"

custom_palette = sns.color_palette(["red", "grey"])
sns.set_palette(custom_palette)
sns.set_style("white")

white_to_red = LinearSegmentedColormap.from_list("white_to_red", ["white", "#d62728"])


data_type_name_dict = {
    'gtex': 'mRNA expression',
    'eprot': 'Protein expr. (steady state)',
    'orthogroup': 'Gene coevolution',
    'depmap': 'CRISPR-KO sensitivity',
    'ubiquitination': 'Protein ubiquitination',
    'humap3': 'Protein complexes',
    'ptmdb': 'Protein phosphorylation',
    'proteomehd': 'Protein expr. (perturbations)',
    'mitchell2023': 'Protein expr. (24h drug perturb.)'
}

data_source_name_dict = {
    'gtex': 'GTEx',
    'eprot': 'PRIDE (tissues)',
    'orthogroup': 'Orthogroups',
    'depmap': 'DepMap',
    'ubiquitination': 'Ubiquitination',
    'humap3': 'hu.MAP 3.0',
    'ptmdb': 'PTMDB',
    'proteomehd': 'ProteomeHD',
    'mitchell2023': 'PRIDE (perturb.)'
}


def plot_missing_feature_cdf(df, na_cumsum_pdf):
    """
    Plots the cumulative distribution of the number of samples missing less than
    a given number of features in a DataFrame.
    Args:
        df: The pandas DataFrame.
    """
    # Count missing features per sample:
    missing_counts = df.isnull().sum(axis=1)
    # Calculate the cumulative distribution:
    n_samples = len(df)
    counts = missing_counts.value_counts().sort_index()  # Sort by number of missing features
    cumulative_counts = counts.cumsum()
    cumulative_probabilities = cumulative_counts / n_samples
    # Plot the CDF:
    plt.clf()
    plt.figure(figsize=(2.5, 2)) 
    plt.plot(cumulative_probabilities.index, cumulative_probabilities.values, marker='.', linestyle='-')
    plt.xlabel("Number of data sources")
    plt.ylabel("Cumulative probability")
    plt.title("Reference genes missing in\nless than X data sources")

    plt.xticks(np.arange(0, df.shape[1] + 1, 1)) # Set x-axis ticks to integers
    plt.xlim(-0.5, df.shape[1] + 0.5) # Set x-axis limits to include 0 and the maximum number of features
    plt.ylim(0, 1.05)
    
    # Make ticks visible on both axes
    plt.tick_params(axis='both', which='major', direction='out', length=4, width=1)
    plt.tick_params(axis='x', which='major', top=False, bottom=True)
    plt.tick_params(axis='y', which='major', left=True, right=False)
    
    sns.despine()
    plt.tight_layout()
    plt.savefig(na_cumsum_pdf)


def plot_msno_bar(df, msno_bar_pdf):
    plt.clf()
    ax = msno.bar(df, color=red_shade, figsize=(2.5, 2.5), fontsize=6, labels=True, sort='descending')
    
    # Set labels on the main axis (ax)
    ax.set_ylabel("Fraction of reference genes")
    ax.set_xlabel("Data sources")
    
    # Clean up other axes
    if len(ax.figure.axes) > 1:
        ax.figure.axes[1].set_yticklabels([])
    if len(ax.figure.axes) > 2:
        ax.figure.axes[2].set_xticklabels([])

    sns.despine()
    plt.tight_layout()
    plt.savefig(msno_bar_pdf)
    

def plot_msno_matrix(df, msno_matrix_pdf):
    plt.clf()
    ax = msno.matrix(df, figsize=(3, 3), fontsize=6)
    sns.despine()
    plt.tight_layout()
    plt.savefig(msno_matrix_pdf)


def plot_msno_heatmap(df, msno_heatmap_pdf):
    plt.clf()
    ax = msno.heatmap(df, figsize=(3, 2.5), fontsize=6, cmap="coolwarm")
    sns.despine()
    plt.tight_layout()
    plt.savefig(msno_heatmap_pdf)


def plot_msno_dendrogram(df, msno_dendro_pdf):
    plt.clf()
    plt.figure(figsize=(3, 2))
    
    ax = msno.dendrogram(df, figsize=(3, 2), fontsize=6)

    # Hide all spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Hide ticks and tick labels
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    plt.tight_layout()
    plt.savefig(msno_dendro_pdf)


def main():
    idmapping_file = sys.argv[1]
    na_cumsum_pdf = sys.argv[2]
    msno_bar_pdf = sys.argv[3]
    msno_matrix_pdf = sys.argv[4]
    msno_heatmap_pdf = sys.argv[5]
    msno_dendro_pdf = sys.argv[6]

    idmapping = pd.read_csv(idmapping_file, sep='\t', index_col=0)
    idmapping = idmapping.drop(columns=['reactome'])

    idmapping = idmapping.rename(columns=data_source_name_dict)

    plot_missing_feature_cdf(idmapping, na_cumsum_pdf)

    plot_msno_bar(idmapping, msno_bar_pdf)

    plot_msno_matrix(idmapping, msno_matrix_pdf)

    plot_msno_heatmap(idmapping, msno_heatmap_pdf)

    plot_msno_dendrogram(idmapping, msno_dendro_pdf)

            
if __name__ == '__main__':
    main()