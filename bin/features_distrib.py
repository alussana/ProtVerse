#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

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

custom_palette = sns.color_palette([grey_shade, red_shade], n_colors=2)
sns.set_palette(custom_palette)


def main():
    """
    examples_file='input/examples.tsv.gz'
    violins_h_out_pdf='features_distrib_h_to_55.pdf'
    violins_h_out_pdf='features_distrib_h_from_55.pdf'
    violins_v_out_pdf='features_distrib_v.pdf'
    """
    examples_file = sys.argv[1]
    violins_h_out_part1_pdf = sys.argv[2]
    violins_h_out_part2_pdf = sys.argv[3]
    violins_v_out_pdf = sys.argv[4]

    examples = pd.read_csv(examples_file, sep='\t')
    examples = examples.drop(labels='index', axis=1)
    class_col = 'label'
    label = examples[class_col].reset_index()
    label.columns = ['Example Index', 'Label']
    label["Label"] = label["Label"].replace({0: "Negative", 1: "Positive"})
    label['Label'] = label['Label'].astype('category')

    # select feature names by discarding the indicators
    features_cols = list(examples.columns)[:-1]
    #pattern = r'.+_indicator$'
    #features_cols = [item for item in features_cols if not re.match(pattern, item)]

    filtered_features = [name for name in features_cols if not re.search(r'Data|Overlap', name)]
    feature_order = sorted(filtered_features, key=str.casefold)
    cols = filtered_features.copy()
    cols.append(class_col)
    df = examples.copy()
    df = df.loc[:, filtered_features]
    
    df_stacked = pd.DataFrame(df.stack()).reset_index()
    
    df_stacked.columns = ['Example Index','Feature','Value']

    df_stacked = pd.merge(df_stacked, label, on='Example Index')

    df_stacked.replace(-1, np.nan, inplace=True)

    df_stacked.dropna(subset=['Value'], inplace=True)


    mid = 55
    halves = [feature_order[:mid], feature_order[mid:]]
    
    # horizontal - part 1    
    plt.clf()
    dims = (4, 7)
    fig, ax = plt.subplots(figsize=dims)
    ax = sns.violinplot(
        ax=ax, data=df_stacked, x='Value', y='Feature',
        hue='Label', split=True, inner=None, linewidth=0.5, orient='h',
        cut=0, bw_adjust=2, scale='width', order=halves[0]
    )
    plt.legend(loc='upper left', fontsize=6, title_fontsize=6, frameon=False)
    ax.legend(
        loc="center left",          
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=6,
        title_fontsize=6
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(violins_h_out_part1_pdf, bbox_inches = "tight")

    # horizontal - part 2    
    plt.clf()
    dims = (4, 7)
    fig, ax = plt.subplots(figsize=dims)
    ax = sns.violinplot(
        ax=ax, data=df_stacked, x='Value', y='Feature',
        hue='Label', split=True, inner=None, linewidth=0.5, orient='h',
        cut=0, bw_adjust=2, scale='width', order=halves[1]
    )
    plt.legend(loc='upper left', fontsize=6, title_fontsize=6, frameon=False)
    ax.legend(
        loc="center left",          
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=6,
        title_fontsize=6
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(violins_h_out_part2_pdf, bbox_inches = "tight")

    
    # vertical
    plt.clf()
    dims = (16, 4)
    fig, ax = plt.subplots(figsize=dims)
    ax = sns.violinplot(
        ax=ax, data=df_stacked, x='Feature', y='Value',
        hue='Label', split=True, inner=None, linewidth=0.5, orient='v',
        cut=0, bw_adjust=2, scale='width', order=feature_order
    )
    plt.xticks(rotation=90)
    ax.legend(
        loc="center left",          
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=6,
        title_fontsize=6
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(violins_v_out_pdf, bbox_inches = "tight")

if __name__ == '__main__':
    main()