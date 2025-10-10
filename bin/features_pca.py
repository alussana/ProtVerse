#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from seaborn import kdeplot
from sklearn.decomposition import PCA
import re


def main():
    """
    examples_file = 'input/examples.tsv.gz'
    pca12_omics_features_png_file = 'pca12_omics_features.png'
    pca12_stringdb_features_png_file = 'pca12_stringdb_features.png'
    pca34_omics_features_png_file = 'pca34_omics_features.png'
    pca34_stringdb_features_png_file = 'pca34_stringdb_features.png'
    pca56_omics_features_png_file = 'pca56_omics_features.png'
    pca56_stringdb_features_png_file = 'pca56_omics_features.png'
    """
    examples_file = sys.argv[1]
    pca12_omics_features_png_file = sys.argv[2]
    pca34_omics_features_png_file = sys.argv[3]
    pca56_omics_features_png_file = sys.argv[4]
    #pca12_stringdb_features_png_file = sys.argv[5]
    #pca34_stringdb_features_png_file = sys.argv[6]
    #pca56_stringdb_features_png_file = sys.argv[7]
    
    data = pd.read_csv(examples_file, sep='\t', index_col=0)

    feature_names = list(data.columns)
    filtered_features = [name for name in feature_names if not re.search(r'Data|Overlap', name)]

    #stringdb_features_cols = list(range(len(data.columns)-9, len(data.columns)-1))
    #d = data.drop(data.columns[stringdb_features_cols], axis=1)
    d = data
    omics_features = d.loc[:, filtered_features]
    omics_features = omics_features.drop('label', axis=1)
    #stringdb_features = data.loc[:, data.columns[stringdb_features_cols]]
    labels = d['label']

    # discard indicator features
    #omics_features = omics_features.filter(regex=r'^(?!.+_indicator$).*$')

    # replace -1 (missing value) with zeroes
    #omics_features.replace(-1, 0, inplace=True)

    # MinMax scale
    #scaler = MinMaxScaler()
    #omics_features = pd.DataFrame(
    #    scaler.fit_transform(omics_features),
    #    columns=omics_features.columns,
    #    index=omics_features.index
    #)
    #stringdb_features = pd.DataFrame(
    #    scaler.fit_transform(stringdb_features),
    #    columns=stringdb_features.columns,
    #    index=stringdb_features.index
    #)

    # compute PCA on omics features
    pca_omics_features = PCA(n_components=6, svd_solver='full')
    pca_omics_features.fit(omics_features)
    pcs_omics_features = pd.DataFrame(pca_omics_features.transform(omics_features))

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

    
    # plot omics features PC1 and PC2
    pcs = pcs_omics_features.iloc[:, :2]
    pc1_name = f'PC1 ({round(pca_omics_features.explained_variance_ratio_[0], 3)})'
    pc2_name = f'PC2 ({round(pca_omics_features.explained_variance_ratio_[1], 3)})'
    pcs.columns = [pc1_name, pc2_name]
    pcs.index = omics_features.index
    pcs = pd.concat([pcs, labels], axis=1)
    pcs["label"] = pcs["label"].astype("category")
    pcs.rename(columns={"label": "Label"}, inplace=True)

    #custom_palette = sns.color_palette(["red", "grey"], n_colors=2)
    custom_palette = sns.color_palette([grey_shade, red_shade], n_colors=2)
    sns.set_palette(custom_palette)

    plt.clf()
    plt.figure(figsize=(2, 2))
    ax = kdeplot(
        data=pcs, x=pc1_name, y=pc2_name, hue="Label",
    )
    ax.set(title='Signalling Interactions')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(loc='upper right', fontsize=6, title_fontsize=6, frameon=False)
    plt.tight_layout()
    plt.savefig(pca12_omics_features_png_file)
    plt.clf()

    # plot omics features PC3 and PC4
    pcs = pcs_omics_features.iloc[:, 2:4]
    pc3_name = f'PC3 ({round(pca_omics_features.explained_variance_ratio_[2], 4)})'
    pc4_name = f'PC4 ({round(pca_omics_features.explained_variance_ratio_[3], 4)})'
    pcs.columns = [pc3_name, pc4_name]
    pcs.index = omics_features.index
    pcs = pd.concat([pcs, labels], axis=1)
    pcs["label"] = pcs["label"].astype("category")
    pcs.rename(columns={"label": "Label"}, inplace=True)
    plt.clf()
    plt.figure(figsize=(2, 2))
    ax = kdeplot(
        data=pcs, x=pc3_name, y=pc4_name, hue="Label"
    )
    ax.set(title='Signalling Interactions')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(loc='upper right', fontsize=6, title_fontsize=6, frameon=False)
    plt.tight_layout()
    plt.savefig(pca34_omics_features_png_file)

    # plot omics features PC5 and PC6
    pcs = pcs_omics_features.iloc[:, 4:6]
    pc5_name = f'PC5 ({round(pca_omics_features.explained_variance_ratio_[4], 4)})'
    pc6_name = f'PC6 ({round(pca_omics_features.explained_variance_ratio_[5], 4)})'
    pcs.columns = [pc5_name, pc6_name]
    pcs.index = omics_features.index
    pcs = pd.concat([pcs, labels], axis=1)
    pcs["label"] = pcs["label"].astype("category")
    pcs.rename(columns={"label": "Label"}, inplace=True)
    plt.clf()
    plt.figure(figsize=(2, 2))
    ax = kdeplot(
        data=pcs, x=pc5_name, y=pc6_name, hue="Label",
    )
    ax.set(title='Signalling Interactions')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(loc='upper right', fontsize=6, title_fontsize=6, frameon=False)
    plt.tight_layout()
    plt.savefig(pca56_omics_features_png_file)

    """
    # compute PCA on stringdb features
    pca_stringdb_features = PCA(n_components=6, svd_solver='full')
    pca_stringdb_features.fit(stringdb_features)
    pcs_stringdb_features = pd.DataFrame(pca_stringdb_features.transform(stringdb_features))
    
    # plot stringdb features PC1 and PC2
    pcs = pcs_stringdb_features.iloc[:, :2]
    pc1_name = f'PC1 ({round(pca_stringdb_features.explained_variance_ratio_[0], 4)})'
    pc2_name = f'PC2 ({round(pca_stringdb_features.explained_variance_ratio_[1], 4)})'
    pcs.columns = [pc1_name, pc2_name]
    pcs.index = stringdb_features.index
    pcs = pd.concat([pcs, labels], axis=1)
    plt.clf()
    ax = kdeplot(
        data=pcs, x=pc1_name, y=pc2_name, hue="label",
    )
    ax.set(title='STRING database-derived features PCA')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig(pca12_stringdb_features_png_file)

    # plot stringdb features PC3 and PC4
    pcs = pcs_stringdb_features.iloc[:, 2:4]
    pc3_name = f'PC3 ({round(pca_stringdb_features.explained_variance_ratio_[2], 4)})'
    pc4_name = f'PC4 ({round(pca_stringdb_features.explained_variance_ratio_[3], 4)})'
    pcs.columns = [pc3_name, pc4_name]
    pcs.index = stringdb_features.index
    pcs = pd.concat([pcs, labels], axis=1)
    plt.clf()
    ax = kdeplot(
        data=pcs, x=pc3_name, y=pc4_name, hue="label",
    )
    ax.set(title='STRING database-derived features PCA')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig(pca34_stringdb_features_png_file)

    # plot stringdb features PC5 and PC6 - not working for some reason
    pcs = pcs_stringdb_features.iloc[:, 4:6]
    pc5_name = f'PC5 ({round(pca_stringdb_features.explained_variance_ratio_[4], 4)})'
    pc6_name = f'PC6 ({round(pca_stringdb_features.explained_variance_ratio_[5], 4)})'
    pcs.columns = [pc5_name, pc6_name]
    pcs.index = stringdb_features.index
    pcs = pd.concat([pcs, labels], axis=1)
    plt.clf()
    ax = kdeplot(
        data=pcs, x=pc5_name, y=pc6_name, hue="label",
    )
    ax.set(title='STRING database-derived features PCA')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig(pca56_stringdb_features_png_file)
    """

if __name__ == '__main__':
    main()