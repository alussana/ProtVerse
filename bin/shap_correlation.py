#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import shap


sns.set_theme(style="ticks")

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


def correlation_heatmap(df: pd.DataFrame, output_path: str, title: str):
    cg = sns.clustermap(df, metric="euclidean", standard_scale=None,
                        method="ward", cmap="magma", robust='TRUE',
                        row_cluster=True, col_cluster=True, yticklabels=True,
                        xticklabels=True)
    cg.ax_col_dendrogram.set_visible(False)
    cg.ax_heatmap.set_ylabel("")
    #title = f'{title} ({len(df.index)} elements)'
    cg.fig.suptitle(title)
    cg.savefig(output_path)

def main():
    """
    examples_X_test_file = 'X_test_nNeg43289.pkl'
    shap_pkl_file = 'shap_values.pkl'

    examples_X_test_file = 'input/examples_X_test.pkl'
    shap_pkl_file = 'input/shap_values.pkl'
    heatmap_png_file = 'test.png'
    """
    shap_pkl_file = sys.argv[1]
    examples_X_test_file = sys.argv[2]
    heatmap_png_file = sys.argv[3]

    with open(shap_pkl_file, 'rb') as shap_pkl:
        shap_values = pickle.load(shap_pkl)
  
    with open(examples_X_test_file, 'rb') as examples_X_test:
        X_test = pickle.load(examples_X_test)

    # take all non stringdb-derived features' names
    feature_names = list(X_test.columns)[:21]

    shap_df = pd.DataFrame(shap_values.values, columns=feature_names)
    shap_pearson = shap_df.corr()
    correlation_heatmap(
        df=shap_pearson,
        output_path=heatmap_png_file,
        title='Shap values Pearson coefficient'
    )

if __name__ == '__main__':
    main()