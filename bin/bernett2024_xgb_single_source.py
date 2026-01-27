#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import xgboost as xgb
from copy import deepcopy
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import matthews_corrcoef, make_scorer, accuracy_score
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

#custom_palette = sns.color_palette([red_shade, grey_shade], n_colors=2)
#sns.set_palette(custom_palette)


def model(n_jobs: int, n_estimators: int, max_depth: int, lr: float, monotonic_constraints: list):   
    model = xgb.XGBClassifier(
        n_jobs=n_jobs,
        missing=-1,
        monotone_constraints=monotonic_constraints,
        n_estimators=n_estimators,
        learning_rate=lr,
        max_depth=max_depth,
        random_state=42,
    )
    return(model)


def main():
    """
    n_jobs = 1
    max_depth = 16
    n_estimators = 50
    examples_X_train_file = 'input/examples_X_train.pkl'
    examples_X_valid_file = 'input/examples_X_valid.pkl'
    examples_y_train_file = 'input/examples_y_train.pkl'
    examples_y_valid_file = 'input/examples_y_valid.pkl'
    """
    n_jobs = int(sys.argv[1])
    max_depth = int(sys.argv[2])
    n_estimators = int(sys.argv[3])
    lr=float(sys.argv[4])
    examples_X_train_file = sys.argv[5]
    examples_X_valid_file = sys.argv[6]
    examples_y_train_file = sys.argv[7]
    examples_y_valid_file = sys.argv[8]
    valid_roc_pdf_file = sys.argv[9]
    valid_restricted_roc_pdf_file = sys.argv[10]
    valid_prc_pdf_file = sys.argv[11]
    valid_restricted_prc_pdf_file = sys.argv[12]

    with open(examples_X_train_file, 'rb') as examples_X_train:
        X_train = pickle.load(examples_X_train)
    with open(examples_y_train_file, 'rb') as examples_y_train:
        y_train = pickle.load(examples_y_train)
    with open(examples_X_valid_file, 'rb') as examples_X_valid:
        X_valid = pickle.load(examples_X_valid)
    with open(examples_y_valid_file, 'rb') as examples_y_valid:
        y_valid = pickle.load(examples_y_valid)
        

    feature_names = list(X_train.columns)

    # filter features, set monotonic contraints
    filtered_features = [name for name in feature_names if not re.search(r'Data|Overlap', name)]
    monotonic_constraints = {name: 1 for name in filtered_features}

    X_train = X_train.loc[:, filtered_features]
    X_valid = X_valid.loc[:, filtered_features]


    # split the features by data source
    features_sources = {
        'PRIDE (tissues)': [i for i in range(0, 13)],
        'ProteomeHD': [i for i in range(13, 26)],
        'PRIDE (perturb.)': [i for i in range(26, 39)],
        'GTEx': [i for i in range(39, 52)],
        'PTMDB': [i for i in range(52, 90)],
        'Ubiquitination': [i for i in range(90, 95)],
        'DepMap': [i for i in range(95, 108)],
        'Orthogroup': [i for i in range(108, 110)],
        'hu.MAP 3.0': [110],
        'Cotranslocation': [i for i in range(111, 114)],
    }


    # train the single-source models (ROC)
    models = {}
    mccs = {}
    aucs = {}
    tprs = {}
    fprs = {}
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    for data_source in features_sources.keys():
        feature_indexes = features_sources[data_source]
        X_train_reduced = X_train.iloc[:, feature_indexes]
        X_valid_reduced = X_valid.iloc[:, feature_indexes]
        features_reduced = [filtered_features[i] for i in feature_indexes]
        monotonic_constraints_reduced = {
            k: v for k, v in monotonic_constraints.items() if k in features_reduced
        }
        forest = model(
            n_jobs=n_jobs,
            n_estimators=n_estimators,
            max_depth=max_depth,
            lr=lr,
            monotonic_constraints=monotonic_constraints_reduced,
        )   
        forest.fit(X_train_reduced, y_train)
        models[data_source] = deepcopy(forest)
        # mcc on valid dataset
        y_pred_valid = forest.predict(X_valid_reduced)
        mcc_valid = matthews_corrcoef(y_valid, y_pred_valid)
        mccs[data_source] = deepcopy(mcc_valid)
        # Get ROC curve on valid data
        viz = RocCurveDisplay.from_estimator(forest,
                X_valid_reduced,
                y_valid,
                name=f'{data_source}',
                alpha=1,
                lw=0.5,
                ax=ax
        )
        tprs[data_source] = deepcopy(viz.tpr)
        aucs[data_source] = deepcopy(viz.roc_auc)
        fprs[data_source] = deepcopy(viz.fpr)

    ax.plot([0, 1], [0, 1], linestyle="--", lw=0.5, color="black", label="Chance", alpha=0.8)

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Physical Interactions (Validation)',
        xlabel="FPR",
        ylabel="TPR"
    )
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(
        loc="center left", 
        bbox_to_anchor=(1.05, 0.5),
        frameon=False
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(valid_roc_pdf_file, bbox_inches="tight")


    # train the single-source models (PRC)
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    for data_source in features_sources.keys():
        feature_indexes = features_sources[data_source]
        X_train_reduced = X_train.iloc[:, feature_indexes]
        X_valid_reduced = X_valid.iloc[:, feature_indexes]
        features_reduced = [filtered_features[i] for i in feature_indexes]
        monotonic_constraints_reduced = {
            k: v for k, v in monotonic_constraints.items() if k in features_reduced
        }
        forest = model(
            n_jobs=n_jobs,
            n_estimators=n_estimators,
            max_depth=max_depth,
            lr=lr,
            monotonic_constraints=monotonic_constraints_reduced,
        )   
        forest.fit(X_train_reduced, y_train)
        models[data_source] = deepcopy(forest)
        # mcc on valid dataset
        y_pred_valid = forest.predict(X_valid_reduced)
        mcc_valid = matthews_corrcoef(y_valid, y_pred_valid)
        mccs[data_source] = deepcopy(mcc_valid)
        # Get PR curve on valid data
        viz = PrecisionRecallDisplay.from_estimator(forest,
                X_valid_reduced,
                y_valid,
                name=f'{data_source}',
                alpha=1,
                lw=0.5,
                ax=ax
        )

    plt.axhline(y=0.5, color='black', linestyle='--', label='Chance', lw=0.5, alpha=0.8)

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Physical Interactions (Validation)',
        xlabel='Recall',
        ylabel='Precision'
    )
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(
        loc="center left", 
        bbox_to_anchor=(1.05, 0.5),
        frameon=False
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(valid_prc_pdf_file, bbox_inches="tight")


    # single-feature version (ROC)
    features_sources = {
        'PRIDE (tissues)': [0],
        'ProteomeHD': [13],
        'PRIDE (perturb.)': [26],
        'GTEx': [39],
        'PTMDB': [52],
        'Ubiquitination': [90],
        'DepMap': [95],
        'Orthogroup': [108],
        'hu.MAP 3.0': [110],
        'Cotranslocation': [111],
    }
    
    plt.clf()
    models = {}
    mccs = {}
    aucs = {}
    tprs = {}
    fprs = {}
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    for data_source in features_sources.keys():
        feature_indexes = features_sources[data_source]
        X_train_reduced = X_train.iloc[:, feature_indexes]
        X_valid_reduced = X_valid.iloc[:, feature_indexes]
        features_reduced = [filtered_features[i] for i in feature_indexes]
        monotonic_constraints_reduced = {
            k: v for k, v in monotonic_constraints.items() if k in features_reduced
        }
        forest = model(
            n_jobs=n_jobs,
            n_estimators=n_estimators,
            max_depth=max_depth,
            lr=lr,
            monotonic_constraints=monotonic_constraints_reduced,
        )   
        forest.fit(X_train_reduced, y_train)
        models[data_source] = deepcopy(forest)
        # mcc on valid dataset
        y_pred_valid = forest.predict(X_valid_reduced)
        mcc_valid = matthews_corrcoef(y_valid, y_pred_valid)
        mccs[data_source] = deepcopy(mcc_valid)
        # Get ROC curve on valid data
        viz = RocCurveDisplay.from_estimator(forest,
                X_valid_reduced,
                y_valid,
                name=f'{data_source}',
                alpha=1,
                lw=0.5,
                ax=ax
        )
        tprs[data_source] = deepcopy(viz.tpr)
        aucs[data_source] = deepcopy(viz.roc_auc)
        fprs[data_source] = deepcopy(viz.fpr)

    ax.plot([0, 1], [0, 1], linestyle="--", lw=0.5, color="black", label="Chance", alpha=0.8)

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Physical Interactions (Validation)',
        xlabel="FPR",
        ylabel="TPR"
    )
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(
        loc="center left", 
        bbox_to_anchor=(1.05, 0.5),
        frameon=False
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(valid_restricted_roc_pdf_file, bbox_inches="tight")


    # single-feature version (PRC)
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    for data_source in features_sources.keys():
        feature_indexes = features_sources[data_source]
        X_train_reduced = X_train.iloc[:, feature_indexes]
        X_valid_reduced = X_valid.iloc[:, feature_indexes]
        features_reduced = [filtered_features[i] for i in feature_indexes]
        monotonic_constraints_reduced = {
            k: v for k, v in monotonic_constraints.items() if k in features_reduced
        }
        forest = model(
            n_jobs=n_jobs,
            n_estimators=n_estimators,
            max_depth=max_depth,
            lr=lr,
            monotonic_constraints=monotonic_constraints_reduced,
        )   
        forest.fit(X_train_reduced, y_train)
        models[data_source] = deepcopy(forest)
        # mcc on valid dataset
        y_pred_valid = forest.predict(X_valid_reduced)
        mcc_valid = matthews_corrcoef(y_valid, y_pred_valid)
        mccs[data_source] = deepcopy(mcc_valid)
        # Get PR curve on valid data
        viz = PrecisionRecallDisplay.from_estimator(forest,
                X_valid_reduced,
                y_valid,
                name=f'{data_source}',
                alpha=1,
                lw=0.5,
                ax=ax
        )

    plt.axhline(y=0.5, color='black', linestyle='--', label='Chance', lw=0.5, alpha=0.8)

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Physical Interactions (Validation)',
        xlabel='Recall',
        ylabel='Precision'
    )
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(
        loc="center left", 
        bbox_to_anchor=(1.05, 0.5),
        frameon=False
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(valid_restricted_prc_pdf_file, bbox_inches="tight")

    
if __name__ == '__main__':
    main()