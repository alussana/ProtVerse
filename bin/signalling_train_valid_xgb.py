#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import xgboost as xgb
from sklearn.metrics import RocCurveDisplay
from sklearn.inspection import permutation_importance
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

custom_palette = sns.color_palette([red_shade, grey_shade], n_colors=2)
sns.set_palette(custom_palette)


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
    n_jobs = 2
    max_depth = 2
    n_estimators = 30
    lr = 0.03
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
    train_roc_pdf_file = sys.argv[9]
    valid_roc_pdf_file = sys.argv[10]
    train_prc_pdf_file = sys.argv[11]
    valid_prc_pdf_file = sys.argv[12]
    model_pkl_file = sys.argv[13]
    valid_mdi_importance_pdf_file = sys.argv[14]
    train_perm_importance_pdf_file = sys.argv[15]
    valid_perm_importance_pdf_file = sys.argv[16]

    with open(examples_X_train_file, 'rb') as examples_X_train:
        X_train = pickle.load(examples_X_train)
    with open(examples_y_train_file, 'rb') as examples_y_train:
        y_train = pickle.load(examples_y_train)
    with open(examples_X_valid_file, 'rb') as examples_X_valid:
        X_valid = pickle.load(examples_X_valid)
    with open(examples_y_valid_file, 'rb') as examples_y_valid:
        y_valid = pickle.load(examples_y_valid)

    feature_names = list(X_train.columns)

    # set monotonic contraints
    filtered_features = [name for name in feature_names if not re.search(r'Data|Overlap', name)]
    monotonic_constraints = {name: 1 for name in filtered_features}

    #X_train = X_train.loc[:, feature_names]
    #X_valid = X_valid.loc[:, feature_names]
    X_train = X_train.loc[:, filtered_features]
    X_valid = X_valid.loc[:, filtered_features]

    #y_train = np.array(y_train)
    #y_valid = np.array(y_valid)

    # train the full model
    forest = model(
        n_jobs=n_jobs,
        n_estimators=n_estimators,
        max_depth=max_depth,
        lr=lr,
        monotonic_constraints=monotonic_constraints,
    )   
    forest.fit(X_train, y_train)

    # serialize the model to disk
    with open(model_pkl_file, 'wb') as model_pkl:
        pickle.dump(forest, model_pkl)
        

    # mcc and accuracy on train and valid datasets
    y_pred_train = forest.predict(X_train)
    mcc_train = round(matthews_corrcoef(y_train, y_pred_train), 2)
    acc_train = round(accuracy_score(y_train, y_pred_train), 2)
    y_pred_valid = forest.predict(X_valid)
    mcc_valid = round(matthews_corrcoef(y_valid, y_pred_valid), 2)
    acc_valid = round(accuracy_score(y_valid, y_pred_valid), 2)
        
    
    # Get ROC curve on training data
    plt.clf()
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    display = RocCurveDisplay.from_estimator(forest,
            X_train,
            y_train,
            name="XGB",
            alpha=1,
            lw=1,
            ax=ax
    )
    ax.plot([0, 1], [0, 1], linestyle="--", lw=0.5, color="black", label="Chance", alpha=0.8)
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Signalling Interactions (Training)\nMCC = {mcc_train}; ACC = {acc_train}',
        xlabel='FPR',
        ylabel='TPR',
    )
    fig, ax = display.figure_, display.ax_
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(loc="lower right", frameon=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.savefig(train_roc_pdf_file)


    # PR curves
    plt.clf()
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    display = PrecisionRecallDisplay.from_estimator(
        forest, X_train, y_train, name="XGB", lw=1
    )
    fig, ax = display.figure_, display.ax_
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Signalling Interactions (Training)',
        xlabel='Recall',
        ylabel='Precision'
    )
    plt.axhline(y=0.5, color='black', linestyle='--', label='Chance', lw=0.5, alpha=0.8)
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(loc="lower right", frameon=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.savefig(train_prc_pdf_file)
    
    plt.clf()
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    display = PrecisionRecallDisplay.from_estimator(
        forest, X_valid, y_valid, name="XGB", lw=1
    )
    fig, ax = display.figure_, display.ax_
    plt.axhline(y=0.5, color='black', linestyle='--', label='Chance', lw=0.5, alpha=0.8)
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Signalling Interactions (Validation)',
        xlabel='Recall',
        ylabel='Precision'
    )
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(loc="lower right", frameon=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.savefig(valid_prc_pdf_file)


    # Get ROC curve on valid data
    plt.clf()
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    display = RocCurveDisplay.from_estimator(forest,
            X_valid,
            y_valid,
            name="XGB",
            alpha=1,
            lw=1,
            ax=ax
    )
    fig, ax = display.figure_, display.ax_
    ax.plot([0, 1], [0, 1], linestyle="--", lw=0.5, color="black", label="Chance", alpha=0.8)
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=f'Signalling Interactions (Validation)\nMCC = {mcc_valid}; ACC = {acc_valid}',
        xlabel='FPR',
        ylabel='TPR',
    )
    ax.set_aspect(1 / ax.get_data_ratio())
    ax.legend(loc="lower right", frameon=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.savefig(valid_roc_pdf_file)


    # feature importance based on feature permutation
    result = permutation_importance(    
        forest, X_train, y_train, n_repeats=64, random_state=42,
        scoring=make_scorer(matthews_corrcoef), n_jobs=n_jobs
    )
    sorted_importances_idx = result.importances_mean.argsort()
    sorted_features = [y for (x,y) in sorted(zip(sorted_importances_idx, feature_names))]
    importances = pd.DataFrame(
        result.importances[sorted_importances_idx].T,
        columns=sorted_features,
    )
    plt.clf()
    ax = importances.plot.box(vert=False, whis=10)
    fig = ax.figure
    fig.set_size_inches(5, 9)
    ax.set_title("Feature Importance (Training)")
    ax.axvline(x=0, color="k", linestyle="--")
    ax.set_xlabel("Decrease in MCC score")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.figure.tight_layout()
    plt.savefig(train_perm_importance_pdf_file)

    result = permutation_importance(
        forest, X_valid, y_valid, n_repeats=64, random_state=42,
        scoring=make_scorer(matthews_corrcoef), n_jobs=n_jobs
    )
    sorted_importances_idx = result.importances_mean.argsort()
    sorted_features = [y for (x,y) in sorted(zip(sorted_importances_idx, feature_names))]
    importances = pd.DataFrame(
        result.importances[sorted_importances_idx].T,
        columns=sorted_features,
    )
    plt.clf()
    ax = importances.plot.box(vert=False, whis=10)
    fig = ax.figure
    fig.set_size_inches(5, 9)
    ax.set_title("Feature Importance (Validation)")
    ax.axvline(x=0, color="k", linestyle="--")
    ax.set_xlabel("Decrease in MCC score")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.figure.tight_layout()
    plt.savefig(valid_perm_importance_pdf_file)


if __name__ == '__main__':
    main()
    print('END', file=sys.stderr)