#!/usr/bin/env python3

# https://scikit-learn.org/stable/modules/calibration.html

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.calibration import CalibrationDisplay
from sklearn.calibration import CalibratedClassifierCV
import pickle

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


def main():
    """
    examples_X_test_file = 'input/X_test.pkl'
    examples_y_test_file = 'input/y_test.pkl'
    model_pkl_file = 'input/model.pkl'
    
    before_calibration_curve_pdf_file = 'training/calibration_before_maxDepth4_nEst64_minSplit2.pdf'
    after_calibration_curve_pdf_file = 'training/calibration_after_maxDepth4_nEst64_minSplit2.pdf'
    calibrated_forest_pkl_file = 'training/forest_calibrated_maxDepth8_nEst128_minSplit2.pkl'
    """
    examples_X_test_file = sys.argv[1]
    examples_y_test_file = sys.argv[2]
    model_pkl_file = sys.argv[3]
    before_calibration_curve_pdf_file = sys.argv[4]
    after_calibration_curve_pdf_file = sys.argv[5]
    calibrated_forest_pkl_file = sys.argv[6]

    with open(examples_X_test_file, 'rb') as examples_X_test:
        X_test = pickle.load(examples_X_test)
    with open(examples_y_test_file, 'rb') as examples_y_test:
        y_test = pickle.load(examples_y_test)
    with open(model_pkl_file, 'rb') as model:
        forest = pickle.load(model)

    # remove String-related features
    #feature_names = list(X_test.columns)[:21]
    
    feature_names = list(X_test.columns)

    # plot calibration curve
    plt.clf()
    calDisp = CalibrationDisplay.from_estimator(forest, X_test, y_test)
    fig, ax = calDisp.figure_, calDisp.ax_
    ax.set_aspect(1 / ax.get_data_ratio())
    fig.tight_layout()
    plt.savefig(before_calibration_curve_pdf_file)

    # perform calibration
    calibrated_forest = CalibratedClassifierCV(
        base_estimator=forest, cv='prefit', method='sigmoid'
    )
    calibrated_forest.fit(X_test, y_test)
    plt.clf()
    calDisp = CalibrationDisplay.from_estimator(calibrated_forest, X_test, y_test)
    fig, ax = calDisp.figure_, calDisp.ax_
    ax.set_aspect(1 / ax.get_data_ratio())
    fig.tight_layout()
    plt.savefig(after_calibration_curve_pdf_file)

    # export calibrated model
    with open(calibrated_forest_pkl_file, 'wb') as model_pkl:
        pickle.dump(calibrated_forest, model_pkl)

if __name__ == '__main__':
    main()