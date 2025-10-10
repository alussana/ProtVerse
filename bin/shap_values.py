#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.ensemble import RandomForestClassifier
import shap

def main():
    """
    examples_X_test_file = 'input/examples_X_test.pkl'
    examples_y_test_file = 'input/examples_y_test.pkl'
    model_file = 'input/model.pkl'
    pkl_out_file = 'training/shap_values.pkl'
    png_out_file = 'training/shap_test_maxDepth10_nEst200_minSplit2.png'
    """
    examples_X_test_file = sys.argv[1]
    examples_y_test_file = sys.argv[2]
    model_file = sys.argv[3]
    pkl_out_file = sys.argv[4]
    png_out_file = sys.argv[5]
    with open(examples_X_test_file, 'rb') as examples_X_test:
        X_test = pickle.load(examples_X_test)
    with open(examples_y_test_file, 'rb') as examples_y_test:
        y_test = pickle.load(examples_y_test)
    with open(model_file, 'rb') as model:
        forest = pickle.load(model)
    
    feature_names = list(X_test.columns)
    X_test = X_test.loc[:, feature_names]

    X_test = np.array(X_test)
    y_test = np.array(y_test)
    explainer = shap.Explainer(forest.predict, X_test)
    shap_values = explainer(X_test)
    shap_values.feature_names = feature_names
    with open(pkl_out_file, 'wb') as pkl_out:
        pickle.dump(shap_values, pkl_out)

    """
    shap.plots.bar(shap_values)
    shap.plots.beeswarm(shap_values) # or shap.summary_plot(shap_values, plot_type='violin')
    shap.plots.bar(shap_values[0])
    shap.plots.waterfall(shap_values[0])
    shap.plots.force(shap_values[0])
    shap.plots.force(shap_values)
    """
    shap.summary_plot(shap_values, plot_type='violin', max_display=30)
    plt.savefig(png_out_file)

if __name__ == '__main__':
    main()