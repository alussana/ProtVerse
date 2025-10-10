#!/usr/bin/env python3

import sys
import pickle
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

def main():
    """
    features_file = 'input/features.tsv'
    model_file = 'input/model.pkl'
    """
    features_file = sys.argv[1]
    model_file = sys.argv[2]

    with open(model_file, 'rb') as model_fh:
        model = pickle.load(model_fh)

    X = pd.read_csv(features_file, sep='\t', index_col=0)

    # get predicted probabilities for class 1 (positive class)
    y_pred = [round(out[1], 3) for out in model.predict_proba(X)]
    edges = pd.DataFrame(y_pred, columns=['score'], index=X.index)

    print(edges.to_csv(sep='\t', header=False))

if __name__ == '__main__':
    main()