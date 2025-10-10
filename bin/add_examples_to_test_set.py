#!/usr/bin/env python3

import sys
import pandas as pd
import pickle

def main():
    """
    new_examples_tsv_file = 'examples_N43289.tsv.gz'
    examples_X_test_file = 'X_test_nNeg43289.pkl'
    examples_y_test_file = 'y_test_nNeg43289.pkl'
    new_X_test_file = 'X_test_nNeg262144.pkl'
    new_y_test_file = 'y_test_nNeg262144.pkl'
    """
    new_examples_tsv_file = sys.argv[1]
    examples_X_test_file = sys.argv[2]
    examples_y_test_file = sys.argv[3]
    new_X_test_file = sys.argv[4]
    new_y_test_file = sys.argv[5]

    with open(examples_X_test_file, 'rb') as examples_X_test:
        X_test = pickle.load(examples_X_test)
    with open(examples_y_test_file, 'rb') as examples_y_test:
        y_test = pickle.load(examples_y_test)
    new_examples = pd.read_csv(new_examples_tsv_file, sep='\t', index_col=0)

    # remove existing negative label examples from the existing test set
    X_test = X_test.loc[y_test==1]
    y_test = y_test.loc[y_test==1]

    # Add the new negative label examples to the test set
    X = new_examples.drop(labels=["label"], axis=1)
    y = new_examples["label"]
    new_X_test = pd.concat([X_test, X], axis=0)
    new_y_test = pd.concat([y_test, y], axis=0)

    with open(new_X_test_file, 'wb') as new_examples_X:
        pickle.dump(new_X_test, new_examples_X)
    with open(new_y_test_file, 'wb') as new_examples_y:
        pickle.dump(new_y_test, new_examples_y)

if __name__ == '__main__':
    main()