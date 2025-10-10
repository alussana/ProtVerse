#!/usr/bin/env python3

import sys
import pandas as pd
import pickle
from sklearn.model_selection import train_test_split

def split_dataset(X, y):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=0)
    return(X_train, X_test, y_train, y_test)

def main():
    examples_file = sys.argv[1]
    examples_X_train_file = sys.argv[2]
    examples_X_test_file = sys.argv[3]
    examples_y_train_file = sys.argv[4]
    examples_y_test_file = sys.argv[5]
    
    examples = pd.read_csv(examples_file, sep='\t', index_col=0)

    X = examples.drop(labels=["label"], axis=1)
    y = examples["label"]
    
    X_train, X_test, y_train, y_test = split_dataset(X, y)

    with open(examples_X_train_file, 'wb') as examples_X_train:
        pickle.dump(X_train, examples_X_train)
    with open(examples_y_train_file, 'wb') as examples_y_train:
        pickle.dump(y_train, examples_y_train)
    with open(examples_X_test_file, 'wb') as examples_X_test:
        pickle.dump(X_test, examples_X_test)
    with open(examples_y_test_file, 'wb') as examples_y_test:
        pickle.dump(y_test, examples_y_test)

if __name__ == '__main__':
    main()