#!/usr/bin/env python3

import sys
import pandas as pd
import pickle


def main():
    """
    examples_tsv_gz = "input/examples.tsv.gz"
    train_ids_txt = "train_ids.txt"
    valid_ids_txt = "valid_ids.txt"
    test_ids_txt = "test_ids.txt"
    X_train_pkl = "training_bernett2024/X_train.pkl"
    X_valid_pkl = "training_bernett2024/X_valid.pkl"
    X_test_pkl = "training_bernett2024/X_test.pkl"
    y_train_pkl = "training_bernett2024/y_train.pkl"
    y_valid_pkl = "training_bernett2024/y_valid.pkl"
    y_test_pkl = "training_bernett2024/y_test.pkl"
    """
    examples_tsv_gz = sys.argv[1]
    train_ids_txt = sys.argv[2]
    valid_ids_txt = sys.argv[3]
    test_ids_txt = sys.argv[4]
    X_train_pkl = sys.argv[5]
    X_valid_pkl = sys.argv[6]
    X_test_pkl = sys.argv[7]
    y_train_pkl = sys.argv[8]
    y_valid_pkl = sys.argv[9]
    y_test_pkl = sys.argv[10]
    
    
    examples = pd.read_csv(examples_tsv_gz, sep='\t', index_col=0)
    
    
    train_ids = []
    with open(train_ids_txt, "r") as train_ids_fh:
        for line in train_ids_fh:
            train_ids.append(line.strip())
            
    valid_ids = []
    with open(valid_ids_txt, "r") as valid_ids_fh:
        for line in valid_ids_fh:
            valid_ids.append(line.strip())
            
    test_ids = []
    with open(test_ids_txt, "r") as test_ids_fh:
        for line in test_ids_fh:
            test_ids.append(line.strip())
    
    
    train_existing_ids = list(set(train_ids) & set(examples.index))
    train_df = examples.loc[train_existing_ids]
    
    valid_existing_ids = list(set(valid_ids) & set(examples.index))
    valid_df = examples.loc[valid_existing_ids]
    
    test_existing_ids = list(set(test_ids) & set(examples.index))
    test_df = examples.loc[test_existing_ids]
    

    X_train = train_df.drop(labels=["label"], axis=1)
    y_train = train_df["label"]
    
    X_valid = valid_df.drop(labels=["label"], axis=1)
    y_valid = valid_df["label"]
    
    X_test = test_df.drop(labels=["label"], axis=1)
    y_test = test_df["label"]


    with open(X_train_pkl, 'wb') as examples_X_train:
        pickle.dump(X_train, examples_X_train)
    with open(y_train_pkl, 'wb') as examples_y_train:
        pickle.dump(y_train, examples_y_train)
    with open(X_valid_pkl, 'wb') as examples_X_valid:
        pickle.dump(X_valid, examples_X_valid)
    with open(y_valid_pkl, 'wb') as examples_y_valid:
        pickle.dump(y_valid, examples_y_valid)
    with open(X_test_pkl, 'wb') as examples_X_test:
        pickle.dump(X_test, examples_X_test)
    with open(y_test_pkl, 'wb') as examples_y_test:
        pickle.dump(y_test, examples_y_test)


if __name__ == '__main__':
    main()