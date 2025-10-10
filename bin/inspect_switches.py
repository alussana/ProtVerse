import sys
import matplotlib.pyplot as plt
import pandas as pd
import pickle

def main():
    examples_X_train_file = sys.argv[1]
    examples_X_test_file = sys.argv[2]
    examples_y_train_file = sys.argv[3]
    examples_y_test_file = sys.argv[4]
    hist_png_file = sys.argv[5]

    with open(examples_X_train_file, 'rb') as examples_X_train:
        X_train = pickle.load(examples_X_train)
    with open(examples_y_train_file, 'rb') as examples_y_train:
        y_train = pickle.load(examples_y_train)
    with open(examples_X_test_file, 'rb') as examples_X_test:
        X_test = pickle.load(examples_X_test)
    with open(examples_y_test_file, 'rb') as examples_y_test:
        y_test = pickle.load(examples_y_test)

    X_train['label'] = y_train
    X_test['label'] = y_test
    X_train_pos = X_train.loc[X_train.label==1, :]
    X_train_neg = X_train.loc[X_train.label==0, :]
    X_test_pos = X_test.loc[X_test.label==1, :]
    X_test_neg = X_test.loc[X_test.label==0, :]

    switches_train_pos = X_train_pos.loc[:, ['eprot_switch','gtex_switch','ptmdb_switch']]
    switches_train_neg = X_train_neg.loc[:, ['eprot_switch','gtex_switch','ptmdb_switch']]
    switches_test_pos = X_test_pos.loc[:, ['eprot_switch','gtex_switch','ptmdb_switch']]
    switches_test_neg = X_test_neg.loc[:, ['eprot_switch','gtex_switch','ptmdb_switch']]
    switches_train_pos_sum = switches_train_pos.apply(sum, axis=1)
    switches_test_pos_sum = switches_test_pos.apply(sum, axis=1)
    switches_train_neg_sum = switches_train_neg.apply(sum, axis=1)
    switches_test_neg_sum = switches_test_neg.apply(sum, axis=1)
    switches_train_pos_count = switches_train_pos.apply(sum, axis=0)
    switches_test_pos_count = switches_test_pos.apply(sum, axis=0)
    switches_train_neg_count = switches_train_neg.apply(sum, axis=0)
    switches_test_neg_count = switches_test_neg.apply(sum, axis=0)

    plt.clf()
    fig, ax = plt.subplots(2,2)
    ax[0][0].hist(switches_train_pos_sum, bins='fd')
    ax[1][0].hist(switches_test_pos_sum, bins='fd')
    ax[0][1].hist(switches_train_neg_sum, bins='fd')
    ax[1][1].hist(switches_test_neg_sum, bins='fd')
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            ax[i][j].spines['right'].set_visible(False)
            ax[i][j].spines['top'].set_visible(False)
    ax[0][0].set_title('Training examples\nPositive label')
    ax[1][0].set_title('Test examples\nPositive label')
    ax[0][1].set_title('Training examples\nNegative label')
    ax[1][1].set_title('Test examples\nNegative label')
    ax[0][0].set_ylabel('# of examples')
    ax[1][0].set_ylabel('# of examples')
    ax[1][0].set_xlabel('# available data sources')
    ax[1][1].set_xlabel('# available data sources')
    fig.tight_layout()
    plt.savefig(hist_png_file)


if __name__ == '__main__':
    main()