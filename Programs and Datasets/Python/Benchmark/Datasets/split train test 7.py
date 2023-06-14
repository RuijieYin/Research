import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit

input_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp4 boxplot/data7.csv')

X = input_data.iloc[:, 0:-2]
y = input_data.iloc[:, -1]

for i in range(0, 21):
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.3, random_state=i)
    split_gen = sss.split(X, y)
    train_index, test_index = next(split_gen)

    X_train_temp=X[X.index.isin(train_index)]
    y_train_temp=y[y.index.isin(train_index)]
    # merge into 1 dataset
    train_temp = pd.concat([X_train_temp, y_train_temp], axis=1)

    X_test_temp=X[X.index.isin(test_index)]
    y_test_temp=y[y.index.isin(test_index)]
    # merge into 1 dataset
    test_temp = pd.concat([X_test_temp, y_test_temp], axis=1)

    # save to csv:
    train_temp.to_csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp4 boxplot/train7_{}.csv".format(i), index=False)
    test_temp.to_csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp4 boxplot/test7_{}.csv".format(i), index=False)

