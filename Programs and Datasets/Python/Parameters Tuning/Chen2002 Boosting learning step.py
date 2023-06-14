# use RFA
import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from numpy import sort
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from shaphypetune import BoostSearch, BoostRFE, BoostRFA, BoostBoruta
import lightgbm as lgb
from sklearn.model_selection import StratifiedShuffleSplit
import re

accuracy = []
exp_id = []


def convertGenepairs(inputData):
    # inputData = np.array(inputData)
    X = inputData
    n = inputData.shape[0]
    d = inputData.shape[1] - 1
    newX = []
    for i in range(d - 1):
        for j in range(i + 1, d):
            newX.append(np.multiply((X[:, i] < X[:, j]), 1))
            # print(' i: {} \n j: {}'.format(i, j))
    newXreshaped = np.reshape(newX, (-1, n)).T

    newX_ = np.concatenate((newXreshaped, np.reshape(X[:, (d)], (n, -1))), axis=1)
    newX_ = pd.DataFrame(newX_)
    # print("newX_ generated")
    # for i in range((newXreshaped.shape[1])):
    #     newX_.rename(columns={i: ("V" + str(i))}, inplace=True)
    newX_.columns = newX_.columns.astype(str)
    # print("column names converted to str")

    newX_.rename(columns={newXreshaped.shape[1]: "subtype"}, inplace=True)

    # saving the dataframe
    # newX_.to_feather(outputData)
    return (newX_)

learning_rate_numbers = [0.005,0.075,0.01,0.125,0.15,0.175,0.20,0.215,0.23,0.25,0.3,0.35,0.4,0.45,0.5]
i=1


train_data = pd.read_csv(
        'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp1/train1.csv')
test_data = pd.read_csv(
        'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp1/test1.csv')
    # train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
    # test_data = test_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
    # convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
test_data_converted = convertGenepairs(np.array(test_data))
# print("training and test data were successfully converted".format(k))

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_test = test_data_converted.iloc[:, 0:-1].astype('int64')
y_test = test_data_converted.iloc[:, -1]

for k in learning_rate_numbers:
    clf = GradientBoostingClassifier(n_estimators=100, learning_rate=k,
                                     max_depth=7, random_state=0).fit(X_train, y_train)
    # print(clf.score(X_test, y_test))
    accuracy_temp = clf.score(X_test, y_test)
    accuracy.append(accuracy_temp)
    print("iteration {} is done".format(i))
    exp_id_temp = i
    exp_id.append(exp_id_temp)
    i = i + 1

result_table = {"exp_id": i, "Accuracy": accuracy}
result_table_df = pd.DataFrame(result_table)
result_table_df.to_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp9 boxplot3/Results/Boosting_learning_rates.csv')
