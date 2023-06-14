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
from sklearn.ensemble import RandomForestClassifier


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

    newX_ = np.concatenate((newXreshaped, np.reshape(X[:, (d)],(n,-1))), axis=1)
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


for k in range(1, 51):
    train_data = pd.read_csv(
        '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data/train_yeoh_{}.csv'.format(
            k))
    test_data = pd.read_csv(
        '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data/test_yeoh_{}.csv'.format(
            k))
    train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
    test_data = test_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
    # convert to gene pairs: X_train_temp, X_validation_temp, X_test
    train_data_converted = convertGenepairs(np.array(train_data))
    test_data_converted = convertGenepairs(np.array(test_data))

    X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
    y_train = train_data_converted.iloc[:, -1]
    X_test = test_data_converted.iloc[:, 0:-1].astype('int64')
    y_test = test_data_converted.iloc[:, -1]

    print("dataset is ready")

    # define the model
    model = RandomForestClassifier()
    # fit the model on the whole dataset
    model.fit(X_train, y_train)
    # make a prediction
    accuracy_temp = model.score(X_test, y_test)

    accuracy.append(accuracy_temp)
    print("iteration {} is done".format(k))
    exp_id_temp = k
    exp_id.append(exp_id_temp)

result_table = {"exp_id": exp_id, "Accuracy": accuracy}
result_table_df = pd.DataFrame(result_table)
result_table_df.to_csv(
    '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Results RF/yeoh.results.RF.csv')

