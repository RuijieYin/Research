import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from numpy import sort
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
import lightgbm as lgb
from sklearn.model_selection import StratifiedShuffleSplit
import re
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV


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


for k in range(1, 21):
    train_data = pd.read_csv(
        '/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp5 boxplot 2/Data/train_khan_{}.csv'.format(k))
    test_data = pd.read_csv(
        '/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp5 boxplot 2/Data/test_khan_{}.csv'.format(k))
    # train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
    # test_data = test_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
    # convert to gene pairs: X_train_temp, X_validation_temp, X_test
    train_data_converted = convertGenepairs(np.array(train_data))
    test_data_converted = convertGenepairs(np.array(test_data))
    print("training and test data {} were successfully converted".format(k))

    X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
    y_train = train_data_converted.iloc[:, -1]
    X_test = test_data_converted.iloc[:, 0:-1].astype('int64')
    y_test = test_data_converted.iloc[:, -1]


    model = RandomForestClassifier()

    # define the parameter grid to search over
    param_grid = {
        'min_samples_leaf': [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5],
        'max_features': ['sqrt', 'log2']
    }

    # create a GridSearchCV object
    grid_search = GridSearchCV(model, param_grid=param_grid)
    # fit the GridSearchCV object to the data
    grid_search.fit(X_train, y_train)

    # get the best parameters found
    best_params = grid_search.best_params_

    # create a new RandomForestClassifier object with the best parameters, and fit the best RandomForestClassifier object to the training data
    best_rfc = RandomForestClassifier(min_samples_leaf=best_params['min_samples_leaf'], max_features=best_params['max_features']).fit(X_train, y_train)

    # make predictions on the test data using the best RandomForestClassifier object
    accuracy_temp = best_rfc.score(X_test, y_test)


    accuracy.append(accuracy_temp)
    print("iteration {} is done".format(k))
    exp_id_temp = k
    exp_id.append(exp_id_temp)

result_table = {"exp_id": exp_id, "Accuracy": accuracy}
result_table_df = pd.DataFrame(result_table)
result_table_df.to_csv(
    '/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp13 boxplot4 multiclass RF/Results/khan_RF_tune_nodesize_mtry.csv')


