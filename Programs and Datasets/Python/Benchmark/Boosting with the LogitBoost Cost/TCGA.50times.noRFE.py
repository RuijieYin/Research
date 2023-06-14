# use RFA
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
import numpy as np


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
        '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data/train_TCGA_{}.csv'.format(
            k))
    test_data = pd.read_csv(
        '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data/test_TCGA_{}.csv'.format(
            k))

    # convert to gene pairs: X_train_temp, X_validation_temp, X_test
    train_data_converted = convertGenepairs(np.array(train_data))
    test_data_converted = convertGenepairs(np.array(test_data))

    X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
    y_train = train_data_converted.iloc[:, -1]
    X_test = test_data_converted.iloc[:, 0:-1].astype('int64')
    y_test = test_data_converted.iloc[:, -1]

    # Set up the k-fold cross-validation
    kfold = KFold(n_splits=3, shuffle=True, random_state=42)
    # Set up the parameter grid to search over
    param_grid = {'max_features': [0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 'sqrt', 'log2']}
    # Set up the GradientBoostingClassifier
    gbm = GradientBoostingClassifier(n_estimators=100, random_state=1, max_depth=3)
    # Set up the GridSearchCV object
    grid_search = GridSearchCV(gbm, param_grid=param_grid, cv=kfold)
    # Fit the GridSearchCV object to the data
    grid_search.fit(X_train, y_train)

    # fit the model with the best parameter:
    model = GradientBoostingClassifier(n_estimators=100, random_state=1, max_depth=3,
                                       max_features=grid_search.best_params_['max_features']).fit(X_train, y_train)

    accuracy_temp = model.score(X_test, y_test)
    accuracy.append(accuracy_temp)
    print("iteration {} is done".format(k))
    exp_id_temp = k
    exp_id.append(exp_id_temp)

result_table = {"exp_id": exp_id, "Accuracy": accuracy}
result_table_df = pd.DataFrame(result_table)
result_table_df.to_csv(
    '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Results 2/TCGA_50_boosting.csv')

