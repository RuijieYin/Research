import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from numpy import sort
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from shaphypetune import BoostRFA
import lightgbm as lgb
import re
import pickle

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

train_data = pd.read_csv(
    '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosting_genepairs_TNBC/Model with RFA/bl1_bl2_train.csv')
validation_data = pd.read_csv(
    '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosting_genepairs_TNBC/Model with RFA/bl1_bl2_validation.csv')
test_data = pd.read_csv(
    '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosting_genepairs_TNBC/Model with RFA/bl1_bl2_test.csv')


train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
validation_data = validation_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
test_data = test_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))

# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
validation_data_converted = convertGenepairs(np.array(validation_data))
test_data_converted = convertGenepairs(np.array(test_data))
print("dataset is ready")

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_validation = validation_data_converted.iloc[:, 0:-1].astype('int64')
y_validation = validation_data_converted.iloc[:, -1]
X_test = test_data_converted.iloc[:, 0:-1].astype('int64')
y_test = test_data_converted.iloc[:, -1]


# fit model on all training data
model = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1,
                                   max_leaf_nodes=6, random_state=0).fit(X_train, y_train)

# retain those gene pairs that have non-zero importance scores:
f_importance_ = sort(model.feature_importances_)
f_importance_greater_than_zero = f_importance_[f_importance_ != 0]
selection = SelectFromModel(model, threshold=min(f_importance_greater_than_zero), prefit=True)
select_X_train = selection.transform(X_train)
select_X_validation = selection.transform(X_validation)
select_X_test = selection.transform(X_test)

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=3, min_features_to_select=1, n_jobs=-1
)
rfa.fit(select_X_train, y_train, eval_set=[(select_X_validation, y_validation)])
# check accuracy on validation data:
# rfa.score(select_X_validation, y_validation)
print("model is ready")

# check accuracy on test data:
rfa.score(select_X_test, y_test)


# accuracy on validation data:
# rfa.score(X_validation, y_validation)
filename = '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosting_genepairs_TNBC/Model with RFA/bl1_vs_bl2.sav'
pickle.dump(rfa, open(filename, 'wb'))



