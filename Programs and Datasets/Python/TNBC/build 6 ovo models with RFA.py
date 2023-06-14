import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from numpy import sort
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from shaphypetune import BoostSearch, BoostRFE, BoostRFA, BoostBoruta
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
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl1_bl2_train.csv')
validation_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl1_bl2_validation.csv')

train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
validation_data = validation_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
validation_data_converted = convertGenepairs(np.array(validation_data))
print("dataset is ready")

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_validation = validation_data_converted.iloc[:, 0:-1].astype('int64')
y_validation = validation_data_converted.iloc[:, -1]

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=3, min_features_to_select=1, n_jobs=-1
)
rfa.fit(X_train, y_train, eval_set=[(X_validation, y_validation)])
print("model is ready")
# accuracy on validation data:
# rfa.score(X_validation, y_validation)
filename = 'bl1_vs_bl2.sav'
pickle.dump(rfa, open(filename, 'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/'))



train_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl1_lar_train.csv')
validation_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl1_lar_validation.csv')

train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
validation_data = validation_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
validation_data_converted = convertGenepairs(np.array(validation_data))
print("dataset is ready")

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_validation = validation_data_converted.iloc[:, 0:-1].astype('int64')
y_validation = validation_data_converted.iloc[:, -1]

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=3, min_features_to_select=1, n_jobs=-1
)
rfa.fit(X_train, y_train, eval_set=[(X_validation, y_validation)])
# accuracy on validation data:
rfa.score(X_validation, y_validation)
filename = 'bl1_vs_lar.sav'
pickle.dump(rfa, open(filename, 'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/'))




train_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl1_M_train.csv')
validation_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl1_M_validation.csv')

train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
validation_data = validation_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
validation_data_converted = convertGenepairs(np.array(validation_data))
print("dataset is ready")

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_validation = validation_data_converted.iloc[:, 0:-1].astype('int64')
y_validation = validation_data_converted.iloc[:, -1]

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=3, min_features_to_select=1, n_jobs=-1
)
rfa.fit(X_train, y_train, eval_set=[(X_validation, y_validation)])
# accuracy on validation data:
rfa.score(X_validation, y_validation)
filename = 'bl1_vs_M.sav'
pickle.dump(rfa, open(filename, 'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/'))



train_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl2_lar_train.csv')
validation_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl2_lar_validation.csv')

train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
validation_data = validation_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
validation_data_converted = convertGenepairs(np.array(validation_data))
print("dataset is ready")

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_validation = validation_data_converted.iloc[:, 0:-1].astype('int64')
y_validation = validation_data_converted.iloc[:, -1]

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=3, min_features_to_select=1, n_jobs=-1
)
rfa.fit(X_train, y_train, eval_set=[(X_validation, y_validation)])
# accuracy on validation data:
rfa.score(X_validation, y_validation)
filename = 'bl2_vs_lar.sav'
pickle.dump(rfa, open(filename, 'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/'))



train_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl2_M_train.csv')
validation_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/bl2_M_validation.csv')

train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
validation_data = validation_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
validation_data_converted = convertGenepairs(np.array(validation_data))
print("dataset is ready")

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_validation = validation_data_converted.iloc[:, 0:-1].astype('int64')
y_validation = validation_data_converted.iloc[:, -1]

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=3, min_features_to_select=1, n_jobs=-1
)
rfa.fit(X_train, y_train, eval_set=[(X_validation, y_validation)])
# accuracy on validation data:
rfa.score(X_validation, y_validation)
filename = 'bl2_vs_M.sav'
pickle.dump(rfa, open(filename, 'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/'))



train_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/lar_M_train.csv')
validation_data = pd.read_csv(
    'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/lar_M_validation.csv')

train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
validation_data = validation_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
validation_data_converted = convertGenepairs(np.array(validation_data))
print("dataset is ready")

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_validation = validation_data_converted.iloc[:, 0:-1].astype('int64')
y_validation = validation_data_converted.iloc[:, -1]

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=3, min_features_to_select=1, n_jobs=-1
)
rfa.fit(X_train, y_train, eval_set=[(X_validation, y_validation)])
# accuracy on validation data:
rfa.score(X_validation, y_validation)
filename = 'lar_vs_M.sav'
pickle.dump(rfa, open(filename, 'F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosted.genepairs TNBC/Model with RFA/'))
