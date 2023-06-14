# use RFA
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from numpy import sort
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from shaphypetune import BoostSearch, BoostRFE, BoostRFA, BoostBoruta
import lightgbm as lgb
from lightgbm import LGBMClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import StratifiedKFold
import math
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score
from sklearn.pipeline import Pipeline



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



train_data = pd.read_csv(
        '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosting_genepairs_TNBC/Data/microarray_TCGA_training.csv')
test_data = pd.read_csv(
        '/Users/jerryruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Boosting_genepairs_TNBC/Data/microarray_TCGA_test.csv')
# train_data = train_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# test_data = test_data.rename(columns=lambda x: re.sub('[^A-Za-z0-9_]+', '', x))
# convert to gene pairs: X_train_temp, X_validation_temp, X_test
train_data_converted = convertGenepairs(np.array(train_data))
test_data_converted = convertGenepairs(np.array(test_data))

X_train = train_data_converted.iloc[:, 0:-1].astype('int64')
y_train = train_data_converted.iloc[:, -1]
X_test = test_data_converted.iloc[:, 0:-1].astype('int64')
y_test = test_data_converted.iloc[:, -1]

# get feature importance scores and eliminate features with 0 importance scores:
model = GradientBoostingClassifier(n_estimators=200, random_state=42, max_depth=3).fit(X_train, y_train)

f_importance_ = sort(model.feature_importances_)
f_importance_greater_than_zero = f_importance_[f_importance_ != 0]
selection = SelectFromModel(model, threshold=min(f_importance_greater_than_zero), prefit=True)

select_X_train = pd.DataFrame(selection.transform(X_train))
select_X_test = pd.DataFrame(selection.transform(X_test))


# step2: fit model with RFA:
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
split_gen = sss.split(select_X_train, y_train)
train_index, validation_index = next(split_gen)

X_train_temp = select_X_train[select_X_train.index.isin(train_index)]
y_train_temp = y_train[y_train.index.isin(train_index)]
X_validation = select_X_train[select_X_train.index.isin(validation_index)]
y_validation = y_train[y_train.index.isin(validation_index)]

rfa = BoostRFA(
    lgb.LGBMClassifier(),
    step=5, n_jobs=-1
)
rfa.fit(X_train_temp, y_train_temp, eval_set=[(X_validation, y_validation)])

accuracy_temp = rfa.score(select_X_test, y_test)
print(accuracy_temp)



