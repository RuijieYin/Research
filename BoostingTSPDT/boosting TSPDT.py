import pandas
from sklearn.ensemble import GradientBoostingClassifier

training_data = pandas.read_csv('/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/boost_training.csv')
test_data = pandas.read_csv('/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/boost_test.csv')
X_train = training_data.iloc[:,0:1224]
y_train = training_data.iloc[:,1225]
X_test = test_data.iloc[:,0:1224]
y_test = test_data.iloc[:,1225]
clf = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1,
      max_leaf_nodes=6, random_state=0).fit(X_train, y_train)
print(clf.score(X_test, y_test))

# some results:
# max_depth=2, 0.8987
# max_leaf_nodes=8, 0.9240
# max_leaf_nodes=32, 0.9240
# max_leaf_nodes=6, 0.9367
# n_estimators=1, max_leaf_nodes=6, 0.7215