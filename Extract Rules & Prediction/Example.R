# load treeRule.Obtain.Judge.R
source("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/treeRule.Obtain.Judge.R")
# load split.extract.predict.R
source("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction Apr4/split.extract.predict.R")



# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))
# # view training data
data.split$train.data
# # view validation data
data.split$validation.data
# # view test data
data.split$test.data

# extract rules:
rules <- extract.rules(data = data.split$train.data, n.trees = 50, mtry = 50, node_depth_max = 3)
# see all rules extracted from the n.tree trees
rules$Rules

# use validation data to select the optimal number of rules for each subclass
selected.rules <- select.rules(validation.data = data.split$validation.data, rule.table = rules$Rules)
# see the plot of number of top rules v.s. accuracy for each subclass
# e.x.
selected.rules$Plot.Num.Accuracy[[1]]
# see the optimal rules selected for each subclass
selected.rules$Selected.Rules


# Predict on the test data (the remaining 95 samples)
pred <- predict.rules(data = data.split$test.data, decision.rules = selected.rules$Selected.Rules)
# see the predicted class labels for the test data:
pred$Prediction.results
# see prediction accuracy and cohen's kappa:
pred$Results.table


# Predict on the test data using weighted voting (the remaining 95 samples)
pred.weighted <- predict.rules.weighted(data = data.split$test.data, decision.rules = selected.rules$Selected.Rules)
# see the predicted class labels for the test data:
pred.weighted$Prediction.results
# see prediction accuracy and cohen's kappa:
pred.weighted$Results.table


