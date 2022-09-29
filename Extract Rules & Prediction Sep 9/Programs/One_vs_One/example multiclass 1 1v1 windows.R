# load treeRule.Obtain.Judge.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction August 18/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction September 9/split_extract_predict_multiclass_balanced_trees_1v1.R")


# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))
# continue to split the training data obtained from function 1 using a 1-vs-rest scheme
training_data <- split_one_vs_one_data_train(data = data.split$train.data)


# extract rules:
rules <- extract.rules(data = training_data, n.trees = 20, mtry = 5, node_depth_max = 5)
# see all rules extracted from the n.tree trees
# rules$Rules
write.csv(rules$Rules, "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction September 9/rules_extracted_1v1_setup1.csv")


gc()
# use validation data to select the optimal number of rules for each subclass
selected.rules <- select.rules(validation.data = data.split$validation.data, rule.table = rules$Rules)
# # see the optimal rules selected for each subclass
write.csv(selected.rules$Selected.Rules, "F:/Dropbox/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction September 9/multiclass1_rules_1v1_setup1.csv")


# Predict on the test data (the remaining 95 samples)
pred <- predict.rules(data = data.split$test.data, decision.rules = selected.rules$Selected.Rules)
# see the predicted class labels for the test data:
# pred$Prediction.results
# see prediction accuracy and cohen's kappa:
print(c("n.trees" = 20, "mtry" = 5,
        "node_depth_max" = 5,
        pred$Results.table))



