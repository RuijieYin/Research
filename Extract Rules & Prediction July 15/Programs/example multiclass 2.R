# load treeRule.Obtain.Judge.R
source("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction July 7/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
source("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction July 7/split_extract_predict_multiclass.R")


# load treeRule.Obtain.Judge.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction July 15/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction July 15/split_extract_predict_multiclass.R")



# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))
data.split <- split.train(data = read.csv(file = "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))


# extract rules:
rules <- extract.rules(data = data.split$train.data, n.trees = 400, mtry = 5, node_depth_max = 4)
# see all rules extracted from the n.tree trees
# rules$Rules

gc()
# use validation data to select the optimal number of rules for each subclass
selected.rules <- select.rules(validation.data = data.split$validation.data, rule.table = rules$Rules)
# # see the optimal rules selected for each subclass
write.csv(selected.rules$Selected.Rules, "/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Results/Decision Rules Results July 7/multiclass1_rules.csv")


# Predict on the test data (the remaining 95 samples)
pred <- predict.rules(data = data.split$test.data, decision.rules = selected.rules$Selected.Rules)
# see the predicted class labels for the test data:
# pred$Prediction.results
# see prediction accuracy and cohen's kappa:
print(c("n.trees" = 400, "mtry" = 5,
        "node_depth_max" = 4,
        pred$Results.table))

