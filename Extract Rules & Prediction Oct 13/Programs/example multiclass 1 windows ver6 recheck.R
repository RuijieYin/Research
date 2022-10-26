# load treeRule.Obtain.Judge.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
# doesn't matter which split_extract_predict_multiclass.R we use
# we will only be using the split.train function to get the test data
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/split_extract_predict_multiclass_ver6_recheck.R")



# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))


# extract rules:
rules <- extract.rules(data = data.split$train.data, n.trees = 200, mtry = 5, node_depth_max = 7)
# see all rules extracted from the n.tree trees
write.csv(rules$Rules, "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Results/Decision Rules Results October 13/rules_extracted_ver6_recheck.csv")


# use validation data to select the optimal number of rules for each subclass
selected.rules <- predict.rules(validation.data = data.split$test.data, rule.table = rules$Rules)
# see the plot of number of top rules v.s. accuracy for each subclass
# e.x.
selected.rules$Plot.Num.Accuracy