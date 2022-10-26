# load treeRule.Obtain.Judge.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
# doesn't matter which split_extract_predict_multiclass.R we use
# we will only be using the split.train function to get the test data
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/split_extract_predict_multiclass.R")



# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))


decision.rules = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Results/Decision Rules Results October 13/rules_extracted_exp4.csv")

# # the top number of rules to be used:
# selected.rules$Optimal.Num

# Predict on the test data (the remaining 95 samples)
pred <- predict.rules(data = data.split$test.data, decision.rules = decision.rules)
# see the predicted class labels for the test data:
# pred$Prediction.results
# see prediction accuracy and cohen's kappa:
print(c("n.trees" = 100, "mtry" = 5,
        "node_depth_max" = 7,
        pred$Results.table))