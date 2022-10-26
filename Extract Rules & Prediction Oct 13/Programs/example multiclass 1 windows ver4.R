# load treeRule.Obtain.Judge.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/split_extract_predict_multiclass ver4.R")



# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))


# see all rules extracted from the n.tree trees
rule.table = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Results/Decision Rules Results October 13/rules_extracted_exp3.csv")

# rules$Rules

# use validation data to select the optimal number of rules for each subclass
selected.rules <- predict.rules(validation.data = data.split$validation.data, rule.table = rule.table)
# see the plot of number of top rules v.s. accuracy for each subclass
# e.x.
selected.rules$Plot.Num.Accuracy
# # see the optimal rules selected for each subclass
# write.csv(selected.rules$Selected.Rules, "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Results/Decision Rules Results October 13/multiclass1_rules_selected_exp5.csv")

# # the top number of rules to be used:
# selected.rules$Optimal.Num

# Predict on the test data (the remaining 95 samples)
pred <- predict.rules(data = data.split$test.data, decision.rules = selected.rules$Selected.Rules)
# see the predicted class labels for the test data:
# pred$Prediction.results
# see prediction accuracy and cohen's kappa:
print(c("n.trees" = 300, "mtry" = 5,
        "node_depth_max" = 7,
        pred$Results.table))

