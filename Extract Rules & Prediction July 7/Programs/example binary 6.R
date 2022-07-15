# binary classification:
# LumA v.s. LumB

# load treeRule.Obtain.Judge.R
source("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction July 7/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
source("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction July 7/split_extract_predict_binary.R")



# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))

# extract rules:
rules <- extract.rules(data = data.split$train.data, n.trees = 200, mtry = 30, node_depth_max = 8)


# use validation data to select the optimal number of rules for each subclass
selected.rules <- select.rules(validation.data = data.split$validation.data, rule.table = rules$Rules)
# see the plot of number of top rules v.s. accuracy for each subclass
# e.x.
selected.rules$Plot.Num.Accuracy
# # see the optimal rules selected for each subclass
write.csv(selected.rules$Selected.Rules, "/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Results/Decision Rules Results July 7/binary6_rules.csv")

# # the top number of rules to be used:
# selected.rules$Optimal.Num

# Predict on the test data (the remaining 95 samples)
pred <- predict.rules(data = data.split$test.data, decision.rules = selected.rules$Selected.Rules)

print(c("n.trees" = 200, "mtry" = 30,
        "node_depth_max" = 8,
        pred$Results.table))
