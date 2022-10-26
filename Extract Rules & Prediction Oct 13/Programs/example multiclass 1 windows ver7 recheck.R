# load treeRule.Obtain.Judge.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/treeRule_Obtain_Judge.R")
# load split.extract.predict.R
source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Extract Rules & Prediction October 13/split_extract_predict_multiclass_ver6_recheck.R")



# split the training data (TCGA BRCA 1095 samples: 800 training, 200 validation, 95 test)
data.split <- split.train(data = read.csv(file = "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))

rule.table = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Results/Decision Rules Results October 13/rules_extracted_exp3.csv")
# remove row numbers
rule.table = rule.table[,-1]

# use validation data to select the optimal number of rules for each subclass
selected.rules <- predict.rules(validation.data = data.split$test.data, rule.table = rule.table)
# see the plot of number of top rules v.s. accuracy for each subclass
# e.x.
selected.rules$Plot.Num.Accuracy
# # see the optimal rules selected for each subclass
