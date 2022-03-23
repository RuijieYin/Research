# example
# load treeRule.Obtain.Judge.R
source(".loc/treeRule.Obtain.Judge.R")
# load extract.rule.func.R
source(".loc/extract.rule.func.R")
# load pred.function.R
source(".loc/pred.func.R")
# load eval.kappa.R
source(".loc/eval.kappa.R")



# load training and test data:
train_data <- read.csv(".loc/train_data_example.csv", row.names = 1)
test_data <- read.csv(".loc/test_data_example.csv", row.names = 1)

# use extract.rule.func to extract rules: set n.trees as small number for faster processing;
# returns a list of top 20 retrieved rules for each of the subtypes
extract.rule.func(data_input = train_data, n.trees = 200, mtry = 35, initial_seed = 1, node_depth_max = 3)

# use pred.func to predict on the new data: 
# returns a table of predicted class labels for samples in test data
pred.func(test_input = test_data)

# evaluation of the performance using Cohen's kappa
eval.kappa(test_input = test_data)




