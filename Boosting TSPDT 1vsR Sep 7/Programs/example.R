load(boosting.TSPDT.multiclass.R)

# load the data, split training and test data
data.model = trans.tcga.brca.data(data = read.csv("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1))
# training data
data.model$training.data
# test data
data.model$test.data

# build models 
models = multiclass.model.build(data = data.model$training.data)

# predict the test data
predictions = model.predict(data = data.model$test.data)
predictions$Confusion.Matrix