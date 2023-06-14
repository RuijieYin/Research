# use 20 folds of data to calculate the accuracy of KTSP
library(caret)
library(data.table)
require(switchBox)
library(lsa)

# function 2: KTSP
# train KTSP classifiers on train_data and test on test_data
# in KTSP, rows are genes and columns are samples
# in train_data above, rows are samples and columns are genes
ktsp_train_test <- function(train = train_data, test = test_data) {
  train_data_no_label <- t(train[, -dim(train)[2]])
  subtypes <- factor(train[, dim(train)[2]])
  classifier <- SWAP.Train.KTSP(train_data_no_label, subtypes)
  num_of_genes_selected <- length(classifier$TSPs)
  test_matrix <- as.matrix(t(test[,-dim(test)[2]]))
  pred_result <- matrix(NA, nrow = dim(test_matrix)[2], ncol = 1)
  for(i in 1:dim(test_matrix)[2]) {
    pred_result[i,1] <- as.character(SWAP.KTSP.Classify(as.matrix(test_matrix[,i], ncol = 1), 
                                                        classifier))
  }
  rownames(pred_result) <- colnames(test_matrix)
  prediction <- matrix(pred_result[,1], nrow = 1)
  
  compare_labels <- rbind(matrix(test_data[,dim(test_data)[2]], nrow = 1), 
                          prediction)
  con_table <- confusionMatrix(data = factor(compare_labels[2,]),
                               reference = factor(compare_labels[1,]))
  
  results.return <- list("ConfusionMatrix" = con_table,
                         "Num.of.Genes.Selected" = num_of_genes_selected)
}


setwd("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp4 boxplot")

# data1
#Save File names to a variable
filenames <- list.files(pattern="train1_.*csv")
filenames2 <- list.files(pattern="test1_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train1_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test1_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy[!accuracy >= 0.92], ncol = 1)
write.csv(accuracy, "data1.results.ktsp.csv", row.names = T)




# data3
#Save File names to a variable
filenames <- list.files(pattern="train3_.*csv")
filenames2 <- list.files(pattern="test3_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train3_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test3_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >= 0.92], ncol = 1)
write.csv(accuracy, "data3.results.ktsp.csv", row.names = T)


# data4
#Save File names to a variable
filenames <- list.files(pattern="train4_.*csv")
filenames2 <- list.files(pattern="test4_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train4_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test4_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy[!accuracy >= 0.75], ncol = 1)
write.csv(accuracy, "data4.results.ktsp.csv", row.names = T)


# data5
#Save File names to a variable
filenames <- list.files(pattern="train5_.*csv")
filenames2 <- list.files(pattern="test5_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train5_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test5_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy[!accuracy >= 0.5], ncol = 1)
write.csv(accuracy, "data5.results.ktsp.csv", row.names = T)


# data6
#Save File names to a variable
filenames <- list.files(pattern="train6_.*csv")
filenames2 <- list.files(pattern="test6_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train6_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test6_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >= 0.5], ncol = 1)
write.csv(accuracy, "data6.results.ktsp.csv", row.names = T)



# data7
#Save File names to a variable
filenames <- list.files(pattern="train7_.*csv")
filenames2 <- list.files(pattern="test7_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train7_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test7_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy[!accuracy == 1], ncol = 1)
write.csv(accuracy, "data7.results.ktsp.csv", row.names = T)



# data11
#Save File names to a variable
filenames <- list.files(pattern="train11_.*csv")
filenames2 <- list.files(pattern="test11_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train11_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test11_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy[!accuracy >=0.86], ncol = 1)
write.csv(accuracy, "data11.results.ktsp.csv", row.names = T)



# data12
#Save File names to a variable
filenames <- list.files(pattern="train12_.*csv")
filenames2 <- list.files(pattern="test12_.*csv")
# Get names without ".CSV" and store in "names"
names <- gsub(".csv","",filenames)
names2 <- gsub(".csv","",filenames2)
# read in all training and test data:
for (i in names) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}
for (i in names2) {
  filepath <- file.path(paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, header=TRUE))
}

accuracy = matrix(NA, 20, 1)
for(i in 1:20) {
  assign("train_data_temp", eval(parse(text = paste(text = "train12_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test12_", i, sep = ""))))
  test_data = get("test_data_temp")
  test_data = test_data[,-1]
  # evaluation
  accuracy_temp = ktsp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy[!accuracy >=0.9], ncol = 1)
write.csv(accuracy, "data12.results.ktsp.csv", row.names = T)

