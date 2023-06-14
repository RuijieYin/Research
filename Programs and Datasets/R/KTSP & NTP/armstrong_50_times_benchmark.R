# perform benchmark testing of KTSP and NTP
require(switchBox)
require(gplots)
library(dplyr)
library(arsenal) #compare two datasets ====> equivalent to proc compare in SAS
library(caret) #confusion matrix
library(epibasix) #Cohen's Kappa
library(data.table) #'%>%'
library(lsa) #cosine

setwd("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data")

#Save File names to a variable
filenames <- list.files(pattern="train_Armstrong_.*csv")
filenames2 <- list.files(pattern="test_Armstrong_.*csv")
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


accuracy = matrix(NA, 50, 1)
for(i in 1:50) {
  assign("train_data_temp", eval(parse(text = paste(text = "train_Armstrong_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  # train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test_Armstrong_", i, sep = ""))))
  test_data = get("test_data_temp")
  # test_data = test_data[,-1]
  # evaluation
  # KTSP
  # 1 vs 1 scheme
  #
  num_class_1 = length(which(train_data$subtype== "ALL"))
  num_class_2 = length(which(train_data$subtype== "MLL"))
  num_class_3 = length(which(train_data$subtype== "AML"))

  
  subtypes1 <- factor(c(rep("ALL", num_class_1), rep("MLL", num_class_2)))
  data1 <- as.matrix(cbind(t(train_data[train_data$subtype== "ALL", -dim(train_data)[2]]),
                           t(train_data[train_data$subtype== "MLL", -dim(train_data)[2]])))
  # training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
  # classifier and a function (SWAP.KTSP.Classify) 
  # for predicting the label of an unseen sample.
  classifier1 <- SWAP.Train.KTSP(data1, subtypes1)
  
  subtypes2 <- factor(c(rep("ALL", num_class_1), rep("AML", num_class_3)))
  data2 <- as.matrix(cbind(t(train_data[train_data$subtype== "ALL", -dim(train_data)[2]]),
                           t(train_data[train_data$subtype== "AML", -dim(train_data)[2]])))
  # training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
  # classifier and a function (SWAP.KTSP.Classify) 
  # for predicting the label of an unseen sample.
  classifier2 <- SWAP.Train.KTSP(data2, subtypes2)
  
  subtypes3 <- factor(c(rep("MLL", num_class_2), rep("AML", num_class_3)))
  data3 <- as.matrix(cbind(t(train_data[train_data$subtype== "MLL", -dim(train_data)[2]]),
                           t(train_data[train_data$subtype== "AML", -dim(train_data)[2]])))
  # training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
  # classifier and a function (SWAP.KTSP.Classify) 
  # for predicting the label of an unseen sample.
  classifier3 <- SWAP.Train.KTSP(data3, subtypes3)

  
  # prediction
  num_of_classifiers = 3
  num_of_classes = 3
  
  # columns represent the samples and the rows represent the features
  test_matrix <- as.matrix(t(test_data[,-dim(test_data)[2]]))
  #initialization
  pred_result <- matrix(NA,nrow = dim(test_matrix)[2],ncol=1)
  inter_result <- matrix(NA,nrow = num_of_classifiers, ncol=1)
  #impossible to have tied predicted class labels in 1vs1 scheme
  
  
  for(k in 1:dim(test_matrix)[2]) {
    for(j in 1:num_of_classifiers) {
      inter_result[j,1] <- as.character(SWAP.KTSP.Classify(as.matrix(test_matrix[,k], ncol = 1), 
                                                           eval(parse(text = paste("classifier", j, sep = "")))))
    }
    pred_result[k,1] <- names(which.max(table(inter_result)))
  }
  rownames(pred_result) <- colnames(test_matrix)
  
  #the prediction results are stored in the pred_result
  #save to local disk:
  #write.table(pred_result,"ktsp.switchBox prediction result.txt")
  
  
  # results 1:
  # Computing prediction accuracy:
  # transpose pred_result 
  pred_result_t <- t(pred_result)
  compare_labels <- rbind(matrix(test_data$subtype, nrow=1), 
                          pred_result_t)
  #overall prediction accuracy:
  count <- 0
  for (j in 1:dim(compare_labels)[2]) {
    if(compare_labels[1,j] == compare_labels[2,j]){
      count <- count + 1
    } else{
      count <- count
    }
  }
  prediction_accuracy <- count / dim(compare_labels)[2]
  accuracy_temp = prediction_accuracy
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy, ncol = 1)
write.csv(accuracy, "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Results Benchmark/Armstrong.results.ktsp.csv", row.names = T)






# NTP
setwd("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data")

#Save File names to a variable
filenames <- list.files(pattern="train_Armstrong_.*csv")
filenames2 <- list.files(pattern="test_Armstrong_.*csv")
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


accuracy = matrix(NA, 50, 1)

for(i in 1:50){
  assign("train_data_temp", eval(parse(text = paste(text = "train_Armstrong_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  # train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test_Armstrong_", i, sep = ""))))
  test_data = get("test_data_temp")
# NTP
t_train_data_nolabel = t(train_data[,-dim(train_data)[2]])


# Transform the expression values into ranks:
train_data_rank <- matrix(NA, nrow(t_train_data_nolabel),
                              ncol(t_train_data_nolabel))
for(k in 1:ncol(train_data_rank)) {
  train_data_rank[,k] <- rank(t_train_data_nolabel[,k])
}
rownames(train_data_rank) <- rownames(t_train_data_nolabel)
colnames(train_data_rank) <- colnames(t_train_data_nolabel)
# range(train_data_rank[,1])
# create a new data frame with class labels
rank_data_labeled <- data.frame(t(train_data_rank), matrix(train_data[,dim(train_data)[2]], ncol=1))
colnames(rank_data_labeled)[dim(rank_data_labeled)[2]] <- "subtype"

# Generate the template for each subtype (using average ranking): 
num_of_classes = 3
template <- matrix(NA, dim(rank_data_labeled)[2]-1, num_of_classes)
num_of_features = dim(rank_data_labeled)[2] - 1
for (j in 1:num_of_features) {
  largest <- which.max(c(mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "ALL")]),
                       mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "AML")]),
                       mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "MLL")])))
  template[j,largest] <- 1
  template[j,-largest] <- c(0,0)
}

template_ALL <- template[,1]
template_AML <- template[,2]
template_MLL <- template[,3]


# Compute the cosine similarity between the new test sample and each of the templates:
t_test_data_nolabel = t(test_data[,-dim(test_data)[2]])
prediction <- matrix(NA, ncol = ncol(t_test_data_nolabel), nrow = 1)
inner_results <- matrix(NA, nrow = ncol(t_test_data_nolabel), ncol = num_of_classes)
index <- c("ALL","AML","MLL")
for (t in 1:ncol(t_test_data_nolabel)) {
  test1 <- cosine(t_test_data_nolabel[,t], template_ALL)
  test2 <- cosine(t_test_data_nolabel[,t], template_AML)
  test3 <- cosine(t_test_data_nolabel[,t], template_MLL)
  inner_results[t,] <- c(test1,test2,test3)
  largest_corr <- which.max(c(test1,test2,test3))
  prediction[1,t] <- index[largest_corr]
}


# Prediction results:

# 1st row: reference labels
# 2nd row: predicted labels
compare_labels <- rbind(test_data[,dim(test_data)[2]], prediction)
#result has been saved to local drive:
#write.csv(compare_labels,"G:/Dropbox/UM Biostatistics/Research/1st Project/
#          TNBCtype_research/compare_labels1.1.csv")

# Concordance: kappa statistic
# will need a matrix where the diagonal elements of the matrix are 
# the agreeing elements; the discordant observations are on the off-diagonal.
# A confusion matrix:
con_table <- confusionMatrix(data = as.factor(compare_labels[2,]),
                             reference = as.factor(compare_labels[1,]))
accuracy_temp = con_table$overall[[1]]
accuracy[i,1] = accuracy_temp
}
accuracy = as.matrix(accuracy, ncol = 1)
write.csv(accuracy, "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Results Benchmark/Armstrong.results.ntp.csv", row.names = T)





