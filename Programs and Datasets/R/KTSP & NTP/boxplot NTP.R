# use 20 folds of data to calculate the accuracy of ntp
library(caret)
library(data.table)
require(switchBox)
library(lsa)


# function 3: NTP
# in train_data above, rows are samples and columns are genes
# need to refactor the class labels
ntp_train_test <- function(train = train_data, test = test_data) {
  train_no_label <- train_data[,-dim(train)[2]]
  train_data_rank <- matrix(NA, nrow(train_no_label),
                            ncol(train_no_label))
  for(i in 1:nrow(train_data_rank)) {
    train_data_rank[i,] <- rank(train_no_label[i,])
  }
  rownames(train_data_rank) <- rownames(train_no_label)
  colnames(train_data_rank) <- colnames(train_no_label)
  # attach class labels:
  rank_data_labeled <- data.frame(train_data_rank, factor(train_data[ ,dim(train)[2]]))
  colnames(rank_data_labeled)[dim(rank_data_labeled)[2]] <- "subtypes"
  
  # Generate the template for each class
  class_labels_1 <- unique(rank_data_labeled[,ncol(rank_data_labeled)])[1]
  class_labels_2 <- unique(rank_data_labeled[,ncol(rank_data_labeled)])[2]
  template <- matrix(NA,dim(rank_data_labeled)[2]-1, 2)
  
  for (i in 1:(ncol(rank_data_labeled) - 1)) {
    largest <- which.max(c(mean(rank_data_labeled[,i][which(rank_data_labeled$subtypes == class_labels_1)]),
                           mean(rank_data_labeled[,i][which(rank_data_labeled$subtypes == class_labels_2)])))
    template[i,largest] <- 1
    template[i,-largest] <- -1
  }
  
  template_class1 <- template[,1]
  template_class2 <- template[,2]
  
  # Compute the cosine similarity between the new test sample and each of the templates:
  test_no_label <- test[,-dim(test)[2]]
  prediction <- matrix(NA, ncol = nrow(test_no_label), nrow = 1)
  corr_collection <- matrix(NA, nrow = ncol(test_no_label), ncol = 2)
  index <- c(class_labels_1, class_labels_2)
  for (i in 1:nrow(test_no_label)) {
    test1 <- as.numeric(cor(as.matrix(t(test_no_label[i,]), nrow = 1), as.matrix(template_class1, nrow = 1), method = "pearson"))
    test2 <- as.numeric(cor(as.matrix(t(test_no_label[i,]), nrow = 1), as.matrix(template_class2, nrow = 1), method = "pearson"))
    corr_collection[i,] <- c(test1,test2)
    inner_result <- c(test1,test2)
    largest_corr <- max(inner_result)
    largest_corr_index <- which.max(inner_result)
    second_corr <- max(inner_result[inner_result != max(inner_result)])
    prediction[1,i] <- index[largest_corr_index]
  }
  
  test_labels <- matrix(NA, nrow = 1, ncol = nrow(test))
  for (j in 1:nrow(test)) {
    if (test[j,dim(test)[2]] == class_labels_1) {
      test_labels[1,j] = 1
    } else{
      test_labels[1,j] = 2
    }
  }
  
  compare_labels <- rbind(test_labels, 
                          prediction)
  # remove <- which(prediction[1,] == "UNS")
  # compare_labels <- compare_labels[,-remove]
  con_table <- confusionMatrix(data = as.factor(compare_labels[2,]),
                               reference = as.factor(compare_labels[1,]))
  results.return <- list("ConfusionMatrix" = con_table)
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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >= 0.92], ncol = 1)
write.csv(accuracy, "data1.results.ntp.csv", row.names = T)


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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >= 0.92], ncol = 1)
write.csv(accuracy, "data3.results.ntp.csv", row.names = T)


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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >= 0.75], ncol = 1)
write.csv(accuracy, "data4.results.ntp.csv", row.names = T)


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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >= 0.5], ncol = 1)
write.csv(accuracy, "data5.results.ntp.csv", row.names = T)


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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >= 0.5], ncol = 1)
write.csv(accuracy, "data6.results.ntp.csv", row.names = T)



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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy == 1], ncol = 1)
write.csv(accuracy, "data7.results.ntp.csv", row.names = T)



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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

# accuracy = as.matrix(accuracy[!accuracy >=0.86], ncol = 1)
write.csv(accuracy, "data11.results.ntp.csv", row.names = T)



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
  accuracy_temp = ntp_train_test(train_data, test_data)$ConfusionMatrix$overall[[1]]
  accuracy[i,1] = accuracy_temp
}

accuracy = as.matrix(accuracy[!accuracy >=0.9], ncol = 1)
write.csv(accuracy, "data12.results.ntp.csv", row.names = T)













