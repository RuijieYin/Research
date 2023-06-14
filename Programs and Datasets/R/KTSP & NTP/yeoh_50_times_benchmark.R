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
filenames <- list.files(pattern="train_yeoh_.*csv")
filenames2 <- list.files(pattern="test_yeoh_.*csv")
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
  assign("train_data_temp", eval(parse(text = paste(text = "train_yeoh_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  # train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test_yeoh_", i, sep = ""))))
  test_data = get("test_data_temp")
  # test_data = test_data[,-1]
  # evaluation
  # KTSP
  # 1 vs 1 scheme
  #
  num_class_1 = length(which(train_data$subtype== "BCR-ABL"))
  num_class_2 = length(which(train_data$subtype== "E2A-PBX1"))
  num_class_3 = length(which(train_data$subtype== "Hyperdiploid>50"))
  num_class_4 = length(which(train_data$subtype== "MLL"))
  num_class_5 = length(which(train_data$subtype== "T-ALL"))
  num_class_6 = length(which(train_data$subtype== "TEL-AML1"))

# KTSP
# 1 vs 1 scheme
#
subtypes1 <- factor(c(rep("BCR-ABL", num_class_1), rep("E2A-PBX1", num_class_2)))
data1 <- as.matrix(cbind(t(train_data[train_data$subtype == "BCR-ABL", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "E2A-PBX1", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier1 <- SWAP.Train.KTSP(data1, subtypes1)

subtypes2 <- factor(c(rep("BCR-ABL", num_class_1), rep("Hyperdiploid>50", num_class_3)))
data2 <- as.matrix(cbind(t(train_data[train_data$subtype == "BCR-ABL", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "Hyperdiploid>50", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier2 <- SWAP.Train.KTSP(data2, subtypes2)

subtypes3 <- factor(c(rep("BCR-ABL", num_class_1), rep("MLL", num_class_4)))
data3 <- as.matrix(cbind(t(train_data[train_data$subtype == "BCR-ABL", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "MLL", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier3 <- SWAP.Train.KTSP(data3, subtypes3)


subtypes4 <- factor(c(rep("BCR-ABL", num_class_1), rep("T-ALL", num_class_5)))
data4 <- as.matrix(cbind(t(train_data[train_data$subtype == "BCR-ABL", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "T-ALL", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier4 <- SWAP.Train.KTSP(data4, subtypes4)


subtypes5 <- factor(c(rep("BCR-ABL", num_class_1), rep("TEL-AML1", num_class_6)))
data5 <- as.matrix(cbind(t(train_data[train_data$subtype == "BCR-ABL", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "TEL-AML1", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier5 <- SWAP.Train.KTSP(data5, subtypes5)


subtypes6 <- factor(c(rep("E2A-PBX1", num_class_2), rep("Hyperdiploid>50", num_class_3)))
data6 <- as.matrix(cbind(t(train_data[train_data$subtype == "E2A-PBX1", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "Hyperdiploid>50", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier6 <- SWAP.Train.KTSP(data6, subtypes6)


subtypes7 <- factor(c(rep("E2A-PBX1", num_class_2), rep("MLL", num_class_4)))
data7 <- as.matrix(cbind(t(train_data[train_data$subtype == "E2A-PBX1", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "MLL", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier7 <- SWAP.Train.KTSP(data7, subtypes7)


subtypes8 <- factor(c(rep("E2A-PBX1", num_class_2), rep("T-ALL", num_class_5)))
data8 <- as.matrix(cbind(t(train_data[train_data$subtype == "E2A-PBX1", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "T-ALL", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier8 <- SWAP.Train.KTSP(data8, subtypes8)

subtypes9 <- factor(c(rep("E2A-PBX1", num_class_2), rep("TEL-AML1", num_class_6)))
data9 <- as.matrix(cbind(t(train_data[train_data$subtype == "E2A-PBX1", -dim(train_data)[2]]),
                         t(train_data[train_data$subtype == "TEL-AML1", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier9 <- SWAP.Train.KTSP(data9, subtypes9)


subtypes10 <- factor(c(rep("Hyperdiploid>50", num_class_3), rep("MLL", num_class_4)))
data10 <- as.matrix(cbind(t(train_data[train_data$subtype == "Hyperdiploid>50", -dim(train_data)[2]]),
                          t(train_data[train_data$subtype == "MLL", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier10 <- SWAP.Train.KTSP(data10, subtypes10)

subtypes11 <- factor(c(rep("Hyperdiploid>50", num_class_3), rep("T-ALL", num_class_5)))
data11 <- as.matrix(cbind(t(train_data[train_data$subtype == "Hyperdiploid>50", -dim(train_data)[2]]),
                          t(train_data[train_data$subtype == "T-ALL", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier11 <- SWAP.Train.KTSP(data11, subtypes11)

subtypes12 <- factor(c(rep("Hyperdiploid>50", num_class_3), rep("TEL-AML1", num_class_6)))
data12 <- as.matrix(cbind(t(train_data[train_data$subtype == "Hyperdiploid>50", -dim(train_data)[2]]),
                          t(train_data[train_data$subtype == "TEL-AML1", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier12 <- SWAP.Train.KTSP(data12, subtypes12)

subtypes13 <- factor(c(rep("MLL", num_class_4), rep("T-ALL", num_class_5)))
data13 <- as.matrix(cbind(t(train_data[train_data$subtype == "MLL", -dim(train_data)[2]]),
                          t(train_data[train_data$subtype == "T-ALL", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier13 <- SWAP.Train.KTSP(data13, subtypes13)

subtypes14 <- factor(c(rep("MLL", num_class_4), rep("TEL-AML1", num_class_6)))
data14 <- as.matrix(cbind(t(train_data[train_data$subtype == "MLL", -dim(train_data)[2]]),
                          t(train_data[train_data$subtype == "TEL-AML1", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier14 <- SWAP.Train.KTSP(data14, subtypes14)

subtypes15 <- factor(c(rep("T-ALL", num_class_5), rep("TEL-AML1", num_class_6)))
data15 <- as.matrix(cbind(t(train_data[train_data$subtype == "T-ALL", -dim(train_data)[2]]),
                          t(train_data[train_data$subtype == "TEL-AML1", -dim(train_data)[2]])))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier15 <- SWAP.Train.KTSP(data15, subtypes15)


# prediction
num_of_classifiers = 15
num_of_classes = 6

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
write.csv(accuracy, "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Results Benchmark/yeoh.results.ktsp.csv", row.names = T)







# NTP
setwd("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data")

#Save File names to a variable
filenames <- list.files(pattern="train_yeoh_.*csv")
filenames2 <- list.files(pattern="test_yeoh_.*csv")
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
  assign("train_data_temp", eval(parse(text = paste(text = "train_yeoh_", i, sep = ""))))
  train_data = get("train_data_temp")
  # remove the first column as row numbers
  # train_data = train_data[,-1]
  assign("test_data_temp", eval(parse(text = paste(text = "test_yeoh_", i, sep = ""))))
  test_data = get("test_data_temp")
  
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
# prediction
num_of_classes = 6
template <- matrix(NA, dim(rank_data_labeled)[2]-1, num_of_classes)
num_of_features = dim(rank_data_labeled)[2] - 1
for (j in 1:num_of_features) {
  # # BCR-ABL        E2A-PBX1 Hyperdiploid>50         MLL           T-ALL        TEL-AML1 
  largest <- which.max(c(mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "BCR-ABL")]),
                         mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "E2A-PBX1")]),
                         mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "Hyperdiploid>50")]),
                         mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "MLL")]),
                         mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "T-ALL")]), 
                         mean(rank_data_labeled[,j][which(rank_data_labeled$subtype == "TEL-AML1")])
                         )
  )
  template[j,largest] <- 1
  template[j,-largest] <- c(0,0,0,0,0)
}

template_BCR <- template[,1]
template_E2A <- template[,2]
template_Hyperdiploid <- template[,3]
template_MLL <- template[,4]
template_TALL <- template[,5]
template_TEL <- template[,6]

# Compute the cosine similarity between the new test sample and each of the templates:
t_test_data_nolabel = t(test_data[,-dim(test_data)[2]])
prediction <- matrix(NA, ncol = ncol(t_test_data_nolabel), nrow = 1)
inner_results <- matrix(NA, nrow = ncol(t_test_data_nolabel), ncol = num_of_classes)
index <- c("BCR-ABL", "E2A-PBX1", "Hyperdiploid>50", "MLL", "T-ALL", "TEL-AML1")
for (u in 1:ncol(t_test_data_nolabel)) {
  test1 <- cosine(t_test_data_nolabel[,u], template_BCR)
  test2 <- cosine(t_test_data_nolabel[,u], template_E2A)
  test3 <- cosine(t_test_data_nolabel[,u], template_Hyperdiploid)
  test4 <- cosine(t_test_data_nolabel[,u], template_MLL)
  test5 <- cosine(t_test_data_nolabel[,u], template_TALL)
  test6 <- cosine(t_test_data_nolabel[,u], template_TEL)
  inner_results[u,] <- c(test1, test2, test3, test4, test5, test6)
  largest_corr <- which.max(c(test1, test2, test3, test4, test5, test6))
  prediction[1,u] <- index[largest_corr]
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
write.csv(accuracy, "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Results Benchmark/yeoh.results.ntp.csv", row.names = T)


