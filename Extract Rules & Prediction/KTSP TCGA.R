# KTSP TCGA

# both training and test data are from the 1095 TCGA BRCA data
# 1000: training; 95: test

#
require(switchBox)
require(gplots)
library(dplyr)
library(arsenal) #compare two datasets ====> equivalent to proc compare in SAS
library(caret) #confusion matrix
library(epibasix) #Cohen's Kappa
library(data.table) #'%>%'
#
data_input <- read.csv(file = "/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1)
#data_input <- read.csv(file = "F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1)

# class type conversion
names(data_input) <- make.names(names(data_input))
data_input$subtype <- factor(data_input$subtype)

# split brca.normal samples with 20:20
normal.sample <- data_input[data_input$subtype == "BRCA.Normal",]
set.seed(1)
normal.sample.train <- normal.sample[sample(nrow(normal.sample), 20), ]
normal.sample.test <- normal.sample[!rownames(normal.sample) %chin% rownames(normal.sample.train),]

# randomly select 980 samples from other subtypes
other.sample <- data_input[data_input$subtype != "BRCA.Normal",]
set.seed(1)
other.sample.train <- other.sample[sample(nrow(other.sample), 980), ]
other.sample.test <- other.sample[!rownames(other.sample) %chin% rownames(other.sample.train),]

# training data (1000 obs)
train_data <- rbind(normal.sample.train, other.sample.train)
# test data (95 obs)
test_data <- rbind(normal.sample.test, other.sample.test)


### 
# note: rows are genes
# 1 vs 1 scheme
# BRCA.Basal BRCA.Her2 BRCA.LumA BRCA.LumB BRCA.Normal
#
train_data_basal <- t(train_data[train_data$subtype == "BRCA.Basal", -dim(train_data)[2]])
# 175 samples
train_data_her2 <- t(train_data[train_data$subtype == "BRCA.Her2", -dim(train_data)[2]])
# 72 samples
train_data_luma <- t(train_data[train_data$subtype == "BRCA.LumA", -dim(train_data)[2]])
# 529 samples
train_data_lumb <- t(train_data[train_data$subtype == "BRCA.LumB", -dim(train_data)[2]])
# 204 samples
train_data_normal <- t(train_data[train_data$subtype == "BRCA.Normal", -dim(train_data)[2]])
# 20 samples

# model
# basal vs her2
subtypes1 <- factor(c(rep("BRCA.Basal",175),rep("BRCA.Her2",72)))
data1 <- as.matrix(cbind(train_data_basal,
                         train_data_her2))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier1 <- SWAP.Train.KTSP(data1, subtypes1)

# basal vs luma
subtypes2 <- factor(c(rep("BRCA.Basal",175),rep("BRCA.LumA",529)))
data2 <- as.matrix(cbind(train_data_basal,
                         train_data_luma))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier2 <- SWAP.Train.KTSP(data2, subtypes2)

# basal vs lumb
subtypes3 <- factor(c(rep("BRCA.Basal",175),rep("BRCA.LumB",204)))
data3 <- as.matrix(cbind(train_data_basal,
                         train_data_lumb))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier3 <- SWAP.Train.KTSP(data3, subtypes3)

# basal vs normal
subtypes4 <- factor(c(rep("BRCA.Basal",175),rep("BRCA.Normal",20)))
data4 <- as.matrix(cbind(train_data_basal,
                         train_data_normal))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier4 <- SWAP.Train.KTSP(data4, subtypes4)

# her2 vs luma
subtypes5 <- factor(c(rep("BRCA.Her2",72), rep("BRCA.LumA",529)))
data5 <- as.matrix(cbind(train_data_her2,
                         train_data_luma))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier5 <- SWAP.Train.KTSP(data5, subtypes5)

# her2 vs lumb
subtypes6 <- factor(c(rep("BRCA.Her2",72), rep("BRCA.LumB",204)))
data6 <- as.matrix(cbind(train_data_her2,
                         train_data_lumb))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier6 <- SWAP.Train.KTSP(data6, subtypes6)

# her2 vs normal
subtypes7 <- factor(c(rep("BRCA.Her2",72), rep("BRCA.Normal",20)))
data7 <- as.matrix(cbind(train_data_her2,
                         train_data_normal))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier7 <- SWAP.Train.KTSP(data7, subtypes7)

# luma vs lumb 
subtypes8 <- factor(c(rep("BRCA.LumA",529), rep("BRCA.LumB",204)))
data8 <- as.matrix(cbind(train_data_luma,
                         train_data_lumb))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier8 <- SWAP.Train.KTSP(data8, subtypes8)

# luma vs normal
subtypes9 <- factor(c(rep("BRCA.LumA",529), rep("BRCA.Normal",20)))
data9 <- as.matrix(cbind(train_data_luma,
                         train_data_normal))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier9 <- SWAP.Train.KTSP(data9, subtypes9)

# lumb vs normal
subtypes10 <- factor(c(rep("BRCA.LumB",204), rep("BRCA.Normal",20)))
data10 <- as.matrix(cbind(train_data_lumb,
                         train_data_normal))
# training function (SWAP.KTSP.Train:deprecated, use SWAP.Train.KTSP instead) for the 
# classifier and a function (SWAP.KTSP.Classify) 
# for predicting the label of an unseen sample.
classifier10 <- SWAP.Train.KTSP(data10, subtypes10)



# Use 'test_data' for prediction: 
#
test_matrix <- as.matrix(t(test_data[,-dim(test_data)[2]]))
#initialization
pred_result <- matrix(NA,nrow = dim(test_matrix)[2], ncol = 1)
inter_result <- matrix(NA,nrow = 10, ncol = 1)
#impossible to have tied predicted class labels in 1vs1 scheme

for(i in 1:dim(test_matrix)[2]) {
  # hard coding: number of classifiers built
  for(j in 1:10) {
    inter_result[j,1] <- as.character(SWAP.KTSP.Classify(as.matrix(test_matrix[,i], ncol = 1), 
                                                       eval(parse(text = paste("classifier", j, sep = "")))))
  }
  pred_result[i,1] <- names(which.max(table(inter_result)))
}
rownames(pred_result) <- colnames(test_matrix)

#the prediction results are stored in the pred_result
#save to local disk:
#write.table(pred_result,"ktsp.switchBox prediction result.txt")

# evaluation
# see performance: Cohen's Kappa
# 1st row: reference labels
# 2nd row: predicted labels
prediction <- matrix(pred_result[,1], nrow = 1)

compare_labels <- rbind(matrix(test_data[,dim(test_data)[2]], nrow = 1), 
                        prediction)

remove <- which(prediction[1,] == "UNS")

if (length(remove) == 0) {
  compare_labels <- compare_labels
} else if (length(remove) != 0) {
  compare_labels <- compare_labels[,-remove]
}

# Concordance: kappa statistic
# will need a matrix where the diagonal elements of the matrix are 
# the agreeing elements; the discordant observations are on the off-diagonal.
# A confusion matrix:
con_table <- confusionMatrix(data = as.factor(compare_labels[2,]),
                             reference = as.factor(compare_labels[1,]))

# results:
# Confusion Matrix and Statistics
# 
# Reference
# Prediction    BRCA.Basal BRCA.Her2 BRCA.LumA BRCA.LumB BRCA.Normal
# BRCA.Basal          15         0         0         0           1
# BRCA.Her2            0        10         0         3           1
# BRCA.LumA            0         0        27         0           2
# BRCA.LumB            0         0         9        10           0
# BRCA.Normal          0         0         1         0          16
# 
# Overall Statistics
# 
# Accuracy : 0.8211         
# 95% CI : (0.729, 0.8922)
# No Information Rate : 0.3895         
# P-Value [Acc > NIR] : < 2.2e-16      
# 
# Kappa : 0.7688         
# 
# Mcnemar's Test P-Value : NA             
# 
# Statistics by Class:
# 
#                      Class: BRCA.Basal Class: BRCA.Her2 Class: BRCA.LumA
# Sensitivity                     1.0000           1.0000           0.7297
# Specificity                     0.9875           0.9529           0.9655
# Pos Pred Value                  0.9375           0.7143           0.9310
# Neg Pred Value                  1.0000           1.0000           0.8485
# Prevalence                      0.1579           0.1053           0.3895
# Detection Rate                  0.1579           0.1053           0.2842
# Detection Prevalence            0.1684           0.1474           0.3053
# Balanced Accuracy               0.9938           0.9765           0.8476
#                      Class: BRCA.LumB Class: BRCA.Normal
# Sensitivity                    0.7692             0.8000
# Specificity                    0.8902             0.9867
# Pos Pred Value                 0.5263             0.9412
# Neg Pred Value                 0.9605             0.9487
# Prevalence                     0.1368             0.2105
# Detection Rate                 0.1053             0.1684
# Detection Prevalence           0.2000             0.1789
# Balanced Accuracy              0.8297             0.8933


