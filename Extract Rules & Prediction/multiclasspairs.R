if (!requireNamespace("multiclassPairs", quietly = TRUE)) {
  install.packages("multiclassPairs")
}

# Install the dependencies from Bioconductor
# BiocManager, Biobase, and switchBox packages from Bioconductor are needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}
if (!requireNamespace("switchBox", quietly = TRUE)) {
  BiocManager::install("switchBox")
}

# load multiclassPairs library
library(multiclassPairs)
library(data.table)
library(caret)

# You have two schemes to train your pair-based classifier:
#   
# First option is a one-vs-rest scheme that assemble one-vs-rest binary classifiers built by ‘switchBox’ package which uses Top-score pairs (TSP) algorithm.
# The second option is a scheme based on a novel implementation of the random forest (RF) algorithm.

# Both begin by filtering the features, then combining the filtered features to make the list of all the possible rules (i.e. rule1: feature1 < feature2, 
# rule2: feature1 < feature3, etc…). Then the list of rules will be filtered and the most important and informative rules will be kept. 
# The informative rules will be assembled in an one-vs-rest model or in an RF model. 


# split data:
data_input <- read.csv(file = "/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/BRCA PAM50/tcga_quantiles_50.csv", header = T, row.names = 1)
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

# use 80% of the training data to create rules and the rest 20% to determine the optimal number of rules
normal.20.sample <- train_data[train_data$subtype == "BRCA.Normal",]
set.seed(1)
normal.10.train.80 <- normal.20.sample[sample(nrow(normal.20.sample), 10), ]
normal.10.train.20 <- normal.20.sample[!rownames(normal.20.sample) %chin% rownames(normal.10.train.80),]

other.train.80 <- train_data[train_data$subtype != "BRCA.Normal",]
set.seed(1)
other.790.train.80 <- other.train.80[sample(nrow(other.train.80), 790),]
other.190.train.80 <- other.train.80[!rownames(other.train.80) %chin% rownames(other.790.train.80),]

# 80% training data for creating rules:
train.80 <- rbind(normal.10.train.80, other.790.train.80)
# 20% training data to determine the optimal number of rules
train.20 <- rbind(normal.10.train.20, other.190.train.80)


# build object: data
# rows are genes and columns are samples
Data <- t.data.frame(train.80[,-51])

# class labels
L1 <- train.80[,51]

# platform/study labels
P1 <- c(rep("Affy", 800))

# create the data object
object <- ReadData(Data = Data,
                   Labels = L1,
                   Platform = P1,
                   verbose = FALSE)


# Gene filtering
# For building a pair-based classifier with a one-vs-rest scheme, we start by selecting top differentially expressed genes using the filter_genes_TSP function. 
# This reduces the number of gene combinations (rules) in the next steps. 
# This function can perform the filtering in different ways and return the top differential expressed genes for each class.

filtered_genes <- filter_genes_TSP(data_object = object,
                                   filter = "one_vs_one",
                                   platform_wise = FALSE,
                                   featureNo = 50,
                                   UpDown = TRUE,
                                   verbose = TRUE)
# > filtered_genes
# sorted genes for One-vs-rest Scheme:
#   Object contains:
#   - filtered_genes 
# - class: BRCA.Normal : 27 genes
# - class: BRCA.LumA : 10 genes
# - class: BRCA.Basal : 37 genes
# - class: BRCA.Her2 : 13 genes
# - class: BRCA.LumB : 13 genes
# 
# - calls 
# filter_genes_TSP(data_object = object,
#                  filter = "one_vs_one",
#                  platform_wise = FALSE,
#                  featureNo = 50,
#                  UpDown = TRUE,
#                  verbose = TRUE)

# train our model
classifier <- train_one_vs_rest_TSP(data_object = object,
                                    filtered_genes = filtered_genes,
                                    k_range = 5:50,
                                    include_pivot = FALSE,
                                    one_vs_one_scores = TRUE,
                                    platform_wise_scores = FALSE,
                                    seed = 1,
                                    verbose = FALSE)
# > classifier
# multiclassPairs - One vs rest Scheme
# *Classifier:
#   contains binary classifiers:
#   - Class: BRCA.Normal ... 10 rules
# - Class: BRCA.LumA ... 5 rules
# - Class: BRCA.Basal ... 5 rules
# - Class: BRCA.Her2 ... 6 rules
# - Class: BRCA.LumB ... 5 rules

# Prediction: predict_one_vs_rest_TSP
# apply on the testing data

test <- as.matrix(t.data.frame(test_data[,-51]))

results_test_1 <- predict_one_vs_rest_TSP(classifier = classifier,
                                        Data = test,
                                        tolerate_missed_genes = TRUE,
                                        weighted_votes = FALSE,
                                       
                                        verbose = TRUE)

# weighted version
results_test_2 <- predict_one_vs_rest_TSP(classifier = classifier,
                                          Data = test,
                                          tolerate_missed_genes = TRUE,
                                          weighted_votes = TRUE,
                                          
                                          verbose = TRUE)


# Confusion Matrix and Statistics on training data
confusionMatrix(data = factor(results_test_1$max_score),
                       reference = factor(test_data[,51]))

# Confusion Matrix and Statistics
# 
# Reference
# Prediction    BRCA.Basal BRCA.Her2 BRCA.LumA BRCA.LumB BRCA.Normal
# BRCA.Basal          14         0         0         0           0
# BRCA.Her2            0        10         2         1           3
# BRCA.LumA            0         0        26         2           6
# BRCA.LumB            0         0         7        10           0
# BRCA.Normal          1         0         2         0          11
# 
# Overall Statistics
# 
# Accuracy : 0.7474          
# 95% CI : (0.6478, 0.8309)
# No Information Rate : 0.3895          
# P-Value [Acc > NIR] : 1.495e-12       
# 
# Kappa : 0.6694          
# 
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: BRCA.Basal Class: BRCA.Her2 Class: BRCA.LumA Class: BRCA.LumB Class: BRCA.Normal
# Sensitivity                     0.9333           1.0000           0.7027           0.7692             0.5500
# Specificity                     1.0000           0.9294           0.8621           0.9146             0.9600
# Pos Pred Value                  1.0000           0.6250           0.7647           0.5882             0.7857
# Neg Pred Value                  0.9877           1.0000           0.8197           0.9615             0.8889
# Prevalence                      0.1579           0.1053           0.3895           0.1368             0.2105
# Detection Rate                  0.1474           0.1053           0.2737           0.1053             0.1158
# Detection Prevalence            0.1474           0.1684           0.3579           0.1789             0.1474
# Balanced Accuracy               0.9667           0.9647           0.7824           0.8419             0.7550


# RF scheme
# In Random Forest (RF) scheme, all steps of gene filtering/sorting, 
# rule filtering/sorting, and final model training are performed using the RF algorithm.

# (500 trees here just for fast example)
genes_RF <- sort_genes_RF(data_object = object,
                          # featureNo_altogether, it is better not to specify a number here
                          # featureNo_one_vs_rest, it is better not to specify a number here
                          rank_data = TRUE,
                          platform_wise = FALSE,
                          num.trees = 10000, # more features, more tress are recommended
                          seed=1, # for reproducibility
                          verbose = TRUE)
genes_RF # sorted genes object

# After we sorted the genes, we need to take the top genes and combine them as binary rules and then sort these rules.
# to get an idea of how many genes we will use
# and how many rules will be generated
summary_genes <- summary_genes_RF(sorted_genes_RF = genes_RF,
                                  genes_altogether = c(10,20,50,100,150,200),
                                  genes_one_vs_rest = c(10,20,50,100,150,200))
knitr::kable(summary_genes)

# Now let's run sort_rules_RF to create the rules and sort them
rules_RF <- sort_rules_RF(data_object = object, 
                          sorted_genes_RF = genes_RF,
                          genes_altogether = 50,
                          genes_one_vs_rest = 50, 
                          num.trees = 10000,# more rules, more tress are recommended 
                          seed=1,
                          verbose = TRUE)
rules_RF # sorted rules object

# Now, we have the rules sorted based on their importance. Now we can train our final RF model.
# go with the default settings in the train_RF function directly

# train the final model
# it is preferred to increase the number of trees and rules in case you have
# large number of samples and features
RF_classifier <- train_RF(data_object = object,
                          sorted_rules_RF = rules_RF,
                          #gene_repetition = 1,
                          #rules_altogether = 10,
                          #rules_one_vs_rest = 10,
                          run_boruta = TRUE, 
                          plot_boruta = FALSE,
                          probability = TRUE,
                          num.trees = 1000,
                          boruta_args = list(),
                          verbose = TRUE)

# apply on test data
results <- predict_RF(classifier = RF_classifier, 
                      Data = test,
                      impute = TRUE) # can handle missed genes by imputation

# get the prediction labels
# if the classifier trained using probability   = FALSE
test_pred <- results$predictions
if (is.factor(test_pred)) {
  x <- as.character(test_pred)
}

# if the classifier trained using probability   = TRUE
if (is.matrix(test_pred)) {
  x <- colnames(test_pred)[max.col(test_pred)]
}

# training accuracy
caret::confusionMatrix(data = factor(x),
                       reference = factor(test_data[,51]))

# > caret::confusionMatrix(data = factor(x),
#                          +                        reference = factor(test_data[,51]))
# Confusion Matrix and Statistics
# 
# Reference
# Prediction    BRCA.Basal BRCA.Her2 BRCA.LumA BRCA.LumB BRCA.Normal
# BRCA.Basal          15         0         0         0           3
# BRCA.Her2            0         9         0         0           3
# BRCA.LumA            0         0        36         3          10
# BRCA.LumB            0         1         1        10           0
# BRCA.Normal          0         0         0         0           4
# 
# Overall Statistics
# 
# Accuracy : 0.7789          
# 95% CI : (0.6822, 0.8577)
# No Information Rate : 0.3895          
# P-Value [Acc > NIR] : 1.163e-14       
# 
# Kappa : 0.6971          
# 
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: BRCA.Basal Class: BRCA.Her2 Class: BRCA.LumA Class: BRCA.LumB Class: BRCA.Normal
# Sensitivity                     1.0000          0.90000           0.9730           0.7692            0.20000
# Specificity                     0.9625          0.96471           0.7759           0.9756            1.00000
# Pos Pred Value                  0.8333          0.75000           0.7347           0.8333            1.00000
# Neg Pred Value                  1.0000          0.98795           0.9783           0.9639            0.82418
# Prevalence                      0.1579          0.10526           0.3895           0.1368            0.21053
# Detection Rate                  0.1579          0.09474           0.3789           0.1053            0.04211
# Detection Prevalence            0.1895          0.12632           0.5158           0.1263            0.04211
# Balanced Accuracy               0.9812          0.93235           0.8744           0.8724            0.60000
# 


