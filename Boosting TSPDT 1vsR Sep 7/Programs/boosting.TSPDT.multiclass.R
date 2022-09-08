library(gbm)
library(caret)
library(data.table)

# Boosting TSPDT functions
# 1-vs-rest scheme, multi-class classification
# function 1:
# transform the training and test data into gene-pair binary data
# the last column should be the class types
trans.tcga.brca.data <- function (data) {
  # class type conversion
  names(data) <- make.names(names(data))
  colnames(data)[dim(data)[2]] <- "subtype"
  data$subtype <- factor(data$subtype)
  
  # split brca.normal samples with 20:20
  normal.sample <- data[data$subtype == "BRCA.Normal",]
  set.seed(1)
  normal.sample.train <- normal.sample[sample(nrow(normal.sample), 20), ]
  normal.sample.test <- normal.sample[!rownames(normal.sample) %chin% rownames(normal.sample.train),]
  
  # randomly select 980 samples other than BRCA.Normal
  other.sample <- data[data$subtype != "BRCA.Normal",]
  set.seed(1)
  other.sample.train <- other.sample[sample(nrow(other.sample), 980), ]
  other.sample.test <- other.sample[!rownames(other.sample) %chin% rownames(other.sample.train),]
  
  # training data (1000 obs)
  train_data <- rbind.data.frame(normal.sample.train, other.sample.train)

  test.data <- rbind.data.frame(normal.sample.test, other.sample.test)
  
  data_luma <- train_data[train_data$subtype == "BRCA.LumA",]
  data_luma[,51] <- factor(data_luma[,51])
  data_not_luma <- train_data[train_data$subtype != "BRCA.LumA",]
  # should re-factorize
  data_not_luma[,51] <- factor(c(rep("other", nrow(data_not_luma))))
  data_luma_train <- rbind.data.frame(data_luma, data_not_luma)
  
  
  data_lumb <- train_data[train_data$subtype == "BRCA.LumB",]
  data_lumb[,51] <- factor(data_lumb[,51])
  data_not_lumb <- train_data[train_data$subtype != "BRCA.LumB",]
  # should re-factorize
  data_not_lumb[,51] <- factor(c(rep("other", nrow(data_not_lumb))))
  data_lumb_train <- rbind.data.frame(data_lumb, data_not_lumb)
  
  
  data_normal <- train_data[train_data$subtype == "BRCA.Normal",]
  data_normal[,51] <- factor(data_normal[,51])
  data_not_normal <- train_data[train_data$subtype != "BRCA.Normal",]
  # should re-factorize
  data_not_normal[,51] <- factor(c(rep("other", nrow(data_not_normal))))
  data_normal_train <- rbind.data.frame(data_normal, data_not_normal)
  
  
  data_her2 <- train_data[train_data$subtype == "BRCA.Her2",]
  data_her2[,51] <- factor(data_her2[,51])
  data_not_her2 <- train_data[train_data$subtype != "BRCA.Her2",]
  # should re-factorize
  data_not_her2[,51] <- factor(c(rep("other", nrow(data_not_her2))))
  data_her2_train <- rbind.data.frame(data_her2, data_not_her2)
  
  
  data_basal <- train_data[train_data$subtype == "BRCA.Basal",]
  data_basal[,51] <- factor(data_basal[,51])
  data_not_basal <- train_data[train_data$subtype != "BRCA.Basal",]
  # should re-factorize
  data_not_basal[,51] <- factor(c(rep("other", nrow(data_not_basal))))
  data_basal_train <- rbind.data.frame(data_basal, data_not_basal)
  
  
  # transform the training data into gene-pairs
  # rows are the samples and columns are genes:
  # note: columns should be genes, rows are samples
  
  # get gene pair data:
  X = data_luma_train
  # X has 100 rows, 20 columns
  n = dim(X)[1] # number of rows
  d = dim(X)[2] - 1 # number of columns
  newX = NULL
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      newX = cbind(newX, as.numeric(X[, i] < X[, j]))
    }
  }
  
  data_luma_train_boosting = cbind.data.frame(newX, data_luma_train[,51])
  
  
  # get gene pair data:
  X = data_lumb_train
  # X has 100 rows, 20 columns
  n = dim(X)[1] # number of rows
  d = dim(X)[2] - 1 # number of columns
  newX = NULL
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      newX = cbind(newX, as.numeric(X[, i] < X[, j]))
    }
  }
  
  data_lumb_train_boosting = cbind.data.frame(newX, data_lumb_train[,51])
  
  
  # get gene pair data:
  X = data_normal_train
  # X has 100 rows, 20 columns
  n = dim(X)[1] # number of rows
  d = dim(X)[2] - 1 # number of columns
  newX = NULL
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      newX = cbind(newX, as.numeric(X[, i] < X[, j]))
    }
  }
  
  data_normal_train_boosting = cbind.data.frame(newX, data_normal_train[,51])
  
  
  # get gene pair data:
  X = data_her2_train
  # X has 100 rows, 20 columns
  n = dim(X)[1] # number of rows
  d = dim(X)[2] - 1 # number of columns
  newX = NULL
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      newX = cbind(newX, as.numeric(X[, i] < X[, j]))
    }
  }
  
  data_her2_train_boosting = cbind.data.frame(newX, data_her2_train[,51])
  
  
  # get gene pair data:
  X = data_basal_train
  # X has 100 rows, 20 columns
  n = dim(X)[1] # number of rows
  d = dim(X)[2] - 1 # number of columns
  newX = NULL
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      newX = cbind(newX, as.numeric(X[, i] < X[, j]))
    }
  }
  
  data_basal_train_boosting = cbind.data.frame(newX, data_basal_train[,51])
  
  
  # transform test data:
  # get gene pair data:
  X = test.data
  # X has 100 rows, 20 columns
  n = dim(X)[1] # number of rows
  d = dim(X)[2] - 1 # number of columns
  newX = NULL
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      newX = cbind(newX, as.numeric(X[, i] < X[, j]))
    }
  }
  
  data_test_boosting = cbind.data.frame(newX, test.data[,51])
  
  # add column name to subtype
  colnames(data_luma_train_boosting)[dim(data_luma_train_boosting)[2]] <- "subtype"
  colnames(data_lumb_train_boosting)[dim(data_lumb_train_boosting)[2]] <- "subtype"
  colnames(data_normal_train_boosting)[dim(data_normal_train_boosting)[2]] <- "subtype"
  colnames(data_her2_train_boosting)[dim(data_her2_train_boosting)[2]] <- "subtype"
  colnames(data_basal_train_boosting)[dim(data_basal_train_boosting)[2]] <- "subtype"
  colnames(data_test_boosting)[dim(data_test_boosting)[2]] <- "subtype"
  
  # re-factorize
  data_luma_train_boosting[,dim(data_luma_train_boosting)[2]] <- factor(data_luma_train_boosting[,dim(data_luma_train_boosting)[2]])
  data_lumb_train_boosting[,dim(data_lumb_train_boosting)[2]] <- factor(data_lumb_train_boosting[,dim(data_lumb_train_boosting)[2]])
  data_normal_train_boosting[,dim(data_normal_train_boosting)[2]] <- factor(data_normal_train_boosting[,dim(data_normal_train_boosting)[2]])
  data_her2_train_boosting[,dim(data_her2_train_boosting)[2]] <- factor(data_her2_train_boosting[,dim(data_her2_train_boosting)[2]])
  data_basal_train_boosting[,dim(data_basal_train_boosting)[2]] <- factor(data_basal_train_boosting[,dim(data_basal_train_boosting)[2]])
  data_test_boosting[,dim(data_test_boosting)[2]] <- factor(data_test_boosting[,dim(data_test_boosting)[2]])
  
  final.training.data = list(data_luma_train_boosting, 
                              data_lumb_train_boosting,
                              data_normal_train_boosting,
                              data_her2_train_boosting,
                              data_basal_train_boosting)
  
  train.test.data = list("training.data" = final.training.data,
                         "test.data" = data_test_boosting)
  
  return(train.test.data)
}



# function 2:
# build model for each class:

multiclass.model.build <- function (data = data.model$training.data, learining_rate = 0.1, max_node_depth = 6, n_trees = 100) {
  # load the training and test data:
  data_luma_train_boosting = data[[1]]
  data_lumb_train_boosting = data[[2]]
  data_normal_train_boosting = data[[3]]
  data_her2_train_boosting = data[[4]]
  data_basal_train_boosting = data[[5]]
  
  # hide warnings: variable has no variation.
  options(warn = -1)
  
  # build models for each class label:
  luma_model = gbm(subtype ~.,
                   data = data_luma_train_boosting,
                   distribution = "multinomial",
                   shrinkage = learining_rate,
                   interaction.depth = max_node_depth,
                   n.trees = n_trees)
  
  # build models for each class label:
  lumb_model = gbm(subtype ~.,
                   data = data_lumb_train_boosting,
                   distribution = "multinomial",
                   shrinkage = learining_rate,
                   interaction.depth = max_node_depth,
                   n.trees = n_trees)
  
  # build models for each class label:
  normal_model = gbm(subtype ~.,
                     data = data_normal_train_boosting,
                     distribution = "multinomial",
                     shrinkage = learining_rate,
                     interaction.depth = max_node_depth,
                     n.trees = n_trees)
  
  # build models for each class label:
  her2_model = gbm(subtype ~.,
                   data = data_her2_train_boosting,
                   distribution = "multinomial",
                   shrinkage = learining_rate,
                   interaction.depth = max_node_depth,
                   n.trees = n_trees)
  
  # build models for each class label:
  basal_model = gbm(subtype ~.,
                    data = data_basal_train_boosting,
                    distribution = "multinomial",
                    shrinkage = learining_rate,
                    interaction.depth = max_node_depth,
                    n.trees = n_trees)
  
  model.list = list(luma_model, lumb_model, normal_model, her2_model, basal_model)
  
  return(model.list)
}
  

# function 3:
model.predict <- function(data = data.model$test.data, models.trained = models, n_trees = 100) {
  # predict the test data:
  luma_pred = predict.gbm(object = models.trained[[1]],
                          newdata = data,
                          n.trees = n_trees,
                          type = "response")
  
  lumb_pred = predict.gbm(object = models.trained[[2]],
                          newdata = data,
                          n.trees = n_trees,
                          type = "response")
  
  normal_pred = predict.gbm(object = models.trained[[3]],
                            newdata = data,
                            n.trees = n_trees,
                            type = "response")
  
  her2_pred = predict.gbm(object = models.trained[[4]],
                          newdata = data,
                          n.trees = n_trees,
                          type = "response")
  
  basal_pred = predict.gbm(object = models.trained[[5]],
                           newdata = data,
                           n.trees = n_trees,
                           type = "response")
  
  # get prediction results:
  labels_1 = colnames(luma_pred)[apply(luma_pred, 1, which.max)]
  labels_2 = colnames(lumb_pred)[apply(lumb_pred, 1, which.max)]
  labels_3 = colnames(normal_pred)[apply(normal_pred, 1, which.max)]
  labels_4 = colnames(her2_pred)[apply(her2_pred, 1, which.max)]
  labels_5 = colnames(basal_pred)[apply(basal_pred, 1, which.max)]
  # merge the predicted labels
  labels_all = cbind.data.frame(labels_1, labels_2, labels_3, labels_4, labels_5)
  
  
  predictions = matrix(NA, nrow = nrow(labels_all), ncol = 1)
  for(i in 1:nrow(predictions)) {
    labels.obtained = unlist(unique(as.vector(labels_all[i,])))
    labels.decision = labels.obtained[!labels.obtained %in% "other"]
    if ( length(labels.decision) == 0 ) {
      predictions[i,1] = "UNS"
    } else if ( length(labels.decision) > 1 ) {
      predictions[i,1] = "UNS"
    } else {
      predictions[i,1] = labels.decision
    }
  }
  
  # check the confusion matrix.
  remove_UNS = which(predictions == "UNS")
  print(paste(length(remove_UNS), "samples are not determined by the model", sep = " "))
  final_prediction = predictions[-remove_UNS,]
  test_labels = as.character(data[-remove_UNS, dim(data)[2]])
  compare_labels = rbind.data.frame(test_labels, final_prediction) 
  
  cm = confusionMatrix(as.factor(compare_labels[1,]), as.factor(compare_labels[2,]))
  
  pred.results = list("Predicted.labels" = predictions,
                 "Confusion.Matrix" = cm)
  
  return(pred.results)
  
}


