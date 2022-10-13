library(data.table)
library(randomForestSRC)
library(caret)
library(ggplot2)
library(dplyr)


# function 1
# split.train
# Split the training data (TCGA BRCA data, 1095 samples) 
# 800 samples are training data (used to extract tree rules and build model)
# 200 samples are validation data, and are used for determining the optimal number of rules
# 95 samples are used as test data

split.train <- function(data = data_input) {
  
  data_input <- data
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
  
  # find the most frequent class label in the training data:
  train_label <- table(as.character(train.80$subtype))
  most_class <- names(train_label[train_label == max(train_label)])
  
  data.split <- list("train.data" = train.80,
                     "validation.data" = train.20,
                     "test.data" = test_data,
                     "most.freq" =  most_class)
  
  return(data.split)
}






# function 2: continue to split the training and validation data obtained from function 1
# using a 1-vs-1 scheme
# example:
# table(data.split$train.data[,51])
# 
# BRCA.Basal   BRCA.Her2   BRCA.LumA   BRCA.LumB BRCA.Normal 
# 141          58         431         160          10 

# split training data:
split_one_vs_one_data_train <- function(data = data.split$train.data) {
  # initialize an indicator: num_count to store each data generated
  num_count = 1
  
  num_classes <- length(unique(data$subtype))
  classes <- unique(data$subtype)
  train_data <- list()
  # for each class, create a one-vs-one dataset, respectively
  for ( i in 1 : (num_classes - 1) ) {
    for ( j in (i + 1): num_classes ) {
      
      data_temp_1 <- data[data$subtype == classes[i],]
      data_temp_2 <- data[data$subtype == classes[j],]
      
      # we have to re-factorize the class labels in data_temp_1 and data_temp_2, to make sure only 2 classes are left
      data_temp_1[ ,ncol(data_temp_1)] <- factor(data_temp_1[ ,ncol(data_temp_1)])
      data_temp_2[ ,ncol(data_temp_2)] <- factor(data_temp_2[ ,ncol(data_temp_2)])
      assign(paste0("data.class.", num_count), rbind(data_temp_1, data_temp_2))
      train_data[[num_count]] <- get(paste0("data.class.", num_count))
      num_count = num_count + 1
    }
  }
  return(train_data)
}

# training_data = split_one_vs_one_data_train(data = data.split$train.data)



# function 3: extract.rules
# model section
# grow a forest for each pair of subclass


# hide warning: replacement element x has x rows to replace x rows, when values(rules) are too long
options(warn = -1)


extract.rules <- function (data = training_data, n.trees = 20, mtry = 5, node_depth_max = 5) {
  # grow n.tree of trees for each class:
  set.seed(1)
  
  # initialize a NULL vector to store all rules obtained
  final_rules_table = list()
  
  # grow forests for a total of length(training_data) times:
  for (c in 1:length(data)) {
    # initialize an indicator: row_counter to store all rules obtained from each tree:
    row_count = 1
    
    data_each_class = as.data.frame(data[[c]])
    tree.collection <- rfsrc(subtype~., data = data_each_class,
                             ntree = n.trees, mtry = mtry, nodedepth = node_depth_max,
                             bootstrap = "by.root", samptype = 'swr',
                             membership = T,
                             # grow class balanced trees
                             case.wt = randomForestSRC:::make.wt(data_each_class$subtype),
                             sampsize = randomForestSRC:::make.size(data_each_class$subtype))
    
    # calculate the total number of rules obtained from the above forest:
    total_num_rules <- length(getTreeRule(tree.collection)$tree.id)
    
    # initialize a data frame to store (rules[[j]], class, else_class, as.numeric(perform.score))
    final_result <- data.frame("Rules" = character(),
                               "Class Label" = character(),
                               "Else Class" = character(),
                               "Performance Score" = numeric(),
                               "Performance Score.1" = numeric(),
                               "Performance Score.2" = numeric())
    
    # extract tree rules from each tree:
    rules <- vector("list", total_num_rules)
    
    for (j in 1:total_num_rules) {
      # extract rule j from tree i:
      rules[[j]] <- parseRule(getTreeRule(tree.collection)$treeRule[[j]])
      # the ith tree where rule j comes from:
      tree.i <- getTreeRule(tree.collection)$tree.id[j]
      
      # index of inbag sample to grow tree i:
      index.inbag <- which(tree.collection$inbag[,tree.i] != 0)
      # find inbag samples suffice the jth rule:
      x <- data_each_class # should tell what x is
      all.sample.fit.rule <- which(eval(parse(text = rules[[j]])))
      inbag.sample.fit.rule <- Reduce(intersect, list(index.inbag, all.sample.fit.rule))
      
      # determine the class label for this rule
      # should consider the tied class: 
      class_label_train <- table(data_each_class[inbag.sample.fit.rule, dim(data_each_class)[2]])
      class <- names(class_label_train[class_label_train == max(class_label_train)])
      
      # see if this rule has an "else": the majority class that does not suffice jth rule in the inbag sample,
      # use majority vote to determine "else":
      inbag.sample.not.fit.rule <- index.inbag[-inbag.sample.fit.rule]
      # find the class labels of samples that does not suffice jth rule;
      # first find and remove the majority class in 'class'label' section:
      majority_class_index <- which(as.character(data_each_class[inbag.sample.not.fit.rule, dim(data_each_class)[2]],class) %in% class)
      else_label_train <- table(as.character(data_each_class[inbag.sample.not.fit.rule, dim(data_each_class)[2]])[-majority_class_index])
      else_class_temp <- names(else_label_train[else_label_train == max(else_label_train)])
      # should consider tied votes:
      if (length(else_class_temp) != 1) {
        else_class <- NA
      } else {
        else_class <- else_class_temp
      }
      
      # calculate the performance score for rule j in tree i using oob data from tree i:
      # get class balanced OOB sample:
      out.of.bag <- which(tree.collection$inbag[,tree.i] == 0)
      oob.sample <- data_each_class[out.of.bag,]
      oob.sample.balanced <- oob.sample[complete.cases(oob.sample), ]
      oob.sample.balanced <- downSample(x = data_each_class[,-dim(data_each_class)[2]],
                                        y = factor(data_each_class[,dim(data_each_class)[2]]))
      colnames(oob.sample.balanced)[dim(oob.sample.balanced)[2]] <- "subtype"
      
      # 2 scenarios: the rule has or does not have an 'else'
      
      # 1st scenario: the rule j does not have an 'else':
      if (is.na(else_class) == T) {
        # store results in the iteration:
        class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                              ncol = 2)
        # 1st element: see if sample suffice the above criteria, 1/0
        # 2nd element: see if sample is indeed the predicted class, 1/0
        
        for (k in 1:dim(oob.sample.balanced)[1]) {
          if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
              eval(parse(text = rules[[j]]))[out.of.bag[k]] == T && 
              oob.sample.balanced[k,]$subtype == class) {
            class_label[k,] <- c(1,1)
          } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                     eval(parse(text = rules[[j]]))[out.of.bag[k]] == T && 
                     oob.sample.balanced[k,]$subtype != class) {
            class_label[k,] <- c(1,0)    
          } else {
            class_label[k,] <- c(0,0)     
          }
        }
        
        # see notes for details on how to calculate performance score
        perform.score.1 <- (sum(class_label[,1])/dim(class_label)[1])*
          (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))
        
        perform.score.2 <- 0
        
        perform.score <- perform.score.1 + perform.score.2
      } 
      # 2nd scenario: the rule j has an 'else':
      else if (is.na(else_class) == F) {
        # part 1:
        # store results in the iteration:
        class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                              ncol = 2)
        # 1st element: see if sample suffice the above criteria, 1/0
        # 2nd element: see if sample is indeed the predicted class, 1/0
        
        for (k in 1:dim(oob.sample.balanced)[1]) {
          if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
              eval(parse(text = rules[[j]]))[out.of.bag[k]] == T && 
              oob.sample.balanced[k,]$subtype == class) {
            class_label[k,] <- c(1,1)
          } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                     eval(parse(text = rules[[j]]))[out.of.bag[k]] == T && 
                     oob.sample.balanced[k,]$subtype != class) {
            class_label[k,] <- c(1,0)    
          } else {
            class_label[k,] <- c(0,0)     
          }
        }
        
        # see notes for details on how to calculate performance score
        perform.score.1 <- (sum(class_label[,1])/dim(class_label)[1])*
          (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))
        
        # part 2:
        # store results in the iteration:
        class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                              ncol = 2)
        # 1st element: see if sample suffice the above criteria, 1/0
        # 2nd element: see if sample is indeed the predicted class, 1/0
        
        for (k in 1:dim(oob.sample.balanced)[1]) {
          if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
              eval(parse(text = rules[[j]]))[out.of.bag[k]] == F && 
              oob.sample.balanced[k,]$subtype == else_class) {
            class_label[k,] <- c(1,1)
          } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                     eval(parse(text = rules[[j]]))[out.of.bag[k]] == F && 
                     oob.sample.balanced[k,]$subtype != else_class) {
            class_label[k,] <- c(1,0)    
          } else {
            class_label[k,] <- c(0,0)     
          }
        }
        
        # see notes for details on how to calculate performance score
        perform.score.2 <- (sum(class_label[,1])/dim(class_label)[1])*
          (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))
        
        perform.score <- perform.score.1 + perform.score.2
      }
      final_result[row_count,] <- list("Rules" = rules[[j]],
                                       "Class Label" = class,
                                       "Else Class" = else_class,
                                       "Performance Score" = perform.score,
                                       "Performance Score.1" = perform.score.1,
                                       "Performance Score.2" = perform.score.2)
      row_count = row_count + 1
      
    }
    final_rules_table[[c]] <- final_result
  }
  
  # first convert final_rules_table ( a list object ) to a data frame
  final_rules_table_df = do.call(rbind, final_rules_table)
  
  # remove rules with NaN performance scores
  final_result_temp_1 <- final_rules_table_df[!is.na(final_rules_table_df[,4]),]
  
  # remove rules with performance scores = 0
  final_result_temp_2 <- final_result_temp_1[final_result_temp_1[,4] != 0,]
  
  # remove duplicated rules by rows (attempt 1)
  final_result_all <- final_result_temp_2[!duplicated(final_result_temp_2),]
  
  # remove duplicated rules by column (attempt 2)
  final_result_all <- final_result_all %>% distinct(final_result_all[,1], .keep_all = TRUE)
  
  # remove the last column: duplicated results
  # the 7th column comes from using %>% in the above step:
  final_result_all <- final_result_all[,-7]
  
  # colnames(final_result_all) <- c("Rules", "Class Label", "Else Class")
  rules.all <- list("Rules" = final_result_all) 
  # function ends here
}


# function 4
# use the remaining 20% training data (validation data) to determine the optimal number of rules:
# whole validation data is used:
# not weighted voting:
# consecutive selection of rules:
library(ggplot2)
library(dplyr)

select.rules <- function(validation.data = data.split$validation.data, rule.table = rules$Rules) {
  
  # select rules for a pair of subclass:
  num_classes <- length(unique(validation.data$subtype))
  classes <- unique(validation.data$subtype)
  
  # initialize a final rule table that will store all rules from all pairs of subclass that are retained:
  rules.retained.list = list()
  
  # initialize a count indicator: to store rules retained in rules.retained.list:
  num_count = 1
  
  for ( i in 1 : (num_classes - 1) ) {
    for ( j in (i + 1): num_classes ) {
      
      rules.pair = subset(rule.table, rule.table$Class.Label == classes[i] & rule.table$Else.Class == classes[j] |
                            rule.table$Class.Label == classes[j] & rule.table$Else.Class == classes[i])
      
      # select validation data for a pair of subclass:
      validation.data.pair = subset(validation.data, validation.data$subtype == classes[i] |
                                      validation.data$subtype == classes[j])
      # re-factorize
      validation.data.pair$subtype = as.factor(validation.data.pair$subtype)
      
      # in validation data: rows are samples and columns are vars
      
      # remove rules with performance score = 0
      # hard coding: the 4th column is the performance scores
      rules.pair.no.0 <- rules.pair[rules.pair[,4] != 0,]
      
      # the probability of predicting correctly by random guessing is 0.5*0.5+0.5*0.5 = 0.5
      # hard coding: the 4th column is the performance scores
      # rules.order.temp <- rules.pair.no.0[rules.pair.no.0[,4] > 0.5,]
      
      # order the rules by their performance scores:
      # hard coding: the 4th column is the performance scores
      rules.order <- rules.pair.no.0[order(rules.pair.no.0[,4], decreasing = T),]
      
      
      # initialize a data frame to store # of rules v.s. accuracy
      accuracy.table <- data.frame("Num.of.Rules" = as.numeric(),
                                   "Accuracy" = as.numeric())
      
      # initialize accuracy
      accuracy <- 0
      
      # initialize rules will be retained:
      rules.retained.new <- data.frame("Rules" = character(), 
                                       "Class Label" = character(), 
                                       "Else Class" = character(), 
                                       "Performance Score" = numeric(),
                                       "Performance Score (Class)" = numeric(), 
                                       "Performance Score (Else Class)" = numeric())
      
      # initialize rules retained in the previous selection step:
      rules.retained.previous <- data.frame("Rules" = character(), 
                                            "Class Label" = character(), 
                                            "Else Class" = character(), 
                                            "Performance Score" = numeric(),
                                            "Performance Score (Class)" = numeric(), 
                                            "Performance Score (Else Class)" = numeric())
      
      
      for (t in 1:nrow(rules.order)) {
        
        if (t > nrow(rules.order)) {
          break
        } else {
          rules.add <- rules.order[t,]
          rules.retained.new <- rbind(rules.retained.previous, rules.add)
          
          pred_class <- matrix(NA, nrow = dim(validation.data.pair)[1], ncol = 1)
          
          for (k in 1:dim(validation.data.pair)[1]) {
            # define x
            x <- validation.data.pair[k,]
            inter.table <- data.frame("label" = as.character(),
                                      "Performance Score" = as.numeric(),
                                      "Indicator" = as.numeric())
            
            for (w in 1:nrow(rules.retained.new))  {
              if (eval(parse(text = rules.retained.new[w,]$Rules)) == T) {
                inter.table[w,] <- list("label" = rules.retained.new[w,][,2], 
                                        "Performance Score" = rules.retained.new[w,][,4], 
                                        "Indicator" = 1)
                # indicator = 1: the sample's label is determined by 'class label'
              } else if (eval(parse(text = rules.retained.new[w,]$Rules)) == F) {
                inter.table[w,] <- list("label" = rules.retained.new[w,][,3], 
                                        "Performance Score" = rules.retained.new[w,][,4], 
                                        "Indicator" = 0)
                # indicator = 0: the sample's label is determined by 'else class'
              }
            }
            
            # determine the class label of the kth sample in validation data:
            label <- table(as.character(inter.table[,1]))
            pred_label <- names(label[label == max(label)])
            if (length(pred_label) != 1) {
              pred_class[k,1] <- as.character(inter.table[,1][which.max(inter.table[,2])])
            } else{
              pred_class[k,1] <- pred_label
            }
          }
          
          
          comp_table <- cbind(pred_class, as.matrix(validation.data.pair$subtype, ncol = 1))
          
          accuracy.temp <- sum(comp_table[,1] == comp_table[,2])/nrow(comp_table)
          
          if (accuracy.temp <= accuracy) {
            accuracy <- accuracy
            rules.retained.previous <- rules.retained.previous
            accuracy.table[t,] <- list("Num.of.Rules" = t,
                                       "Accuracy" = accuracy)
          } else if (accuracy.temp > accuracy) {
            accuracy <- accuracy.temp
            accuracy.table[t,] <- list("Num.of.Rules" = t,
                                       "Accuracy" = accuracy)
            rules.retained.previous <- rules.retained.new
          }
          
        } 
      } 
      rules.retained.list[[num_count]] = rules.retained.previous
      num_count = num_count + 1
    }
  }
  
  # df <- data.frame(x = accuracy.table[,1], y = accuracy.table[,2])
  # plot.list <- ggplot(data = df, aes(x = x, y = y, group=1)) +
  #   geom_line() +
  #   geom_point() +
  #   xlab("ith Decision Rule") +
  #   ylab("Accuracy") +
  #   ggtitle("Accuracy Figure of Consecutive Selection of Rules") +
  #   theme_classic()
  
  
  # # output the final selected rules:
  final.decision.rules = do.call(rbind, rules.retained.list) 
  
  # final number of rules:
  final_num <- nrow(final.decision.rules)
  
  optimal.rules <- list(# "Plot.Num.Accuracy" = plot.list,
    "Selected.Rules" = final.decision.rules,
    "Final.Num" = final_num)
  return(optimal.rules)
} 



# function 5: Prediction 
# test data
# use majority voting (non-weighted) to predict the class labels of new data
#

predict.rules <- function(data = data.split$test.data, decision.rules) {

      pred_class <- matrix(NA, nrow = dim(data)[1], ncol = 1)
      
      for (k in 1:dim(data)[1]) {
        # define x
        x <- data[k,]
        inter.table <- data.frame("label" = as.character(),
                                  "Performance Score" = as.numeric(),
                                  "Indicator" = as.numeric())
        
        for (w in 1:nrow(decision.rules))  {
          if (eval(parse(text = decision.rules[w,]$Rules)) == T) {
            inter.table[w,] <- list("label" = decision.rules[w,][,2], 
                                    "Performance Score" = decision.rules[w,][,4], 
                                    "Indicator" = 1)
            # indicator = 1: the sample's label is determined by 'class label'
          } else if (eval(parse(text = decision.rules[w,]$Rules)) == F) {
            inter.table[w,] <- list("label" = decision.rules[w,][,3], 
                                    "Performance Score" = decision.rules[w,][,4], 
                                    "Indicator" = 0)
            # indicator = 0: the sample's label is determined by 'else class'
          }
        }
        
        # determine the class label of the kth sample in validation data:
        label <- table(as.character(inter.table[,1]))
        pred_label <- names(label[label == max(label)])
        if (length(pred_label) != 1) {
          pred_class[k,1] <- as.character(inter.table[,1][which.max(inter.table[,2])])
        } else{
          pred_class[k,1] <- pred_label
        }
      }

  
  prediction.results <- data.frame("sample no." = 1:nrow(data),
                                   "predicted class labels" = pred_class[,1])
  
  # evaluation
  # see performance: Cohen's Kappa
  # 1st row: reference labels
  # 2nd row: predicted labels
  prediction <- matrix(prediction.results[,2], nrow = 1)
  
  compare_labels <- rbind(matrix(data[,dim(data)[2]], nrow = 1), 
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
  
  prediction.labels.kappa <- list("Prediction.results" = prediction.results,
                                  "Results.table" = con_table)
  return(prediction.labels.kappa)
}