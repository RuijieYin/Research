# function 1
# split.train
# Split the training data (considering only BRCA.LumA and BRCA.LumB)  

library(data.table)


split.train <- function(data = data_input) {
    
    data_input <- data

    # use the most common 2 classes as an example: BRCA.LumA   BRCA.LumB
    # BRCA.Basal   BRCA.Her2   BRCA.LumA   BRCA.LumB BRCA.Normal 
    # 190          82         566         217          40 
    
    # 90%: training data (80% training; 20% validation); the rest 10%: test data
    luma.data <- data_input[data_input$subtype == "BRCA.LumA",]
    set.seed(1)
    luma.data.train.validation <- luma.data[sample(nrow(luma.data), 509), ]
    set.seed(1)
    luma.data.training <- luma.data.train.validation[sample(nrow(luma.data.train.validation), 407), ]
    luma.data.validation <- luma.data.train.validation[!rownames(luma.data.train.validation) %chin% rownames(luma.data.training),]
    luma.data.test <- luma.data[!rownames(luma.data) %chin% rownames(luma.data.train.validation),]
    
    lumb.data <- data_input[data_input$subtype == "BRCA.LumB",]
    set.seed(1)
    lumb.data.train.validation <- lumb.data[sample(nrow(lumb.data), 195), ]
    set.seed(1)
    lumb.data.training <- lumb.data.train.validation[sample(nrow(lumb.data.train.validation), 156), ]
    lumb.data.validation <- lumb.data.train.validation[!rownames(lumb.data.train.validation) %chin% rownames(lumb.data.training),]
    lumb.data.test <- lumb.data[!rownames(lumb.data) %chin% rownames(lumb.data.train.validation),]
    
    train.80 <- rbind(luma.data.training, lumb.data.training)
    train.20 <- rbind(luma.data.validation, lumb.data.validation)
    test_data <- rbind(luma.data.test,lumb.data.test)
    
    # class type conversion
    # re-factorize
    names(train.80) <- make.names(names(train.80))
    train.80$subtype <- factor(train.80$subtype)
    names(train.20) <- make.names(names(train.20))
    train.20$subtype <- factor(train.20$subtype)
    names(test_data) <- make.names(names(test_data))
    test_data$subtype <- factor(test_data$subtype)
    
    
    
    data.split <- list("train.data" = train.80,
                       "validation.data" = train.20,
                       "test.data" = test_data)
    
    return(data.split)
  }



# function 2: extract.rules
# model section
library(data.table)

# hide warning: replacement element x has x rows to replace x rows, when values(rules) are too long
options(warn = -1)


extract.rules <- function (data = data.split$train.data, n.trees = 500, mtry = 50, node_depth_max = 3) {
  # grow n.tree of trees
  set.seed(1)
  tree.collection <- rfsrc(subtype~., data = data,
                           ntree = n.trees, mtry = mtry, nodedepth = node_depth_max,
                           bootstrap = "by.root", samptype = 'swr',
                           membership = T,
                           # grow class balanced trees
                           # case.wt = randomForestSRC:::make.wt(data$subtype),
                           sampsize = randomForestSRC:::make.size(data$subtype))
  
  # summary(tree.collection)
  
  # calculate the total number of rules obtained from all trees
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
    x <- data # should tell what x is
    all.sample.fit.rule <- which(eval(parse(text = rules[[j]])))
    inbag.sample.fit.rule <- Reduce(intersect, list(index.inbag, all.sample.fit.rule))
    
    # determine the class label for this rule
    # should consider the tied class: 
    class_label_train <- table(data[inbag.sample.fit.rule, dim(data)[2]])
    class <- names(class_label_train[class_label_train == max(class_label_train)])
    
    # see if this rule has an "else": the majority class that does not suffice jth rule in the inbag sample,
    # use majority vote to determine "else":
    inbag.sample.not.fit.rule <- index.inbag[-inbag.sample.fit.rule]
    # find the class labels of samples that does not suffice jth rule;
    # first find and remove the majority class in 'class'label' section:
    majority_class_index <- which(as.character(data[inbag.sample.not.fit.rule, dim(data)[2]],class) %in% class)
    else_label_train <- table(as.character(data[inbag.sample.not.fit.rule, dim(data)[2]])[-majority_class_index])
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
    oob.sample <- data[out.of.bag,]
    oob.sample.balanced <- oob.sample[complete.cases(oob.sample), ]
    #oob.sample.balanced <- downSample(x = data[,-dim(data)[2]],
    # y = factor(data[,dim(data)[2]]))
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

final_result[j,] <- list("Rules" = rules[[j]],
                         "Class Label" = class,
                         "Else Class" = else_class,
                         "Performance Score" = perform.score,
                         "Performance Score.1" = perform.score.1,
                         "Performance Score.2" = perform.score.2)

  }
  # remove rules with NaN performance scores
  final_result_temp_1 <- final_result[!is.na(final_result[,4]),]
  
  # remove rules with performance scores = 0
  final_result_temp_2 <- final_result_temp_1[final_result_temp_1[,4] != 0,]
  
  # remove duplicated rules
  final_result_all <- final_result_temp_2[!duplicated(final_result_temp_2),]
  
  # remove duplicated rules by column (attempt 2)
  final_result_all <- final_result_all %>% distinct(final_result_all[,1], .keep_all = TRUE)
  
  colnames(final_result_all) <- c("Rules", "Class Label", "Else Class", "Performance Score",
                                  "Performance Score (Class)", "Performance Score (Else Class)")
  rules.all <- list("Rules" = final_result_all) 
  # function ends here
}


# function 3
# use the remaining 20% training data (validation data) to determine the optimal number of rules:
# whole validation data is used:
# not weighted voting:
library(ggplot2)
library(dplyr)

select.rules <- function(validation.data = data.split$validation.data, rule.table=rules$Rules) {
  
  # in validation data: rows are samples and columns are vars
  
  # remove rules with performance score = 0
  rule.table <- rule.table[rule.table$`Performance Score` != 0,]
  
  # the probability of predicting correctly by random guessing is 0.5*0.5+0.5*0.5 = 0.5
  rules.order.temp <- rule.table[rule.table$`Performance Score` > 0.5,]
  
  # order the rules by their performance scores:
  rules.order <- rules.order.temp[order(rules.order.temp$`Performance Score`, decreasing = T),]
  
  # initialize a data frame to store # of rules v.s. accuracy
  accuracy.table <- data.frame("Num.of.Rules" = as.numeric(),
                               "Accuracy" = as.numeric())
  
  # initialize accuracy
  accuracy <- 0
  
  for (t in 1:nrow(rules.order)) {
    # i: top i rules
    # i has to be odd numbers to break ties
    if (t == 1) {
      i <- t
    } else {
      i <- 2*t-1
    }
    
    if (i > nrow(rules.order)) {
      break
    } else {
      rules.top.i <- rules.order[1:i,]
      
      pred_class <- matrix(NA, nrow = dim(validation.data)[1], ncol = 1)
      # 1st element: the number of top t rules used
      # 2nd element: the accuracy for the number of top t rules used
      
      for (k in 1:dim(validation.data)[1]) {
        # define x
        x <- validation.data[k,]
        inter.table <- data.frame("label" = as.character(),
                                  "Performance Score" = as.numeric(),
                                  "Indicator" = as.numeric())
        
        for (w in 1:nrow(rules.top.i))  {
          if (eval(parse(text = rules.top.i[w,]$Rules)) == T) {
            inter.table[w,] <- list("label" = rules.top.i[w,][,2], 
                                    "Performance Score" = rules.top.i[w,][,4], 
                                    "Indicator" = 1)
            # indicator = 1: the sample's label is determined by 'class label'
          } else if (eval(parse(text = rules.top.i[w,]$Rules)) == F) {
            inter.table[w,] <- list("label" = rules.top.i[w,][,3], 
                                    "Performance Score" = rules.top.i[w,][,4], 
                                    "Indicator" = 0)
            # indicator = 0: the sample's label is determined by 'else class'
          }
        }
        
        # determine the class label of the kth sample in validation data:
        if (sum(inter.table[,3]) != 0) {
          # means the sample's label is determined by 'class label'
          label <- table(as.character(inter.table[,1][which(inter.table[,3] == 1)]))
          if (length(label) != 1) {
            pred_class[k,1] <- as.character(inter.table[,1][which.max(inter.table[,3])])
          } else{
            pred_class[k,1] <- names(label[label == max(label)])
          }
        } # means the sample's label is determined by 'else class label'
        else if(sum(inter.table[,3]) == 0) {
          label <- table(as.character(inter.table[,1][which(inter.table[,3] == 0)]))
          if (length(label) != 1) {
            pred_class[k,1] <- as.character(inter.table[,1][which.max(inter.table[,3])])
          } else{
            pred_class[k,1] <- names(label[label == max(label)])
          }
        }
      }
      
      
      comp_table <- cbind(pred_class, as.matrix(validation.data$subtype, ncol = 1))
      
      accuracy.temp <- sum(comp_table[,1] == comp_table[,2])/nrow(comp_table)
      
      if (accuracy.temp <= accuracy) {
        break
      } else if (accuracy.temp > accuracy) {
        accuracy <- accuracy.temp
        accuracy.table[(i+1)/2,] <- list("Num.of.Rules" = i,
                                         "Accuracy" = accuracy)
      }
      
    } 
  }
  
  df <- data.frame(x = accuracy.table[,1], y = accuracy.table[,2])
  plot.list <- ggplot(data = df, aes(x = x, y = y, group=1)) + 
    geom_line() + 
    geom_point() +
    xlab("Number of Rules") +
    ylab("Accuracy") +
    ggtitle("Number of Rules v.s. Accuracy") +
    theme_classic()
  
  
  # # determine the optimal # of rules:
  optimal_num <- accuracy.table[,1][which.max(accuracy.table[,2])]
  final.decision.rules <- rules.order[1:optimal_num,]
  
  optimal.rules <- list("Plot.Num.Accuracy" = plot.list,
                        "Optimal.Num" = optimal_num,
                        "Selected.Rules" = final.decision.rules)
  return(optimal.rules)
  
} 



# function 4: Prediction 
# test data
# use majority voting (non-wighted) to predict the class labels of new data
#
library(caret)

predict.rules <- function(data = data.split$test.data, decision.rules) {
  # store predicted labels from each rule in a list
  decisions <- vector("list", length = nrow(decision.rules)) 
  # make a final decision for each sample in the test data
  final.decisions <- matrix(NA, nrow = nrow(data), ncol = 1)
  
  
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(decision.rules)) {
      # define x:
      x <- data[i,]
      if (eval(parse(text = decision.rules[j,1])) == T) {
        decisions[[j]] <- as.character(decision.rules[j,2])
      } else if (eval(parse(text = decision.rules[j,1])) == F &&
                 is.na(decision.rules[j,3]) == F) {
        decisions[[j]] <- as.character(decision.rules[j,3])
      } else{
        decisions[[j]] <- NA
      }
    }
    
    # make decisions for ith sample in test_input
    # first remove NAs:
    decisions <- decisions[!is.na(decisions)]
    if (all(is.na(decisions)) == F) {
      decisions <- as.matrix(decisions)
      uniqx <- unique(decisions[,1])
      # consider the case with tied decisions: in that case, a sample will be assigned as "UNS"
      if (length(unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])) == 0) {
        final.decisions[i,1] <- "UNS"
      } else if (length(unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])) == 1) {
        final.decisions[i,1] <- unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])
      } else if (length(unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])) != 1 && length(unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])) != 0) {
        final.decisions[i,1] <- "UNS"
      }
    } else if(all(is.na(decisions)) == T) {
      final.decisions[i,1] <- NA
    }
  }
  
  prediction.results <- data.frame("sample no." = 1:nrow(data),
                                   "predicted class labels" = final.decisions[,1])
  
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


