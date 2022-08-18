library(data.table)
library(randomForestSRC)
library(caret)
library(ggplot2)
library(dplyr)


# function 1
# split.train
# Split the training data (TCGA BRCA data, 1095 samples)  
# 1000 samples are used as training data (of which 800 samples are used extracting tree rules and the)
# rest 200 are used for determining the optimal number of rules)
# the remaining 95 samples are used as test data

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


extract.rules <- function (data = data.split$train.data, n.trees = 50, mtry = 5, node_depth_max = 4) {
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
                             "Class Label" = character())
  
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
    
    final_result[j,] <- list("Rules" = rules[[j]],
                             "Class Label" = class)
    
  }
  
  # remove duplicated rules by rows (attempt 1)
  final_result_all <- final_result[!duplicated(final_result),]
  
  # remove duplicated rules by column (attempt 2)
  final_result_all <- final_result_all %>% distinct(final_result_all[,1], .keep_all = TRUE)
  
  # remove the last column: duplicated results
  final_result_all <- final_result_all[,-3]
  
  colnames(final_result_all) <- c("Rules", "Class Label")
  rules.all <- list("Rules" = final_result_all) 
  # function ends here
}


# function 3
# use the remaining 20% training data (validation data) to determine the optimal number of rules:
# whole validation data is used:
# not weighted voting:
# consecutive selection of rules:
library(ggplot2)
library(dplyr)

select.rules <- function(validation.data = data.split$validation.data, rule.table=rules$Rules) {
  
  # in validation data: rows are samples and columns are vars
  
  # continue to use the old name: rules.order; but rules are NOT ordered by their performance scores.
  rules.order <- rule.table
  
  
  # part 1:
  # initialize a data frame to store # of rules v.s. accuracy
  accuracy.table <- data.frame("Num.of.Rules" = as.numeric(),
                               "Accuracy" = as.numeric())
  
  # initialize accuracy
  accuracy <- 0
  
  # kick start: select the first rule to be included:
  classes <- unique(rules.order[,2])
  
  # divide the rules.order by class labels:
  for(i in 1:length(classes)) {
    assign(paste0("rules.order.class",i), rules.order[rules.order[,2] == classes[i],])
  }
  
  # use the 1st class as a starter; it doesn't matter which class we use.
  rules.class.kick.start <- rules.order[rules.order[,2] == classes[1],]
  
  # inter.table: whether a sample in validation data suffice the rule & has the correct class label predicted
  inter.table <- matrix(NA, nrow = nrow(rules.class.kick.start), ncol = nrow(validation.data))
  for (w in 1:nrow(rules.class.kick.start))  {
    for (k in 1:nrow(validation.data)) {
      # define x
      x <- validation.data[k,]
      if (eval(parse(text = rules.class.kick.start[w,]$Rules)) == T &
          x$subtype == classes[1]) {
        inter.table[w,k] <- 1
      } else {
        inter.table[w,k] <- 0
      }
    }
  } 
  
  # calculate the performance score for each rule (each row in inter.table):
  performance.table <- matrix(NA, nrow = nrow(inter.table), ncol = 1)
  for (i in 1:nrow(performance.table)) {
    performance.table[i,1] <- sum(inter.table[i,])/ncol(inter.table)
  }
  
  # initialize rules will be retained:
  rules.retained.new <- data.frame("Rules" = character(), 
                                   "Class Label" = character())
  
  # initialize rules retained in the previous selection step:
  rules.retained.previous <- data.frame("Rules" = character(), 
                                        "Class Label" = character())
  
  # initialize rules retained in the very last step:
  rules.retained.final <- data.frame("Rules" = character(), 
                                     "Class Label" = character())
  
  # rules.retained.new is the first rule that is retained in the final rules table:
  rules.retained.previous <- rules.class.kick.start[which.max(performance.table),]
  
  # select the first set of rules:
  for (j in 2:length(classes)) {
    for (i in 1:nrow(get(paste0("rules.order.class",j)))) {
      rules.add <- get(paste0("rules.order.class",j))[i,]
      rules.retained.new <- rbind(rules.retained.previous, rules.add)
      
      pred_class <- matrix(NA, nrow = dim(validation.data)[1], ncol = 1)
      
      for (k in 1:dim(validation.data)[1]) {
        # define x
        x <- validation.data[k,]
        inter.table <- data.frame("label" = as.character())
        
        for (w in 1:nrow(rules.retained.new))  {
          if (eval(parse(text = rules.retained.new[w,]$Rules)) == T) {
            inter.table[w,] <- list("label" = rules.retained.new[w,2])
          } else if (eval(parse(text = rules.retained.new[w,]$Rules)) == F) {
            inter.table[w,] <- NA
          }
        }
        
        # determine the class label of the kth sample in validation data:
        # note: table() ignores NAs
        label <- table(as.character(inter.table[,1]))
        if (length(label) != 1) {
          pred_class[k,1] <- NA
        } else{
          pred_class[k,1] <- names(label[label == max(label)])
        }
        
        comp_table <- cbind(pred_class, as.matrix(validation.data$subtype, ncol = 1))
        
        accuracy.temp <- sum(comp_table[,1] == comp_table[,2], na.rm = T)/nrow(comp_table)
        
        if (accuracy.temp <= accuracy) {
          accuracy <- accuracy
          rules.retained.previous <- rules.retained.previous
        } else if (accuracy.temp > accuracy) {
          accuracy <- accuracy.temp
          rules.retained.previous <- rules.retained.new
        }
      } 
      if (all(all.equal(rules.retained.previous, rules.retained.new) == T)) {
        break
      } 
    }
  }
  
  
  
  
  # part 2:
  # remove the selected rules so far from rules pool, and start to select next 
  # set of rules:
  # Starting point: rules.retained.new
  # Note: can use either rules.retained.new or rules.retained.previous as starting point,
  # since they are exactly the same in the last step.
  remaining.rules <- rules.order[!rules.order[,1] %chin% rules.retained.new[,1],]
  # the maximum number of rules to be considered at this point should be the 
  # number of rules of a least populated class - 1, because we just obtained a rule from the least
  # populated class
  max.num.of.rules <- min(table(rules.order[,2]))
  
  # 2 scenarios: v = 1 and v > 1, since rules.retained.previous are different
  # for v = 1, rules.retained.previous is the one obtained from part 1:
  # for both v = 1 and v > 1, we always use rules.retained.previous as the starting variable.
  v = 1
  
  rules.retained.previous <- rules.retained.new
  # divide the remaining.rules by class labels:
  for (i in 1:length(classes)) {
    assign(paste0("rules.order.class",i), remaining.rules[remaining.rules[,2] == classes[i],])
  }
  # select the next set of rules:
  for (j in 1:length(classes)) {
    for (i in 1:nrow(get(paste0("rules.order.class",j)))) {
      rules.add <- get(paste0("rules.order.class",j))[i,]
      rules.retained.new <- rbind(rules.retained.previous, rules.add)
      
      pred_class <- matrix(NA, nrow = dim(validation.data)[1], ncol = 1)
      
      for (k in 1:dim(validation.data)[1]) {
        # define x
        x <- validation.data[k,]
        inter.table <- data.frame("label" = as.character())
        
        for (w in 1:nrow(rules.retained.new))  {
          if (eval(parse(text = rules.retained.new[w,]$Rules)) == T) {
            inter.table[w,] <- list("label" = rules.retained.new[w,2])
          } else if (eval(parse(text = rules.retained.new[w,]$Rules)) == F) {
            inter.table[w,] <- NA
          }
        }
        
        # predict the class label of the kth sample in validation data:
        # note: table() ignores NAs
        label <- table(as.character(inter.table[,1]))
        if (length(label) != 1) {
          pred_class[k,1] <- NA
        } else{
          pred_class[k,1] <- names(label[label == max(label)])
        }
        
        comp_table <- cbind(pred_class, as.matrix(validation.data$subtype, ncol = 1))
        
        accuracy.temp <- sum(comp_table[,1] == comp_table[,2], na.rm = T)/nrow(comp_table)
        
        if (accuracy.temp <= accuracy) {
          accuracy <- accuracy
          rules.retained.previous <- rules.retained.previous
          
        } else if (accuracy.temp > accuracy) {
          accuracy <- accuracy.temp
          rules.retained.previous <- rules.retained.new
        }
      } 
      if (all(all.equal(rules.retained.previous, rules.retained.new) == T)) {
        break
      } 
    }
  }
  # remove the selected rules so far from rules pool, and start to select next 
  # set of rules:
  remaining.rules <- remaining.rules[!remaining.rules[,1] %chin% rules.retained.new[,1],]
  
  # v > 1
  # note: should check if the rules obtained so far (kick starter + v = 1) are class-balanced, if not --> abort
  if (length(rules.retained.previous[,2]) %% length(classes) != 0) {
    # retain only balanced rules:
    upper <- dim(rules.retained.previous)[1]%/%length(classes) * length(classes)
    rules.retained.final <- rules.retained.previous[1:upper,]
    # this means the rules are not balanced --> should abort 
  } else {
    for (v in 2:(max.num.of.rules - 1)) {
      # this means the rules are balanced so far (kick starter + v = 1) --> continue to find the next sets of rules
      # break condition:
      # check if all class labels are balanced: use %% to check the remainder of a division
      if (length(rules.retained.previous[,2]) %% length(classes) != 0) {
        # retain only balanced rules:
        upper <- dim(rules.retained.previous)[1]%/%length(classes) * length(classes)
        rules.retained.final <- rules.retained.previous[1:upper,]
        break
      } else {
        # divide the remaining.rules by class labels:
        for (i in 1:length(classes)) {
          assign(paste0("rules.order.class",i), remaining.rules[remaining.rules[,2] == classes[i],])
        }
        
        # select the next set of rules:
        for (j in 1:length(classes)) {
          for (i in 1:nrow(get(paste0("rules.order.class",j)))) {
            rules.add <- get(paste0("rules.order.class",j))[i,]
            rules.retained.new <- rbind(rules.retained.previous, rules.add)
            
            pred_class <- matrix(NA, nrow = dim(validation.data)[1], ncol = 1)
            
            for (k in 1:dim(validation.data)[1]) {
              # define x
              x <- validation.data[k,]
              inter.table <- data.frame("label" = as.character())
              
              for (w in 1:nrow(rules.retained.new))  {
                if (eval(parse(text = rules.retained.new[w,]$Rules)) == T) {
                  inter.table[w,] <- list("label" = rules.retained.new[w,2])
                } else if (eval(parse(text = rules.retained.new[w,]$Rules)) == F) {
                  inter.table[w,] <- NA
                }
              }
              
              # predict the class label of the kth sample in validation data:
              # note: table() ignores NAs
              label <- table(as.character(inter.table[,1]))
              if (length(label) != 1) {
                pred_class[k,1] <- NA
              } else{
                pred_class[k,1] <- names(label[label == max(label)])
              }
              
              comp_table <- cbind(pred_class, as.matrix(validation.data$subtype, ncol = 1))
              
              accuracy.temp <- sum(comp_table[,1] == comp_table[,2], na.rm = T)/nrow(comp_table)
              
              if (accuracy.temp <= accuracy) {
                accuracy <- accuracy
                rules.retained.previous <- rules.retained.previous
                
              } else if (accuracy.temp > accuracy) {
                accuracy <- accuracy.temp
                rules.retained.previous <- rules.retained.new
              }
            } 
            # when the above 'else if' suffices, should go to next v:
            if (all(all.equal(rules.retained.previous, rules.retained.new) == T)) {
              break
            } 
          }
        }
        # remove the selected rules so far from rules pool, and start to select next 
        # set of rules:
        remaining.rules <- remaining.rules[!remaining.rules[,1] %chin% rules.retained.new[,1],]
      }
    }
  }
  
  # # output the final selected rules:
  final.decision.rules <- rules.retained.final
  
  # final number of rules:
  final_num <- nrow(final.decision.rules)
  
  optimal.rules <- list("Selected.Rules" = final.decision.rules,
                        "Final.Num" = final_num)
  return(optimal.rules)
  
} 



# function 4: Prediction 
# test data
# use majority voting (non-weighted) to predict the class labels of new data
#

predict.rules <- function(data = data.split$test.data, decision.rules) {
  for (i in 1:nrow(decision.rules)) {
    if (i > nrow(decision.rules)) {
      break
    } else {
      pred_class <- matrix(NA, nrow = dim(data)[1], ncol = 1)
      
      for (k in 1:dim(data)[1]) {
        # define x
        x <- data[k,]
        inter.table <- data.frame("label" = as.character())
        
        for (w in 1:nrow(decision.rules))  {
          if (eval(parse(text = decision.rules[w,]$Rules)) == T) {
            inter.table[w,] <- list("label" = decision.rules[w,][,2])
          } else if (eval(parse(text = decision.rules[w,]$Rules)) == F) {
            inter.table[w,] <- list("label" = "UNS")
          }
        }
        
        # determine the class label of the kth sample in test data:
        label <- table(as.character(inter.table[,1]))
        if (length(label) != 1) {
          pred_class[k,1] <- "UNS"
        } else{
          pred_class[k,1] <- names(label[label == max(label)])
        }
      }
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
  
  # Concordance: kappa statistic
  # will need a matrix where the diagonal elements of the matrix are 
  # the agreeing elements; the discordant observations are on the off-diagonal.
  # A confusion matrix:
  con_table <- confusionMatrix(data = as.factor(compare_labels[2,]),
                               reference = as.factor(compare_labels[1,]))
  # 8.16.2022
  # if error occurs, it means the predicted values are UNS
  
  prediction.labels.kappa <- list("Prediction.results" = prediction.results,
                                  "Results.table" = con_table)
  return(prediction.labels.kappa)
}





