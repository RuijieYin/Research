# function 1
# split.train
# Split the training data (TCGA BRCA data, 1095 samples)  
# 1000 samples are used as training data (of which 800 samples are used extracting tree rules and the)
# rest 200 are used for determining the optimal number of rules)
# the remaining 95 samples are used as test data

library(data.table)
library(randomForestSRC)
library(caret)


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

extract.rules <- function (data = data.split$train.data, n.trees = 500, mtry = 50, node_depth_max = 3) {
# grow n.tree of trees
set.seed(1)
tree.collection <- rfsrc(subtype~., data = data,
                         ntree = n.trees, mtry = mtry, nodedepth = node_depth_max,
                         bootstrap = "by.root", samptype = 'swr',
                         membership = T,
                         # grow class balanced trees
                         case.wt = randomForestSRC:::make.wt(data$subtype),
                         sampsize = randomForestSRC:::make.size(data$subtype))

# calculate the total number of rules obtained from all trees
total_num_rules <- length(getTreeRule(tree.collection)$tree.id)


# initialize a data frame to store c("rules", "class label", "performance score")
final_result <- data.frame()

# extract tree rules from each tree:
rules <- vector("list", total_num_rules)

for (j in 1:total_num_rules) {
  # calculate how many rules are generated from the ith tree:
  #number.of.rules <- length(which(getTreeRule(tree.collection)$tree.id == i))
  #rules.tree.i <- vector("list", number.of.rules)
  
  # extract rule j from tree i:
  rules[[j]] <- parseRule(getTreeRule(tree.collection)$treeRule[[j]])
  # the ith tree where rule j comes from:
  tree.i <- getTreeRule(tree.collection)$tree.id[j]
  
  # index of inbag sample to grow tree i:
  index.inbag <- which(tree.collection$inbag[,tree.i] != 0)
  # find inbag samples suffice the jth rule:
  x <- data # should tell what x is
  all.sample.fit.rule <- which(eval(parse(text = rules[[j]])))
  inbag.sample.fit.rule <- Reduce(intersect, list(index.inbag,all.sample.fit.rule))
  
  # determine the class label for this rule
  # should consider the tied class: 
  class_label_train <- table(data[inbag.sample.fit.rule, dim(data)[2]])
  class <- names(class_label_train[class_label_train == max(class_label_train)])
  
  # see if this rule has an "else": the majority class that does not suffice jth rule in the inbag sample,
  # and has a proportion over 50%:
  inbag.sample.not.fit.rule <- index.inbag[-inbag.sample.fit.rule]
  # find the class labels of samples that does not suffice jth rule:
  else_data <- as.character(data[inbag.sample.not.fit.rule, dim(data)[2]])
  else_label_train <- table(data[inbag.sample.not.fit.rule, dim(data)[2]])
  else_class_temp <- names(class_label_train[class_label_train == max(class_label_train)])
  # calculate proportion:
  prop <- sum(else_data == else_class_temp, na.rm = TRUE)/length(else_data == else_class_temp)
  if (is.nan(prop) == T) {
    else_class <- NA
  } else if (prop > 0.5) {
    else_class <- else_class_temp
  } else if (prop <= 0.5) {
    else_class <- NA
  }
  
  # calculate the performance score for rule j in tree i using oob data from tree i:
  # get class balanced OOB sample:
  out.of.bag <- which(tree.collection$inbag[,tree.i] == 0)
  oob.sample <- data[out.of.bag,]
  oob.sample <- oob.sample[complete.cases(oob.sample), ]
  oob.sample.balanced <- downSample(x = data[,-dim(data)[2]],
                                    y = factor(data[,dim(data)[2]]))
  colnames(oob.sample.balanced)[dim(oob.sample.balanced)[2]] <- "subtype"
  
  # 2 scenarios: the rules has or does not have an 'else'
  
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
    perform.score <- (sum(class_label[,1])/dim(class_label)[1])*
      (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))
  } else if (is.na(else_class) == F) {
    # part 1:
    # store results in the iteration:
    class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                          ncol = 2)
    # 1st element: see if sample suffice the above criteria, 1/0
    # 2nd element: see if sample is indeed the predicted class, 1/0
    
    for (k in 1:dim(oob.sample.balanced)[1]) {
      if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
          eval(parse(text = rules[[j]]))[out.of.bag[k]] == T && 
          oob.sample.balanced[i,]$subtype == class) {
        class_label[k,] <- c(1,1)
      } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                 eval(parse(text = rules[[j]]))[out.of.bag[k]] == T && 
                 oob.sample.balanced[i,]$subtype != class) {
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
          oob.sample.balanced[i,]$subtype == else_class) {
        class_label[k,] <- c(1,1)
      } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                 eval(parse(text = rules[[j]]))[out.of.bag[k]] == F && 
                 oob.sample.balanced[i,]$subtype != else_class) {
        class_label[k,] <- c(1,0)    
      } else {
        class_label[k,] <- c(0,0)     
      }
    }
    
    # see notes for details on how to calculate performance score
    perform.score.2 <- (sum(class_label[,1])/dim(class_label)[1])*
      (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))
    
    perform.score = perform.score.1 + perform.score.2
  }
  
  # retain rules with performance scores and is not NaN
  if (is.nan(perform.score) == F && is.numeric(perform.score) == T) {
    add.row <- data.frame(rules[[j]], class, else_class, as.numeric(perform.score))
    final_result <- rbind.data.frame(final_result, add.row)
  } else {
    final_result <- final_result
  }
  
  
  final_result[,4] <- as.numeric(final_result[,4])
  # remove NAs if there's any
  # final_result <- final_result[complete.cases(final_result),]
  
  rules_table <- final_result
  colnames(rules_table) <- c("Rules", "Class Label", "Else Class", "Performance Score")
}

  rules.all <- list("Rules" = rules_table) 
  return(rules.all)
# function ends here
}


# function 3
# use the remaining 20% training data (validation data) to determine the optimal number of rules for each class:
# not weighted voting:

library(ggplot2)
library(tidyr)

select.rules <- function(validation.data = data.split$validation.data, rule.table) {
# calculate the total number of class:
num_class <- length(levels(validation.data$subtype))

# store number of rules v.s. accuracy for each class
# the list is also used for plotting
plot.table.t <- list()

for (i in 1:length(levels(validation.data$subtype))) {
  
  class.data <- validation.data[validation.data$subtype == levels(validation.data$subtype)[i],]
  
  rules.class <- rule.table[rule.table$`Class Label` == levels(validation.data$subtype)[i],]
  # remove rules with performance score = 0
  rules.class <- rules.class[rules.class$`Performance Score` != 0,]
  # order the rules by their performance scores:
  rules.class <- rules.class[order(rules.class$`Performance Score`, decreasing = T),]
  
  
  # store results in the iteration:
  accuracy.table <- matrix(NA, nrow = dim(rules.class)[1],
                           ncol = 2)
  # 1st element: the number of top t rules used
  # 2nd element: the accuracy for the number of top t rules used
  
  # if (is.na(rules.class$`Else Class`) == T ) {
  for (k in 1:dim(class.data)[1]) {
    # define x
    x <- class.data[k,]
    
    # t: top t rules for each class
    for (t in 1:dim(rules.class)[1]) {
      
      rules.top.t <- rules.class[1:t,]
      inter.table <- matrix(NA, nrow = t, ncol = 1)
      for (w in 1:t) {
        if (eval(parse(text = rules.top.t[w,]$Rules)) == T &&
            x$subtype == levels(class.data$subtype)[i]) {
          inter.table[w,1] <- 1
        } else {
          inter.table[w,1] <- 0
        }
      }
      accuracy.t <- sum(inter.table[,1])/(dim(inter.table)[1])
      accuracy.table[t,] <- c(t, accuracy.t)
    }
    
  }
  #}     
  plot.table.t[[i]] <- accuracy.table
}     

# loop ends here

# plot top # of rules v.s. accuracy:
# if((num_class %% 2) == 0) {
#   par(mfrow = c(num_class/2, 2))
# } else {
#   par(mfrow = c((num_class+1)/2, 2))
# }

plot.list <- list()
for ( i in 1:num_class) {
  df <- data.frame(x = plot.table.t[[i]][,1], y = plot.table.t[[i]][,2])
  plot.list[[i]] <- ggplot(data = df, aes(x = x, y = y, group=1)) + 
    geom_line() + 
    geom_point() +
    xlab("Number of Rules") +
    ylab("Accuracy") +
    ggtitle(levels(validation.data$subtype)[i]) +
    theme_classic()
}

# determine the optimal # of rules:
optimal_num <- matrix(NA, nrow = num_class, ncol = 1)
final.decision.rules <- data.frame()
decision.rules <- list()
for(k in 1:num_class) {
  rules.class <- rule.table[rule.table$`Class Label` == levels(validation.data$subtype)[k],]
  # remove rules with performance score = 0
  rules.class <- rules.class[rules.class$`Performance Score` != 0,]
  # order the rules by their performance scores:
  rules.class <- rules.class[order(rules.class$`Performance Score`, decreasing = T),]
  
  rule.data <- plot.table.t[[k]]
  optimal_num[k,1] <- which(rule.data[,2] == max(rule.data[,2]))[length(which(rule.data[,2] == max(rule.data[,2])))]
  
  decision.rules[[k]] <- rules.class[1:optimal_num[k,1],]
  
  final.decision.rules <- rbind(final.decision.rules, rules.class[1:optimal_num[k,1],])
}
# function ends here
  optimal.rules <- list("Plot.Num.Accuracy" = plot.list,
                        "Optimal.Num" = optimal_num,
                        "Selected.Rules" = final.decision.rules)
  return(optimal.rules)

}



# function 4: Prediction 
# test data
# use majority voting (non-wighted) to predict the class labels of new data
#

predict.rules <- function(data = data.split$test.data, decision.rules) {
# store predicted labels from each rule in a list
decisions <- vector("list", length = nrow(decision.rules)) 
# make a final decision for each sample in the test data
final.decisions <- matrix(NA, nrow = nrow(data), ncol = 1)


for (i in 1:nrow(data)) {
  for (j in 1:nrow(decision.rules)) {
    # define x:
    x <- data[i,]
    if (is.na(eval(parse(text = decision.rules[j,1]))) == T) {
      decisions[[j]] <- NA
    } else if (eval(parse(text = decision.rules[j,1])) == T) {
      decisions[[j]] <- as.character(decision.rules[j,2])
    } else if (eval(parse(text = decision.rules[j,1])) == F) {
      decisions[[j]] <- NA
    }
  }
  
  # make decisions for ith sample in test_input
  # first remove NAs:
  decisions <- decisions[!is.na(decisions)]
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


# function 5: weighted voting prediction
# Prediction 
# test data
# use majority voting (wighted) to predict the class labels of new data
#
# store predicted labels from each rule in a list
predict.rules.weighted <- function(data = data.split$test.data, decision.rules) {
decisions.df <- data.frame( Predicted = character(), 
                            Scores = numeric())
# make a final decision for each sample in the test data
final.decisions <- matrix(NA, nrow = nrow(data), ncol = 1)


for (i in 1:nrow(data)) {
  for (j in 1:nrow(decision.rules)) {
    # define x:
    x <- data[i,]
    # weighted voting: use performance scores of each rule as weights
    if (is.na(eval(parse(text = decision.rules[j,1]))) == T) {
      decisions.df[j,] <- c(NA, NA)
    } else if (eval(parse(text = decision.rules[j,1])) == T) {
      decisions.df[j,] <- c(decision.rules[j,2], decision.rules[j,4])
    } else if (eval(parse(text = decision.rules[j,1])) == F) {
      decisions.df[j,] <- c(NA, NA)
    }
  }
  
  # make decisions for ith sample in test_input
  # first remove NAs:
  decisions.df <- decisions.df[complete.cases(decisions.df),]
  # weighted voting: use performance scores as weights
  weighted.df <- data.frame( Predicted = character(), 
                             Scores = numeric())
  # see number of unique predicted class labels
  num_pred <- length(unique(decisions.df$Predicted))
  for (k in 1:num_pred) {
    weighted.df[k,] <- c(unique(decisions.df$Predicted)[k],
                         sum(as.numeric(decisions.df[decisions.df$Predicted == unique(decisions.df$Predicted)[k], 2])))
  }
  
  if(is.na(weighted.df[which.max(weighted.df$Scores),1]) == T) {
    final.decisions[i,1] <- "UNS"
  } # consider tied votes
  else if(length(weighted.df[which.max(weighted.df$Scores),1]) > 1) {
    final.decisions[i,1] <- "UNS"
  } else {
    final.decisions[i,1] <- weighted.df[which.max(weighted.df$Scores),1]
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






