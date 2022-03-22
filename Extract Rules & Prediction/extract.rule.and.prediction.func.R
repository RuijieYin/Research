source("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/treeRule.Obtain.Judge.R")

train_data <- read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/train_data_example.csv", row.names = 1)
test_data <- read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/test_data_example.csv", row.names = 1)

# extract.rule.and.prediction.func(data_input = train_data, n.trees = 1000, mtry = 35, initial_seed = 1, node_depth_max = 3, test_input = test_data)

data_input = train_data
n.trees = 1000
mtry = 35
initial_seed = 1
node_depth_max = 3
test_input = test_data

#####################
# parameters: 
# n.trees: number of trees to be grown
# initial_seed: seed number to be used
# node_depth_max: sets the value of nodedepth in rfsrc()

# ignore warnings when coercing:
options(warn=-1)

# initialize a data frame to store c("rules", "class label", "performance score")
final_result <- data.frame()

# initialize a list to store trees
tree.collection <- vector("list", n.trees)

# class type conversion
names(data_input) <- make.names(names(data_input))
data_input$subtype <- factor(data_input$subtype)

# grow n.tree of trees
for (i in 1:n.trees) {
  set.seed(initial_seed)
  tree.collection[[i]] <- rfsrc(subtype~., data = data_input,
                                ntree = 2, mtry = mtry, nodedepth = node_depth_max,
                                bootstrap = "by.root", samptype = 'swr',
                                membership = T,
                                # grow class balanced trees
                                case.wt = randomForestSRC:::make.wt(data_input$subtype),
                                sampsize = randomForestSRC:::make.size(data_input$subtype))
  initial_seed <- initial_seed + 1
} 

# calculate the total number of rules obtained from all trees
total_num_rules <- 0
for (i in 1:n.trees) {
  total_num_rules <- total_num_rules + length(getTreeRule(tree.collection[[i]])$treeRule)
}

# extract tree rules from each tree:

for (i in 1:n.trees) {
  number.of.rules <- length(getTreeRule(tree.collection[[i]])$treeRule)
  rules.tree.i <- vector("list", number.of.rules)
  for (j in 1:number.of.rules) {
    # extract rule j from tree i:
    rules.tree.i[[j]] <- parseRule(getTreeRule(tree.collection[[i]])$treeRule[[j]])
    # index of inbag sample to grow tree i:
    index.inbag <- which(tree.collection[[i]]$inbag != 0)
    # find inbag samples suffice this rule:
    x <- data_input # should tell what x is
    all.sample.fit.rule <- which(eval(parse(text = rules.tree.i[[j]])))
    inbag.sample.fit.rule <- Reduce(intersect, list(index.inbag,all.sample.fit.rule))
    
    # determine the class label for this rule
    # should consider the tied class: 
    class_label_train <- table(data_input[inbag.sample.fit.rule, dim(data_input)[2]])
    class <- names(class_label_train[class_label_train == max(class_label_train)])
    
    # calculate the performance score for rule j in tree i using oob data from tree i:
    # get class balanced OOB sample:
    out.of.bag <- which(tree.collection[[i]]$inbag == 0)
    oob.sample <- data_input[out.of.bag,]
    oob.sample <- oob.sample[complete.cases(oob.sample), ]
    oob.sample.balanced <- downSample(x = data_input[,-dim(data_input)[2]],
                                      y = factor(data_input[,dim(data_input)[2]]))
    colnames(oob.sample.balanced)[dim(oob.sample.balanced)[2]] <- "subtype"
    
    # store results in the iteration:
    class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                          ncol = 2)
    # 1st element: see if sample suffice the above criteria, 1/0
    # 2nd element: see if sample is indeed the predicted class, 1/0
    
    for (k in 1:dim(oob.sample.balanced)[1]) {
      if (is.na(eval(parse(text = rules.tree.i[[j]]))[out.of.bag[k]]) == F &&
          eval(parse(text = rules.tree.i[[j]]))[out.of.bag[k]] == T && 
          oob.sample.balanced[i,]$subtype == class) {
        class_label[k,] <- c(1,1)
      } else if (is.na(eval(parse(text = rules.tree.i[[j]]))[out.of.bag[k]]) == F &&
                 eval(parse(text = rules.tree.i[[j]]))[out.of.bag[k]] == T && 
                 oob.sample.balanced[i,]$subtype != class) {
        class_label[k,] <- c(1,0)    
      } else {
        class_label[k,] <- c(0,0)     
      }
    }
    
    # see notes for details on how to calculate performance score
    perform.score <- (sum(class_label[,1])/dim(class_label)[1])*
      (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))
    
    # retain rules with performance scores and is not NaN
    if (is.nan(perform.score) == F && is.numeric(perform.score) == T) {
      final_result <- rbind.data.frame(final_result, c(rules.tree.i[[j]], class, as.numeric(perform.score)))
    } else {
      final_result <- final_result
    }
    
  }
  final_result[,3] <- as.numeric(final_result[,3])
  # remove NAs if there's any
  final_result <- final_result[complete.cases(final_result),]
  
  
  colnames(final_result) <- c("Rules", "Class Label", "Performance Score")
  rules_table <<- final_result
} 
# print(rules_table)

# predict on the test data:
# count the number of classes:
num.of.classes <- length(levels(data_input$subtype))

# retain the top 20 rules for each class:
# rules are stored in data frame final.decision.rules

final.decision.rules <- matrix(NA, nrow = length(levels(data_input$subtype))*10, ncol = 3)
# create colnames to use rbind
colnames(final.decision.rules) <- c("Rules", "Class Label", "Performance Score")

for (i in 1:length(levels(data_input$subtype))) {
  class.i.rules <- rules_table[rules_table$"Class Label" == levels(data_input$subtype)[i],]
  class.i.rules.20 <- class.i.rules[order(class.i.rules$"Performance Score", decreasing = T),][1:20,]
  final.decision.rules <- rbind(final.decision.rules, class.i.rules.20)
}

# remove NAs and the 3rd column "Performance Score"
final.decision.rules <- final.decision.rules[complete.cases(final.decision.rules), -3]

# use majority vote to predict the class labels of new data
decisions <- vector("list", nrow(final.decision.rules)) 
final.decisions <- matrix(NA, nrow = nrow(test_input), ncol = 1)
for (i in 1:nrow(test_input)) {
  for (j in 1:nrow(final.decision.rules)) {
    x <- test_input[i,]
    if (is.na(eval(parse(text = final.decision.rules[j,1]))) == T) {
      decisions[[j]] <- NA
    } else if (eval(parse(text = final.decision.rules[j,1])) == T) {
      decisions[[j]] <- as.character(final.decision.rules[j,2])
    } else if (eval(parse(text = final.decision.rules[j,1])) == F) {
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
  }
  if (length(unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])) == 1) {
    final.decisions[i,1] <- unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])
  } else if (length(unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])) != 1 && length(unlist(uniqx[which.max(tabulate(match(decisions[,1], uniqx)))])) != 0) {
    final.decisions[i,1] <- "UNS"
  }
} 

prediction.results <- data.frame("sample no." = 1:nrow(test_input),
                                 "predicted class labels" = final.decisions[,1])
#return(prediction.results)


