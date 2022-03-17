# a function to extract rules
# note 1: performance score of a rule is also reported
# note 2: sanity check of input is not included

# preferred structure of training data: 
# columns are variables and the last column should be class labels for each sample
# rows are samples; should name the last column as 'subtype'
library(randomForestSRC)
library(caret)

source("/Users/ruijieyin/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/treeRule.Obtain.Judge.R")

extract_rules <- function (data, n.trees, mtry, initial_seed, node_depth_max) {
  # parameters: 
  # n.trees: number of trees to be grown
  # initial_seed: seed number to be used
  # node_depth_max: sets the value of nodedepth in rfsrc()
  
  
  # initialize a data frame to store c("rules", "class label", "performance score")
  final_result <- data.frame()
  
  # initialize a list to store trees
  tree.collection <- vector("list", n.trees)
  
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
  # initialize a null list to store all rules obtained:
  rules.collection <- vector("list", total_num_rules)
  
  for (i in 1:n.trees) {
    number.of.rules <- length(getTreeRule(tree.collection[[i]])$treeRule)
    rules.tree.i <- vector("list", number.of.rules)
    for (j in 1:number.of.rules) {
      # extract rule j from tree i:
      rules.tree.i[[j]] <- parseRule(getTreeRule(tree.collection[[i]])$treeRule[[j]])
      # index of inbag sample to grow tree i:
      index.inbag <- which(tree.collection[[i]]$inbag != 0)
      # find inbag samples suffice this rule:
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
      final_result <- rbind.data.frame(final_result, c(rules.tree.i[[j]], class, perform.score))
      
    }
    colnames(final_result) <- c("Rules", "Class Label", "Performace Score")
    results_table <<- final_result
  } 
  print(results_table)
}




