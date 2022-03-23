
pred.func <- function(test_input) {
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

prediction.results <<- data.frame("sample no." = 1:nrow(test_input),
                                 "predicted class labels" = final.decisions[,1])
return(prediction.results)
}
