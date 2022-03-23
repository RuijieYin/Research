eval.kappa <- function (test_input){
  # see performance: Cohen's Kappa
  # 1st row: reference labels
  # 2nd row: predicted labels
  prediction <- matrix(prediction.results[,2], nrow = 1)
  
  compare_labels <- rbind(matrix(test_input[,dim(test_input)[2]], nrow = 1), 
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
  # return(con_table$table)
  return(con_table$overall)
}