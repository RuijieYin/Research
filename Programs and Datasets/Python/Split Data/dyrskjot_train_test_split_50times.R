# function 1: divide.data
# divide the datasets into 70% training and 30% test
# class labels in both training and test data are balanced;
divide.data <- function(input_data, seed_number) {
  data_pre_division <- input_data
  # transpose data_pre_division, so that rows are samples, columns are genes
  # and merge with its class labels
  data_pre_division_df <- data_pre_division
  # count how many classes are there in the data, for multiclass datasets, this number should be greater than 2
  # the last column are the class labels
  num_class <- length(unique(data_pre_division_df[, ncol(data_pre_division_df)]))
  # find out the class labels
  class_labels_collection <- unique(data_pre_division_df[, ncol(data_pre_division_df)])
  for (i in 1:num_class) {
    nam <- paste("class_labels", i, sep = "_")
    assign(nam, class_labels_collection[i])
  }
  
  # subset data_pre_division_df into multiple groups, by class labels
  for (i in 1:num_class) {
    data_nam <- paste("class_data_temp", i, sep = "_")
    nam <- paste("class_labels", i, sep = "_")
    assign(data_nam, data_pre_division_df[data_pre_division_df[,ncol(data_pre_division_df)] == eval(parse(text = nam)),])
  }
  
  # add a row number, to tell us which samples are selected as training samples:
  for (i in 1:num_class) {
    data_nam <- paste("class_data_temp", i, sep = "_")
    ready_to_split_data_nam <- paste("class_data", i, sep = "_")
    assign(ready_to_split_data_nam, cbind.data.frame(eval(parse(text = data_nam)), matrix(c(1:nrow(eval(parse(text = data_nam)))), ncol = 1)))
    # rename as row_number:
    class_data = eval(parse(text = ready_to_split_data_nam))
    colnames(class_data)[dim(class_data)[2]] <- "row_number"
    # randomly sample 70% data as training data
    set.seed(seed_number)
    train_nam <- paste("class_train", i, sep = "_")
    test_nam <- paste("class_test", i, sep = "_")
    
    assign(train_nam, class_data[sample(nrow(class_data), round(0.7*nrow(class_data))), ])
    # train_data is used only for train_data $row_number, to identify the row numbers that are indicators of training data
    # train_data is the final training data for class i 
    train_data <- eval(parse(text = train_nam))
    
    assign(test_nam, class_data[!class_data$row_number %in% train_data$row_number,])
    
    # test_data is the final test data for class i 
    test_data <- eval(parse(text = test_nam))
    # drop row_number
    # use the same name train_nam and test_data to overwrite the data
    assign(train_nam, train_data[,-dim(train_data)[2], ])
    assign(test_nam, test_data[,-dim(test_data)[2],])
  }
  
  # merge the training and two test data
  # rows are samples and columns are genes
  train_data_final <- data.frame()
  test_data_final <- data.frame()
  for(i in 1:num_class) {
    train_nam <- paste("class_train", i, sep = "_")
    test_nam <- paste("class_test", i, sep = "_")
    train_data_final <- rbind.data.frame(train_data_final, eval(parse(text = train_nam)))
    test_data_final <- rbind.data.frame(test_data_final, eval(parse(text = test_nam)))
  }
  
  data.list <- list("train.data" = train_data_final,
                    "test.data" = test_data_final)
  
  return(data.list)
}


# divide data into 20 training and 20 test data
setwd("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp17 50 replica/Data")
input_data = read.csv("dyrskjot_all_with_label.csv")
for(i in 1:50) {
  train = divide.data(input_data, i)$train.data
  test = divide.data(input_data, i)$test.data
  file_name <- paste0("train_dyrskjot_", i, ".csv")
  write.csv(train, file_name, row.names = F)
  file_name2 <- paste0("test_dyrskjot_", i, ".csv")
  write.csv(test, file_name2, row.names = F)
}









