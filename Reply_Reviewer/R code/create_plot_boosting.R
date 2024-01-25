library(ggplot2)
library(tidyverse)
library(dplyr)

# data 1:
results_data1 = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/data1_results_boosting.csv")
results_data1[,1]=c(rep("Liver",nrow(results_data1)))
colnames(results_data1) = c("dataset", "step1", "step2")

results_data3 = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/data3_results_boosting.csv")
results_data3[,1]=c(rep("CNS",nrow(results_data3)))
colnames(results_data3) = c("dataset", "step1", "step2")

results_data4 = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/data4_results_boosting.csv")
results_data4[,1]=c(rep("Glioblastoma",nrow(results_data4)))
colnames(results_data4) = c("dataset", "step1", "step2")

results_data6 = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/data6_results_boosting.csv")
results_data6[,1]=c(rep("Prostate",nrow(results_data6)))
colnames(results_data6) = c("dataset", "step1", "step2")

results_data7 = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/data7_results_boosting.csv")
results_data7[,1]=c(rep("NHL",nrow(results_data7)))
colnames(results_data7) = c("dataset", "step1", "step2")

results_data11 = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/data11_results_boosting.csv")
results_data11[,1]=c(rep("Breast",nrow(results_data11)))
colnames(results_data11) = c("dataset", "step1", "step2")

results_data_khan = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/khan_results_boosting.csv")
results_data_khan[,1]=c(rep("SRBCTs",nrow(results_data_khan)))
colnames(results_data_khan) = c("Khan", "step1", "step2")

results_data_armstrong = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/armstrong_results_boosting.csv")
results_data_armstrong[,1]=c(rep("Leukemia",nrow(results_data_armstrong)))
colnames(results_data_armstrong) = c("Khan", "step1", "step2")

results_data_bhattacharjee = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/bhattacharjee_results_boosting.csv")
results_data_bhattacharjee[,1]=c(rep("Lung",nrow(results_data_bhattacharjee)))
colnames(results_data_bhattacharjee) = c("Khan", "step1", "step2")

results_data_dyrskjot = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/dyrskjot_results_boosting.csv")
results_data_dyrskjot[,1]=c(rep("Bladder",nrow(results_data_dyrskjot)))
colnames(results_data_dyrskjot) = c("Khan", "step1", "step2")

results_data_yeoh = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/yeoh_results_boosting.csv")
results_data_yeoh[,1]=c(rep("ALL",nrow(results_data_yeoh)))
colnames(results_data_yeoh) = c("Khan", "step1", "step2")

results_data_tcga = read.csv("F:/Dropbox/UM Biostatistics/Research/1st Project/TNBCtype_research/Programs/Benchmarking/exp22_reviewer/experiment3/Results/tcga_results_boosting.csv")
results_data_tcga[,1]=c(rep("TNBC",nrow(results_data_tcga)))
colnames(results_data_tcga) = c("Khan", "step1", "step2")



results_final = rbind.data.frame(results_data1,results_data3, results_data4,
                                 results_data6,results_data7,results_data11,
                                 results_data_khan,results_data_armstrong,results_data_bhattacharjee,
                                 results_data_dyrskjot,results_data_yeoh,results_data_tcga,
                                 
)

results_final[,1] = as.character(results_final[,1])

# grouped boxplot
# show labels on x-axis not in non-alphabetical order
# results_final$depth <- factor(results_final$depth, levels = unique(results_final$depth), ordered = T)
results_final$dataset <- factor(results_final$dataset, levels = c("Liver", "CNS","Glioblastoma","Prostate",
                                                                  "NHL",
                                                                  "Breast",
                                                                  "SRBCTs",
                                                                  "Leukemia",
                                                                  "Lung",
                                                                  "Bladder",
                                                                  "ALL",
                                                                  "TNBC"), ordered = T)
melted_df <- reshape2::melt(results_final, id.vars = "dataset")


#
boosting_plot <- ggplot(melted_df, aes(x = dataset, y = value, fill = variable)) +
  geom_boxplot() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5,7.5,8.5,9.5), linetype = "longdash") +
  xlab("") + ylab("Number of Variables") +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "right") +
  theme_linedraw() +
  theme_light() 
boosting_plot +  scale_fill_manual(values = c("dodgerblue4", "deepskyblue1")) + labs(fill = NULL)





