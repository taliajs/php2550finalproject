# libraries
library(randomForest)
library(rsample)
library(tidyverse)
library(pROC)
library(rfUtilities)
library(kableExtra)

# load in dataset
isolates <- read.csv("~/Downloads/php2550finalproject-main 3/data/isolates2.csv")


# select interest independent variables for learning its association to the new outbreak variable
outbreak_df <- isolates %>% select(c("new_outbreak","Host.category","Host.disease","Virulence.genotypes.count",
                                     "AMR.genotypes.count","isolate_source_type","SNP.cluster.rank",
                                     "Stress.genotypes.count","Year","region","Month","Week","Serovar"))

# modifying the categorical variable into numeric variable 
# converting variables with NA into unknown
outbreak_df$isolate_source_type <- unclass(factor(outbreak_df$isolate_source_type))
outbreak_df$Host.category[is.na(outbreak_df$Host.category)] <- "unknown"
outbreak_df$Host.disease[is.na(outbreak_df$Host.disease)] <- "unknown"
outbreak_df$Host.disease <- unclass(factor(outbreak_df$Host.disease))
outbreak_df$region[is.na(outbreak_df$region)] <-"unknown"
outbreak_df$region <- as.factor(outbreak_df$region)

# month model dataset: compute the number of outbreaks by interest groups and time variables
rf_model_df_m <- outbreak_df %>% select(-c(Week)) %>% group_by(Serovar,Year,Month,region,isolate_source_type) %>%
  # count number of isolate by groups, modify year variables starting 2017
  mutate(Outbreak_num=sum(new_outbreak)) %>%
  select(-c("new_outbreak"))

# week model dataset: compute the number of outbreaks by interest groups and time variables
rf_model_df_w <- outbreak_df %>% select(-c(Month)) %>%group_by(Serovar,Year,Week,region,isolate_source_type) %>%
  # count number of isolate by groups, modify year variables starting 2017
  mutate(Outbreak_num=sum(new_outbreak)) %>%
  select(-c("new_outbreak"))

# Simple Training/Test Set Splitting
set.seed(123)

# split data used for month model
partition_m <- rf_model_df_m %>% 
  initial_split(strata = Year, prop = 0.8)

train_df_m<- training(partition_m)
test_df_m <- testing(partition_m)

# split data used for week model
partition_w <- rf_model_df_w%>% 
  initial_split(strata = Year, prop = 0.8)

train_df_w<- training(partition_w)
test_df_w <- testing(partition_w)

# fit the random forest month model
bag.outbreak_m <- randomForest(Outbreak_num~.,data=train_df_m,mtry=6,importance=TRUE)
bag.outbreak_m

# fit the random forest week model
bag.outbreak_w <- randomForest(Outbreak_num~.,data=train_df_w,mtry=6,importance=TRUE)
bag.outbreak_w

# predicting on the test dataset
predict_m <- predict(bag.outbreak_m,test_df_m)
predict_w <- predict(bag.outbreak_w,test_df_w)

# visualize the observes vs. predicted for knowing the predictive power on test data
par(mfrow=c(1,2))
plot(predict_m, test_df_m$Outbreak_num,xlab="Predicted",ylab="Observed")
title(main="Random Forest Model \nPredicting in Months")
abline(0,1)

plot(predict_w, test_df_w$Outbreak_num,xlab="Predicted",ylab="Observed")
title(main="Random Forest Model \nPredicting in Weeks")
abline(0,1)

# checking the importance of each variable
varImpPlot(bag.outbreak_m,main="Variable Importance of Month Model",cex=0.8)
varImpPlot(bag.outbreak_w,main="Variable Importance of Week Model",cex=0.8)

train_ROC_m <- roc(train_df_m$Outbreak_num,predict(bag.outbreak_m,train_df_m))
test_ROC_m <- roc(test_df_m$Outbreak_num,predict_m)
plot.roc(train_ROC_m,col="red",main="ROC for RF Model \non Months")
plot.roc(test_ROC_m,add=TRUE,col="green")
legend("bottomright",c("Train", "Test"), lty=1, 
       col = c("red", "green"),bty="n",cex=0.8)

train_ROC_w <- roc(train_df_w$Outbreak_num,predict(bag.outbreak_w,train_df_w))
test_ROC_w <- roc(test_df_w$Outbreak_num,predict_w)
plot.roc(train_ROC_w,col="red",main="ROC for RF Model \non Weeks")
plot.roc(test_ROC_w,add=TRUE,col="green")
legend("bottomright",c("Train", "Test"), lty=1, 
       col = c("red", "green"),bty="n",cex=0.8)


# MSE
train_MSE_m <- round(mean((predict(bag.outbreak_m,train_df_m)-train_df_m$Outbreak_num)^2),2)
train_MSE_w <- round(mean((predict(bag.outbreak_w,train_df_w)-train_df_w$Outbreak_num)^2),2)
test_MSE_m <- round(mean((predict_m-test_df_m$Outbreak_num)^2),2)
test_MSE_w <- round(mean((predict_w-test_df_w$Outbreak_num)^2),2)

# Measuring metrics

rf_metrics_df <- data.frame(AUC=c(round(auc(train_ROC_m),2),round(auc(test_ROC_m),2)),
                            MSE=c(train_MSE_m,test_MSE_m),
                            auc=c(round(auc(train_ROC_w),2),round(auc(test_ROC_w),2)),
                            mse=c(train_MSE_w,test_MSE_w))
rf_metrics_df<- rbind(rf_metrics_df, c("Month","Month","Week","Week"))
rf_metrics_df<-t(rf_metrics_df)
colnames(rf_metrics_df) <-c("Training","Testing","Model")
rf_metrics_df<-data.frame(rf_metrics_df)

# comparison table
rf_metrics_df[,-3] %>% kbl(caption = 'Comparison Table of 2 RF Models') %>% 
  pack_rows(index = table(rf_metrics_df$Model)) %>%
  kable_classic(full_width = F, html_font = "Cambria")

