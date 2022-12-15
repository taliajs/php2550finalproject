# libraries
library(tidyverse)
library(ggplot2)
library(lmtest)
library(rpart.plot)
library(kableExtra)
library(caret)
library(pROC)
library(DescTools)

# load in dataset
isolates <- read.csv("~/Downloads/php2550finalproject-main 3/data/isolates2.csv")


# select interest independent variables for learning its association to the new outbreak variable
outbreak_df <- isolates %>% select(c("new_outbreak","Host.category","Host.disease",
                                   "Virulence.genotypes.count","AMR.genotypes.count",
                                   "isolate_source_type","SNP.cluster.rank","Stress.genotypes.count",
                                   "Year","region","Month","Week","Serovar"))

# modifying the categorical variable into numeric variable 
# converting variables with NA into unknown
outbreak_df$isolate_source_type <- unclass(factor(outbreak_df$isolate_source_type))
outbreak_df$Host.category[is.na(outbreak_df$Host.category)] <- "unknown"
outbreak_df$Host.disease[is.na(outbreak_df$Host.disease)] <- "unknown"
outbreak_df$Host.disease <- unclass(factor(outbreak_df$Host.disease))
outbreak_df$region[is.na(outbreak_df$region)] <-"unknown"
outbreak_df$region <- as.factor(outbreak_df$region)

# exploring the key parameters relating to the new outbreak value
positiveWeight <- 1.0 / (nrow(subset(outbreak_df, new_outbreak== 1)) / nrow(outbreak_df))
negativeWeight <- 1.0 / (nrow(subset(outbreak_df, new_outbreak== 0)) / nrow(outbreak_df))
modelWeights <- ifelse(outbreak_df$new_outbreak== 0,negativeWeight, positiveWeight)
tree_model <- rpart(factor(new_outbreak)~., data=outbreak_df[,-c(10:13)],weights = modelWeights,minbucket=3,minsplit=10,cp=0.003)
rpart.plot(tree_model)
plotcp(tree_model)
printcp(tree_model)

# modify the key parameters from the decision tree
outbreak_df$Stress.genotypes.count <- ifelse(outbreak_df$Stress.genotypes.count>=3,1,0)
outbreak_df$AMR.genotypes.count <- ifelse(outbreak_df$AMR.genotypes.count>=6,1,0)
outbreak_df$SNP.cluster.rank <- ifelse(outbreak_df$SNP.cluster.rank<=21,1,2)

# use those key parameters for further building predictive model on isolates numbers over time
model_df <- outbreak_df %>% select(c("AMR.genotypes.count","Stress.genotypes.count",
                                      "isolate_source_type","Host.category","SNP.cluster.rank",
                                      "region","Year","Month","Week","Serovar","new_outbreak")) %>%
  group_by(Serovar,Year,Month,Week,region) %>%
  # count number of isolate by groups, modify year variables starting 2017
  mutate(Outbreak_num=sum(new_outbreak),Year=Year-2016) %>%
  select(-c("new_outbreak"))

model_df$isolate_source_type <- levels(model_df$isolate_source_type)[model_df$isolate_source_type]

########################### Gaussian Linear Regression Model ###########################

###### predicting in months ######

# fit the model and checking the significance of variables
modelm_1 <-glm(Outbreak_num~.,model_df[,-9],family="gaussian")
#summary(modelm_1) 

# need some transformation on month variable
model_df <- model_df %>% mutate(Month=log(Month))
modelm_2 <- glm(Outbreak_num~.,model_df[,-9],family="gaussian")
# summary(modelm_2)

# relevel variables for model improving
model_df$region <- relevel(model_df$region, ref = "unknown")
modelm_3 <- glm(Outbreak_num~.,model_df[,-9],family="gaussian")
# summary(modelm_3)

# modify isolate source type of pork, beef into the category of other meat
model_df$isolate_source_type[model_df$isolate_source_type=="pork"] <- "other meat"
model_df$isolate_source_type[model_df$isolate_source_type=="beef"] <- "other meat"
model_df$isolate_source_type <- factor(model_df$isolate_source_type)
model_df$isolate_source_type <- relevel(model_df$isolate_source_type,ref="other meat")

modelm_4 <- glm(Outbreak_num~.,model_df[,-9],family="gaussian")
# summary(modelm_4)

# modify host category of Human into mammal and reptile into other
model_df$Host.category[model_df$Host.category=="Human"] <- "Mammal"
model_df$Host.category[model_df$Host.category=="Reptile"] <- "Other"
model_df$Host.category <- factor(model_df$Host.category)
model_df$Host.category <- relevel(model_df$Host.category,ref="unknown")

modelm_5 <- glm(Outbreak_num~.,model_df[,-9],family="gaussian")
# summary(modelm_5)

# transform the year variable
model_df <- model_df %>% mutate(Year=exp(Year))
modelm_6 <- glm(Outbreak_num~.,model_df[,-9],family="gaussian")
# summary(modelm_6) #final month predictive model


###### predicting in weeks ######

# fit the model and checking the significance of variables
modelw_1 <-glm(Outbreak_num~.,model_df[,-8],family="gaussian")
#summary(modelw_1) 

# transform week variable
model_df <- model_df %>% mutate(Week=log(Week))

modelw_2 <-glm(Outbreak_num~.,model_df[,-8],family="gaussian")
#summary(modelw_2) #final week predictive model

########################### data validation ###########################
# cross-validation: split into train-80% and test-20% using year variable
# dividing data into five sets, each with a training set and a test set.
train_index <- createDataPartition(model_df$Year, p=0.8, list = FALSE, times = 5)
train_df1 <- model_df[train_index[,1],]
train_df2 <- model_df[train_index[,2],]
train_df3 <- model_df[train_index[,3],]
train_df4 <- model_df[train_index[,4],]
train_df5 <- model_df[train_index[,5],]

test_df1 <- model_df[-train_index[,1],]
test_df2 <- model_df[-train_index[,2],]
test_df3 <- model_df[-train_index[,3],]
test_df4 <- model_df[-train_index[,4],]
test_df5 <- model_df[-train_index[,5],]

# fitting the model on all training data and average the coefficients
model1_coef <- cbind(summary(glm(Outbreak_num~.,train_df1[,-9],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df2[,-9],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df3[,-9],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df4[,-9],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df5[,-9],family="gaussian"))$coefficients[,1])

avg_model1_coef <- apply(model1_coef,1,mean)


# fitting the model on all training data and average the coefficients
model2_coef <- cbind(summary(glm(Outbreak_num~.,train_df1[,-8],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df2[,-8],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df3[,-8],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df4[,-8],family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df5[,-8],family="gaussian"))$coefficients[,1])

avg_model2_coef <- apply(model2_coef,1,mean)


########################### Evaluation ###########################

test_df1$predict <- predict(glm(Outbreak_num~.,train_df1[,-9],family="gaussian"),test_df1,type="response")
test_df2$predict <- predict(glm(Outbreak_num~.,train_df2[,-9],family="gaussian"),test_df2,type="response")
test_df3$predict <- predict(glm(Outbreak_num~.,train_df3[,-9],family="gaussian"),test_df3,type="response")
test_df4$predict <- predict(glm(Outbreak_num~.,train_df4[,-9],family="gaussian"),test_df4,type="response")
test_df5$predict <- predict(glm(Outbreak_num~.,train_df5[,-9],family="gaussian"),test_df5,type="response")

test_df1$predict2 <- predict(glm(Outbreak_num~.,train_df1[,-8],family="gaussian"),test_df1,type="response")
test_df2$predict2 <- predict(glm(Outbreak_num~.,train_df2[,-8],family="gaussian"),test_df2,type="response")
test_df3$predict2 <- predict(glm(Outbreak_num~.,train_df3[,-8],family="gaussian"),test_df3,type="response")
test_df4$predict2 <- predict(glm(Outbreak_num~.,train_df4[,-8],family="gaussian"),test_df4,type="response")
test_df5$predict2 <- predict(glm(Outbreak_num~.,train_df5[,-8],family="gaussian"),test_df5,type="response")

roc.model1_test1 <- roc(test_df1$Outbreak_num,test_df1$predict)
roc.model1_test2 <- roc(test_df2$Outbreak_num,test_df2$predict)
roc.model1_test3 <- roc(test_df3$Outbreak_num,test_df3$predict)
roc.model1_test4 <- roc(test_df4$Outbreak_num,test_df4$predict)
roc.model1_test5 <- roc(test_df5$Outbreak_num,test_df5$predict)

roc.model2_test1 <- roc(test_df1$Outbreak_num,test_df1$predict2)
roc.model2_test2 <- roc(test_df2$Outbreak_num,test_df2$predict2)
roc.model2_test3 <- roc(test_df3$Outbreak_num,test_df3$predict2)
roc.model2_test4 <- roc(test_df4$Outbreak_num,test_df4$predict2)
roc.model2_test5 <- roc(test_df5$Outbreak_num,test_df5$predict2)

plot(roc.model1_test1,main="ROC of 2 Models on One Training \nDataset",col="red")
lines(roc.model2_test1,col="blue")
legend("bottomright",c("Month Pred", "Week Pred"), lty=1, 
       col = c("red", "blue"),bty="n",cex=0.8)

AUC_df <- data.frame(test1=c(auc(roc.model1_test1),auc(roc.model2_test1)),
                     test2=c(auc(roc.model1_test2),auc(roc.model2_test1)),
                     test3=c(auc(roc.model1_test3),auc(roc.model2_test1)),
                     test4=c(auc(roc.model1_test4),auc(roc.model2_test1)),
                     test5=c(auc(roc.model1_test5),auc(roc.model2_test1)))

r_squared_df <- data.frame(test1=c(summary(lm(Outbreak_num~.,train_df1[,-8]))$adj.r.squared,summary(lm(Outbreak_num~.,train_df1[,-7]))$adj.r.squared),
test2=c(summary(lm(Outbreak_num~.,train_df2[,-8]))$adj.r.squared,summary(lm(Outbreak_num~.,train_df2[,-7]))$adj.r.squared),
test3=c(summary(lm(Outbreak_num~.,train_df3[,-8]))$adj.r.squared,summary(lm(Outbreak_num~.,train_df3[,-7]))$adj.r.squared),
test4=c(summary(lm(Outbreak_num~.,train_df4[,-8]))$adj.r.squared,summary(lm(Outbreak_num~.,train_df4[,-7]))$adj.r.squared),
test5=c(summary(lm(Outbreak_num~.,train_df5[,-8]))$adj.r.squared,summary(lm(Outbreak_num~.,train_df5[,-7]))$adj.r.squared))

AIC_df <- data.frame(test1=c(AIC(lm(Outbreak_num~.,train_df1[,-8])),AIC(lm(Outbreak_num~.,train_df1[,-7]))),
                     test2=c(AIC(lm(Outbreak_num~.,train_df2[,-8])),AIC(lm(Outbreak_num~.,train_df2[,-7]))),
                     test3=c(AIC(lm(Outbreak_num~.,train_df3[,-8])),AIC(lm(Outbreak_num~.,train_df3[,-7]))),
                     test4=c(AIC(lm(Outbreak_num~.,train_df4[,-8])),AIC(lm(Outbreak_num~.,train_df4[,-7]))),
                     test5=c(AIC(lm(Outbreak_num~.,train_df5[,-8])),AIC(lm(Outbreak_num~.,train_df5[,-7]))))

BIC_df <- data.frame(test1=c(BIC(lm(Outbreak_num~.,train_df1[,-8])),BIC(lm(Outbreak_num~.,train_df1[,-7]))),
                     test2=c(BIC(lm(Outbreak_num~.,train_df2[,-8])),BIC(lm(Outbreak_num~.,train_df2[,-7]))),
                     test3=c(BIC(lm(Outbreak_num~.,train_df3[,-8])),BIC(lm(Outbreak_num~.,train_df3[,-7]))),
                     test4=c(BIC(lm(Outbreak_num~.,train_df4[,-8])),BIC(lm(Outbreak_num~.,train_df4[,-7]))),
                     test5=c(BIC(lm(Outbreak_num~.,train_df5[,-8])),BIC(lm(Outbreak_num~.,train_df5[,-7]))))

Comparison_df <- data.frame(AUC=round(apply(AUC_df,1,mean),2),
                            R_squared=round(apply(r_squared_df,1,mean),2),
                            AIC=round(apply(AIC_df,1,mean),0),
                            BIC=round(apply(BIC_df,1,mean),0))

names(Comparison_df)<- c("AUC","Ajusted R Squared","AIC","BIC")
eval <- t(Comparison_df)
colnames(eval) <- c("Model of Months","Model of Weeks")

# evaluation of Gaussian linear regression model
eval %>%
  kbl() %>%
  kable_paper("hover", full_width = F)



########################### Example ###########################
# observed data for a specific group: Typhimurium Serovar in egg, chicken, turkey and poultry
x_vars_month <- model.matrix(Outbreak_num~.,model_df[,-9])
predict_month <- x_vars_month %*% avg_model1_coef
x_vars_week <- model.matrix(Outbreak_num~.,model_df[,-8])
predict_week <- x_vars_week %*% avg_model2_coef
example <- cbind(model_df,predict_month=predict_month,predict_week=predict_week) %>% 
  filter(Serovar=="Typhimurium",isolate_source_type %in% c("egg","chicken","turkey","poultry")) %>%
  mutate(Year=log(Year),Month=exp(Month),Week=exp(Week)) %>%
  group_by(isolate_source_type,Year,region) %>%
  summarize(obs_outbreak = sum(Outbreak_num),pred_m_outbreak=abs(sum(predict_month)),pred_w_outbreak=abs(sum(predict_week)))

ggplot(example,aes(Year+2016,obs_outbreak))+
  geom_smooth(aes(col=isolate_source_type))+ 
  labs(title = "Observed the Number of Outbreak of Typhimurium",
       x = "Year of observation",
       y = "Number of Observed Outbreak") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top")+
  facet_wrap(~region)

ggplot(example)+
  geom_smooth(aes(Year+2016,pred_m_outbreak,col=isolate_source_type))+ 
  geom_smooth(aes(Year+2016,pred_w_outbreak,col=isolate_source_type),linetype="dashed")+
  labs(title = "Predicted Number of Outbreak of Typhimurium",
       x = "Year of observation",
       y = "Number of Observed Outbreak") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top")+
  facet_wrap(~region)

