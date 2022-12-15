# libraries
library(tidyverse)
library(ggplot2)
library(rpart.plot)
library(kableExtra)
library(caret)
library(pROC)
library(DescTools)

# load in dataset
isolates <- read.csv("data/isolates2.csv")


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
model_df_m <- outbreak_df %>% select(c("AMR.genotypes.count","Stress.genotypes.count",
                                      "isolate_source_type","Host.category","SNP.cluster.rank",
                                      "region","Year","Month","Serovar","new_outbreak")) %>%
  group_by(Serovar,Year,Month,region,isolate_source_type) %>%
  # count number of isolate by groups, modify year variables starting 2017
  mutate(Outbreak_num=sum(new_outbreak),Year=Year-2016) %>%
  select(-c("new_outbreak"))

model_df_m$isolate_source_type <- levels(model_df_m$isolate_source_type)[model_df_m$isolate_source_type]

model_df_w <- outbreak_df %>% select(c("AMR.genotypes.count","Stress.genotypes.count",
                                       "isolate_source_type","Host.category","SNP.cluster.rank",
                                       "region","Year","Week","Serovar","new_outbreak")) %>%
  group_by(Serovar,Year,Week,region,isolate_source_type) %>%
  # count number of isolate by groups, modify year variables starting 2017
  mutate(Outbreak_num=sum(new_outbreak),Year=Year-2016) %>%
  select(-c("new_outbreak"))

model_df_w$isolate_source_type <- levels(model_df_w$isolate_source_type)[model_df_w$isolate_source_type]


########################### Gaussian Linear Regression Model ###########################

###### predicting in months ######

# fit the model and checking the significance of variables
modelm_1 <-glm(Outbreak_num~.,model_df_m,family="gaussian")
#summary(modelm_1) 

# modify isolate source type of pork, beef into the category of other meat
model_df_m$isolate_source_type <- factor(model_df_m$isolate_source_type)
model_df_m$isolate_source_type <- relevel(model_df_m$isolate_source_type,ref="eggs")

modelm_2 <- glm(Outbreak_num~.,model_df_m,family="gaussian")
#summary(modelm_2)

# modify host category of Human into mammal and reptile into other
model_df_m$Host.category[model_df_m$Host.category=="Human"] <- "Mammal"
model_df_m$Host.category[model_df_m$Host.category=="Reptile"] <- "Other"
model_df_m$Host.category <- factor(model_df_m$Host.category)
model_df_m$Host.category <- relevel(model_df_m$Host.category,ref="Mammal")

modelm_3 <- glm(Outbreak_num~.,model_df_m,family="gaussian")
# summary(modelm_3) #final predictive model in month


###### predicting in weeks ######

# fit the model and checking the significance of variables
modelw_1 <-glm(Outbreak_num~.,model_df_w,family="gaussian")
#summary(modelw_1) 

# relevel variables for model improving
model_df_w$region <- relevel(model_df_w$region, ref ="4")
modelw_2 <- glm(Outbreak_num~.,model_df_w,family="gaussian")
#summary(modelw_2)

# modify isolate source type of pork, beef into the category of other meat
model_df_w$isolate_source_type <- factor(model_df_w$isolate_source_type)
model_df_w$isolate_source_type <- relevel(model_df_w$isolate_source_type,ref="feces")
modelw_3 <- glm(Outbreak_num~.,model_df_w,family="gaussian")
#summary(modelw_3)

# modify host category of Human into mammal and reptile into other
model_df_w$Host.category[model_df_w$Host.category=="Human"] <- "Mammal"
model_df_w$Host.category[model_df_w$Host.category=="Reptile"] <- "Other"
model_df_w$Host.category <- factor(model_df_w$Host.category)
model_df_w$Host.category <- relevel(model_df_w$Host.category,ref="Other")

modelw_4 <- glm(Outbreak_num~.,model_df_w,family="gaussian")
#summary(modelw_4)

# remove Stress.genotype.count for its insignificance
model_df_w <- model_df_w %>% select(-c("Stress.genotypes.count"))

modelw_5 <-glm(Outbreak_num~.,model_df_w,family="gaussian")
#summary(modelw_5) #final predictive model in weeks


########################### data validation ###########################
# cross-validation: split into train-80% and test-20% using year variable
# dividing data into five sets, each with a training set and a test set.
set.seed(123)
train_index <- createDataPartition(model_df_m$Year, p=0.8, list = FALSE, times = 5)
train_df1_m <- model_df_m[train_index[,1],]
train_df2_m <- model_df_m[train_index[,2],]
train_df3_m <- model_df_m[train_index[,3],]
train_df4_m <- model_df_m[train_index[,4],]
train_df5_m <- model_df_m[train_index[,5],]

test_df1_m <- model_df_m[-train_index[,1],]
test_df2_m <- model_df_m[-train_index[,2],]
test_df3_m <- model_df_m[-train_index[,3],]
test_df4_m <- model_df_m[-train_index[,4],]
test_df5_m <- model_df_m[-train_index[,5],]

train_df1_w <- model_df_w[train_index[,1],]
train_df2_w <- model_df_w[train_index[,2],]
train_df3_w <- model_df_w[train_index[,3],]
train_df4_w <- model_df_w[train_index[,4],]
train_df5_w <- model_df_w[train_index[,5],]

test_df1_w <- model_df_w[-train_index[,1],]
test_df2_w <- model_df_w[-train_index[,2],]
test_df3_w <- model_df_w[-train_index[,3],]
test_df4_w <- model_df_w[-train_index[,4],]
test_df5_w <- model_df_w[-train_index[,5],]

# fitting the model on all training data and average the coefficients
model1_coef <- cbind(summary(glm(Outbreak_num~.,train_df1_m,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df2_m,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df3_m,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df4_m,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df5_m,family="gaussian"))$coefficients[,1])

avg_model1_coef <- apply(model1_coef,1,mean)


# fitting the model on all training data and average the coefficients
model2_coef <- cbind(summary(glm(Outbreak_num~.,train_df1_w,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df2_w,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df3_w,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df4_w,family="gaussian"))$coefficients[,1],
                     summary(glm(Outbreak_num~.,train_df5_w,family="gaussian"))$coefficients[,1])

avg_model2_coef <- apply(model2_coef,1,mean)


########################### Evaluation ###########################

test_df1_m$predict <- predict(glm(Outbreak_num~.,train_df1_m,family="gaussian"),test_df1_m,type="response")
test_df2_m$predict <- predict(glm(Outbreak_num~.,train_df2_m,family="gaussian"),test_df2_m,type="response")
test_df3_m$predict <- predict(glm(Outbreak_num~.,train_df3_m,family="gaussian"),test_df3_m,type="response")
test_df4_m$predict <- predict(glm(Outbreak_num~.,train_df4_m,family="gaussian"),test_df4_m,type="response")
test_df5_m$predict <- predict(glm(Outbreak_num~.,train_df5_m,family="gaussian"),test_df5_m,type="response")

test_df1_w$predict2 <- predict(glm(Outbreak_num~.,train_df1_w,family="gaussian"),test_df1_w,type="response")
test_df2_w$predict2 <- predict(glm(Outbreak_num~.,train_df2_w,family="gaussian"),test_df2_w,type="response")
test_df3_w$predict2 <- predict(glm(Outbreak_num~.,train_df3_w,family="gaussian"),test_df3_w,type="response")
test_df4_w$predict2 <- predict(glm(Outbreak_num~.,train_df4_w,family="gaussian"),test_df4_w,type="response")
test_df5_w$predict2 <- predict(glm(Outbreak_num~.,train_df5_w,family="gaussian"),test_df5_w,type="response")

roc.model1_test1 <- roc(test_df1_m$Outbreak_num,test_df1_m$predict)
roc.model1_test2 <- roc(test_df2_m$Outbreak_num,test_df2_m$predict)
roc.model1_test3 <- roc(test_df3_m$Outbreak_num,test_df3_m$predict)
roc.model1_test4 <- roc(test_df4_m$Outbreak_num,test_df4_m$predict)
roc.model1_test5 <- roc(test_df5_m$Outbreak_num,test_df5_m$predict)

roc.model2_test1 <- roc(test_df1_w$Outbreak_num,test_df1_w$predict2)
roc.model2_test2 <- roc(test_df2_w$Outbreak_num,test_df2_w$predict2)
roc.model2_test3 <- roc(test_df3_w$Outbreak_num,test_df3_w$predict2)
roc.model2_test4 <- roc(test_df4_w$Outbreak_num,test_df4_w$predict2)
roc.model2_test5 <- roc(test_df5_w$Outbreak_num,test_df5_w$predict2)

plot(roc.model1_test1,main="ROC of 2 Models on One Training \nDataset",col="red")
lines(roc.model2_test1,col="blue")
legend("bottomright",c("Month Pred", "Week Pred"), lty=1, 
       col = c("red", "blue"),bty="n",cex=0.8)

AUC_df <- data.frame(test1=c(auc(roc.model1_test1),auc(roc.model2_test1)),
                     test2=c(auc(roc.model1_test2),auc(roc.model2_test1)),
                     test3=c(auc(roc.model1_test3),auc(roc.model2_test1)),
                     test4=c(auc(roc.model1_test4),auc(roc.model2_test1)),
                     test5=c(auc(roc.model1_test5),auc(roc.model2_test1)))

r_squared_df <- data.frame(test1=c(summary(lm(Outbreak_num~.,train_df1_m))$adj.r.squared,summary(lm(Outbreak_num~.,train_df1_w))$adj.r.squared),
test2=c(summary(lm(Outbreak_num~.,train_df2_m))$adj.r.squared,summary(lm(Outbreak_num~.,train_df2_w))$adj.r.squared),
test3=c(summary(lm(Outbreak_num~.,train_df3_m))$adj.r.squared,summary(lm(Outbreak_num~.,train_df3_w))$adj.r.squared),
test4=c(summary(lm(Outbreak_num~.,train_df4_m))$adj.r.squared,summary(lm(Outbreak_num~.,train_df4_w))$adj.r.squared),
test5=c(summary(lm(Outbreak_num~.,train_df5_m))$adj.r.squared,summary(lm(Outbreak_num~.,train_df5_w))$adj.r.squared))

AIC_df <- data.frame(test1=c(AIC(lm(Outbreak_num~.,train_df1_m)),AIC(lm(Outbreak_num~.,train_df1_w))),
                     test2=c(AIC(lm(Outbreak_num~.,train_df2_m)),AIC(lm(Outbreak_num~.,train_df2_w))),
                     test3=c(AIC(lm(Outbreak_num~.,train_df3_m)),AIC(lm(Outbreak_num~.,train_df3_w))),
                     test4=c(AIC(lm(Outbreak_num~.,train_df4_m)),AIC(lm(Outbreak_num~.,train_df4_w))),
                     test5=c(AIC(lm(Outbreak_num~.,train_df5_m)),AIC(lm(Outbreak_num~.,train_df5_w))))

BIC_df <- data.frame(test1=c(BIC(lm(Outbreak_num~.,train_df1_m)),BIC(lm(Outbreak_num~.,train_df1_w))),
                     test2=c(BIC(lm(Outbreak_num~.,train_df2_m)),BIC(lm(Outbreak_num~.,train_df2_w))),
                     test3=c(BIC(lm(Outbreak_num~.,train_df3_m)),BIC(lm(Outbreak_num~.,train_df3_w))),
                     test4=c(BIC(lm(Outbreak_num~.,train_df4_m)),BIC(lm(Outbreak_num~.,train_df4_w))),
                     test5=c(BIC(lm(Outbreak_num~.,train_df5_m)),BIC(lm(Outbreak_num~.,train_df5_w))))

Comparison_df <- data.frame(AUC=round(apply(AUC_df,1,mean),2),
                            R_squared=round(apply(r_squared_df,1,mean),2),
                            AIC=round(apply(AIC_df,1,mean),0),
                            BIC=round(apply(BIC_df,1,mean),0))

names(Comparison_df)<- c("AUC","Ajusted R Squared","AIC","BIC")
eval <- t(Comparison_df)
colnames(eval) <- c("Model of Months","Model of Weeks")

# evaluation of Gaussian linear regression model
eval %>%
  kbl(caption="Comparison Table of Measuring Metrics of 2 Models") %>%
  kable_classic(full_width = F, html_font = "Cambria")



########################### Example ###########################
# observed data for a specific group: Typhimurium Serovar in egg, chicken, turkey and poultry
x_vars_month <- model.matrix(Outbreak_num~.,model_df_m)
predict_month <- x_vars_month %*% avg_model1_coef
x_vars_week <- model.matrix(Outbreak_num~.,model_df_w)
predict_week <- x_vars_week %*% avg_model2_coef
outbreak_df$isolate_source_type <- levels(outbreak_df$isolate_source_type)[outbreak_df$isolate_source_type]
example <- cbind(outbreak_df,predict_month=predict_month,predict_week=predict_week) %>% 
  filter(Serovar=="Typhimurium" & isolate_source_type %in% c("egg","chicken","turkey","poultry")) %>%
  group_by(isolate_source_type,Year,region) %>%
  summarize(obs_outbreak = sum(new_outbreak),pred_m_outbreak=abs(sum(predict_month)),pred_w_outbreak=abs(sum(predict_week)))

ggplot(example,aes(Year,obs_outbreak))+
  geom_smooth(aes(col=isolate_source_type))+ 
  labs(title = "Observed the Number of Outbreak of Typhimurium",
       x = "Year of observation",
       y = "Number of Observed Outbreak") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top")+
  facet_wrap(~region)

ggplot(example)+
  geom_smooth(aes(Year,pred_m_outbreak,col=isolate_source_type))+ 
  geom_smooth(aes(Year,pred_w_outbreak,col=isolate_source_type),linetype="dashed")+
  labs(title = "Predicted Number of Outbreak of Typhimurium",
       x = "Year of observation",
       y = "Number of Observed Outbreak") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top")+
  facet_wrap(~region)

