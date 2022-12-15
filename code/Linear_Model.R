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
isolates <- read.csv("data/isolates2.csv")

# select interest independent variables
isolates2 <- isolates %>% select(c("new_outbreak","Host.category","Host.disease",
                                   "Host","AST.phenotypes","Virulence.genotypes",
                                   "isolate_source_type","N50","Length","SNP.cluster","Contigs",
                                   "Computed.types","TaxID","Stress.genotypes","WGS.prefix",
                                   "region","Year","Month","Week","Serovar"))

# modifying the categorical variable into numeric variable 
# converting variables with NA into unknown
isolates2$isolate_source_type <- unclass(factor(isolates2$isolate_source_type))
isolates2$Host.category[is.na(isolates2$Host.category)] <- "unknown"
isolates2$Host.category <- unclass(factor(isolates2$Host.category))
isolates2$Host.disease[is.na(isolates2$Host.disease)] <- "unknown"
isolates2$Host.disease <- unclass(factor(isolates2$Host.disease))
isolates2$Host[is.na(isolates2$Host)] <- "unknown"
isolates2$Host <- unclass(factor(isolates2$Host))
isolates2$AST.phenotypes[is.na(isolates2$AST.phenotypes)] <- "unknown"
isolates2$AST.phenotypes<- unclass(factor(isolates2$AST.phenotypes))
isolates2$Virulence.genotypes[is.na(isolates2$Virulence.genotypes)] <-"unknown"
isolates2$Virulence.genotypes<- unclass(factor(isolates2$Virulence.genotypes))
isolates2$SNP.cluster[is.na(isolates2$SNP.cluster)] <- "unknown"
isolates2$SNP.cluster<- unclass(factor(isolates2$SNP.cluster))
isolates2$Computed.types[is.na(isolates2$Computed.types)] <- "unknown"
isolates2$Computed.types <- unclass(factor(isolates2$Computed.types))
isolates2$Stress.genotypes[is.na(isolates2$Stress.genotypes)]<- "unknown"
isolates2$Stress.genotypes<- unclass(factor(isolates2$Stress.genotypes))
isolates2$WGS.prefix[is.na(isolates2$WGS.prefix)] <- "unknown"
isolates2$WGS.prefix<- unclass(factor(isolates2$WGS.prefix))
isolates2$region[is.na(isolates2$region)] <-"unknown"
isolates2$region <- as.factor(isolates2$region)

# exploring the key parameters relating to the new outbreak value
positiveWeight <- 1.0 / (nrow(subset(isolates2, new_outbreak== 1)) / nrow(isolates2))
negativeWeight <- 1.0 / (nrow(subset(isolates2, new_outbreak== 0)) / nrow(isolates2))
modelWeights <- ifelse(isolates2$new_outbreak== 0,negativeWeight, positiveWeight)
tree_model <- rpart(factor(new_outbreak)~., data=isolates2[,-c(16:21)],weights = modelWeights,minbucket=3,minsplit=10,cp=0.003)
rpart.plot(tree_model)
plotcp(tree_model)
printcp(tree_model)

isolates2$isolate_source_type <- ifelse(isolates2$isolate_source_type>=3,1,0)
isolates2$Computed.types <- ifelse(isolates2$Computed.types<48,1,0)
isolates2$Host.category <- ifelse(isolates2$Host.category>=5,1,0)
isolates2$SNP.cluster <- case_when((isolates2$SNP.cluster<479)~1,(isolates2$SNP.cluster>=868)~2,
                                   (isolates2$SNP.cluster<868&isolates2$SNP.cluster>=479)~0)

# use those key parameters for further building predictive model on isolates numbers over time
isolates_df <- isolates2 %>% select(c("Computed.types","isolate_source_type","Host.category",
           "SNP.cluster","region","Year","Month","Week","Serovar","new_outbreak")) %>%
  group_by(Serovar,Year,Month,Week,region) %>%
  # count number of isolate by groups, modify year variables starting 2017 as year 1
  mutate(Outbreak_num=sum(new_outbreak)) %>%
  select(-c("new_outbreak"))

isolates_df$Host.category[isolates_df$Host.category==1] <- levels(isolates_df$Host.category)[1]
isolates_df$Host.category[isolates_df$Host.category==2] <- levels(isolates_df$Host.category)[2]
isolates_df$Host.category[isolates_df$Host.category==3] <- levels(isolates_df$Host.category)[3]
isolates_df$Host.category[isolates_df$Host.category==4] <- levels(isolates_df$Host.category)[4]
isolates_df$Host.category[isolates_df$Host.category==5] <- levels(isolates_df$Host.category)[5]
isolates_df$Host.category[isolates_df$Host.category==6] <- levels(isolates_df$Host.category)[6]
isolates_df$Host.category[isolates_df$Host.category==7] <- levels(isolates_df$Host.category)[7]

isolates_df$isolate_source_type[isolates_df$isolate_source_type==1] <- levels(isolates_df$isolate_source_type)[1]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==2] <- levels(isolates_df$isolate_source_type)[2]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==3] <- levels(isolates_df$isolate_source_type)[3]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==4] <- levels(isolates_df$isolate_source_type)[4]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==5] <- levels(isolates_df$isolate_source_type)[5]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==6] <- levels(isolates_df$isolate_source_type)[6]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==7] <- levels(isolates_df$isolate_source_type)[7]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==8] <- levels(isolates_df$isolate_source_type)[8]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==9] <- levels(isolates_df$isolate_source_type)[9]
isolates_df$isolate_source_type[isolates_df$isolate_source_type==10] <- levels(isolates_df$isolate_source_type)[10]



########################### data validation ###########################
# cross-validation: split into train-80% and test-20%
# dividing data into five sets, each with a training set and a test set.
train_index <- createDataPartition(isolates_df$Serovar, p=0.8, list = FALSE, times = 5)
train_df1 <- isolates_df[train_index[,1],]
train_df2 <- isolates_df[train_index[,2],]
train_df3 <- isolates_df[train_index[,3],]
train_df4 <- isolates_df[train_index[,4],]
train_df5 <- isolates_df[train_index[,5],]

test_df1 <- isolates_df[-train_index[,1],]
test_df2 <- isolates_df[-train_index[,2],]
test_df3 <- isolates_df[-train_index[,3],]
test_df4 <- isolates_df[-train_index[,4],]
test_df5 <- isolates_df[-train_index[,5],]

########################### predicting in months ###########################

# fit the model on one training dataset and checking the significance of variables
model1 <-lm(Outbreak_num~.,train_df1[,-8])
#summary(model1) 

# fitting the model on all training data and average the coefficients
model1_coef <- cbind(summary(lm(Outbreak_num~.,train_df1[,-8]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df2[,-8]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df3[,-8]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df4[,-8]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df5[,-8]))$coefficients[,1])

avg_model1_coef <- apply(model1_coef,1,mean)



########################### predicting in weeks ###########################

# fit the model on one training dataset and checking the significance of variables
model2 <-lm(Outbreak_num~.,train_df1[,-7])
#summary(model2) # checking the significance of variables


# fitting the model on all training data and average the coefficients
model2_coef <- cbind(summary(lm(Outbreak_num~.,train_df1[,-7]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df2[,-7]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df3[,-7]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df4[,-7]))$coefficients[,1],
                     summary(lm(Outbreak_num~.,train_df5[,-7]))$coefficients[,1])

avg_model2_coef <- apply(model2_coef,1,mean)


########################### Example ###########################
# observed data for a specific group: Region 8, Typhimurium Serovar
x_vars <- model.matrix(Outbreak_num~.,isolates_df[,-8])
predict_month <- x_vars %*% avg_model1_coef
observed <- cbind(isolates_df,predict=predict_month) %>% filter(Serovar=="Typhimurium" & region=="8"& Computed.types==31 
                                   &(isolate_source_type %in% c("beef","chicken","eggs","other meat","poultry","turkey")))

ggplot(observed,aes(Year,Outbreak_num))+
  geom_smooth(aes(col=isolate_source_type))+ 
  labs(title = "Observed Typhimurium Outbreak in Region 8 across years",
                                                x = "Year of observation",
                                                y = "Number of Observed Outbreak") +
  ylim(c(0,10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top")


#prediction for specific groups 6 years later: region 8, Typhimurium Serovar, isolate source of poultry, chicken, and eggs, host from environment
ggplot(observed,aes(Year+6,predict))+
  geom_smooth(aes(col=isolate_source_type))+ 
  labs(title = "Predicted Typhimurium Outbreak in Region 8 across years",
       x = "Year of observation",
       y = "Number of Observed Outbreak") +
  ylim(c(0,10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top")

########################### Evaluation ###########################

test_df1$predict <- predict(model1,test_df1,type="response")
test_df2$predict <- predict(model1,test_df2,type="response")
test_df3$predict <- predict(model1,test_df3,type="response")
test_df4$predict <- predict(model1,test_df4,type="response")
test_df5$predict <- predict(model1,test_df5,type="response")

test_df1$predict2 <- predict(model2,test_df1,type="response")
test_df2$predict2 <- predict(model2,test_df2,type="response")
test_df3$predict2 <- predict(model2,test_df3,type="response")
test_df4$predict2 <- predict(model2,test_df4,type="response")
test_df5$predict2 <- predict(model2,test_df5,type="response")

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

#plot(roc.model1,main="ROC of Model1 on \n the Training Dataset")

AUC_df <- data.frame(test1=c(auc(roc.model1_test1),auc(roc.model2_test1)),
                     test2=c(auc(roc.model1_test2),auc(roc.model2_test1)),
                     test3=c(auc(roc.model1_test3),auc(roc.model2_test1)),
                     test4=c(auc(roc.model1_test4),auc(roc.model2_test1)),
                     test5=c(auc(roc.model1_test5),auc(roc.model2_test1)))

bs.model1_test1 <- BrierScore(test_df1$Outbreak_num,test_df1$predict)
bs.model1_test2 <- BrierScore(test_df2$Outbreak_num,test_df2$predict)
bs.model1_test3 <- BrierScore(test_df3$Outbreak_num,test_df3$predict)
bs.model1_test4 <- BrierScore(test_df4$Outbreak_num,test_df4$predict)
bs.model1_test5 <- BrierScore(test_df5$Outbreak_num,test_df5$predict)

bs.model2_test1 <- BrierScore(test_df1$Outbreak_num,test_df1$predict2)
bs.model2_test2 <- BrierScore(test_df2$Outbreak_num,test_df2$predict2)
bs.model2_test3 <- BrierScore(test_df3$Outbreak_num,test_df3$predict2)
bs.model2_test4 <- BrierScore(test_df4$Outbreak_num,test_df4$predict2)
bs.model2_test5 <- BrierScore(test_df5$Outbreak_num,test_df5$predict2)

Brier_score_df <- data.frame(test1=c(bs.model1_test1,bs.model2_test1),
                             test2=c(bs.model1_test2,bs.model2_test2),
                             test3=c(bs.model1_test3,bs.model2_test3),
                             test4=c(bs.model1_test4,bs.model2_test4),
                             test5=c(bs.model1_test5,bs.model2_test5))


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

Comparison_df <- data.frame(AUC=round(apply(AUC_df,1,mean),2),Brier_score=round(apply(Brier_score_df,1,mean),2),
                            R_squared=round(apply(r_squared_df,1,mean),2),AIC=round(apply(AIC_df,1,mean),0),
                            BIC=round(apply(BIC_df,1,mean),0))

names(Comparison_df)<- c("AUC","Brier score","Ajusted R Squared","AIC","BIC")
eval <- t(Comparison_df)
colnames(eval) <- c("Model of Months","Model of Weeks")

eval %>%
  kbl() %>%
  kable_paper("hover", full_width = F)
