# libraries
library(tidyverse)
library(ggplot2)
library(lmtest)
library(rpart.plot)
library(kableExtra)

# load in dataset
isolates <- read.csv("~/Downloads/php2550finalproject-main/data/isolates2.csv")

# select interest independent variables
isolates2 <- isolates %>% select(c("new_outbreak","Host.category","Host.disease",
                                   "Host","AST.phenotypes","Virulence.genotypes",
                                   "isolate_source_type","N50","Length","SNP.cluster","Contigs",
                                   "Computed.types","TaxID","Stress.genotypes","WGS.prefix",
                                   "region","Year","Month","Week","Serovar","Isolate"))

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
           "SNP.cluster","region","Year","Month","Week","Serovar","Isolate")) %>%
  group_by(Serovar,Year,Month,Week,region) %>%
  # count number of isolate by groups, modify year variables starting 2017 as year 1
  mutate(Isolate_num = n(),Year=Year-2016) %>%
  select(-c("Isolate"))

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
model1 <-lm(Isolate_num~.,train_df1[,-8])
#summary(model1) 
 
# removing variable of Computed.types as it is insignificant
model1_new <- lm(Isolate_num~.,train_df1[,-c(1,8)])
#summary(model1_new) #recheck the p-values

# fitting the model on all training data and average the coefficients
model1_coef <- cbind(summary(lm(Isolate_num~.,train_df1[,-c(1,8)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df2[,-c(1,8)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df3[,-c(1,8)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df4[,-c(1,8)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df5[,-c(1,8)]))$coefficients[,1])

avg_model1_coef <- apply(model1_coef,1,mean)

########################### predicting in weeks ###########################

# fit the model on one training dataset and checking the significance of variables
model2 <-lm(Isolate_num~.,train_df1[,-7])
#summary(model2) # checking the significance of variables

# removing variable of Computed.types as it is insignificant
model2_new <- lm(Isolate_num~.,train_df1[,-c(1,7)])
#summary(model2_new) #recheck the p-values

# fitting the model on all training data and average the coefficients
model2_coef <- cbind(summary(lm(Isolate_num~.,train_df1[,-c(1,7)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df2[,-c(1,7)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df3[,-c(1,7)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df4[,-c(1,7)]))$coefficients[,1],
                     summary(lm(Isolate_num~.,train_df5[,-c(1,7)]))$coefficients[,1])

avg_model2_coef <- apply(model2_coef,1,mean)


#prediction for specific groups
#newdata <- with(model3, data.frame(region = c(7), Serovar=c("Enteritidis"), complete.year=c(2025),season=c(1,2,3,4),Isolation.type=c("environmental/other")))
#predict(model3, newdata)



