# Model Building
  # predict probability for a region/state of having a specific Serovar Samonella

# libraries
library(tidyverse)
library(ggplot2)
library(lubridate)
library(mice)
library(modelr)
library(lmtest)

# load in dataset
isolates <- read.csv("isolates2.csv")


########################### imputations before fitting the model ###########################
# 45.31% missing data of regions
# 2.24% missing data of collection year)
# 54.1% missing data of collection month

# subset the data for imputations as the large dataset would run long time
month_df <- isolates %>% select(c("Serovar","Host.category","region","Collection.year","Collection.month"))%>% group_by(Serovar,Host.category,Collection.year)

# using predictive mean matching for multiple imputations
imputed_Data_month <- mice(month_df, m=5, maxit = 50, method = 'pmm', seed = 500)

# Add the data back to original data using one of the iterations
completeData_month <- complete(imputed_Data_month, 3)

# check if full replacement of missing data in Collection.month
# sum(is.na(completeData_month$Collection.month))/nrow(isolates)

# add on complete data of month and year onto the big dataset
isolates$complete.month <- completeData_month$Collection.month
isolates$complete.year <- completeData_month$Collection.year

# compute the frequency of each Serovar in each region in each year
isolates_clean <- isolates %>%
  #count the number for each Serovar for each isolation type in a region for a month and year
  group_by(Serovar, region, complete.year,Isolation.type,complete.month)%>%
  mutate(cat_count=n()) %>%
  ungroup() %>%
  #count the number for total Serovars for each isolation type in a region for a month and year
  group_by(region, complete.year,Isolation.type,complete.month) %>%
  mutate(serovar_count=n(),freq=cat_count/serovar_count)
isolates_clean$id <- c(1:nrow(isolates_clean))

# fit the model selecting variables in overall
model1 <-lm(freq~region+complete.year+Isolation.type+complete.month,isolates_clean)
summary(model1)

# month not significance and so thinking about transformation on the variable
# group into seasons instead of by months
for (i in 1:nrow(isolates_clean)){
  # if the month is Mar, Apr, May
  if (isolates_clean$complete.month[i] %in% c(3,4,5)){
    # then new variable season has the value 1
    isolates_clean$season[i] <- 1
  }# if the month is Jun, Jul, Aug
  else if (isolates_clean$complete.month[i] %in% c(6,7,8)){
    # then new variable season has the value 2
    isolates_clean$season[i] <- 2
  }# if the month is Sep,Oct,Nov
  else if(isolates_clean$complete.month[i] %in% c(9,10,11)){
    # then new variable season has the value 3
    isolates_clean$season[i] <- 3
  }# if the month is Dec, Jan, Feb
  else if (isolates_clean$complete.month[i] %in% c(12,1,2)){
    # then new variable season has the value 4
    isolates_clean$season[i] <- 4
  }
}

# compute the frequency but with each season instead of each month
isolates_clean <- isolates_clean %>%
  select(-c(freq)) %>%
  group_by(Serovar, region, complete.year,Isolation.type,season)%>%
  mutate(cat_count=n()) %>%
  ungroup() %>%
  group_by(region, complete.year,Isolation.type,season) %>%
  mutate(serovar_count=n(),freq=cat_count/serovar_count)


# fit the model with new variable season on overall dataset: season variable shows significant as p-value <0.05
model2 <-lm(freq~region+complete.year+Isolation.type+season,isolates_clean)
summary(model2)

# p-value of interaction term of year and season shows significance
model3 <-lm(freq~region+log(complete.year*season)+Isolation.type,isolates_clean)
summary(model3)

# comparing if there is the significant difference between model2 and model3
lrest(model2,model3)

########################### data validation ###########################
# cross-validation: split into train-80% and test-20%
# dividing data into five sets, each with a training set and a test set.
train_ind <- createDataPartition(isolates_clean$id, p=0.8, list = FALSE, times = 5)
train_df1 <- isolates_clean[train[,1],]
train_df2 <- isolates_clean[train[,2],]
train_df3 <- isolates_clean[train[,3],]
train_df4 <- isolates_clean[train[,4],]
train_df5 <- isolates_clean[train[,5],]

test_df1 <- isolates_clean[-train[,1],]
test_df2 <- isolates_clean[-train[,2],]
test_df3 <- isolates_clean[-train[,3],]
test_df4 <- isolates_clean[-train[,4],]
test_df5 <- isolates_clean[-train[,5],]

#prediction for specific groups
newdata <- with(model3, data.frame(region = c(7), Serovar=c("Enteritidis"), complete.year=c(2025),season=c(1,2,3,4),Isolation.type=c("environmental/other")))
predict(model3, newdata)





