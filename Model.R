# Model Building
  # predict probability for a region/state of having a specific Serovar Samonella

# libraries
library(tidyverse)
library(ggplot2)
library(lubridate)
library(mice)
library(modelr)

# load in dataset
isolates <- read.csv("isolates2.csv")

########################### data processing ###########################

# replacing the NA
isolates$Collection.date[is.na(isolates$Collection.date)] <- ""

# process the year and month
for (i in 1:length(isolates$Collection.date)){
  # some string has the date in the version of "%mm/%dd/%YYYY", filtering those strings out
  if(str_detect(isolates$Collection.date[i],"/")==TRUE){
    # character length of the string for "%m/%d/%YYYY"
    if(nchar(isolates$Collection.date[i])==8){
      # extract the year at the position 5 to 8
      isolates$Collection.year[i] <- str_sub(isolates$Collection.date[i],5,8)
      # extract the month at the position 1
      isolates$Collection.month[i] <- str_sub(isolates$Collection.date[i],1,1)
    }
    # character length of the string for "%mm/%d/%YYYY" or "%m/%dd/%YYYY"
    else if(nchar(isolates$Collection.date[i])==9){
      # extract the year at the position 6 to 9
      isolates$Collection.year[i] <- str_sub(isolates$Collection.date[i],6,9)
      ≈
      isolates$Collection.month[i] <- str_sub(isolates$Collection.date[i],1,2)
      # when the string is in the format of "%m/%dd/%YYYY", it will include "/" in the month, so we need further extraction
      if(str_detect(isolates$Collection.month[i],"/")==TRUE){
        # removing the "/" in the month from the string "%m/%dd/%YYYY"
        isolates$Collection.month[i] <- str_extract(isolates$Collection.month[i],"(\\d)+")
      }
    }
    # character length of the string for "%mm/%dd/%YYYY"
    else if(nchar(isolates$Collection.date[i])==10){
      # extract the year at the position 7 to 10
      isolates$Collection.year[i] <- str_sub(isolates$Collection.date[i],7,10)
      # extract the year at the position 1 to 2
      isolates$Collection.month[i] <- str_sub(isolates$Collection.date[i],1,2)
    }
  }
  # other dates are either missing or in the format as "%Y-%m-%d"
  else if(str_detect(isolates$Collection.date[i],"/")==FALSE){
    # extract the year at position 1 to 4
    isolates$Collection.year[i] <- str_sub(isolates$Collection.date[i],1,4)
    # extract the month at position 6 to 7
    isolates$Collection.month[i] <- str_sub(isolates$Collection.date[i],6,7)
  }
}
# convert the class of year and month to numeric from characters
isolates$Collection.month <- as.numeric(isolates$Collection.month)
isolates$Collection.year <- as.numeric(isolates$Collection.year)

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
model1 <-glm(freq~region+complete.year+Isolation.type+complete.month,family=binomial(),isolates_clean)
summary(model1)

# month not significance and so thinking about transformation on the variable
# group into seasons instead of by months
for (i in 1:nrow(isolates_clean)){
  # if the month is Jan, Feb or March
  if (isolates_clean$complete.month[i] %in% c(1,2,3)){
    # then new variable season has the value 1
    isolates_clean$season[i] <- 1
  }# if the month is April, May, Jun
  else if (isolates_clean$complete.month[i] %in% c(4,5,6)){
    # then new variable season has the value 2
    isolates_clean$season[i] <- 2
  }# if the month is Jul, Aug, Sep
  else if(isolates_clean$complete.month[i] %in% c(7,8,9)){
    # then new variable season has the value 3
    isolates_clean$season[i] <- 3
  }# if the month is Oct, Nov, Dec
  else if (isolates_clean$complete.month[i] %in% c(10,11,12)){
    # then new variable season has the value 4
    isolates_clean$season[i] <- 4
  }
}

# fit the model with new variable season on overall dataset
model2 <-glm(freq~region+complete.year+Isolation.type+season,family=binomial(),isolates_clean)
summary(model2)

# p-value of season gets smaller and closer to 0.05
model3 <-glm(freq~region+complete.year*season+Isolation.type,family=binomial(),isolates_clean)
summary(model3)
# both season and the interaction term of season and year shows significance

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

# prediction for specific groups
# newdata <- with(model3, data.frame(region = c(7), Season=c(1,2,3,4))





