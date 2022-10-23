library (readr)
library(naniar)
library(ggplot2)
library(tidyverse)
library(VIM)
library(mice)
#Read the dataset from the Github website of raw material
urlfile="https://raw.githubusercontent.com/taliajs/php2550finalproject/main/isolates.csv"
isolates<-read_csv(url(urlfile))

# Visualization for percent missing
missing_val <- isolates %>%
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, isna) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total * 100)
levels <- (missing_val %>% filter(isna == T) %>%
             arrange(desc(pct)))$key
percentage.plot <- missing_val %>%
  ggplot() +
  geom_bar(aes(x = reorder(key, desc(pct)),
               y = pct, fill=isna),
           stat = 'identity', alpha=0.8) +
  scale_x_discrete(limits = levels) +
  scale_fill_manual(name = "",
                    values = c('steelblue', 'tomato3'),
                    labels = c("Present", "Missing")) +
  coord_flip() +
  labs(title = "Percentage of missing values",
       x = 'Variable', y = "% of missing values")+
  theme(plot.title = element_text(size = 10),axis.text.y = element_text(size=6))
percentage.plot

# Visualization for each row
row.plot <- isolates %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  ggplot(aes(key, id, fill = isna)) +
  geom_raster(alpha=0.8) +
  scale_fill_manual(name = "",
                    values = c('steelblue', 'tomato3'),
                    labels = c("Present", "Missing")) +
  scale_x_discrete(limits = levels) +
  labs(x = "Variable",
       y = "Row Number", title = "Missing values in rows") +
  coord_flip()+
  theme(axis.text.y = element_text(size=6))
row.plot

#calculate what percentage of data is missing from each variable
#p_miss <- function(x) {round(sum(is.na(x))/length(x)*100,2)}
#apply(isolates, 2, p_miss)

# Get number of missings per variable (n and %)
miss_var_summary(isolates)
miss_var_table(isolates)

#explore if there is the missing relation between isolation source and location
#Blue values are observed values and red ones are missing values
#both categorical character variables transform into categorical numeric variables
marginplot(isolates[,c('Location', 'Isolation source')])

#explore if there is the missing relation between isolation type and location
marginplot(isolates[,c('Location', 'Isolation type')])

#corresponding location character values with the numeric version
#unique(isolates$Location)

# check if a mixture of variables has missing values simultaneously
gg_miss_upset(isolates)
#There is a substantial number of cases in which some missings happen to occur across certain variables. 
#We thought that the missing data type is not missing completely at random (MCAR). 

#the impact of the presence of missing values in isolation type on year
isolates_test <- isolates %>%
  mutate(iso_type_na = is.na(`Isolation type`))

na_year_2000 <- isolates_test %>%
  filter(`Collection date` <= "2000") %>%
  pull(iso_type_na)

na_year_2020 <- isolates_test %>%
  filter(`Collection date` > "2000") %>%
  pull(iso_type_na)

#check whether the percentage of missings in isolation type differ per level of collection year
t.test(na_year_2000, na_year_2020)
#The percentage of missing sin isolation type has p-value smaller than 0.05.
# We can reject null hypothesis as there is significant differnce between isolation type missingness for collecting before 2000 and after 2000.
#Then it is expected for some reasons that the percentage of missing values in isolation type differs depending on the level of collection years.
# In the future analysis, we may consider to separate the analysis based on levels of collection year to minimize the difference among groups that consists of isolation type variable.


