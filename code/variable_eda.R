# exploring variables in model
library(tidyverse)

# read in data
isolates2 <- read.csv("data/isolates2.csv")

# ISOLATE SOURCES

#unique isolate sources and types
unique(isolates2$isolate_source) #6
unique(isolates2$isolate_source_type) #10

# frequency table for sources
isolate_source_freq <- isolates2 %>% 
  group_by(isolate_source) %>%
  summarise(freq = n()) %>%
  mutate(sum = sum(freq)) %>%
  mutate(percent = (freq/sum)*100)

isolate_source_freq <- isolate_source_freq[order(isolate_source_freq$freq, decreasing = TRUE),]

# STATE
state_freq <- isolates2 %>% 
  group_by(state) %>% 
  summarise(freq = n()) %>%
  mutate(sum = sum(freq)) %>%
  mutate(percent = (freq/sum)*100)

# REGIONS
region_freq <- isolates2 %>%
  group_by(region) %>%
  summarise(freq = n()) %>%
  mutate(sum = sum(freq)) %>%
  mutate(percent = (freq/sum)*100)

regions <- isolates2 %>% select(X, state, location2, region) %>%
  distinct(state, .keep_all = TRUE)

#############################
### EDA Figures for Paper ###
#############################

# Create table 1
table1(~ Serovar + Host.category + isolate_source + isolate_source_type + region + 
         Min.same + Min.diff + AMR.genotypes.count +
         SNP.cluster.rank + Stress.genotypes.count, data = isolates2)

# Summarize serovars over time
isolates2$region <- factor(isolates2$region)
isolates2$Year <- factor(isolates2$Year)

label(isolates2$region) <- "Region"

table1(~ region + Year | Serovar, data = isolates2)







