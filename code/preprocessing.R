
# Data Pre-Processing
  # recoding variables (NA's and Serovar variable)
  # creating host categories
  # standardizing location
  # creating regions
  # processing year and month and time
  # seasons
  # genetic info
  # isolate source categories/groups (TBD/in progress)

# Libraries
library(tidyverse)
library(dplyr)
library(readr)
library(stringr)

# read in isolates dataset
isolates <- read.csv("isolates.csv")
# rename first column 
colnames(isolates)[colnames(isolates) == "X.Organism.group"] <- "Organism.group"

## Recoding Variables ##

# Recode blank spaces as NA's
isolates2 <- isolates %>% 
  mutate_all(~replace(., . == "", NA))

# Recode Serovar
isolates2$Serovar[isolates2$Serovar=="Enteriditis"] <- "Enteritidis"
isolates2$Serovar[isolates2$Serovar=="enteritidis"] <- "Enteritidis"
isolates2$Serovar[isolates2$Serovar=="Enteritidis (predicted)"] <- "Enteritidis"
isolates2$Serovar[isolates2$Serovar=="Enteritidis (Predicted)"] <- "Enteritidis"
isolates2$Serovar[isolates2$Serovar=="Enteritidus"] <- "Enteritidis"
isolates2$Serovar[isolates2$Serovar=="Salmonella enterica subsp. enterica serovar Enteritidis"] <- "Enteritidis"
isolates2$Serovar[isolates2$Serovar=="potential monophasic variant of Typhimurium"] <- "Typhimurium (Monophasic)"
isolates2$Serovar[isolates2$Serovar=="Salmonella enterica subsp. enterica serovar Typhimurium"] <- "Typhimurium"
isolates2$Serovar[isolates2$Serovar=="Typhimurium (Predicted)"] <- "Typhimurium"
isolates2$Serovar[isolates2$Serovar=="Typhimurium Copenhagen"] <- "Typhimurium (Copenhagen)"
isolates2$Serovar[isolates2$Serovar=="Typhimurium var. 5-"] <- "Typhimurium (Copenhagen)"
isolates2$Serovar[isolates2$Serovar=="Typhimurium var. Copenhagen"] <- "Typhimurium (Copenhagen)"
isolates2$Serovar[isolates2$Serovar=="Typhimurium var. O 5 - (Copenhagen)"] <- "Typhimurium (Copenhagen)"
isolates2$Serovar[isolates2$Serovar=="Typhimurium var. O 5-(Copenhagen)"] <- "Typhimurium (Copenhagen)"
isolates2$Serovar[isolates2$Serovar=="Typhimurium var. O:5-"] <- "Typhimurium (Copenhagen)"
isolates2$Serovar[isolates2$Serovar=="Typhimurium* (Cerro)"] <- "Typhimurium (Cerro)"


## Creating categories for host ##

# Create host categories
Human <- "Homo sapiens"
Environment <- c("Dust","Environment","Environmental","soil","Soil")
Bird <- c("Bird","chicken","Chicken","Columba livia","Eudocimus albus","Gallus gallus","Gallus gallus domesticus",
          "Meleagris gallopavo","Parrot","Pelecanus occidentalis","Poultry","Spinus pinus","turkey","Turkey")
Reptile <- c("Lepidochyelys olivacea","Malaclemys terrapin","Pogona vitticeps","SNAKE","Terrepene carolina")
Animal <- "animal"
Food <- "Raw almond"
Mammal <- c("Alces alces","Alpaca","Bos taurus","bovine","Bovine","Canis sp.","Cat","cattle","cow","Deer","Dog",
            "Enhydra","Equine","Equus caballus","Equus ferus caballus","Horse","Lama glama","lamb","mouse",
            "Mus musculus","Neogale vison","opposum","Ovis aries","Pig","porcine","Rabbit","rodent","Sus scrofa",
            "Sus sp.","swine")

isolates2$Host.category <- vector(length = length(isolates2$Host))

for (i in 1:length(isolates2$Host)) {
  if (isolates2$Host[i]%in%Human) {
    isolates2$Host.category[i] <- "Human"
    
  } else if (isolates2$Host[i]%in%Environment) {
    isolates2$Host.category[i] <- "Environment"
    
  } else if (isolates2$Host[i]%in%Bird) {
    isolates2$Host.category[i] <- "Bird"
    
  } else if (isolates2$Host[i]%in%Reptile) {
    isolates2$Host.category[i] <- "Reptile"
    
  } else if (isolates2$Host[i]%in%Animal) {
    isolates2$Host.category[i] <- "Animal"
    
  } else if (isolates2$Host[i]%in%Food) {
    isolates2$Host.category[i] <- "Food"
    
  } else if (isolates2$Host[i]%in%Mammal) {
    isolates2$Host.category[i] <- "Mammal"
    
  } else {
    isolates2$Host.category[i] <- NA
  }
}




############################
## STANDARDIZING LOCATION ##
############################

# GOAL: Create a function that goes through the location column and standardizes! (create a new column for state)

standardize_location <- function(df) {
  # standardizes location variable
  #'@param data frame
  #'@return standardized location (state abbreviation)
  
  # create standardized location variable (state abbreviation)
  isolates2$state <- vector(length = length(isolates2$Location)) # this will be the standardized location 
  
  # create variable that represents the initial extracted location values 
  isolates2$location2 <- vector(length = length(isolates2$Location)) 
  
  for (i in 1:nrow(isolates2)) {
    
    # getting the location values
    location <- isolates2$Location[i]
    
    # extract everything after the colon (USA:[.....])
    location_names <- sapply(strsplit(location, ":"), `[`,2) #source 1
    
    # add the split location names to dataframe
    isolates2$location2[i] <- location_names
    
    # trim white space in front: source #2
    isolates2$location2[i] <- trimws(location_names)
    
    # if only USA (aka NA or no state) --> not specified
    if (is.na(isolates2$location2[i])){
      isolates2$state[i] <- "Not Specified"
    }
    
    # if state is in abbreviation form, keep it like that
    else if (nchar(isolates2$location2[i], keepNA = FALSE) == 2) {
      isolates2$state[i] <- isolates2$location2[i]
    }
    
    # if state name, get abbreviated form
    else if (isolates2$location2[i] %in% state.name) {
      # match state name to state abbreviation
      isolates2$state[i] <- state.abb[match(isolates2$location2[i], state.name)] 
    }
    
    # splitting on dash 
    else if (grepl("-", isolates2$location2[i])) {
      split_dash <- sapply(strsplit(isolates2$location2[i], "-"), `[`,1)
      isolates2$state[i] <- split_dash
    }
    
    # splitting on commas
    else if (grepl(",", isolates2$location2[i])) {
      split_comma <- strsplit(isolates2$location2[i], ",")

      # getting first argument from split
      split_comma_1 <- sapply(strsplit(isolates2$location2[i], ","), `[`,1)
      
      # if extracted value is in state name list, get state abbreviation
      states_comma_1 <- ifelse(split_comma_1 %in% state.name, isolates2$state[i] <- state.abb[match(split_comma_1, state.name)], "no")

      # getting second argument from list 
      split_comma_2 <- sapply(strsplit(isolates2$location2[i], ","), `[`,2)
      split_comma_2 <- trimws(split_comma_2) # get rid of white space
      
      # if extracted value is in state name list, get abbreviation 
      states_comma_2 <- ifelse(split_comma_2 %in% state.name, isolates2$state[i] <- state.abb[match(split_comma_2, state.name)], "no")
      #print(states_comma_2)
      
      #if argument is state abbreviation, add it
      states_comma_2 <- ifelse(split_comma_2 %in% state.abb, isolates2$state[i] <- state.abb[match(split_comma_2, state.abb)], "no")
      
      # dealing with the North Carolina University --> it's only a single entry (row 9883)
      split_space1 <- sapply(strsplit(split_comma_1, " "), `[`,1)
      split_space2 <- sapply(strsplit(split_comma_1, " "), `[`,2)
      loc <- paste(split_space1, split_space2)
      state_loc <- ifelse(loc %in% state.name, isolates2$state[i] <- state.abb[match(loc, state.name)], "no")
    }
    
    #BIFSCo regions, or other regions/areas (e.g. Midwest, Western Region) --> assign to "OTHER"
    else {
      isolates2$state[i] <- "Other"
    }
  }
  
  return(isolates2)
}
standardize_location(isolates2)
isolates2 <- standardize_location(isolates2)


######################
## CREATING REGIONS ##
######################

# Create regions based on the BIFSco Regions (Source #4)
  # 48 states in total
  # Note: missing Alaska (AK) and Hawaii (HI)
    # ** WHAT REGIONS WILL THESE STATES BE ADDED TO?? --> Hawaii will be added to the West (Region 2), and Alaska will be added to Northwest (Region 1). 
  # referred to the CDC NCCDPHP regions, which are based loosely off the Dept. of Health and Human Services regions (source: https://www.cdc.gov/coordinatedchronic/docs/nccdphp-regions-map.pdf)

# Region 1: Northwest (3 states)
region1 <- c("WA", "OR", "ID") #AK
# Region 2: West (2 states)
region2 <- c("CA", "NV") #HI
# Region 3: Southwest (3 states)
region3 <- c("AZ", "NM", "TX")
# Region 4 (4 states)
region4 <- c("MT", "CO", "WY", "UT")
# Region 5: Upper Midwest (5 states)
region5 <- c("NE", "ND", "SD", "MN", "WI")
# Region 6: Central (3 states)
region6 <- c("IA", "KS", "MO")
# Region 7: Southeast (10 states)
region7 <- c("OK", "AR", "LA", "NC", "SC", "FL", "AL", "MS", "GA", "TN")
# Region 8: Northeast (18 states)
region8 <- c("IL", "IN", "KY", "MA", "ME", "MD", "MI", "NJ", "NY", "NH", "CT", "RI", "OH", "WV", "VA", "VT", "PA", "DE")


regions <- function(df) {
  # creates regions for the states
  #'@param data frame with the states variable
  #'@return regions
  
  df <- df
  
  # make a new variable for region (and add it to dataframe with states)
  df$region <- vector(length = length(df$state))
  
  for (i in 1:nrow(df)) {
    if (df$state[i] %in% region1) {
      df$region[i] <- 1
    } else if (df$state[i] %in% region2) {
      df$region[i] <- 2 
    } else if (df$state[i] %in% region3) {
      df$region[i] <- 3
    } else if (df$state[i] %in% region4) {
      df$region[i] <- 4 
    } else if (df$state[i] %in% region5) {
      df$region[i] <- 5
    } else if (df$state[i] %in% region6) {
      df$region[i] <- 6 
    } else if (df$state[i] %in% region7) {
      df$region[i] <- 7
    } else if (df$state[i] %in% region8) {
      df$region[i] <- 8
    } 
    
    # some location had the BIFsco Regions mentioned --> so get those!
    else if (grepl(" ", df$location2[i])) {
      split_bifsco <- sapply(strsplit(df$location2[i], " "), `[`, 3) 
      #print(split_bifsco) 
      df$region[i] <- split_bifsco
      
      # for the longer BIFsco Regions 
      split_bifsco <- sapply(strsplit(df$location2[i], " "), `[`, 5) 
      df$region[i] <- split_bifsco
    }
    
    # other regions 
    else if (grepl("Western Region", df$location2[i])) {
      # assume Western Region is region 2
      df$region[i] <- 2
    }
    
    else if (grepl("Midwest", df$location2[i])) {
      # Midwest states (source 5)
      #region 5: 5/5
      #region 6: 3/3
      #region 8: 4/18
      # most of the Midwest states fall in between BIFSCo regions 5 and 6 (all 8 states listed in regions 5 and 6 are included as the "Midwest). So randomly choose either region 5 or 6
      df$region[i] <- sample(5:6, 1)
    }
    
    # if state not specified --> NA 
    else {
      df$region[i] <- NA
    }
  }
  return(df)
}
df <- isolates2 
isolates2 <- regions(isolates2)


#################################
## DATE, TIME, YEAR, and MONTH ##
#################################

# Extracting the date and time from `Create.date` variable, and creating `Date` and `Time` variables
isolates2$Date <- str_sub(isolates2$Create.date, 1, -str_locate(isolates2$Create.date,"T")[1])

isolates2$Time <- str_sub(isolates2$Create.date,
                          (str_locate(isolates2$Create.date,"T")[1]+1),
                          (str_locate(isolates2$Create.date,"Z")[1]-1))

# Create `Year`, `Month`, and `Date` variables using lubridate functions
isolates2$Year <- year(isolates2$Date)

isolates2$Month <- month(isolates2$Date)

isolates2$Week <- week(isolates2$Date)

# Restrict data to 2017 to present
isolates2 <- isolates2 %>% filter(Year == 2017 | Year == 2018 | Year == 2019 |
                                    Year == 2020 | Year == 2021 | Year == 2022)


#############
## SEASONS ##
#############

# grouping months into seasons 

#for (i in 1:nrow(isolates2)){
#  # if the month is Mar, Apr, May --> spring
#  if (isolates2$Collection.month[i] %in% c(3,4,5)){
#    # then new variable season has the value 1
#    isolates2$season[i] <- 1
#  }# if the month is Jun, Jul, Aug --> summer
#  else if (isolates2$Collection.month[i] %in% c(6,7,8)){
#    # then new variable season has the value 2
#    isolates2$season[i] <- 2
#  }# if the month is Sep,Oct,Nov --> fall
#  else if(isolates2$Collection.month[i] %in% c(9,10,11)){
#    # then new variable season has the value 3
#    isolates2$season[i] <- 3
#  }# if the month is Dec, Jan, Feb --> winter
#  else if (isolates2$Collection.month[i] %in% c(12,1,2)){
#    # then new variable season has the value 4
#    isolates2$season[i] <- 4
#  } else { #if month not specified, NA 
#    isolates2$season[i] <- NA
#  }
#}



####################
## GENETIC INFO ##
####################

# Creating outbreak variables

# finding associations between Min.same and Min.diff on the provided Outbreak info
outbreak_df <- isolates %>% filter(Outbreak!="") %>% 
  select(c(Min.same,Min.diff,Outbreak))
outbreak_df$Outbreak <- as.factor(outbreak_df$Outbreak)

outbreak_logit <- glm(Outbreak ~sqrt(Min.diff), 
                      data = outbreak_df, family = "binomial")

exp_outbreak <- predict(outbreak_logit, newdata=isolates2, type="response")
exp_outbreak[is.na(exp_outbreak)] <- 0

# threshold of predicted outbreak 0.9813, representing 75% of data with Outbreak variable pressent
# outbreak <- cbind(outbreak,pred)
# summary(outbreak$pred)

new_outbreak <- ifelse(exp_outbreak >= 0.9813, 1, 0)
isolates2 <- cbind(isolates2,new_outbreak)

####################
## ISOLATE SOURCE ##
####################

#isolates2 <- read.csv("isolates2.csv")


# Creating groups 
#main category (e.g meat, poultry)
isolates2$isolate_source <- vector(length = length(isolates2$Isolation.source)) 
# sub categories (eg. for poultry --> chicken and turkey, for meat: beef, pork)
isolates2$isolate_source_type <- vector(length = length(isolates2$Isolation.source)) 

for (i in 1:nrow(isolates2)) {
  
  # -- POULTRY --
  if (grepl("chicken|turkey|Chicken|CHICKEN|Chicekn|Cicken|Turkey|poultry|Poultry", isolates2$Isolation.source[i])) {
    #print(isolates2$Isolation.source[i])
    isolates2$isolate_source[i] <- "poultry"
  }
  # poultry categories
  if (grepl("chicken|Chicken|CHICKEN|Chicekn|Cicken", isolates2$Isolation.source[i])) {
    #testing <- c(testing, isolates2$Isolation.source[i])
    #print(isolates2$Isolation.source[i])
    isolates2$isolate_source_type[i] <- "chicken"
  }
  if (grepl("turkey|Turkey", isolates2$Isolation.source[i])) {
    #testing <- c(testing, isolates2$Isolation.source[i])
    #print(isolates2$Isolation.source[i])
    isolates2$isolate_source_type[i] <- "turkey"
  } else if (grepl("poultry|Poultry", isolates2$Isolation.source[i])) {
    isolates2$isolate_source_type[i] <- "poultry"
  } 
  
  #unknown poultry
  # else if (grepl("Thighs|Breast|Wings", isolates2$Isolation.source)) {
  #   isolates2$isolate_source_type[i] <- "poultry"
  #   isolates2$isolate_source_type[i] <- "unknown poultry"
  # }
  
  # -- MEAT --
  if (grepl("cattle|Cattle|beef|Beef|Trimmings|swine|Swine|pork|Pork|goat|Goat|sheep|Sheep|lamb|Lamb", isolates2$Isolation.source[i])) {
    #testing <- c(testing, isolates2$Isolation.source[i])
    #print(isolates2$Isolation.source[i])
    isolates2$isolate_source[i] <- "meat"
  } 
  # meat categories
  if (grepl("cattle|Cattle|beef|Beef|Trimmings", isolates2$Isolation.source[i])) {
    isolates2$isolate_source_type[i] <- "beef"
  } 
  if (grepl("swine|Swine|pork|Pork", isolates2$Isolation.source[i])){
    isolates2$isolate_source_type[i] <- "pork" 
  } 
  
  if (grepl("goat|Goat|sheep|Sheep|lamb|Lamb", isolates2$Isolation.source[i])) {
    isolates2$isolate_source_type[i] <- "other meat"
  } 
  
  # -- EGGS --
  if (grepl("egg", isolates2$Isolation.source[i]) | 
      grepl("Egg", isolates2$Isolation.source[i]) |
      grepl("Egg yolks", isolates2$Isolation.source[i]) |
      grepl("Egg", isolates2$Isolation.source[i]) & grepl("yolks", isolates2$Isolation.source[i]) |
      grepl("Egg", isolates2$Isolation.source[i]) & grepl("whites", isolates2$Isolation.source[i])
  ) {
    isolates2$isolate_source[i] <- "eggs"
    isolates2$isolate_source_type[i] <- "eggs"
  }
  
  # -- FECES (feces + stool) -- 
  if (grepl("stool|Stool|feces|Feces|fecal", isolates2$Isolation.source[i])) {
    isolates2$isolate_source[i] <- "feces"
    isolates2$isolate_source_type[i] <- "feces"
  }
  
  # # -- CLINICAL --
  # else if (grepl("swab|Swab", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "clinical"
  #   isolates2$isolate_source_type[i] <- "swab"
  # } 
  # # organs
  # else if (grepl("Abdomen|aorta|Aorta|lung|Lung|colon|Colon|intestine|Intestine|intestinal|liver|Liver|lymph node|Lymph Node|brain|kidney|Kidney|jejunum|pericardium|Pericardium|brain|breast|heart|Heart|spleen|joints|trachea|Trachea|organ|small|Small", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "clinical"
  #   isolates2$isolate_source_type[i] <- "organs"
  # } 
  # #anatomy
  # else if (grepl("joint|Joint|hip|LEGS|Bone|foot|yok sac|Huck|yolks|Yolk sac|Rectal|nasal|tonsil|pluera|Groin|Hip|
  #           SCROTUM|bone", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "clinical"
  #   isolates2$isolate_source_type[i] <- "anatomical"
  # } 
  # # other clinical 
  # else if (grepl("Vaginal|wound|Abscess|abscess|necropsy|carcass|lab|sputum|perionitis|aspirate|Cyst|CSF|autopsy|
  #           auxotroph|ISOLATE", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "clinical"
  #   isolates2$isolate_source_type[i] <- "other clinical"
  # } 
  # 
  # -- BIOLOGICAL (biological fluids, tissues) -- 
  # if (grepl("biological fluid|fluid|Fluid|tissue|TISSUE", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "biological"
  #   isolates2$isolate_source_type[i] <- "biological"
  # }
  # 
  # # -- BODILY FLUIDS -- 
  # if (grepl("blood|Blood|body fluid|urine|Urine|URINE", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "bodily fluids"
  #   isolates2$isolate_source_type[i] <- "bodily fluids"
  # }
  
  # -- ENVIRONMENT -- 
  # else if (grepl("environment|water|lagoon|lake|river|River|Water|sewage|soil|stream|Creek|Foam|Enteric pool", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "environmental"
  #   isolates2$isolate_source_type[i] <- "environmental"
  # }
  
  # -- FOOD --
  # if (grepl("spinach", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "food"
  #   isolates2$isolate_source_type[i] <- "leafy greens"
  # } else if (grepl("Animal feed|Pet food|pet food|Cat food|treat|Kibble|kibble", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "food"
  #   isolates2$isolate_source_type[i] <- "animal food"
  # } else if (grepl("apple|cantaloupe|pudding|olive|Food|Hazelnut|hazelnut|Meal|Peanut butter|Pecans|pecans|almonds|Almonds|Tomato|tomato|salad|celery|wheat|meal|Celery|corn|Ready to eat|hamburger|peanut butter", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "food"
  #   isolates2$isolate_source_type[i] <- "food"
  # }
  
  # # -- ANIMAL --
  # if (grepl("human|ratite|sea lion|Reptile|Siluriformes|feline|equine|Hock|hock|turtle||goose|dog|snake|peafowl|Wild|mammal|Deer|parrot|mouse|horse|caprine|opossum|Canis|rabbit", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "animal"
  #   isolates2$isolate_source_type[i] <- "animal"
  # }
  # 
  # if (grepl("bovine|Bos taurus|Veal", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "animal"
  #   isolates2$isolate_source_type[i] <- "cattle"
  # }
  # 
  # if (grepl("porcine", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "animal"
  #   isolates2$isolate_source_type[i] <- "swine"
  # }
  # 
  # if (grepl("Gallus|gallus|avian clocae|Columba", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "animal"
  #   isolates2$isolate_source_type[i] <- "bird"
  # }
  # 
  
  # -- NOT AVAILABLE -- 
  if (grepl("Not collected|not known|Unknown", isolates2$Isolation.source[i])) {
    isolates2$isolate_source[i] <- "not available"
    isolates2$isolate_source_type[i] <- "not available"
  } 
  
  # -- NA -- 
  if (is.na(isolates2$Isolation.source[i])) {
    isolates2$isolate_source[i] <- "NA"
    isolates2$isolate_source_type[i] <- "NA"
  }
  
  # -- OTHER
  # else {
  #   isolates2$isolate_source[i] <- "other"
  #   isolates2$isolate_source_type[i] <- "other"
  # }
  # 
  
  #-- OTHER
  # if (grepl("treat|soil|node|feed|lagoon|food", isolates2$Isolation.source[i])) {
  #   isolates2$isolate_source[i] <- "other"
  #   isolates2$isolate_source_type[i] <- "other"
  # }
  
  # -- OTHER -- 
  if (isolates2$isolate_source[i] == "FALSE") {
    isolates2$isolate_source[i] <- "other"
    isolates2$isolate_source_type[i] <- "other"
  }
}

# export isolates2 csv
write.csv(isolates2, file = "isolates2.csv")



# SOURCES
# 1. sapply, string split -- https://stackoverflow.com/questions/33683862/first-entry-from-string-split
# 2. https://stat.ethz.ch/R-manual/R-devel/library/base/html/trimws.html
# 3. grepl: https://stackoverflow.com/questions/64035151/using-grepl-within-an-if-else-statement-within-a-for-loop-in-r
#4. Which states are in each Bifsco Regions: https://journals.asm.org/doi/10.1128/aem.02833-10 (under "materials and methods" section")
#5. Midwest States (according to the Census Bureau): https://en.wikipedia.org/wiki/Midwestern_United_States
