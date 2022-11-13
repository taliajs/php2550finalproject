
# Data Pre-Processing
  # recoding variables (NA's and Serovar variable)
  # creating host categories
  # standardizing location
  # creating regions
  # processing year and month 
  # seasons
  # categories/groups for isolate source (TBD)

# Libraries
library(tidyverse)
library(dplyr)
library(readr)
library(lubridate)

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
  isolates2$state <- vector(length = length(isolates2$Location)) # this will be the standardized location, aka 
  
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

# if want to select the original location columns, lat long, and the new columns from function: 
  #isolates2 <- isolates2 %>% select(7,32,51,52)


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


####################
## YEAR AND MONTH ##
####################

# Extracting the year and month from `collection.date` variable

# replacing the NA
isolates2$Collection.date[is.na(isolates2$Collection.date)] <- ""

# process the year and month
for (i in 1:length(isolates2$Collection.date)){
  # some string has the date in the version of "%mm/%dd/%YYYY", filtering those strings out
  if(str_detect(isolates2$Collection.date[i],"/")==TRUE){
    # character length of the string for "%m/%d/%YYYY"
    if(nchar(isolates2$Collection.date[i])==8){
      # extract the year at the position 5 to 8
      isolates2$Collection.year[i] <- str_sub(isolates2$Collection.date[i],5,8)
      # extract the month at the position 1
      isolates2$Collection.month[i] <- str_sub(isolates2$Collection.date[i],1,1)
    }
    # character length of the string for "%mm/%d/%YYYY" or "%m/%dd/%YYYY"
    else if(nchar(isolates2$Collection.date[i])==9){
      # extract the year at the position 6 to 9
      isolates2$Collection.year[i] <- str_sub(isolates2$Collection.date[i],6,9)
      
      isolates2$Collection.month[i] <- str_sub(isolates2$Collection.date[i],1,2)
      # when the string is in the format of "%m/%dd/%YYYY", it will include "/" in the month, so we need further extraction
      if(str_detect(isolates2$Collection.month[i],"/")==TRUE){
        # removing the "/" in the month from the string "%m/%dd/%YYYY"
        isolates2$Collection.month[i] <- str_extract(isolates2$Collection.month[i],"(\\d)+")
      }
    }
    # character length of the string for "%mm/%dd/%YYYY"
    else if(nchar(isolates2$Collection.date[i])==10){
      # extract the year at the position 7 to 10
      isolates2$Collection.year[i] <- str_sub(isolates2$Collection.date[i],7,10)
      # extract the year at the position 1 to 2
      isolates2$Collection.month[i] <- str_sub(isolates2$Collection.date[i],1,2)
    }
  }
  # other dates are either missing or in the format as "%Y-%m-%d"
  else if(str_detect(isolates2$Collection.date[i],"/")==FALSE){
    # extract the year at position 1 to 4
    isolates2$Collection.year[i] <- str_sub(isolates2$Collection.date[i],1,4)
    # extract the month at position 6 to 7
    isolates2$Collection.month[i] <- str_sub(isolates2$Collection.date[i],6,7)
  }
}
# convert the class of year and month to numeric from characters
isolates2$Collection.month <- as.numeric(isolates2$Collection.month)
isolates2$Collection.year <- as.numeric(isolates2$Collection.year)


#############
## SEASONS ##
#############

# grouping months into seasons 


## checking the processing for year, month and season

  #isolates3 <- isolates2 %>% select(31,55,56)
  


####################
## ISOLATE SOURCE ##
####################

# Creating groups 




# export isolates2 csv
write.csv(isolates2, file = "isolates2.csv")


# SOURCES
# 1. sapply, string split -- https://stackoverflow.com/questions/33683862/first-entry-from-string-split
# 2. https://stat.ethz.ch/R-manual/R-devel/library/base/html/trimws.html
# 3. grepl: https://stackoverflow.com/questions/64035151/using-grepl-within-an-if-else-statement-within-a-for-loop-in-r
#4. Which states are in each Bifsco Regions: https://journals.asm.org/doi/10.1128/aem.02833-10 (under "materials and methods" section")
#5. Midwest States (according to the Census Bureau): https://en.wikipedia.org/wiki/Midwestern_United_States
