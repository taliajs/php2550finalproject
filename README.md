# Identifying Targets for Poultry Vaccination by Predicting Salmonella Outbreaks

### Sirui Cui, Tova Ibbotson, Talia Seshaiah


## About 
This project was done as part of the final project for PHP 2550: Practical Data Analysis class. Our project is focused on identifying targets for Salmonella vaccines for poultry.

Our goal for this project is to identify strains/regions/time of year, etc. that should be prioritized for vaccination campaigns. For instance, Are there certain/specific regions and time of year that should be focused on when thinking about vaccination? Are there certain strains of Salmonella that are more prevalant than others? We fit a prediction model to predict outbreaks of Salmonella enterica in the United States.

## Setup 

R was used for this project.

R packages required: 

(Preprocessing and EDA)
- tidyverse
- dplyr
- readr
- stringr 
- lubridate
- ggplot2
- table1
- kableExtra

(Specific R packages related to model/analysis)
- rpart.plot
- caret
- pROC
- DescTools
- randomForest
- rsample
- rfUtilities 


## Files

- `code`: This folder contains the code from the project. 
  - `EDA.Rmd`: This file contains the code for the exploratory analysis.
  - `preprocessing.R`: Contains the code used in our data cleaning/preprocessing
  - `Linear_Model.R`: Linear model.
  - `RF_Model.R`: Random forest model. 
  - `geographic_methods.Rmd` and `geographic_methods.pdf`: Exploring methods for geographic visualizations.
  - `variable_eda.R`: Exploring variables in model (frequency tables). 

- `data`: This folder contains the data used for this project:
  - `isolates.csv`: This is the original data (from the NCBI Isolates Browser) 
  - `isolates2.csv`: This is the data that was used for our analysis (the data after pre-processing). 

- `lit-review`: This folder contains the literature review (more background information on salmonella poultry vaccines, and methods for identifying/predicting outbreaks). 

- `supplementary-material`: This folder contains the supplementary material for project.
