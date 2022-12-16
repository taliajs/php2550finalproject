# Identifying Targets for Poultry Vaccination by Predicting Salmonella Outbreaks

### Sirui Cui, Tova Ibbotson, Talia Seshaiah


## About 
This project was done as part of the final project for PHP 2550: Practical Data Analysis class. Our project is focused on identifying targets for Salmonella vaccines for poultry.

Our goal for this project is to predict outbreaks of Salmonella enterica in the United States, and recommend when and where poultry vaccination campaigns should be introduced to reduce outbreaks. 

## Setup 

R was used for this project.

R packages required/needed:

| R package | Use |
| ------------- | ------------- |
| tidyverse  | Overall  |
| table1 | Overall (displaying tables) |
| kableExtra | Overall (displaying tables) | 
| dpylr | Data pre-processing, and EDA |
| readr | Data pre-processing |
| stringr | Data pre-processing |
| lubridate | Data pre-processing, and EDA | 
| ggplot2 | Visualizations |
| rpart.plot | Linear Model |
| caret | Linear Model |
| randomForest | Random Forest Model |
| rsample | Random Forest Model |
| rfUtilities | Random Forest Model |
| pROC | Analyzing models |
| DescTools | Analyzing models |


## Files

- `code`: This folder contains the code from the project. 
  -  `EDA.Rmd`: This file contains the code for the exploratory analysis.
  -  `geographic_methods.Rmd` and `geographic_methods.pdf`: Exploring methods for geographic visualizations.
  - `preprocessing.R`: Contains the code used in our data cleaning/preprocessing.
  - `variable_eda.R`: Exploring variables in model (frequency tables). 
  - `Linear_Model.R`: Linear model (code for creating and analyzing model)
  - `RF_Model.R`: Random forest model (code for creating and analyzing model)

- `data`: This folder contains the data used for this project:
  - `isolates.csv`: This is the original data (from the NCBI Isolates Browser) 
  - `isolates2.csv`: This is the data that was used for our analysis (the data after pre-processing). 

- `lit-review`: This folder contains the literature review (more background information on salmonella poultry vaccines, and methods for identifying/predicting outbreaks). 

- `supplementary-material`: This folder contains the supplementary material for project.
