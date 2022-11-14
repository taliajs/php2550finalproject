# PHP 2550 Final Project

### Sirui Cui, Tova Ibbotson, Talia Seshaiah


## About 
This repository contains the code associated with our final project, which is focused on identifying "targets" for Salmonella vaccines. Our goal for this project is to identify strains/regions/time of year, etc. that should be prioritized for vaccination campaigns. For instance, Are there certain/specific regions and time of year that should be focused on when thinking about vaccination? Are there certain strains of Salmonella that are more prevalant than others?

## Files

- `data`: This folder contains the data used. The original data (from the NCBI Isolates Browser) can be found in the folder titled "original_data". The data that was used for our analysis is `isolates2.csv` (this is our data after pre-processing). 

- `code`: This folder contains the code from the project. 
  - `missingdata.R`: Missing data analysis. This file contains the code used for our initial exploratory analysis of missing data.
  - `EDA.Rmd`: This file contains the code for the exploratory analysis.
  - `preprocessing.R`: Contains the code used in our data cleaning/preprocessing
  - `Model.R`: Initial model building for a prediction model that predicts the probability of a specific serovar for Salmonella. 
  - `geographic_methods.Rmd`: Exploring methods for geographic visualizations.
