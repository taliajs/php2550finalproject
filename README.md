# PHP 2550 Final Project

### Sirui Cui, Tova Ibbotson, Talia Seshaiah


## About 
This repository contains the code associated with our final project, which is focused on Salmonella vaccines for poultry. Our goal for this project is to identify strains/flocks/regions, etc. that should be prioritized for vaccination campaigns (e.g specific strains, certain flocks and/or certain regions). For instance, looking at where, when and how are outbreaks happening generally, and looking at Salmonella outbreaks that cause illness. 


## Files

- `data`: This folder contains the data used. The original data (from the NCBI Isolates Browser) can be found in the folder titled "original_data". The data that was used for our analysis is `isolates2.csv` (this is our data after pre-processing). 

- `code`: This folder contains the code from the project. 
  - `missingdata.R`: Missing data analysis. This file contains the code used for our initial exploratory analysis of missing data.
  - `EDA.Rmd`: This file contains the code for the exploratory analysis.
  - `preprocessing.R`: Contains the code used in our data cleaning/preprocessing
  - `Model.R`: Initial model building for a prediction model that predicts the probability of a specific serovar for Salmonella. 
  - `geographic_methods.Rmd`: Exploring methods for geographic visualizations.
