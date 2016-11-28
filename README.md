[![DOI](https://zenodo.org/badge/74072078.svg)](https://zenodo.org/badge/latestdoi/74072078)

#Overview

This repository contains the scripts used to run analyses and generate figures in Hampton et al. 2016. (Hampton, S.E., et al. (2016) Ecology under lake ice. Ecology Letters, doi: 10.1111/ele.12699)

The repository and scripts are set up so that necessary functions are pulled from the Scripts folder, data are pulled from the Data folder, and outputs are placed in the appropriate folders (e.g. Output/Figures, Output/Files).

Each script begins with an overview of (a) what the script does and (b) what portion of Hampton et al. 2016 the script corresponds to. 

People who contributed code to each script are noted at the beginning of the script.

#Using this repository

To use the scripts here, please download the publicly available dataset from the [Knowledge Network for Biocomplexity](https://knb.ecoinformatics.org/#), DOI: 10.5063/F12V2D1V. This is the "base" data that will be used in all subsequent data reformatting, analysis, and figure generation.

**Please carefully read the metadata of the publicly available dataset for information about how the data were collected and data use caveats and restrictions. We are making these scripts available to complement the public data and be transparent about when and how data were filtered and subsetted prior to analyses.** We hope this will help future users of the public data make informed decisions about re-use of the dataset.

Please note that results output using the publicly available dataset will be different than those in Hampton et al. 2016 due to the absence of certain lakes in the publicly available dataset. 
                   
###System requirements and packages needed

Scripts included in this repository were developed on a variety of systems. See REQS.md for more detailed information about R version used. REQS.md also has information about all packages needed to run the scripts in this repository, as well as details about package versions used.

###Set working directory
Once packages are installed, open R/RStudio and use setwd() to set the working directory to the main folder that contains "Data", "Outputs", and "Scripts" folders (e.g. folder that README.md is in).

###Running scripts

**The data formatting scripts must be run first**, otherwise subsequent scripts will not run. All other scripts can be run in any order.

####1) Reformat data 
These three scripts prepare and reformat data into "long" and "long time series" formats.
Please note: the scripts call functions from the Functions folder within the script.

+ `Scripts/FormatData/pre_edits.R` 
+ `Scripts/FormatData/create.data.long.R` 
+ `Scripts/FormatData/focal.long.timeSeries.R`  
<br>

####2) Calculate seasonal means
This script calculates iceon and iceoff means and standard deviations.
It also calculates total number of observations and total number of iceon observations.

+ `Scripts/Analysis-Figures/seasonal_means.R`   
<br>

####3) Regression trees
These scripts construct regression trees for water chemistry variables and 
zooplankton and phytoplankton seasonal community differences (multivariate regression trees).

+ `Scripts/Analysis-Figures/regression_trees_nutrients.R`   
+ `Scripts/Analysis-Figures/regression_trees_zoop_phyto.R`  
<br>


####4) Analysis/Figure: Are ice-on and ice-off seasons different across lakes?
These scripts run analyses (PERMANOVA for proportional data) and create figures 
(scatterplots and boxplots for non-proportional data and stacked bar charts for proportional data). 

+ `Scripts/Analysis-Figures/seasonal_compare_scatters_boxplots.R` 
+ `Scripts/Analysis-Figures/composition_compare_zoop_phyto.R`   
<br>


####5) Analysis: Seasonal patterns within lakes
The first script detrends and deseasons data (by lake and variable).
The subsequent scripts create linear mixed effects models. 
Please note: the LME scripts calls `lme_iteratethis.R` (also in the `Analysis-Figures` folder)

+ `Scripts/Analysis-Figures/timeseries.R` 
+ `Scripts/Analysis-Figures/lme_style.R`  
<br>

####6) Figure: Examples of time series
This script creates a figure of five example time series.
 
+ `Scripts/Analysis-Figures/timeseries_example.R`    
<br>

#### This project is licensed under the MIT license (see LICENSE.MD) Please see CITATION.md for citation instructions.
