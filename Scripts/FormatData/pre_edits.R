########################################################################
# This script is a pre-edit to ready data for subsequent formatting    #
# and analyses from Hampton et al. 2016.                               #
# This script does the following:                                      #
#    1) Remove ice-free winter data. This is included in the public    #
#        data so other researchers may use it, but was not included    #
#        in our analyses.                                              #
#    2) Removes extremely high DOC value from one lake.                #
########################################################################

#################### Code by Stephanie Labou ###########################

#set working directory - setwd() - to "core" folder (e.g. the name of this entire repository)

#### clear workspace ####
rm(list = ls()) 
graphics.off()


#### load packages ####
library(dplyr)

#### read in data ####

#Download the publicly available dataset from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V
#Save the dataset to the "Data" folder of this repository.
#The file should be named "ecology_under_lake_ice_public.csv" - if not, update the file name below

data.raw <- read.csv("Data/ecology_under_lake_ice_public.csv", stringsAsFactors = FALSE)

data.fix <- data.raw

#### remove ice-free winter data ####

# Remove Lake Erie 2012 iceon data. 
# Lake Erie was >90% uncovered by ice this winter (per communication with research who submitted data).
# The data was included in the main dataset which is made public,
# but is removed here from these synthesis analyses, to focus solely on ice-covered lakes during iceon/winter.

data.fix <- filter(data.fix, !(lakename == "Lake Erie" & researcher == "Michael.Twiss" & 
                                             year == 2012 & season == "iceon"))

#### remove suspect high DOC value ####
# Remove Lake Peipsi extremely high ice DOC value (suspect), replace with NA
# Note: because value is unpaired (in fact, only DOC from Lake Peipsi) and not part of a longer time series, only impact of removal is on overall DOC mean across lakes

data.fix <- data.fix %>% 
  mutate(avetotdoc = ifelse(lakename== "Lake Peipsi" & !is.na(avetotdoc), NA, avetotdoc)) 

#### write to csv ####
# Note: this is the data that will be read into create.data.long.R.

write.csv(data.fix, "Data/public_data_prepped.csv", row.names = FALSE)

