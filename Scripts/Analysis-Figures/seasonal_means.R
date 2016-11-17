###################################################################################
# This script calculates seasonal means (e.g. iceon and iceoff means) and         #
# standard deviations and errors for each variable and unique lake/station        #
# combination. Output also includes the total number of (non-NA) observations     #
# as well as the number of iceon (non-NA) observations. This is Table S2 in       #
# Hampton et al. 2016.                                                            #
###################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
################################################################################### 

########################### Code by Stephanie Labou ###############################

#### clear workspace ####
rm(list = ls()) 
graphics.off()

######### load packages ############
library(dplyr)
library(reshape2)
library(plotrix)

###############################################################
################## Load and process data ######################
###############################################################

# The data .csv file is created by the script create.data.long.R
data.orig <- read.csv("Data/data.long.csv", stringsAsFactors = FALSE)

#select columns of interest
data <- data.orig %>% select(lakename, stationname, stationlat, 
                             stationlong, lakecountry, lakeregloc,
                             year, season, variable, stat, value)

#keep only average values
data <- filter(data, stat == "ave")

#keep only variables of interest
vars.keep <- c("chla","propotherphyto","propcyano",
               "propchloro","propcrypto","propdiatom",
               "phytocount","phytomass", "phytoppr", 
               "zoopcount","propdaphnia","propcalanoid",
               "propcyclopoid", "propothercladoc",
               "propotherzoop","proprotifer","zoopmass",
               "propdino", "totphos","totnitro","totdissphos",
               "totdissnitro", "totdoc", "TDN.TDP", "TN.TP", 
               "C.TDN", "C.TDP", "watertemp")

data <- filter(data, variable %in% vars.keep)

## find non-rotifer zoop counts ##
data_cast <- dcast(data, lakename + stationname + stationlat + stationlong +lakecountry + lakeregloc + year + season ~ variable, value.var = "value")


#create new zoop count with excluding rotifers and "other" zooplankton
data_new <- data_cast %>% 
  #calculate zoop counts
  mutate(n.rotifer = proprotifer*zoopcount,
         n.otherzoop = propotherzoop*zoopcount,
         #fix NAs to 0 so subtraction works
         n.rotifer = ifelse(is.na(n.rotifer), 0, n.rotifer),
         n.otherzoop = ifelse(is.na(n.otherzoop), 0, n.otherzoop),
        #new total
        new.zoop.count = zoopcount - n.rotifer - n.otherzoop) %>% 
  select(-n.rotifer, n.otherzoop)

#melt back
data_use <- melt(data_new, id.vars = c("lakename", "stationname", "stationlat", "stationlong",
                                        "lakecountry", "lakeregloc", "year", "season"))

## add number of total observations (non-NA) and number of ice observations (non-NA)
data_use.tot.yrs <- data_use %>% 
  group_by(lakename, stationname, variable, lakecountry, lakeregloc, stationlat, stationlong) %>% 
  filter(!is.na(value)) %>% 
  summarize(n.obvs = length(value))

data_use.ice.yrs <- data_use %>% 
  filter(season == "iceon") %>% 
  group_by(lakename, stationname, variable, lakecountry, lakeregloc, stationlat, stationlong) %>% 
  filter(!is.na(value)) %>% 
  summarize(n.ice = length(value))

#merge
data_use.meta <- merge(data_use.tot.yrs, data_use.ice.yrs, by = c("lakename", "stationname", "variable", 
                                                   "lakecountry", "lakeregloc", "stationlat", "stationlong"))

#calculate averages and standard deviations (one seasonal average and stdev per lake/station/variable/season)
data_use.sm <- data_use %>% select(lakename, stationname, variable, season, value)

data_use.summary <- data_use.sm %>% 
                group_by(lakename, stationname, variable, season) %>% 
                summarize(mean.value = mean(value, na.rm = TRUE),
                          stdev = sd(value, na.rm = TRUE),
                          sterror = std.error(value, na.rm = TRUE)) %>% 
                as.data.frame() %>% 
                #replace NaN values with NA
                mutate(mean.value = ifelse(is.nan(mean.value), NA, mean.value),
                       stdev = ifelse(is.nan(stdev), NA, stdev),
                       sterror = ifelse(is.nan(sterror), NA, sterror))

#reshape
data_use.cast <-reshape(data_use.summary, direction  = "wide", idvar=c("lakename","stationname", "variable"), timevar="season")

data_use.cast <- data_use.cast %>% 
              arrange(variable, lakename, stationname) 

#merge with metadata_use
data_use.write <- merge(data_use.cast, data_use.meta, by = c("lakename", "stationname", "variable"))

data_use.write <- select(data_use.write, lakename, stationname, lakecountry, lakeregloc, stationlat, stationlong, 
                     n.obvs, n.ice, variable, mean.value.iceoff, stdev.iceoff, sterror.iceoff, 
                     mean.value.iceon, stdev.iceon, sterror.iceon)

#prep for write to csv
#no scientific notation
options(scipen=999)
data_use.write <- data_use.write %>% 
              arrange(variable, lakename, stationname)

#round to 3 significant figures
data_use.write$mean.value.iceon <- signif(data_use.write$mean.value.iceon, 3)
data_use.write$stdev.iceon <- signif(data_use.write$stdev.iceon, 3)
data_use.write$sterror.iceon <- signif(data_use.write$sterror.iceon, 3)

data_use.write$mean.value.iceoff <- signif(data_use.write$mean.value.iceoff, 3)
data_use.write$stdev.iceoff <- signif(data_use.write$stdev.iceoff, 3)
data_use.write$sterror.iceoff <- signif(data_use.write$sterror.iceoff, 3)

#write to csv
write.csv(data_use.write, "Outputs/Files/lake_seasonal_means.csv", row.names = FALSE)
