################################################################################################
## This script uses a regression tree approach with zooplankton and phytoplankton community   ## 
## composition.                                                                               ##
##                                                                                            ##
## This analysis requires the mvpart package, which can be downloaded from the CRAN archive:  ##
##      https://cran.r-project.org/src/contrib/Archive/mvpart/                                ##
##                                                                                            ##
## Once downloaded, the package can be installed using:                                       ##
##      install.packages(".../Downloads/mvpart_1.6-2.tar.gz", repos = NULL, type = "source")  ##
##                                                                                            ##
## For help downloading/installing mvpart, see the stackoverflow page:                        ##
## http://stackoverflow.com/questions/29656320/r-mvpart-package-any-option-to-use-in-r-3-1-x  ##
################################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
################################################################################### 

#### clear workspace ####
rm(list = ls()) 
graphics.off()

#### load libraries ####
library(dplyr)
library(reshape2)
library(vegan)
library(tidyr)
library(mvpart)

#### read in data ####
dat <- read.csv("Data/data.long.csv", stringsAsFactors = FALSE)

head(dat)

####### PREDICTORS ########

#### physical characteristics of lakes and summer chla (explanatory for reg trees) ####
meta <- select(dat, lakename, stationname, year, season, 
               stationlat, lakearea, lakemaxdepth, lakeelevation,
               variable, stat, value) %>% 
        filter(stat == "ave" & variable == "chla" & season == "iceoff")

#summer chla
meta <- meta %>% 
  group_by(lakename, stationname) %>% 
  #find average summer chla
  #and average across physical 
  #(since some people report different values for ex: max depth by year)
  summarize(lat = mean(stationlat, na.rm = TRUE),
            area = mean(lakearea, na.rm = TRUE),
            maxdepth = mean(lakemaxdepth, na.rm = TRUE),
            elev = mean(lakeelevation, na.rm = TRUE),
            #note: using chla as explanatory will knock out some Finland lakes
            #which have zoop counts but no chla data
            summer.chla = mean(value, na.rm = TRUE)) %>% 
  as.data.frame()

#### change in water temp ######

water <- select(dat, lakename, stationname, year, season, variable, stat, value) %>% 
  filter(variable == "watertemp" & stat == "ave")

water_avg <- water %>% 
  group_by(lakename, stationname, season, variable) %>% 
  summarize(avg_value = mean(value, na.rm = TRUE)) %>% 
  as.data.frame()

water_avg_wide <- dcast(water_avg, lakename + stationname + variable ~ season, value.var = "avg_value")

water_avg_diff <- mutate(water_avg_wide, SW_diff_watemp = iceoff - iceon) %>% 
  select(-variable, -iceoff, -iceon)

#### merge #####
dat.full <- merge(meta, water_avg_diff, by = c("lakename", "stationname"), all = TRUE)


######################### ZOOPLANKTON REGRESSION TREE ################################

#### only zoop info ####
dat_zoop <- filter(dat, variable %in% c("propotherzoop", "propothercladoc", "proprotifer", 
                                        "propdaphnia", "propcyclopoid", "propcalanoid", "zoopcount")) %>% 
  select(lakename, stationname, year, season, variable, stat, 
         value) %>% 
  filter(stat == "ave") %>% 
  select(-stat)

#reshape data 
zoop_counts <- dcast(dat_zoop, lakename + stationname + year + season ~ variable, value.var = "value")

#find non-rotifer and non-other zoop proportions
zoop_counts <- zoop_counts %>% 
  #calculate zoop counts
  mutate(n.calanoid = propcalanoid*zoopcount,
         n.cyclopoid= propcyclopoid*zoopcount,
         n.daphnia = propdaphnia*zoopcount,
         n.othercladoc = propothercladoc*zoopcount,
         
         n.rotifer = proprotifer*zoopcount,
         n.otherzoop = propotherzoop*zoopcount,
         
         #replace NAs with 0
         n.rotifer = ifelse(is.na(n.rotifer), 0, n.rotifer),
         n.otherzoop = ifelse(is.na(n.otherzoop), 0, n.otherzoop)) %>% 
           
  select(lakename, stationname, year, season, zoopcount,
         n.calanoid, n.cyclopoid, n.daphnia, n.rotifer, n.othercladoc, n.otherzoop) %>% 
    
  #new total
  mutate(new.total.count = zoopcount - n.rotifer - n.otherzoop,
         #new proportions
         prop.calanoid = n.calanoid/new.total.count,
         prop.cyclopoid = n.cyclopoid/new.total.count,
         prop.daphnia = n.daphnia/new.total.count,
         prop.othecladoc = n.othercladoc/new.total.count) %>% 
  select(-zoopcount, -new.total.count, -n.rotifer, -n.otherzoop) 

#melt back
dat_zoop_counts <- melt(zoop_counts, id.vars = c("lakename", "stationname", "year", "season"))

#find average winter and summer communitites across years
dat_zoop_avg <- dat_zoop_counts %>% 
  #find average winter and average summer for lake/station/variable
  group_by(lakename, stationname, season, variable) %>% 
  summarize(avg_value = mean(value, na.rm = TRUE)) %>% 
  as.data.frame()

#find seasonal difference in average values
diff_zoop <- dcast(dat_zoop_avg, lakename + stationname + variable ~ season, value.var = "avg_value")

diff_zoop <- diff_zoop %>% 
                  mutate(WS_diff = iceon-iceoff) %>% 
                  select(-iceoff, -iceon)

#make wide
diff_zoop_wide <- dcast(diff_zoop, lakename + stationname ~ variable, value.var = "WS_diff")
#here the values are the W-S difference

#merge in metadata on physical characteristics
diff_zoop_wide <- merge(dat.full, diff_zoop_wide, by = c("lakename", "stationname"), all = TRUE)

diff_zoop_wide_mv <- na.omit(diff_zoop_wide) #27 lakes (lost the Finnish lakes)

head(diff_zoop_wide_mv)


#### reg tree for W-S difference in zoop proportions ####

prop_zoop_reg_dat <- diff_zoop_wide_mv[,13:16] #diff in proportion winter-sumer (not transformed)

prop_diff_tree <- mvpart(data.matrix(prop_zoop_reg_dat) ~ lat + area + maxdepth + elev + 
                                                          SW_diff_watemp + summer.chla, 
                         data = diff_zoop_wide_mv,
                         xv = "min", xval=nrow(diff_zoop_wide_mv), xvmult=100,
                         which=4, legend = TRUE)

#find specific lakes that split out based on reg tree
filter(diff_zoop_wide_mv, maxdepth < 2.15) 

filter(diff_zoop_wide_mv, maxdepth >= 2.15 & lat < 42.89)

bio_tree <- mvpart(data.matrix(prop_zoop_reg_dat) ~ SW_diff_watemp + summer.chla, 
                         data = diff_zoop_wide_mv,
                         xv = "min", xval=nrow(diff_zoop_wide_mv), xvmult=100,
                         which=4, legend = TRUE)

filter(diff_zoop_wide_mv, SW_diff_watemp >= 20.372) 


#adding summer chla added depth split
#but not as helpful, since also meant we lost the Finnish lakes
summary(prop_diff_tree)
printcp(prop_diff_tree)

#prune
pfit.prop.diff <- prune(prop_diff_tree, cp= prop_diff_tree$cptable[which.min(prop_diff_tree$cptable[,"xerror"]),"CP"])
plot(pfit.prop.diff, uniform = TRUE, compress = T)
text(pfit.prop.diff) #pruning doesn't change tree


######################### PHYTOPLANKTON REGRESSION TREE ##############################

#### only phyto info ####
dat_phyto <- filter(dat, variable %in% c("propchloro", "propcrypto", "propcyano",
                                         "propdiatom", "propdino", "propotherphyto")) %>% 
  select(lakename, stationname, year, season, variable, stat, 
         value) %>% 
  filter(stat == "ave") %>% 
  select(-stat)

#find average winter and summer communitites across years
dat_phyto_avg <- dat_phyto %>% 
  #find average winter and average summer for lake/station/variable
  group_by(lakename, stationname, season, variable) %>% 
  summarize(avg_value = mean(value, na.rm = TRUE)) %>% 
  as.data.frame()

#find seasonal difference
diff_phyto <- dcast(dat_phyto_avg, lakename + stationname + variable ~ season, value.var = "avg_value")

diff_phyto <- diff_phyto %>% 
  mutate(WS_diff = iceon-iceoff) %>% 
  select(-iceoff, -iceon)

#make wide
diff_phyto_wide <- dcast(diff_phyto, lakename + stationname ~ variable, value.var = "WS_diff")
#here the values are the W-S difference

#merge in metadata on physical characteristics
diff_phyto_wide <- merge(dat.full, diff_phyto_wide, by = c("lakename", "stationname"), all = TRUE)

diff_phyto_wide_mv <- na.omit(diff_phyto_wide) #17 lakes

head(diff_phyto_wide_mv)

#### regression tree for W-S difference in phyto proportions ####

prop_phyto_reg_dat <- diff_phyto_wide_mv[,9:14] #diff in proportion winter-sumer (not transformed)

#test_dat <- test[,12:16]

prop_diff_tree <- mvpart(data.matrix(prop_phyto_reg_dat) ~ lat + area + maxdepth + elev + 
                                                            SW_diff_watemp + summer.chla, 
                         data = diff_phyto_wide_mv,
                         xv = "min", xval=nrow(diff_phyto_wide_mv), xvmult=100,
                         which=4, legend = TRUE)

summary(prop_diff_tree)
printcp(prop_diff_tree)

