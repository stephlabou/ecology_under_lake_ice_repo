###################################################################################
# This script makes regression trees for under ice variables of interest,         #
# such as chla and DOC, using lake physical characteristics as predictors.        #
# This script also outputs the base for Figure S7 in Hampton et al. 2016.         #
#                                                                                 #
# Note that we are using only complete cases for explanatory variables,           #
# so the final counts for each variable will be slightly off from the values in   #
# Hampton et al. 2016. For instance, Fish lake chla won't get included because    #
# the maxdepth is NA.                                                             # 
###################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
################################################################################### 

############################ Code by Stephanie Labou ##############################

#### clear workspace ####
rm(list = ls()) 
graphics.off()

#### load packages ####
library(dplyr)
library(rpart)
library(rpart.plot)
library(reshape2)
library(plotrix)

#################### read in data #############################
# The data .csv file is created by the script create.data.long.R
data.orig <-read.csv("Data/data.long.csv", stringsAsFactors = FALSE)

#################### reorganize data #######################
#we want absolute latitude
data <- data.orig 
data <- data %>% 
              mutate(stationlat = abs(stationlat))

#metadata about lake physical properties
meta <- data %>% 
            select(lakename, stationname, stationlat, 
                    stationlong, lakearea, lakemaxdepth,
                    lakeelevation, year, season) %>% 
            unique()

meta <- meta %>% 
        group_by(lakename, stationname) %>% 
        #because some eople varied by year for lake max depth
        summarize(lat = mean(stationlat, na.rm = TRUE),
                  area = mean(lakearea, na.rm = TRUE),
                  maxdepth = mean(lakemaxdepth, na.rm = TRUE),
                  elev = mean(lakeelevation, na.rm = TRUE)) %>% 
        as.data.frame()

#data
data_vals <-data %>% 
            select(lakename, stationname, variable, 
                    stat, year, season, value) %>% 
            filter(stat == "ave") %>% 
            select(-stat)

#Replace zero values for nutrients (likely actually just below detection limit)
#TDN: 1, TP: 0.5, TDP: 0.5
#(note that there are only zero values for TP, TDN, and TDP, so TN replace isn't necessary)
#CAVEAT: Data provided in the Morpho data package is provided as is. 
#Researchers are encouraged to use their expert judgement as to values potentially below detectable limits
data_vals$value[which(data_vals$variable=="totphos" & data_vals$value==0)]<-0.5
data_vals$value[which(data_vals$variable=="totdissphos" & data_vals$value==0)]<-0.5
data_vals$value[which(data_vals$variable=="totdissnitro" & data_vals$value==0)]<-1


#### find average winter-summer difference ####

#reshape
dat_seas <- dcast(data_vals, lakename +stationname + year + variable ~ season, value.var = "value")

#find average winter-summer diff (diff of winter avg and summer avg across years)
dat_avg <- dat_seas %>% 
          group_by(lakename, stationname, variable) %>% 
          summarize(winter.avg = mean(iceon, na.rm = TRUE),
                    summer.avg = mean(iceoff, na.rm = TRUE),
                    #finding avg winter and summer and finding diff by avg winter - avg summer
                    winter_summer_diff_avg = winter.avg - summer.avg) %>% 
          as.data.frame()

#merge with metadata
dat.full <- merge(meta, dat_avg, by = c("lakename", "stationname"), all = TRUE) 
str(dat.full) #6210

#get rid of NAs
dat_regtree <- filter(dat.full, !is.na(winter_summer_diff_avg))
str(dat_regtree) #1688
summary(dat_regtree)

#remove rows with no max depth
dat_regtree <- filter(dat_regtree, !is.na(maxdepth))
str(dat_regtree) #1680
summary(dat_regtree)

#using lat/area/maxdepth/elev as predictors (simple tree)
dat_tree_basic <- dat_regtree %>% select(-winter.avg, -summer.avg) 
str(dat_tree_basic) #1680

#rename variables
dat_tree_basic <- rename(dat_tree_basic, Latitude = lat, Area = area, 
                         MaxDepth = maxdepth, Elevation = elev)

######################################################################
##################### Make regression trees ##########################
######################################################################

####################### chla #######################

dat.chla <- filter(dat_tree_basic, variable == "chla")

#difference in chla (winter minus summer) is response (Y)
fit.chla <- rpart(winter_summer_diff_avg ~ Latitude + Area + MaxDepth + Elevation, xval =100, 
                  data = dat.chla)

printcp(fit.chla)
plotcp(fit.chla)
summary(fit.chla)

plot(fit.chla)
text(fit.chla)

#prune tree using best practices
pfit.chla <- prune(fit.chla, cp= fit.chla$cptable[which.min(fit.chla$cptable[,"xerror"]),"CP"])
summary(pfit.chla) 
printcp(pfit.chla)

#write to png
png("Outputs/Figures/regtree_chla.png", width = 10, height = 8, units = "in", res = 1200)
rpart.plot(pfit.chla, type = 4, extra=1, box.palette = 0, #branch.lty=3,
           main = "Chl a winter-summer difference")
dev.off()

#tag chla df to find se for means of each branch

head(dat.chla)

dat.chla.tag <- dat.chla %>% 
                mutate(branch = ifelse(MaxDepth < 2.9, "-52.3 branch", "-2.28 branch")) %>% 
                group_by(branch) %>% 
                summarize(mean = mean(winter_summer_diff_avg, na.rm = TRUE),
                          SD = sd(winter_summer_diff_avg, na.rm = TRUE),
                          SE = std.error(winter_summer_diff_avg, na.rm = TRUE)) %>% 
                as.data.frame()
dat.chla.tag

########################### TDN ############################

dat.totdissnitro <- filter(dat_tree_basic, variable == "totdissnitro")

#difference in totdissnitro (winter minus summer) is response (Y)
fit.totdissnitro <- rpart(winter_summer_diff_avg ~ Latitude + Area + MaxDepth + Elevation, xval =100, 
                          data = dat.totdissnitro)

printcp(fit.totdissnitro)
plotcp(fit.totdissnitro)
summary(fit.totdissnitro)

plot(fit.totdissnitro)
text(fit.totdissnitro)

#prune tree using best practices
pfit.totdissnitro <- prune(fit.totdissnitro, cp= fit.totdissnitro$cptable[which.min(fit.totdissnitro$cptable[,"xerror"]),"CP"])
summary(pfit.totdissnitro) 
printcp(pfit.totdissnitro)

#write to png
png("Outputs/Figures/regtree_TDN.png", width = 10, height = 8, units = "in", res = 1200)
rpart.plot(pfit.totdissnitro, type = 4, extra=1, box.palette = 0, #branch.lty=3,
           main = "TDN winter-summer difference")
dev.off()

#tag dat.TDN for SE and SD calcs

dat.TDN.tag <- dat.totdissnitro %>% 
                mutate(branch = ifelse(MaxDepth < 2.1, "2070 group", "other"),
                branch = ifelse(MaxDepth >= 2.1 & MaxDepth < 5.2, "758 group", branch),
                branch = ifelse(MaxDepth >= 2.1 & MaxDepth >= 5.2, "123 group", branch)) %>% 
                group_by(branch)%>% 
                summarize(mean = mean(winter_summer_diff_avg, na.rm = TRUE),
                          SD = sd(winter_summer_diff_avg, na.rm = TRUE),
                          SE = std.error(winter_summer_diff_avg, na.rm = TRUE)) %>% 
                as.data.frame()

dat.TDN.tag

######################### TP #########################

dat.totphos <- filter(dat_tree_basic, variable == "totphos")

#difference in totphos (winter minus summer) is response (Y)
fit.totphos <- rpart(winter_summer_diff_avg ~ Latitude + Area + MaxDepth + Elevation, xval =100, 
                     data = dat.totphos)

printcp(fit.totphos)
plotcp(fit.totphos)
summary(fit.totphos)

plot(fit.totphos)
text(fit.totphos)

#prune tree using best practices
pfit.totphos <- prune(fit.totphos, cp= fit.totphos$cptable[which.min(fit.totphos$cptable[,"xerror"]),"CP"])
summary(pfit.totphos) 

#write to png
png("Outputs/Figures/regtree_TP.png", width = 10, height = 8, units = "in", res = 1200)
rpart.plot(pfit.totphos, type = 4, extra=1, box.palette = 0, #branch.lty=3,
           main = "TP winter-summer difference")
dev.off()

###################### TDP ######################

dat.totdissphos <- filter(dat_tree_basic, variable == "totdissphos")

#difference in totdissphos (winter minus summer) is response (Y)
fit.totdissphos <- rpart(winter_summer_diff_avg ~ Latitude + Area + MaxDepth + Elevation, xval =100, 
                         data = dat.totdissphos)

printcp(fit.totdissphos)
plotcp(fit.totdissphos)
summary(fit.totdissphos)

plot(fit.totdissphos)
text(fit.totdissphos)

#prune tree using best practices
pfit.totdissphos <- prune(fit.totdissphos, cp= fit.totdissphos$cptable[which.min(fit.totdissphos$cptable[,"xerror"]),"CP"])
summary(pfit.totdissphos) 

#write to png
png("Outputs/Figures/regtree_totdissphos.png", width = 10, height = 8, units = "in", res = 1200)
rpart.plot(pfit.totdissphos, type = 4, extra=1, box.palette = 0, #branch.lty=3,
           main = "TDP winter-summer difference")
dev.off()


######################## DOC ######################

dat.totdoc <- filter(dat_tree_basic, variable == "totdoc")

#difference in totdoc (winter minus summer) is response (Y)
fit.totdoc <- rpart(winter_summer_diff_avg ~ Latitude + Area + MaxDepth + Elevation, xval =100, 
                    data = dat.totdoc)

printcp(fit.totdoc)
plotcp(fit.totdoc)
summary(fit.totdoc)

plot(fit.totdoc)
text(fit.totdoc)

#prune tree using best practices
pfit.totdoc <- prune(fit.totdoc, cp= fit.totdoc$cptable[which.min(fit.totdoc$cptable[,"xerror"]),"CP"])
summary(pfit.totdoc) 
printcp(pfit.totdoc)

#write to png
png("Outputs/Figures/regtree_DOC.png", width = 10, height = 8, units = "in", res = 1200)
rpart.plot(pfit.totdoc, type = 4, extra=1, box.palette = 0, #branch.lty=3,
           main = "DOC winter-summer difference")
dev.off()


#tag DOC data
dat.doc.tag <- dat.totdoc %>% 
                mutate(branch = ifelse(Area >= 0.373, "-0.145 group", "other"),
                       branch = ifelse(Area < 0.373 & Elevation >= 366, "0.810 group", branch),
                       branch = ifelse(Area < 0.373 & Elevation < 366, "6.69 group", branch)) %>% 
                group_by(branch) %>% 
                summarize(mean = mean(winter_summer_diff_avg, na.rm = TRUE),
                          SD = sd(winter_summer_diff_avg, na.rm = TRUE),
                          SE = std.error(winter_summer_diff_avg, na.rm = TRUE)) %>% 
                as.data.frame()

dat.doc.tag


######################## water temp ######################

dat.watertemp <- filter(dat_tree_basic, variable == "watertemp")

#difference in watertemp (winter minus summer) is response (Y)
fit.watertemp <- rpart(winter_summer_diff_avg ~ Latitude + Area + MaxDepth + Elevation, xval =100, 
                    data = dat.watertemp)

printcp(fit.watertemp)
plotcp(fit.watertemp)
summary(fit.watertemp)

plot(fit.watertemp)
text(fit.watertemp)

#prune tree using best practices
pfit.watertemp <- prune(fit.watertemp, cp= fit.watertemp$cptable[which.min(fit.watertemp$cptable[,"xerror"]),"CP"])
summary(pfit.watertemp) 

#write to png
png("Outputs/Figures/regtree_watertemp.png", width = 10, height = 8, units = "in", res = 1200)
rpart.plot(pfit.watertemp, type = 4, extra=1, box.palette = 0, #branch.lty=3,
           main = "Water temp winter-summer difference")
dev.off()

