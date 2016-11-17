###################################################################################
# This script runs linear mixed models on each variable (chla, DOC, etc.) to      #
# estimate the "bice" fixed effect (difference between winter and summer),        #
# with year as a random effect. The summary lme stats are then compiled, variable # 
# by variable, into an output file named "lme.txt."                               #                             
#                                                                                 #                                             
# Note: This script depends on another script, lme_iteratethis.R, which is called # 
# many times within.                                                              #
###################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
################################################################################### 

############################ Code by Steve Powers #################################

#### clear workspace ####
rm(list = ls()) 
graphics.off()

######### load packages ############
library(grid)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(data.table)
library(exactRankTests)
library(RColorBrewer)
library(tidyr)
library(nlme)
library(MuMIn)
library(pls)

###############################################################
################## Load and process data ######################
###############################################################

# The data .csv file is created by the script create.data.long.R
Data.orig <-read.csv("Data/data.long.csv", stringsAsFactors = FALSE)

data<-select(Data.orig, lakename, stationname, variable, stat, 
             stationlat, stationlong, year, season, value, 
             old.variable, poolstation, value.before, value.after, 
             nYears, n.season, nFullYears, lakemaxdepth,lakearea,lakeelevation)


#restrict to variable averages (not max, min, or CV)
data<-subset(data,data$stat=="ave")

#replace zero values for nutrients
#TDN: 1, TP: 0.5, TDP: 0.5
#(note that there are only zero values for TP, TDN, and TDP, so TN replace isn't necessary)
#CAVEAT: Data provided in the Morpho data package is provided as is. 
#Researchers are encouraged to use their expert judgement as to values potentially below detectable limits
data$value[which(data$variable=="totphos" & data$value==0)]<-0.5
data$value[which(data$variable=="totdissphos" & data$value==0)]<-0.5
data$value[which(data$variable=="totdissnitro" & data$value==0)]<-1

#drop NA values and select subset of columns
data<-data[which(is.na(data$value)==FALSE),]
data<-select(data, lakename, stationname, stationlat, stationlong, 
             year, season, stat, variable,value,
             lakemaxdepth,lakearea,lakeelevation)


############################################################################
###################### Find TN/TP and TDN/TDP ##############################
############################################################################

#include only average values, case to wide format 
np.data <- filter(data, stat == "ave")
np.data <- select(np.data, lakename, stationname, stationlat, stationlong,
                  year, season, stat, variable, value,
                  lakemaxdepth,lakearea,lakeelevation)
np.data.cast <- dcast(np.data, lakename + stationname + stationlat + stationlong +
                        year + season + stat + lakemaxdepth + lakearea + lakeelevation ~ variable, value.var = "value")
np.data.cast <- select(np.data.cast, lakename, stationname, stationlat, stationlong, year, 
                       season, stat, totnitro, totphos, totdissnitro, totdissphos,totdoc,
                       lakemaxdepth,lakearea,lakeelevation)


#calculate TN/TP and TDN/TDP
#use atomic N:P, N=14 g/mmol, P=31 g/mol
np.calc <- np.data.cast %>% 
  mutate(TN.TP = (totnitro/14)/(totphos/31),
         TDN.TDP = (totdissnitro/14)/(totdissphos/31))

#here, if N or P is NA, then the ratio is also NA
np.final <- np.calc %>% select(lakename, stationname, stationlat, stationlong, year, 
                               season, stat, TN.TP, TDN.TDP,
                               lakemaxdepth,lakearea,lakeelevation)

#melt back into long format
np.melt <- melt(np.final, id.var = c("lakename", "stationname", "stationlat", "stationlong", "year", 
                                     "season", "stat", "lakemaxdepth","lakearea","lakeelevation"))

#merge np.data and bulk data 
data.nutrient<-rbind(data,np.melt)
data.nutrient<-arrange(data.nutrient, lakename, stationname, 
                       stationlat, stationlong, year, season, variable,
                       lakemaxdepth,lakearea,lakeelevation)

data <- data.nutrient

#######################################################################
######## find non-rotifer, non-other zoop zooplankton count ###########
#######################################################################

#find zoopcount when non-rotifer, non-otherzoop
data.zoop1<-subset(data,data$variable=="zoopcount")
data.zoop2<-subset(data,data$variable=="proprotifer")
data.zoop3<-subset(data,data$variable=="propotherzoop")


data.zoop.merge<-merge(data.zoop1,data.zoop2,by=c("lakename","stationname","year","season","stationlat","stationlong",
                                                  "lakemaxdepth","lakearea","lakeelevation","stat"),all=TRUE)

data.zoop.merge<-merge(data.zoop.merge,data.zoop3,by=c("lakename","stationname","year","season","stationlat","stationlong",
                                      "lakemaxdepth","lakearea","lakeelevation","stat"),all=TRUE)

data.zoop.merge$value.y[which(is.na(data.zoop.merge$value.y)==TRUE)]<-0
data.zoop.merge$value[which(is.na(data.zoop.merge$value)==TRUE)]<-0
data.zoop.merge<-data.zoop.merge[-which(is.na(data.zoop.merge$value.x)==TRUE),]
data.zoop.merge$zoopcount_nonrotifer<-data.zoop.merge$value.x-(data.zoop.merge$value.x*data.zoop.merge$value.y)-(data.zoop.merge$value.x*data.zoop.merge$value)
data.zoop.merge$variable.x<-"zoopcount_nonrotifer"

data.zoop.nonrotifer<-data.zoop.merge %>% select(lakename,stationname,stationlat,stationlong,year,
                                                 season,stat,variable=variable.x,value=zoopcount_nonrotifer,
                                                 lakemaxdepth,lakearea,lakeelevation)

###########################################################
############# Merge all, clean data frame #################
###########################################################

#bind new zoop data to rest of data
data<-rbind(data,data.zoop.nonrotifer)

data.byyear<-data
data.byyear$lakestation<-paste(data.byyear$lakename,data.byyear$stationname)

#reshape and clean
data.byyear.season<-dcast(data.byyear, lakename + stationname + stationlat + stationlong +
                        year + stat + lakemaxdepth + lakearea + lakeelevation+variable~season, value.var = "value")
data.byyear.season$lakestation<-paste(data.byyear.season$lakename,data.byyear.season$stationname)
data.byyear.season$stationlat.abs<-abs(data.byyear.season$stationlat)
data.byyear.season$diff<-data.byyear.season$iceon-data.byyear.season$iceoff

###############################################################################################
############## Check for long-term linear relationships in moving-averaged data ###############
###############################################################################################

data.covars <- 
  data.byyear %>% group_by(lakename,lakestation) %>%
  dplyr::summarize(lakemaxdepth=mean(lakemaxdepth,na.rm=TRUE),
                   lakearea=mean(lakearea,na.rm=TRUE),
                   lakeelevation=mean(lakeelevation,na.rm=TRUE)) %>% as.data.frame()

data.byyear<-merge(data.byyear,data.covars,by=c("lakename","lakestation"),all=TRUE)
data.byyear<-data.byyear %>% select(-c(lakemaxdepth.x,lakearea.x,lakeelevation.x))#%>% head()
data.byyear<-data.byyear %>% mutate(lakemaxdepth=lakemaxdepth.y,lakearea=lakearea.y,lakeelevation=lakeelevation.y)# %>% head()
data.byyear<-data.byyear %>% select(-c(lakemaxdepth.y,lakearea.y,lakeelevation.y))

#find median values
lakearea.median<-median(data.covars$lakearea)
lakemaxdepth.median<-median(na.omit(data.covars$lakemaxdepth))
lakeelevation.median<-median(data.covars$lakeelevation)
stationlat.median<-median(unique(abs(data.byyear$stationlat)))

#tag data by type (i.e. deeper, larger, polar, etc.)
data.byyear$deeperTF<-0
data.byyear$deeperTF[which(data.byyear$lakemaxdepth>=lakemaxdepth.median)]<-1
data.byyear$largerTF<-0
data.byyear$largerTF[which(data.byyear$lakearea>=lakearea.median)]<-1
data.byyear$higherTF<-0
data.byyear$higherTF[which(data.byyear$lakeelevation>=lakeelevation.median)]<-1
data.byyear$polarTF<-0
data.byyear$polarTF[which(abs(data.byyear$stationlat)>=stationlat.median)]<-1

data.byyear$winterTF<-0
data.byyear$winterTF[which(data.byyear$season=="iceon")]<-1


########################## chla ################################
data.lme<-subset(data.byyear,data.byyear$variable=="chla" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="chla" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt")
print("#####################################")
print("chla")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)

print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

#################### old zoopcount (includes all categories) #########################
data.lme<-subset(data.byyear,data.byyear$variable=="zoopcount" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="zoopcount" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("zoopcount")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)

print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

######################## zoopcount (non-rotifer, non-other zoop ###########################
data.lme<-subset(data.byyear,data.byyear$variable=="zoopcount_nonrotifer" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="zoopcount_nonrotifer" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("zoopcount_nonrotifer")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)

print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

############################### TP #######################################
data.lme<-subset(data.byyear,data.byyear$variable=="totphos" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="totphos" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("TP")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)

print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

############################## TDN ####################################
data.lme<-subset(data.byyear,data.byyear$variable=="totdissnitro" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="totdissnitro" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("TDN")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

############################### DOC ###################################
data.lme<-subset(data.byyear,data.byyear$variable=="totdoc" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="totdoc" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("DOC")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

########################## phyto mass #################################
data.lme<-subset(data.byyear,data.byyear$variable=="phytomass" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="phytomass" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("phytomass")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

############################## TN #################################
data.lme<-subset(data.byyear,data.byyear$variable=="totnitro" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="totnitro" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("TN")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

############################### TDP ################################
data.lme<-subset(data.byyear,data.byyear$variable=="totdissphos" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="totdissphos" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("TDP")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

############################### TN:TP ###################################
data.lme<-subset(data.byyear,data.byyear$variable=="TN.TP" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="TN.TP" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

#WARNING, TN.TP stats in the covariate analysis are very sensitive to a huge outlier
# let's drop the outlier
# # 749 and 750
data.lme<-data.lme[-c(749,750),]

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("TN:TP")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

########################### TDN:TDP #############################
data.lme<-subset(data.byyear,data.byyear$variable=="TDN.TDP" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="TDN.TDP" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

#WARNING, TDN.TDP stats in the covariate analysis are very sensitive to a huge outlier
# let's drop the outlier
# data.lme[150,]
data.lme<-data.lme[-150,]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("TDN.TDP")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

########################### water temp ################################
data.lme<-subset(data.byyear,data.byyear$variable=="watertemp" & data.byyear$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$value)==FALSE),]
main.lme<-lme(value ~ season,random=~1|year, data=data.lme)

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="watertemp" & data.byyear.season$stat=="ave")
data.lme<-data.lme[which(is.na(data.lme$diff)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakearea)==FALSE),]
data.lme<-data.lme[which(is.na(data.lme$lakemaxdepth)==FALSE),]

#WARNING, watertemp stats in the covariate analysis are very sensitive to a huge outlier
# let's drop the outlier
# data.lme[183,]
data.lme<-data.lme[-183,]

data.lme$lakearea<-log10(data.lme$lakearea)
data.lme$lakemaxdepth<-log10(data.lme$lakemaxdepth)
data.lme$lakeelevation<-log10(data.lme$lakeelevation)

source("Scripts/Analysis-Figures/lme_iteratethis.R")

sink("Outputs/Files/lme.txt",append=TRUE)
print("#####################################")
print("watertemp")
print("main.lme")
summary(main.lme)
print("#####################################")
print("r2 main")
r.squaredGLMM(main.lme)
print("#####################################")
print("anova")
#anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme)
print(df)
print("#####################################")
print("covars.lme")
summary(covars.lme)$tTable
print("#####################################")
print("depth.lme")
summary(depth.lme)$tTable
print("area.lme")
summary(area.lme)$tTable
print("elev.lme")
summary(elev.lme)$tTable
print("lat.lme")
summary(lat.lme)$tTable
print("#####################################")
print("#####################################")
print("#####################################")
sink()

#######################################################################
######################### water temp example ##########################
#######################################################################

data.lme.agg<-data.lme %>% group_by(lakename,stationname) %>%
  dplyr::summarize(mean.diff=mean(diff),
                   lakemaxdepth=mean(lakemaxdepth,na.rm=TRUE),
                   lakearea=mean(lakearea,na.rm=TRUE),
                   lakeelevation=mean(lakeelevation,na.rm=TRUE),
                   stationlat.abs=mean(stationlat.abs,na.rm=TRUE)) %>% as.data.frame()

name<-paste("Outputs/Figures/Example.covar.",data.lme$variable[1],".png",sep="")
png(name,height=8,width=8,units="in",res=300)
par(mfrow=c(2,2))

plot((data.lme.agg$lakemaxdepth),data.lme.agg$mean.diff, 
     xlab="log10 max depth (m)", ylab="winter minus summer")#,ylim=c(-1000,5000))

plot((data.lme.agg$lakearea),data.lme.agg$mean.diff, 
     xlab="log10 lake area (km2)", ylab="winter minus summer")#,ylim=c(-1000,5000))

plot((data.lme.agg$lakeelevation),data.lme.agg$mean.diff, 
     xlab="log10 lake elevation (m)", ylab="winter minus summer")#,ylim=c(-1000,5000))

plot((data.lme.agg$stationlat.abs),data.lme.agg$mean.diff, 
     xlab="latitude (absolute)", ylab="winter minus summer")#,ylim=c(-1000,5000))

mtext(text=data.lme$variable[1],side=3,outer=TRUE,padj=3)#,padj=-1)
dev.off()

###############################################
############ covariate crossplots #############
###############################################

data.lme<-subset(data.byyear.season,data.byyear.season$variable=="chla" & data.byyear.season$stat=="ave")

png("Outputs/Figures/covariate.crossplot.png",height=6,width=8,units="in",res=300)
par(mfrow=c(2,3))
plot(data.lme.agg$stationlat.abs,data.lme.agg$lakemaxdepth, 
     xlab="latitude", ylab="max depth (m)")#,ylim=c(-1000,5000))
plot(data.lme.agg$stationlat.abs,data.lme.agg$lakearea, 
     xlab="latitude", ylab="surface area (km2)")#,ylim=c(-1000,5000))
plot(data.lme.agg$stationlat.abs,data.lme.agg$lakeelevation, 
     xlab="latitude", ylab="elevation (m)")#,ylim=c(-1000,5000))
plot(data.lme.agg$lakemaxdepth, data.lme.agg$lakearea, 
     xlab="max depth (m)", ylab="surface area (km2)")#,ylim=c(-1000,5000))
plot(data.lme.agg$lakeelevation,data.lme.agg$lakedepth, 
     xlab="elevation (m)", ylab="max depth (m)")#,ylim=c(-1000,5000))
plot(data.lme.agg$lakeelevation,data.lme.agg$lakearea, 
     xlab="elevation (m)", ylab="surface area (km2)")#,ylim=c(-1000,5000))
dev.off()
