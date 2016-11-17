###################################################################################
# This script produces:                                                           #                                          
# i.  A multipanel boxplot + scatterplot figure of winter-summer average          #
# conditions for multiple variables                                               #
#       across all lakes                                                          #                                                       
#        -The variables (one per panel) are: chl a, zoop abundance, TP, TDN,      #
#            DOC, and phyto biomass                                               #
#        -This is Fig. 2 in Hampton et al. 2016                                   #                                     
# ii. Statistical summaries for "primary"  limnological variables, and a few      #
#     additional variables                                                        #
#        -This is Table S3 in Hampton et al. 2016                                 #
#                                                                                 #                             
# This script also produces a table of summary "coverage" statistics of each      #
# variable (e.g. min, max, etc.). This is Table S1 in Hampton et al. 2016.        #
###################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
################################################################################### 

################## Code by Steve Powers and Stephanie Labou #######################

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
library(plotrix)

###############################################################
################## Load and process data ######################
###############################################################

# The data .csv file is created by the script create.data.long.R
Data.orig <-read.csv("Data/data.long.csv", stringsAsFactors = FALSE)

data<-select(Data.orig, lakename, stationname, variable, stat, 
             stationlat, stationlong, year, season, value, 
             old.variable, poolstation, value.before, value.after, 
             nYears, n.season, nFullYears)

#restrict to variable averages (not max, min, or CV)
data<-subset(data,data$stat=="ave")

## find non-rotifer zoop counts ##
data_cast <- dcast(data, lakename + stationname + stationlat + stationlong +
                     year + season ~ variable, value.var = "value")

#non-rotifer and non-other zoop new zoop counts
data_new <- data_cast %>% 
  #calculate zoop counts
  mutate(n.rotifer = proprotifer*zoopcount,
         n.otherzoop = propotherzoop*zoopcount, 
         n.rotifer = ifelse(is.na(n.rotifer), 0, n.rotifer),
         n.otherzoop = ifelse(is.na(n.otherzoop), 0, n.otherzoop),
         #new total
         new.zoop.count = zoopcount - n.rotifer - n.otherzoop) %>% 
  select(-n.rotifer, n.otherzoop)

#melt back
data_use <- melt(data_new, id.vars = c("lakename", "stationname", "stationlat", "stationlong",
                                        "year", "season"))

data <- data_use

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
data<-select(data, lakename, stationname, stationlat, stationlong, year, season, variable,value)


############################################################################
###################### Find TN/TP and TDN/TDP ##############################
############################################################################
np.data <- data

np.data <- select(np.data, lakename, stationname, stationlat, stationlong,
                  year, season, variable, value)
np.data.cast <- dcast(np.data, lakename + stationname + stationlat + stationlong +
                        year + season ~ variable, value.var = "value")
np.data.cast <- select(np.data.cast, lakename, stationname, stationlat, stationlong, year, 
                       season, totnitro, totphos, totdissnitro, totdissphos,totdoc)

#calculate TN/TP and TDN/TDP
#use atomic N:P, N=14 g/mmol, P=31 g/mol
np.calc <- np.data.cast %>% 
  mutate(TN.TP = (totnitro/14)/(totphos/31),
         TDN.TDP = (totdissnitro/14)/(totdissphos/31))

#here, if N or P is NA, then the ratio is also NA
np.final <- np.calc %>% select(lakename, stationname, stationlat, stationlong, year, 
                               season, TN.TP, TDN.TDP)

#melt back into long format
np.melt <- melt(np.final, id.var = c("lakename", "stationname", "stationlat", "stationlong", "year", "season"))

#merge np.data and bulk data 
data.nutrient<-rbind(data,np.melt)
data.nutrient<-arrange(data.nutrient, lakename, stationname, 
                       stationlat, stationlong, year, season, variable)

data <- data.nutrient

##############################################################################
##################### Iceon, iceoff correlations #############################
##############################################################################

#calculate lake/station/season variable averages (across years)
data<-data %>% 
  group_by(lakename, stationname, variable,season,stationlat, stationlong) %>% 
  dplyr::summarize(value=mean(value,na.rm = TRUE)) %>% 
  as.data.frame()

#reshape data
data<-reshape(data, direction  = "wide", idvar=c("lakename","stationname", "stationlat", "stationlong", "variable"), timevar="season")

#calculate correlations
data.cor<-na.omit(data) %>%
  group_by(variable) %>%
  dplyr::summarize(corr=cor(value.iceon,value.iceoff),n=length(value.iceon)) %>% 
  as.data.frame()

#check correlations by decreasing n
data<-merge(data,data.cor,by="variable")
print(data.cor[order(data.cor$n,decreasing=TRUE),])

#prepare data for wilcox calculations and scatterplots
data.all.for.wilcox<-data

##############################################################################
######################## Scatterplots and boxplots ###########################
##############################################################################

#Note: Wilcox results not used in final manuscript.
#Wilcox results as a hold over from prior analyses; included here in case
#of future need/use.
#Important!: p-values in plots are from LME models, NOT from Wilcox

#set latitude limits for color bar - used for graphical coloration later
#latitudes < abs(40) set to 40 and later labeled "<40"
#likewise, latitudes > abs(65) set to 65 and later labeled ">65"
data$stationlat[which(data$stationlat<40)]<-40
data$stationlat[which(data$stationlat>65)]<-65
colour.lat.lims<-c(40,65)

fontsize<-12

plot.scatter<-function(dataplot,label,lims){
  plotit<-ggplot(data=na.omit(dataplot),aes(x=value.iceoff, y=value.iceon))+
    geom_point(aes(size=2, color=abs(stationlat)), show_guide = FALSE)+
    xlab("Summer")+ylab("Winter")+
    scale_x_log10()+scale_y_log10()+
    expand_limits(y=c(xylims), x=c(xylims))+
    geom_abline(intercept=0,slope=1, linetype=2) +
    scale_colour_gradient(limits=colour.lat.lims, low="green",high="blue",
                          guide=FALSE,breaks=c(40,45,50,55,60,65))+
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=element_text(vjust=1),
          axis.title.x=element_text(vjust=-0.5))+
    annotation_logticks()+
    #text
    annotation_custom(label)
  }

#### run wilcox test and create scatterplot for each variable of interest ####

#note: p-values are from lme models
#n and r are from data.cor

#### totphos ####
dataplotTP<-subset(data,data$variable=="totphos")
# wilcox.TP<-wilcox.exact(dataplotTP$value.iceon,dataplotTP$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.TP)
TP.text <- grobTree(textGrob("TP (µg/L) \nn=106, p=0.49 \nr=0.85", 
                             hjust=0, x=0.08,  y=0.85,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotTP$value.iceon,dataplotTP$value.iceoff))),
          1.25*max(na.omit(c(dataplotTP$value.iceon,dataplotTP$value.iceoff))))
scat.plotTP<-plot.scatter(dataplotTP,label=TP.text,lims=xylims)
scat.plotTP


#### totnitro ####
dataplotTN<-subset(data,data$variable=="totnitro")
# wilcox.TN<-wilcox.exact(dataplotTN$value.iceon,dataplotTN$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.TN)
TN.text <- grobTree(textGrob(paste("TN (µg/L) \nn=75, p<0.001** \nr=0.97",sep=""), 
                             hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotTN$value.iceon,dataplotTN$value.iceoff))),
          1.25*max(na.omit(c(dataplotTN$value.iceon,dataplotTN$value.iceoff))))
scat.plotTN<-plot.scatter(dataplotTN,label=TN.text,lims=xylims)
scat.plotTN


#### totdissphos ####
dataplotTDP<-subset(data,data$variable=="totdissphos")
# wilcox.TDP<-wilcox.exact(dataplotTDP$value.iceon,dataplotTDP$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.TDP)
TDP.text <- grobTree(textGrob("TDP (µg/L) \nn=72, p=0.21 \nr=0.96", 
                              hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotTDP$value.iceon,dataplotTDP$value.iceoff))),
          1.25*max(na.omit(c(dataplotTDP$value.iceon,dataplotTDP$value.iceoff))))
scat.plotTDP<-plot.scatter(dataplotTDP,label=TDP.text,lims=xylims)
scat.plotTDP


#### totdissnitro ####
dataplotTDN<-subset(data,data$variable=="totdissnitro")
# wilcox.TDN<-wilcox.exact(dataplotTDN$value.iceon,dataplotTDN$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.TDN)
TDN.text <- grobTree(textGrob("TDN (µg/L) \nn=73, p<0.001** \nr=0.93", 
                              hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotTDN$value.iceon,dataplotTDN$value.iceoff))),
          1.25*max(na.omit(c(dataplotTDN$value.iceon,dataplotTDN$value.iceoff))))
scat.plotTDN<-plot.scatter(dataplotTDN,label=TDN.text,lims=xylims)
scat.plotTDN


#### TN/TP ####
dataplotTNTP<-subset(data,data$variable=="TN.TP")
# wilcox.TN.TP<-wilcox.exact(dataplotTNTP$value.iceon,dataplotTNTP$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.TN.TP)
TNTP.text <- grobTree(textGrob("TN:TP \nn=74, p<0.001** \nr=0.39", 
                               hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotTNTP$value.iceon,dataplotTNTP$value.iceoff))),
          1.25*max(na.omit(c(dataplotTNTP$value.iceon,dataplotTNTP$value.iceoff))))
scat.plotTNTP<-plot.scatter(dataplotTNTP,label=TNTP.text,lims=xylims)
scat.plotTNTP


#### TDN/TDP ####
dataplotTDNTDP<-subset(data,data$variable=="TDN.TDP")
# wilcox.TDN.TDP<-wilcox.exact(dataplotTDNTDP$value.iceon,dataplotTDNTDP$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.TDN.TDP)
TDNTDP.text <- grobTree(textGrob("TDN:TDP \nn=66, p=0.50 \nr=0.03 (0.47 with outlier removed)", 
                                 hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
#Note: this very low correlation is due to an outlier at Lake Erie station Erie-central-basin 
#(very low TDP one year gives a TDN.TDP of over 5000)
xylims<-c(0.75*min(na.omit(c(dataplotTDNTDP$value.iceon,dataplotTDNTDP$value.iceoff))),
          1.25*max(na.omit(c(dataplotTDNTDP$value.iceon,dataplotTDNTDP$value.iceoff))))
scat.plotTDNTDP<-plot.scatter(dataplotTDNTDP,label=TDNTDP.text,lims=xylims)
scat.plotTDNTDP


#### chla ####
dataplotchla<-subset(data,data$variable=="chla")
# wilcox.chla<-wilcox.exact(dataplotchla$value.iceon,dataplotchla$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.chla)
chla.text <- grobTree(textGrob("Chl a (µg/L) \nn=118, p<0.001** \nr=0.52", 
                               hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotchla$value.iceon,dataplotchla$value.iceoff))),
          1.25*max(na.omit(c(dataplotchla$value.iceon,dataplotchla$value.iceoff))))
scat.plotchla<-plot.scatter(dataplotchla,label=chla.text,lims=xylims)
scat.plotchla


#### totdoc ####
dataplotDOC<-subset(data,data$variable=="totdoc")
# wilcox.DOC<-wilcox.exact(dataplotDOC$value.iceon,dataplotDOC$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.DOC)
DOC.text <- grobTree(textGrob("DOC (mg/L) \nn=81, p=0.86 \nr=0.98", 
                              hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotDOC$value.iceon,dataplotDOC$value.iceoff))),
          1.25*max(na.omit(c(dataplotDOC$value.iceon,dataplotDOC$value.iceoff))))
scat.plotDOC<-plot.scatter(dataplotDOC,label=DOC.text,lims=xylims)
scat.plotDOC


#### phytomass ####
dataplotphyto<-subset(data,data$variable=="phytomass")
# wilcox.phyto<-wilcox.exact(dataplotphyto$value.iceon,dataplotphyto$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.phyto)
phyto.text <- grobTree(grid.text((bquote(atop(atop(paste("Phyto biovolume (mm" ^3, "/L)"), "n=17, p<0.001**               "), 
                                              atop("r=0.66                             ", "")))),
                                              just = "left", hjust=0, x=0.06, y=0.83, gp=gpar(fontsize=18)))
xylims<-c(0.75*min(na.omit(c(dataplotphyto$value.iceon,dataplotphyto$value.iceoff))),
          1.25*max(na.omit(c(dataplotphyto$value.iceon,dataplotphyto$value.iceoff))))
scat.plotphyto<-plot.scatter(dataplotphyto,label=phyto.text,lims=xylims)
scat.plotphyto


#### zoopcount (non-rotifer, non-other zoop) ####
dataplotzoop<-subset(data,data$variable=="new.zoop.count")
# wilcox.zoop<-wilcox.exact(dataplotzoop$value.iceon,dataplotzoop$value.iceoff, paired = TRUE, alternative = "two.sided")
# print(wilcox.zoop)
zoop.text <- grobTree(textGrob("Crustacean zoop density (#/L) \nn=36, p<<0.001**\nr=0.51", 
                               hjust=0, x=0.06,  y=0.83,gp=gpar(fontsize=fontsize)))
xylims<-c(0.75*min(na.omit(c(dataplotzoop$value.iceon,dataplotzoop$value.iceoff))),
          1.25*max(na.omit(c(dataplotzoop$value.iceon,dataplotzoop$value.iceoff))))
scat.plotzoop<-plot.scatter(dataplotzoop,label=zoop.text,lims=xylims)
scat.plotzoop

#note 1:1 line on zoop graph
one.to.one <- grobTree(textGrob("1:1", hjust=0, x=0.17,  y=0.3,gp=gpar(fontsize=fontsize)))
scat.plotzoop<-scat.plotzoop+annotation_custom(one.to.one)+
  scale_colour_gradient(name="Abs(lat)",limits=colour.lat.lims, low="green",high="blue",
                        breaks=c(40,45,50,55,60,65),  
                        labels = c("<40","45","50","55","60",">65"))+
  theme(legend.position = c(0.85, 0.28),legend.key.size=unit(0.55, "cm"))
scat.plotzoop


######## Summary stats results for all variables ########

#Note: Wilcox analyses were not included in the Hampton et al. 2016 manuscript
#but remain in this script as a hold over from previous analysis attempts.

#drop pointprob column so can be stacked with other wilcox outputs
wilcox.phyto<-wilcox.phyto[-2] 
wilcox.zoop<-wilcox.zoop[-2] 

doc.test <- data.all.for.wilcox %>% filter(variable == "totdoc")
t.test(doc.test$value.iceon, doc.test$value.iceoff)

vars.all<-unique(data.all.for.wilcox$variable)
var<-c();pval<-c();stat<-c();null<-c();alternative<-c();method<-c();
  mean.iceon<-c();mean.iceoff<-c();
  median.iceon<-c();median.iceoff<-c();
  se.iceon<-c();se.iceoff<-c();
  mean.on.minus.off<-c();mean.on.over.off<-c();
  median.on.minus.off<-c();median.on.over.off<-c();
  n.ice<-c();n.paired<-c()

#create dataframe of wilcox values
for(i in 1:length(vars.all)){
  
  #set up per variable
  vari<-as.character(vars.all[i])
  datai<-subset(data.all.for.wilcox,as.character(data.all.for.wilcox$variable)==vari)
  wilcox.i<-wilcox.exact(datai$value.iceon,datai$value.iceoff, paired = TRUE, alternative = "two.sided")
  
  #calculate means (across lakes)
  mean.iceon.i<-mean(datai$value.iceon,na.rm=TRUE)
  mean.iceoff.i<-mean(datai$value.iceoff,na.rm=TRUE)
  
  #calculate medians (across lakes)
  median.iceon.i<-median(datai$value.iceon,na.rm=TRUE)
  median.iceoff.i<-median(datai$value.iceoff,na.rm=TRUE)
  
  #standard deviations of seasonal (across lakes) values
  se.iceon.i<-std.error(datai$value.iceon,na.rm=TRUE)
  se.iceoff.i<-std.error(datai$value.iceoff, na.rm=TRUE)
  
  #number of iceon values and number of paired iceon/iceoff values
  n.ice.i<-length(na.omit(datai$value.iceon))
  n.paired.i<-length(na.omit(data.frame(datai$value.iceon,datai$value.iceoff))[,1])
  
  #means for difference (iceon-iceoff)
  mean.on.minus.off.i<-mean(datai$value.iceon-datai$value.iceoff,na.rm=TRUE)
  
  #median for difference (iceon-iceoff)
  median.on.minus.off.i<-median(datai$value.iceon-datai$value.iceoff,na.rm=TRUE)
  
  #calculating ratios - for paired data, and considering special cases
  #remove instances where iceoff is zero but iceon is non-zero (can't use ratio in those cases)
  sub1 <- subset(datai,(datai$value.iceon==0 & datai$value.iceoff==0)|abs(datai$value.iceoff)>0)
  #just in cases where both iceon and iceoff are zero
  set1 <- which(sub1$value.iceon==0 & sub1$value.iceoff==0)
  sub2 <- sub1
  #in those cases, ratio should be 1
  sub2$value.iceon [set1] <- 1
  sub2$value.iceoff [set1] <- 1
  #take care of instances of NAs (e.g. more complex way of na.rm = TRUE)
  set2 <-which(is.na(sub2$value.iceon)==FALSE & is.na(sub2$value.iceoff)==FALSE)
  #calculate ratios - find mean
  mean.on.over.off.i <- mean(sub2$value.iceon[set2]/sub2$value.iceoff[set2])
  #calculate ratios - find median
  median.on.over.off.i <- median(sub2$value.iceon[set2]/sub2$value.iceoff[set2])
  
  #get stats
  var<-c(var,vari)
  pval<-c(pval,wilcox.i$p.value)
  stat<-c(stat,wilcox.i$statistic)
  null<-c(null,wilcox.i$null)
  alternative<-c(alternative,wilcox.i$alternative)
  method<-c(method,wilcox.i$method)
  
  #append values
  mean.iceon<-c(mean.iceon, mean.iceon.i)
  mean.iceoff<-c(mean.iceoff, mean.iceoff.i)
  
  median.iceon<-c(median.iceon, median.iceon.i)
  median.iceoff<-c(median.iceoff, median.iceoff.i)
  
  se.iceon<-c(se.iceon, se.iceon.i)
  se.iceoff<-c(se.iceoff, se.iceoff.i)
  
  n.ice<-c(n.ice,n.ice.i)
  n.paired<-c(n.paired,n.paired.i)
  
  mean.on.minus.off<-c(mean.on.minus.off,mean.on.minus.off.i)
  median.on.minus.off<-c(median.on.minus.off, median.on.minus.off.i)
  
  mean.on.over.off<-c(mean.on.over.off, mean.on.over.off.i)
  median.on.over.off<-c(median.on.over.off,median.on.over.off.i)  
}


######## final wilcox dataframe ########

wilcox.df<-data.frame(var, n.ice, n.paired, 
                      mean.iceon, mean.iceoff, 
                      median.iceon, median.iceoff,
                      se.iceon, se.iceoff,
                      mean.on.over.off.paired=mean.on.over.off, 
                      median.on.over.off.paired=median.on.over.off, 
                      mean.on.minus.off.paired=mean.on.minus.off, 
                      median.on.minus.off.paired=median.on.minus.off,
                      pval, stat, null, alternative, method)

wilcox.df.clean <- wilcox.df
wilcox.df.clean$var <- as.character(wilcox.df.clean$var)

########################################################
################ summary statistics ####################
########################################################

#data coverage summary statistics for variables

n.counts <- wilcox.df.clean %>% select(var, n.ice, n.paired)

#find number of total observations (non-NA) and number of ice observations (non-NA) 
#per lake/station/variable
Data.yrs <- data.nutrient %>% 
  group_by(lakename, stationname, variable) %>% 
  filter(!is.na(value)) %>% 
  summarize(n.obvs = length(value))

Data.ice.yrs <- data.nutrient %>% 
  filter(season == "iceon") %>% 
  group_by(lakename, stationname, variable) %>% 
  filter(!is.na(value)) %>% 
  summarize(n.ice = length(value))

Data.summary <- merge(Data.yrs, Data.ice.yrs, by = c("lakename", "stationname", "variable"))

#calculate summary statistics
#for each variable, output of summary() for n.ice and n.obvs
n.ice.summary <- tapply(Data.summary$n.ice, Data.summary$variable, summary)
n.ice.summary <- do.call("rbind", n.ice.summary) %>% as.data.frame()
n.ice.summary$n.type <- "n.ice"
n.ice.names <- row.names(n.ice.summary)
n.ice.summary <- cbind(n.ice.names, n.ice.summary)
row.names(n.ice.summary) <- NULL
n.ice.summary <- rename(n.ice.summary, variable = n.ice.names)

n.obvs.summary <- tapply(Data.summary$n.obvs, Data.summary$variable, summary)
n.obvs.summary <- do.call("rbind", n.obvs.summary) %>% as.data.frame()
n.obvs.summary$n.type <- "n.obvs"
n.obvs.names <- row.names(n.obvs.summary)
n.obvs.summary <- cbind(n.obvs.names, n.obvs.summary)
row.names(n.obvs.summary) <- NULL
n.obvs.summary <- rename(n.obvs.summary, variable = n.obvs.names)

n.all.summary <- rbind(n.ice.summary, n.obvs.summary)
n.all.summary <- arrange(n.all.summary, variable, desc(n.type))

#merge with counts (n.ice and n.paired)
n.final <- merge(n.all.summary, n.counts, by.x = "variable", by.y = "var")

n.final <- n.final %>%select(variable, n.ice, n.paired, n.type, Min., `1st Qu.`, Median, Mean, `3rd Qu.`, Max.)

#set groups
n.final$group<-""
n.final$group[which(n.final$variable %in% c("chla","propotherphyto","propcyano",
                                                    "propchloro","propcrypto","propdiatom",
                                                    "phytocount","phytomass", "phytoppr"))]<-"1-phyto"

n.final$group[which(n.final$variable %in% c("zoopcount","propdaphnia","propcalanoid",
                                                    "propcyclopoid", "propothercladoc",
                                                    "propotherzoop","proprotifer","zoopmass",
                                                    "propdino", "zoopcount", "new.zoop.count"))]<-"2-zoop"

n.final$group[which(n.final$variable %in% c("totphos","totnitro","totdissphos","totdissnitro",
                                                    "totdoc", "TDN.TDP", "TN.TP", 
                                                    "C.TDN", "C.TDP"))]<-"3-chem"

n.final$group[which(n.final$variable %in% c("airtemp", "watertemp", "color", "icedepth", 
                                                    "photicdepth", "radiation", "sampledepth",
                                                    "suva", "secchidepth", "snowdepth"))]<-"4-phys"

n.final$group[which(n.final$variable %in% c("bactcount", "bactprod", "benamphdens", 
                                                    "benbivalvedens", "bengastrodens", "beninsectdens",
                                                    "benoligodens", "ciliacount", "ciliamass",
                                                    "hnfcount", "hnfmass"))]<-"5-otherbio"

#rename variables
n.final$variable <- as.character(n.final$variable)
n.final$variable[which(n.final$variable=="chla")]<-"chl a"
n.final$variable[which(n.final$variable=="totdissnitro")]<-"TDN"
n.final$variable[which(n.final$variable=="totdissphos")]<-"TDP"
n.final$variable[which(n.final$variable=="totdoc")]<-"DOC"
n.final$variable[which(n.final$variable=="totnitro")]<-"TN"
n.final$variable[which(n.final$variable=="totphos")]<-"TP"
n.final$variable[which(n.final$variable=="zoopcount")]<-"zoop density"
n.final$variable[which(n.final$variable=="new.zoop.count")]<-"crustacean zoop density"
n.final$variable[which(n.final$variable=="phytocount")]<-"phyto density"
n.final$variable[which(n.final$variable=="phytomass")]<-"phyto biovolume"

n.final <- n.final %>% select(group, variable, n.ice, n.paired, n.type,
                              Min., `1st Qu.`, Median, Mean, `3rd Qu.`, Max.)

n.final <- arrange(n.final, group, variable)

#prep for write to csv
n.final$Mean <- signif(n.final$Mean, 3)

#write to csv
write.csv(n.final, "Outputs/Files/n_final_stats.csv", row.names = FALSE)


########################################################
################## final data frame ####################
########################################################

#keep only variables of interest
vars.keep <- c("chla","propotherphyto","propcyano",
               "propchloro","propcrypto","propdiatom",
               "phytocount","phytomass", "phytoppr", 
               "zoopcount","propdaphnia","propcalanoid",
               "propcyclopoid", "propothercladoc",
               "zoopmass", 
               "new.zoop.count",
               "propdino", "totphos","totnitro","totdissphos",
               "totdissnitro", "totdoc", "TDN.TDP", "TN.TP", 
               "C.TDN", "C.TDP", "watertemp")
  
wilcox.df.clean <- filter(wilcox.df.clean, var %in% vars.keep)


#set groups
wilcox.df.clean$group<-""
wilcox.df.clean$group[which(wilcox.df.clean$var %in% c("chla","propotherphyto","propcyano",
                                                          "propchloro","propcrypto","propdiatom",
                                                          "phytocount","phytomass", "phytoppr"))]<-"1-phyto"

wilcox.df.clean$group[which(wilcox.df.clean$var %in% c("zoopcount","propdaphnia","propcalanoid",
                                                          "propcyclopoid", "propothercladoc",
                                                          "zoopmass",
                                                       "new.zoop.count",
                                                       "propdino"))]<-"2-zoop"

wilcox.df.clean$group[which(wilcox.df.clean$var %in% c("totphos","totnitro","totdissphos","totdissnitro",
                                                          "totdoc", "TDN.TDP", "TN.TP", 
                                                       "C.TDN", "C.TDP"))]<-"3-chem"

wilcox.df.clean$group[which(wilcox.df.clean$var %in% c("watertemp"))]<-"4-phys"

#rename vars
wilcox.df.clean$var[which(wilcox.df.clean$var=="chla")]<-"chl a"
wilcox.df.clean$var[which(wilcox.df.clean$var=="totdissnitro")]<-"TDN"
wilcox.df.clean$var[which(wilcox.df.clean$var=="totdissphos")]<-"TDP"
wilcox.df.clean$var[which(wilcox.df.clean$var=="totdoc")]<-"DOC"
wilcox.df.clean$var[which(wilcox.df.clean$var=="totnitro")]<-"TN"
wilcox.df.clean$var[which(wilcox.df.clean$var=="totphos")]<-"TP"

wilcox.df.clean <- wilcox.df.clean %>% select(group, var, n.ice, n.paired, mean.iceon, 
                                              mean.iceoff, median.iceon, median.iceoff,
                                              se.iceon, se.iceoff, mean.on.over.off.paired, 
                                              median.on.over.off.paired, 
                                              mean.on.minus.off.paired, median.on.minus.off.paired,
                                              pval, stat, null, alternative, method) %>% 
                                       arrange(group, var)

#round to 3 significant digits for output
options(scipen=999)

wilcox.df.clean$mean.iceon <- signif(wilcox.df.clean$mean.iceon, digits=3)
wilcox.df.clean$mean.iceoff <- signif(wilcox.df.clean$mean.iceoff, digits=3)

wilcox.df.clean$median.iceon <- signif(wilcox.df.clean$median.iceon, digits=3)
wilcox.df.clean$median.iceoff <- signif(wilcox.df.clean$median.iceoff, digits=3)

wilcox.df.clean$se.iceon <- signif(wilcox.df.clean$se.iceon, digits=3)
wilcox.df.clean$se.iceoff <- signif(wilcox.df.clean$se.iceoff, digits=3)

wilcox.df.clean$mean.on.over.off.paired <- signif(wilcox.df.clean$mean.on.over.off.paired, digits = 3)
wilcox.df.clean$median.on.over.off.paired <- signif(wilcox.df.clean$median.on.over.off.paired, digits = 3)

wilcox.df.clean$mean.on.minus.off.paired <- signif(wilcox.df.clean$mean.on.minus.off.paired, digits = 3)
wilcox.df.clean$median.on.minus.off.paired <- signif(wilcox.df.clean$median.on.minus.off.paired, digits = 3)

wilcox.df.clean$pval <- signif(wilcox.df.clean$pval, digits = 3)

#write to .csv
write.csv(wilcox.df.clean, "Outputs/Files/variable_summaries.csv", row.names = FALSE)


#smaller dataframe with main variables of interest
wilcox.small <- wilcox.df.clean %>% 
                filter(var %in% c("chl a", "phytomass", "new.zoop.count",
                                  "DOC", "TN", "TP", "TN.TP", "TDN",
                                  "TDP", "TDN.TDP", "watertemp"))

write.csv(wilcox.small, "Outputs/Files/variable_summaries_small.csv", row.names = FALSE)


###########################################################################
############################# Boxplots ####################################
###########################################################################

#format data for use in boxplots
data.bp<-select(data,lakename,stationname,stationlat,variable,value.iceon,value.iceoff)
data.bp<-melt(data.bp, id.var = c("lakename", "stationname", "stationlat","variable"))
names(data.bp)[5]<-"season"

data.bp$season<-as.character(data.bp$season)
data.bp$season[which(data.bp$season=="value.iceon")]<-"Winter"
data.bp$season[which(data.bp$season=="value.iceoff")]<-"Summer"
data.bp$season<-factor(data.bp$season,levels=c("Winter","Summer"))

#remove carbon variables
data.bp<-data.bp[which(!data.bp$variable %in% c("C.TDN","C.")),]


#boxplot code as a function
plot.bp<-function(dataplot,label){
  plot.bp<-ggplot(dataplot,aes(factor(season),log10(value+0.01),colour=abs(stationlat)))+
    geom_boxplot(outlier.size=0.5)+
    scale_y_continuous(breaks=c(-5:5),labels=10^c(-5:5))+
    geom_point(size = 0.1)+
    theme_bw()+
    theme(strip.text.x=element_text(size=10))+
    scale_colour_gradient(colour.lat.lims,low="green",high="blue",guide=FALSE)+
    geom_jitter(position = position_jitter(height=0,width = 0.1))+
    #ylab("value")+
    facet_wrap(~ variable,scales="free")+
    xlab("")+
    ylab("")
}


#individual variables, run through boxplot plot function

dataplotTP<-subset(data.bp,data.bp$variable=="totphos")
dataplotTP$variable<-"TP"
plot.bpTP<-plot.bp(dataplotTP,label=dataplotTP$variable[1])

dataplotTN<-subset(data.bp,data.bp$variable=="totnitro")
dataplotTN$variable<-"TN"
plot.bpTN<-plot.bp(dataplotTN,label=dataplotTN$variable[1])

dataplotTDP<-subset(data.bp,data.bp$variable=="totdissphos")
dataplotTDP$variable<-"TDP"
plot.bpTDP<-plot.bp(dataplotTDP,label=dataplotTDP$variable[1])

dataplotTDN<-subset(data.bp,data.bp$variable=="totdissnitro")
dataplotTDN$variable<-"TDN"
plot.bpTDN<-plot.bp(dataplotTDN,label=dataplotTDN$variable[1])

dataplotTNTP<-subset(data.bp,data.bp$variable=="TN.TP")
dataplotTNTP$variable<-"TN:TP ratio"
plot.bpTNTP<-plot.bp(dataplotTNTP,label=dataplotTNTP$variable[1])

dataplotTDNTDP<-subset(data.bp,data.bp$variable=="TDN.TDP")
dataplotTDNTDP$variable<-"TDN:TDP ratio"
plot.bpTDNTDP<-plot.bp(dataplotTDNTDP,label=dataplotTDNTDP$variable[1])

dataplotchla<-subset(data.bp,data.bp$variable=="chla")
dataplotchla$variable<-"Chl a"
plot.bpchla<-plot.bp(dataplotchla,label=dataplotchla$variable[1])

dataplotDOC<-subset(data.bp,data.bp$variable=="totdoc")
dataplotDOC$variable<-"DOC"
plot.bpDOC<-plot.bp(dataplotDOC,label=dataplotDOC$variable[1])

dataplotphyto<-subset(data.bp,data.bp$variable=="phytomass")
dataplotphyto$variable<-"Phyto biovolume"
plot.bpphyto<-plot.bp(dataplotphyto,label=dataplotphyto$variable[1])

dataplotzoop<-subset(data.bp,data.bp$variable=="new.zoop.count")
dataplotzoop$variable<-"Crustacean zoop"
plot.bpzoop<-plot.bp(dataplotzoop,label=dataplotzoop$variable[1])


####################################################################
#################### Layout final figure ###########################
####################################################################

#laying out scatterplots and boxplots for final image - 6 variables
png("Outputs/Figures/seasonal_variable_boxplot_scatter.png",height=11,width=13,units="in",res=1200)
grid.newpage()
#3 rows, 6 cols (scatters over two cols)
pushViewport(viewport(layout = grid.layout(3, 6)))

print(plot.bpchla,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(scat.plotchla,vp=viewport(layout.pos.row = 1, layout.pos.col = 2:3))

print(plot.bpzoop,vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(scat.plotzoop,vp=viewport(layout.pos.row = 1, layout.pos.col = 5:6))

print(plot.bpTP,vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(scat.plotTP,vp=viewport(layout.pos.row = 2, layout.pos.col = 2:3))

print(plot.bpTDN,vp=viewport(layout.pos.row = 2, layout.pos.col = 4))
print(scat.plotTDN,vp=viewport(layout.pos.row = 2, layout.pos.col = 5:6))

print(plot.bpDOC,vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(scat.plotDOC,vp=viewport(layout.pos.row = 3, layout.pos.col = 2:3))

print(plot.bpphyto,vp=viewport(layout.pos.row = 3, layout.pos.col = 4))
print(scat.plotphyto,vp=viewport(layout.pos.row = 3, layout.pos.col = 5:6))

dev.off()
