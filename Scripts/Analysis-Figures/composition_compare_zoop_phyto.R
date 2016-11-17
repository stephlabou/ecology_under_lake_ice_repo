###################################################################################
# This script uses permutational analysis of variance (PERMANOVA) to test         #
# significance of winter vs summer zooplankton and phytoplankton communities      #
# (based on proportional data). This script also creates Fig. 3 in                #
# Hampton et al. 2016.                                                            #
#                                                                                 #
# Note that the zooplankon data used in PERMANOVA is non-rotifer, non-other zoop  #
# proportional composition data, extrapolated from original proportions and total #
# ooplankton density provided by researchers.                                     #
###################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
################################################################################### 

################### Stephanie Hampton and Stephanie Labou #########################


######## clear workspace ########
rm(list = ls()) 
graphics.off()


######## load packages #########
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(grid)
library(gridExtra)
library(vegan)

###################################################################
#################### load and process data ########################
###################################################################

#### read in data (long form - already pooled aggregates) ####
dat <- read.csv("Data/data.long.csv", stringsAsFactors = FALSE)

head(dat)

#### only interested in zoop & phyto info ####
dat_zoop_phyto <- filter(dat, variable %in% c( 
                                        #zoop proportions
                                        "propotherzoop", "propothercladoc", "proprotifer", 
                                        "propdaphnia", "propcyclopoid", "propcalanoid", "zoopcount",
                                         #phyto proportions
                                        "propchloro", "propcrypto", "propcyano", 
                                        "propdiatom", "propdino", "propotherphyto", "phytocount")) %>% 
  select(lakename, stationname, year, season, variable, stat, 
         value) %>% 
  filter(stat == "ave") %>% 
  select(-stat)

#Find actual zoop counts (i.e. prop*zoopcount) after removing rotifers and "other zoop" from zoops
#Rationale: otifer zero counts potentially confounded due to variable net size 
#(i.e. net too large to capture rotifer count accurately)
#and otherzoop is potentially too broad a category to make meaningful interpretations.

zoop_phyto <- dcast(dat_zoop_phyto, lakename + stationname + year + season ~ variable, value.var = "value")

#Extrapolated density calcs

zoop_phyto_counts <- zoop_phyto %>% 
  #calculate taxa zoop counts (excluding rotifers from later analyses) as #organisms/L
  mutate(n.calanoid = propcalanoid*zoopcount,
         n.cyclopoid= propcyclopoid*zoopcount,
         n.daphnia = propdaphnia*zoopcount,
         n.rotifer = proprotifer*zoopcount,
         n.othercladoc = propothercladoc*zoopcount,
         n.otherzoop = propotherzoop*zoopcount,
    
         #change NAs to zero for rotifer and otherzoop
         #so subtraction works
         n.rotifer = ifelse(is.na(n.rotifer), 0, n.rotifer),
         n.otherzoop = ifelse(is.na(n.otherzoop), 0, n.otherzoop),
    
    #new proportions for zoops (proportions based on non-rotifer zoops)
         new.total.count = zoopcount - n.rotifer - n.otherzoop,
         prop.calanoid = n.calanoid/new.total.count,
         prop.cyclopoid = n.cyclopoid/new.total.count,
         prop.daphnia = n.daphnia/new.total.count,
         prop.othercladoc = n.othercladoc/new.total.count) %>% 
  
    #final result: new proportions for zoop, original proportions for phyto
    select(lakename, stationname, year, season,
           prop.calanoid, prop.cyclopoid, prop.daphnia, prop.othercladoc,
           propchloro, propcrypto, propcyano, propdiatom, propdino, propotherphyto)

#Melt back to desired shape
zoop_phyto_counts_long <- melt(zoop_phyto_counts, id.vars = c("lakename", "stationname", "year", "season"))

#Find seasonal averages
dat.ave <- zoop_phyto_counts_long %>% 
  #find average winter and average summer for lake/station/variable
  group_by(lakename, stationname, season, variable) %>% 
  summarize(avg_value = mean(value, na.rm = TRUE)) %>% 
  as.data.frame()


####################################################################
######### Make stacked bar graphs for zoops and phytos #############
####################################################################

#set colors for bar graphs - all same chroma and luminance
zoop_col <- rainbow_hcl(6, c=80, l=80, start = 0, end = 100)
phyto_col <- rainbow_hcl(6, c=80, l=80, start = 150, end = 300)

phyto_col2 <- phyto_col[c(1, 4, 6, 3, 5, 2)]
zoop_col2 <- zoop_col[c(1, 4, 6, 3, 5, 2)]

#### phyto data ####
dat.phyto <- subset(dat.ave, dat.ave$variable %in% c("propchloro", "propcrypto", "propcyano", 
                                                     "propdiatom", "propdino", "propotherphyto"))

dat.phyto.bar <- rename(dat.phyto, Group = variable)

dat.phyto.bar$season <- factor(dat.phyto.bar$season, levels = c("iceon", "iceoff"))

dat.phyto.bar$Group <- factor(dat.phyto.bar$Group, levels = c("propchloro", "propcrypto", "propcyano",
                                                              "propdiatom", "propdino", "propotherphyto"))

#summarize by group, season
data.phyto.plot <- dat.phyto.bar %>% 
  group_by(Group,season) %>%
  dplyr::summarize(value=mean(avg_value,na.rm=TRUE)) %>% 
  as.data.frame()

#make stacked bar chart
dataplot1<-data.phyto.plot

names(phyto_col2)<-rev(levels(dataplot1$Group))

phyto.prop.plot<-ggplot(dataplot1, aes(x=season, y=value, width=0.9,fill=Group))+
  geom_bar(stat="identity",width=0.75)+
  theme_bw()+
  ggtitle(quote(Phytoplankton))+
  xlab("")+
  scale_fill_manual(values=phyto_col2,labels=c("Chlorophyta","Cryptophyta","Cyano","Diatom","Dinoflag","Other"),
                    guide = guide_legend(reverse=TRUE))+
  scale_x_discrete(breaks=c("iceon", "iceoff"),
                   labels=c("Winter", "Summer")) +
  ylab("Proportion")+
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(vjust=1),
        panel.grid.minor = element_blank(),
        text=element_text(size=18))

#phyto.prop.plot


#### zoop data ####
#non-rotifer, non-other zoop proportions
dat.zoop <- subset(dat.ave, dat.ave$variable %in% c("prop.daphnia", "prop.othercladoc", "prop.cyclopoid", 
                                                    "prop.calanoid"))

dat.zoop.bar <- rename(dat.zoop, Group = variable)

dat.zoop.bar$season <- factor(dat.zoop.bar$season, levels = c("iceon", "iceoff"))

dat.zoop.bar$Group <- factor(dat.zoop.bar$Group, levels = c("prop.calanoid","prop.cyclopoid","prop.daphnia",
                                                            "prop.othercladoc"))

#summarize by group, season              
data.zoop.plot <- dat.zoop.bar %>% 
  group_by(Group,season) %>%
  dplyr::summarize(value=mean(avg_value,na.rm=TRUE)) %>% 
  as.data.frame()


#make stacked bar chart
dataplot2<-data.zoop.plot

names(zoop_col2)<-rev(levels(dataplot2$Group))

zoop.prop.plot<-ggplot(dataplot2, aes(x=season, y=value, width=0.9, fill=Group))+
  geom_bar(stat="identity",width=0.75)+
  theme_bw()+
  ggtitle(quote(Zooplankton))+
  xlab("")+
  scale_fill_manual(values=zoop_col2,labels=c("Calanoida","Cyclopoida","Daphnia", "Other cladoc"),
                    guide = guide_legend(reverse=TRUE))+
  scale_x_discrete(breaks=c("iceon", "iceoff"),
                   labels=c("Winter", "Summer")) +
  ylab("Proportion")+
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(vjust=1),
        panel.grid.minor = element_blank(),
        text=element_text(size=18))
#zoop.prop.plot


#### combine plots and print to png ####
png("Outputs/Figures/phyto_zoop_barplot.png",height=4,width=10,units="in",res=480)
grid.arrange(phyto.prop.plot, zoop.prop.plot, ncol = 2)
dev.off()


########################################################
#################### PERMANOVA #########################
########################################################

#reshape to "wide" format
dat.wide <- dcast(dat.ave, season + lakename + stationname ~ variable, value.var = "avg_value")

#### phyto data ####

dat.phyto.prop <- dat.wide %>% 
  select(lakename, stationname, season, #metadata
         propchloro, propcrypto, propcyano, propdiatom, propdino, propotherphyto) #phyto

dat.phyto.prop.filt <- na.omit(dat.phyto.prop) #42 obvs

dat.phyto.prop.filt <- mutate(dat.phyto.prop.filt, lakesta = paste(lakename, stationname, sep = "_"))

#want only matching lakes (i.e. matching winter and summer rows)
winter.phyto.lakes <- select(dat.phyto.prop.filt, season, lakesta) %>% 
                          filter(season == "iceon")
winter.phyto.lakes <- winter.phyto.lakes$lakesta

dat.phyto.pairs <- filter(dat.phyto.prop.filt, lakesta %in% winter.phyto.lakes) #34 obvs (17 pairs)

#### zoop data ####

dat.zoop.prop <- dat.wide %>% 
  select(lakename, stationname, season, #metadata
         prop.daphnia, prop.othercladoc, prop.cyclopoid, prop.calanoid) #zoop

dat.zoop.prop.filt <- na.omit(dat.zoop.prop) #72 obvs

dat.zoop.prop.filt <- mutate(dat.zoop.prop.filt, lakesta = paste(lakename, stationname, sep = "_"))

#want only pairs
winter.zoop.lakes <- select(dat.zoop.prop.filt, season, lakesta) %>% filter(season == "iceon")
winter.zoop.lakes <- winter.zoop.lakes$lakesta

dat.zoop.pairs <- filter(dat.zoop.prop.filt, lakesta %in% winter.zoop.lakes) #36 pairs (72 obvs)

#### data summaries ####

#using only complete cases (no NAs) for PERMANVOA analysis

## phyto data ##

phyto <- dat.phyto.pairs %>% 
  select(propchloro, propcrypto, propcyano, propdiatom, propdino, propotherphyto)

head(phyto)
summary(phyto)

phyto.env <- dat.phyto.pairs %>% 
  select(lakename, stationname, season)

head(phyto.env)
summary(phyto.env)

## zoop data ##

zoop <- dat.zoop.pairs %>% 
  select(prop.daphnia, prop.othercladoc, prop.cyclopoid, prop.calanoid)

head(zoop)
summary(zoop)

zoop.env <- dat.zoop.pairs %>% 
  select(lakename, stationname, season)

head(zoop.env)
summary(zoop.env)


## PERMANOVA analysis using the vegan package ##

#default method is bray
phyto.adonis<-adonis(phyto~season, data=phyto.env, permutations = 99, na.rm=TRUE)
phyto.adonis
summary(phyto.adonis)

#default method is bray
zoop.adonis<-adonis(zoop~season, data=zoop.env, permutations = 99, na.rm=TRUE)
zoop.adonis
summary(zoop.adonis)
