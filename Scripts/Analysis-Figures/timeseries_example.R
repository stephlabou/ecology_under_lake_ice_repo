#################################################################
# This script makes a plot of example time series               #
# (Fig. 4 in Hampton et al. 2016):                              #
# 1) first order autoregressive                                 #
# 2) first order autoregressive and moving average              #
# 3) seasonal difference                                        #
# 4) seasonal difference with moving average                    #
# 5) seasonal difference with first order autocorrelation       #
#    and moving average                                         #
#################################################################

################## Code by Steve Powers #########################

#### clear workspace ####
rm(list = ls()) 
graphics.off()

#### load packages ####
library("dplyr")
library("ggplot2")
library("reshape2")
library("scales")
library("grid")
library("gridExtra")

###############################################################
################## Load and process data ######################
###############################################################

#read in data (data.longTS is created by focal.longtimeSeries.R)
data.longTS <- read.csv("Data/data.longTS.csv",sep=",")

#prepare data for subsetting
#aggregate by lake and variable and season
#calculate nyears and nvalues, also calculate stdev of the values
countrecord.data.longTS<- data.longTS %>% 
  group_by(lakename, stationname, variable,iceonTF=season=="iceon",is.na=is.na(value)) %>% 
  dplyr::summarize(nyear=length(unique(substr(year,1,4))),nval = length(value),stdev=sd(value,na.rm = TRUE))

#subset the aggregated data
#use stdev!=0 to get rid of a few time series that have non-varying (uniform) iceon values, 
#which is nonsensical and causes auto.arima function to fail
countrecord.data.longTS<-subset(countrecord.data.longTS,countrecord.data.longTS$iceonTF==TRUE & countrecord.data.longTS$is.na==FALSE 
                                & countrecord.data.longTS$stdev!=0)

#subset to time series that have a minimum number of values
#note that time series with 25 or fewer data points are prone to causing failure of the auto-arima function 
use.lakestavar<-countrecord.data.longTS[which(countrecord.data.longTS$nyear>=15 & countrecord.data.longTS$nval>=10),]
use.lakestavar<-subset(use.lakestavar,!use.lakestavar$variable %in% c("icedepth","snowdepth"))

#get unique combinations for lakename, stationname, and variable that meet the above criteria
unique.lakestavar<-unique(data.frame(lakename=use.lakestavar$lakename,stationname=use.lakestavar$stationname,variable=use.lakestavar$variable))

#subset data to lakename, stationname, and variable combos that are "good" for timeseries 
#(and used for making example time series figs)
long.df.use<-subset(data.longTS,data.longTS$lakename %in% unique.lakestavar$lakename & 
                      data.longTS$stationname %in% unique.lakestavar$stationname & 
                      data.longTS$variable %in% unique.lakestavar$variable)

################################################################
################# Plot time series examples ####################
################################################################

############### Buffalo Pound Lake, DOC #################

#example of first order autoregressive structure

lakei<-"Buffalo Pound Lake"
stationi<-"OfftakePipe"
variablei<-"totdoc"

datai<-subset(long.df.use,long.df.use$lakename==lakei & long.df.use$stationname==stationi & long.df.use$variable==variablei)
datai<-datai[-which(is.na(datai$value)),]

plot<-ggplot(datai,aes(x=year,y=value,colour=season))+
  geom_point(size=2)+ ylab(quote(DOC~(mg~L^-1)))+xlab("Year") +
  scale_color_discrete(name = "Season", breaks = c("iceoff", "iceon"), labels = c("Summer", "Winter"))+
  theme_bw()+
  theme(panel.border=element_blank(),
        axis.line.x=element_line(color = "black"),
        axis.line.y=element_line(color = "black"),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank())+
  ggtitle("Buffalo Pound Lake") #AR+MA

one.plot<-plot
one.plot

############# Big Muskellunge Lake, chl a #############

#example of seasonal difference

lakei<-"Big Muskellunge Lake"
stationi<-"unnamed"
variablei<-"chla"

datai<-subset(long.df.use,long.df.use$lakename==lakei & long.df.use$stationname==stationi & long.df.use$variable==variablei)
datai<-datai[-which(is.na(datai$value)),]

plot<-ggplot(datai,aes(x=year,y=value,colour=season))+
  geom_point(size=2)+ylab(quote(Chl~a~(ug~L^-1)))+xlab("Year")+
  scale_color_discrete(name = "Season", breaks = c("iceoff", "iceon"), labels = c("Summer", "Winter"))+
  theme_bw()+
  theme(panel.border=element_blank(),
        axis.line.x=element_line(color = "black"),
        axis.line.y=element_line(color = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ggtitle("Big Muskellunge Lake") #iceon

two.plot<-plot
two.plot

############## Allequash Lake, TP ##############

#example of seasonal difference with moving average

lakei<-"Allequash Lake"
stationi<-"unnamed"
variablei<-"totphos"

datai<-subset(long.df.use,long.df.use$lakename==lakei & long.df.use$stationname==stationi & long.df.use$variable==variablei)
datai<-datai[-which(is.na(datai$value)),]

plot<-ggplot(datai,aes(x=year,y=value,colour=season))+
  geom_point(size=2)+ylab(quote(TP~(ug~L^-1)))+xlab("Year")+
  scale_color_discrete(name = "Season", breaks = c("iceoff", "iceon"), labels = c("Summer", "Winter"))+
  theme_bw()+
  theme(panel.border=element_blank(),
        axis.line.x=element_line(color = "black"),
        axis.line.y=element_line(color = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ggtitle("Allequash Lake") #trend+iceon

three.plot<-plot
three.plot

####### Lake Superior (Thunder Bay), TN #######

#example of seasonal difference with first order autocorrelation structure and moving average

lakei<-"Lake Superior"
stationi<-"ThunderBay"
variablei<-"totnitro"

datai<-subset(long.df.use,long.df.use$lakename==lakei & long.df.use$stationname==stationi & long.df.use$variable==variablei)
datai<-datai[-which(is.na(datai$value)),]

plot<-ggplot(datai,aes(x=year,y=value,colour=season))+
  geom_point(size=2)+ylab(quote(TN~(ug~L^-1)))+xlab("Year")+
  scale_color_discrete(name = "Season", breaks = c("iceoff", "iceon"), labels = c("Summer", "Winter"))+
  theme_bw()+
  theme(panel.border=element_blank(),
        axis.line.x=element_line(color = "black"),
        axis.line.y=element_line(color = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ggtitle("Lake Superior at Thunder Bay") #AR+MA+iceon

four.plot<-plot
four.plot

############### Sparkling Lake, suva ################

#example of first order autocorrelation structure

lakei<-"Sparkling Lake"
stationi<-"unnamed"
variablei<-"suva"

datai<-subset(long.df.use,long.df.use$lakename==lakei & long.df.use$stationname==stationi & long.df.use$variable==variablei)
datai<-datai[-which(is.na(datai$value)),]

plot<-ggplot(datai,aes(x=year,y=value,colour=season))+
geom_point(size=2)+ylab(quote(SUVA~(L~mgC^-1~m^-1)))+xlab("Year")+
  scale_color_discrete(name = "Season", breaks = c("iceoff", "iceon"), labels = c("Summer", "Winter"))+
  theme_bw()+
  theme(panel.border=element_blank(),
        axis.line.x=element_line(color = "black"),
        axis.line.y=element_line(color = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ggtitle("Sparkling Lake") #AR

five.plot<-plot
five.plot

########### arrange plots and print out ##############

#format so axes line up
g1.plot <- ggplotGrob(one.plot)
g2.plot <- ggplotGrob(two.plot)
g3.plot <- ggplotGrob(three.plot)
g4.plot <- ggplotGrob(four.plot)
g5.plot <- ggplotGrob(five.plot)

#print to png
png(filename="Outputs/Figures/example_timeseries.png",units="in",res=180,width=8,height=11)
grid::grid.draw(rbind(g5.plot, g1.plot, g2.plot, g3.plot, g4.plot, size = "last"))
dev.off()

