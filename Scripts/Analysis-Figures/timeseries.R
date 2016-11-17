#################################################################################################                               
# The code helps to answer these questions...                                                   #
# A) How many lakes show winter-summer differences, variable-by-variable? In which direction?   #
# B) How many lakes show seasonal connections,                                                  #
#    as expressed by the AR1 coef, or residuals in the detrended/deseasoned data?               #
#################################################################################################
#This script does the following:                                                                #
#0. Set directories, load packages, load data                                                   #
#   (This script requires "data.longTS.csv", which is created by focal.long.timeSeries.R)       #
#1. Identify "good" time series                                                                 #
#2. Do the analysis...                                                                          #
#     -Detrended the data (moving average)                                                      #
#     -Determine sign of seasonal difference, if present, using AIC                             #
#     -Determine sign of seasonal correlations (AB,BC), if present, using AIC                   #
#     -Redo fits for seasonal differences and seasonal correlations when residuals are          #
#       autocorrelated, adding AR1 error structure                                              #
#     -Save diagnostic plots (incl. Figures S5 and S6 in Hampton et al. 2016)                   #
#3. Stack results outputs into data frame (Table S5 in Hampton et al. 2016)                     #
#4. Summarize and aggregate the results (Table 2 in Hampton et al. 2016)                        #
#5. Send model summaries to file(s)                                                             #
################################################################################################# 

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
################################################################################### 

########################## Code by Steve Powers ###################################

######## clear workspace ########
rm(list = ls()) 
graphics.off()

############ load packages ##########
library(zoo) 
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(AICcmodavg)
library(nlme)

#########################################
###### Load Data and Define Values ######
#########################################

#This .csv file is created by the script focal.long.timeSeries.R
data.longTS <- read.csv("Data/data.longTS.csv",sep=",")

############ functions and values ##############

#set AIC thresholds
AIC_thresh<-2
AIC_thresh_connect<-2
ma.use<-7

#AIC function - calculate AIC manually
AICc_manual<-function(k,RSS,n){
  AICc<-n*log(RSS/n)+2*k+(2*k*(k+1))/(n-k-1)
  return(AICc)
}

#consecutive NA's function
consecna <- function(x, n=4) {
  # function to identify elements with n or more consecutive NA values
  # returns TRUE if consecutive NA (runs of 3 consecutive NA or more)
  y <- rle(is.na(x))
  y$values <- y$lengths > (n - 0.5) & y$values
  inverse.rle(y)
}

############## ready data ################

#add nonrotifer, non-other zoop new zoopcount

data.zoop1<-subset(data.longTS,data.longTS$variable=="zoopcount")
data.zoop2<-subset(data.longTS,data.longTS$variable=="proprotifer")
data.zoop3<-subset(data.longTS,data.longTS$variable=="propotherzoop")

data.zoop.merge<-merge(data.zoop1,data.zoop2,by=c("lakename","stationname","year","season","stationlat","stationlong"),all=TRUE)
data.zoop.merge<-merge(data.zoop.merge,data.zoop3,by=c("lakename","stationname","year","season","stationlat","stationlong"
),all=TRUE)

data.zoop.merge$value.y[which(is.na(data.zoop.merge$value.y)==TRUE)]<-0
data.zoop.merge$value[which(is.na(data.zoop.merge$value)==TRUE)]<-0
data.zoop.merge<-data.zoop.merge[-which(is.na(data.zoop.merge$value.x)==TRUE),]
data.zoop.merge$zoopcount_nonrotifer<-data.zoop.merge$value.x-(data.zoop.merge$value.x*data.zoop.merge$value.y)-(data.zoop.merge$value.x*data.zoop.merge$value)
data.zoop.merge$variable.x<-"zoopcount_nonrotifer"

data.zoop.nonrotifer<-data.zoop.merge %>% 
                        select(lakename,stationname,stationlat,stationlong,year,
                               season,variable=variable.x,value=zoopcount_nonrotifer)

merge2<-merge(data.zoop1,data.zoop.nonrotifer,by=c("lakename","stationname","year","season","stationlat","stationlong"),ALL=TRUE)

data.zoop.nonrotifer<-merge2 %>% select(lakename, stationname,   year,  variable=variable.y, season, 
                                        poolstation, stat, sampletype, stationlat, 
                                        stationlong,researcher, lakeregloc, lakecountry, lakearea, 
                                        lakemeandepth, lakemaxdepth, lakeelevation, watershedarea, 
                                        h2oresidence, lakefetch, stationdistance, stationdepth, multiplestations, 
                                        startday, startmonth, startyear, endday, endmonth, endyear, 
                                        iceduration, periodn, nYears, n.season, nFullYears, old.variable, 
                                        value.before, value=value.y, value.after, stationname2, n.stationname2)

data.longTS<-rbind(data.longTS,data.zoop.nonrotifer)

######################################################
#### Identify the "good" time series for modeling ####
######################################################

#prepare data for subsetting
#aggregate by lake and variable and season
#calculate nyears and nvalues, also calculate stdev of the values
countrecord.data.longTS<- data.longTS %>% 
  group_by(lakename, stationname, variable,iceonTF=season=="iceon",is.na=is.na(value)) %>% 
  dplyr::summarize(nyear=length(unique(substr(year,1,4))),nval = length(value),stdev=sd(value,na.rm = TRUE))

#subset the aggregated data
#use stdev!=0 to get rid of a few time series that have non-varying (uniform) iceon values
countrecord.data.longTS<-subset(countrecord.data.longTS,countrecord.data.longTS$iceonTF==TRUE & 
                                  countrecord.data.longTS$is.na==FALSE & 
                                  countrecord.data.longTS$stdev!=0)

#subset to time series that have a minimum number of values
use.lakestavar<-countrecord.data.longTS[which(countrecord.data.longTS$nval>=10),]


#drop sampledepth - not needed here
use.lakestavar<-subset(use.lakestavar,!use.lakestavar$variable %in% c("icedepth","snowdepth","sampledepth"))#,"airtemp","radiation"))

#get unique combinations for lakename, stationname, and variable that meet the above criteria
unique.lakestavar<-unique(data.frame(lakename=use.lakestavar$lakename,stationname=use.lakestavar$stationname,variable=use.lakestavar$variable))

#subset data to lakename, stationname, and variable combos that are "good" for time series
data.TS.run<-subset(data.longTS,data.longTS$lakename %in% unique.lakestavar$lakename & 
                      data.longTS$stationname %in% unique.lakestavar$stationname & 
                      data.longTS$variable %in% unique.lakestavar$variable)

#######################################################################
################ Iterate and analyze each time series #################
#######################################################################

######################## start for loop ###############################

#loops through each lake/station/variable unique combination

for (i in 1:length(unique.lakestavar[,1])){
  
  print(i)
  lakestavari<-unique.lakestavar[i,]
  
  ########## get the ith time series ############
  datai<-subset(data.TS.run,
                data.TS.run$lakename==lakestavari$lakename & 
                  data.TS.run$stationname==lakestavari$stationname & 
                  data.TS.run$variable==lakestavari$variable)
  
  #replace zoop proportions with zoop counts
  data.zoopcount<-subset(data.TS.run,
                         data.TS.run$lakename==lakestavari$lakename & 
                           data.TS.run$stationname==lakestavari$stationname & data.TS.run$variable=="zoopcount") %>% 
    select(lakename, stationname, year, season, variable, value) %>% 
    rename(value.zoop = value)
  
  if(lakestavari$variable %in% c("propcalanoid","propcyclopoid","propdaphnia",
                                 "propothercladoc","propotherzoop","proprotifer")) {
    
    if(length(data.zoopcount[,1])==0){next()}
    
    mergei<-merge(datai,data.zoopcount,by=c("lakename", "stationname", "year","season"),all = FALSE)
    #calculate count by proportion * zoopcount
    valuei<-mergei$value*mergei$value.zoop
    datai$value<-valuei
  }
  
  
  #find the longest contiguous run of years with non-NA values
  yr_run <- na.contiguous(data.frame(year=datai$year,value=datai$value))
  
  # subset to years NOT in contiguous run
  # looking if there is a second not-quite-as-long contiguous run
  # (this takes care not to split long time series and throw away half if there is one missing in the middle)
  data2i<-subset(datai,!datai$year %in% yr_run$year)
  
  #if there is no second longest run, then the final year run is just the original longest year run
  if (length(data2i$year) == 0) {
    yr_run12<-yr_run$year
  }
  
  #if second longest is all NAs, use only first longest run
  if (length(data2i$year) != 0 & sum(is.na(data2i$value)==FALSE)==0) {
    yr_run12<-yr_run$year
  }
  
  #find second longest run, if it exists (for non NA runs)
  if (length(data2i$year) != 0 & sum(is.na(data2i$value)==FALSE)>0) {
    yr_run2 <- na.contiguous(data.frame(year=data2i$year,value=data2i$value))
    yr_run12<-sort(c(yr_run$year,yr_run2$year))
  }
  
  #subset data to longest run(s)
  datai_in<- subset(datai, datai$year %in% yr_run12)
  datai<-datai_in
  
  if (length(datai$year) < 10 | sum(na.omit(datai$value))==0) 
  {next()}
  
  #drop rows with value=NA!
  #datai<-datai[which(is.na(datai$value)==FALSE),]
  
  #checks for any NAs
  na.count<-sum(is.na(datai$value))
  
  #find consecutive NA's (4) and drop those rows, only applies to a few time series (i.e. Lake Peipsi)
  consec.nai<-consecna(datai$value,n=4)
  if(sum(consec.nai==TRUE)>0){
    datai<-datai[which(consec.nai==FALSE),]
  }
  
  # find total number of observations
  sub<-datai[which(is.na(datai$value)==FALSE),]
  n.obvs<-length(sub$value)
  
  #checks for any NAs - after removing consecutive NAs
  na.count2<-sum(is.na(datai$value))
  
  # find number of iceon observations
  sub1<-subset(datai,datai$season=="iceon")
  sub2<-sub1[which(is.na(sub1$value)==FALSE),]
  n.ice <- length(sub2$value)
  lakecountry <- datai$lakecountry[1]
  lakeregloc <- datai$lakeregloc[1]
  stationlat<-datai$stationlat[1]
  stationlong<-datai$stationlong[1]
  
  
  ########### detrending and deseasoning ###########
  ty<-datai$value 
  
  # for time series with huge values, convert units
  if(max(na.omit(ty))>10^6){ty<-ty/10^6}
  
  # compute moving average (set above, near AIC threshold)
  ty.ma<-rollapply(ty,ma.use,mean,na.rm=TRUE)
  
  # since moving average cannot be calculated for edges of time series, 
  # next few lines handle the edges (rather than discard edge values)
  # for the edges we substitute the simple average of the first few points, 
  # and simple average of last few points
  
  if(ma.use==3){  
    ty.first <- mean(ty[1:2])
    ty.last <- mean(ty[(length(ty)-1):length(ty)])
    ty.ma <- c(ty.first, ty.ma, ty.last)
  }
  if(ma.use==5){  
    ty.first <- mean(ty[1:4])
    ty.last <- mean(ty[(length(ty)-3):length(ty)])
    ty.ma <- c(ty.first, ty.first, ty.ma, ty.last, ty.last)
  }
  if(ma.use==7){  
    ty.first <- mean(ty[1:5])
    ty.last <- mean(ty[(length(ty)-4):(length(ty))])
    ty.ma <- c(ty.first, ty.first, ty.first, ty.ma, ty.last, ty.last, ty.last)
  }
  
  #add moving average vectors
  datai$value_ma<-ty.ma
  
  #detrend (value - moving average)
  datai$value_detrended<-datai$value-datai$value_ma
  
  #compute means/medians, raw and detrended 
  mean_iceoni<- mean(subset(datai,datai$season=="iceon")$value,na.rm=TRUE)
  mean_iceoffi<- mean(subset(datai,datai$season=="iceoff")$value,na.rm=TRUE)
  
  meandiffi<-mean_iceoni-mean_iceoffi
  meanratioi<-mean_iceoni/mean_iceoffi
  
  median_iceoni<- median(subset(datai,datai$season=="iceon")$value,na.rm=TRUE)
  median_iceoffi<- median(subset(datai,datai$season=="iceoff")$value,na.rm=TRUE)
  
  mediandiffi<-median_iceoni-median_iceoffi
  medianratioi<-median_iceoni/median_iceoffi
  
  mean_iceon_detrendedi<- mean(subset(datai,datai$season=="iceon")$value_detrended,na.rm=TRUE)
  mean_iceoff_detrendedi<- mean(subset(datai,datai$season=="iceoff")$value_detrended,na.rm=TRUE)
  
  meandiff_detrendedi<-mean_iceon_detrendedi-mean_iceoff_detrendedi
  meanratio_detrendedi<-mean_iceon_detrendedi/mean_iceoff_detrendedi
  
  median_iceon_detrendedi<- median(subset(datai,datai$season=="iceon")$value_detrended,na.rm=TRUE)
  median_iceoff_detrendedi<- median(subset(datai,datai$season=="iceoff")$value_detrended,na.rm=TRUE)
  
  mediandiff_detrendedi<-median_iceon_detrendedi-median_iceoff_detrendedi
  medianratio_detrendedi<-median_iceon_detrendedi/median_iceoff_detrendedi
  
  #detrended+deseasoned, calculated even if the seasonal difference is not meaningful
  datai$value_dtds<-NA
  datai$value_dtds[which(datai$season=="iceon")]<-datai$value_detrended[which(datai$season=="iceon")]-mean_iceon_detrendedi
  datai$value_dtds[which(datai$season=="iceoff")]<-datai$value_detrended[which(datai$season=="iceoff")]-mean_iceoff_detrendedi
  
  
  #do AICc calcs for seasonal model (lmi) and intercept model (lm0i)
  #testing detrended vs simple intercept model
  lm0i<-lm(value_detrended~1,data=datai)
  lmi<-lm(value_detrended~season,data=datai)
  
  box.test.i<-Box.test(lmi$residuals,lag=1)
  box.test.i
  box.test.AR1sigTFi<-0
  if(box.test.i$p.value<=0.05){
    box.test.AR1sigTFi<-1
  }
  
  if(box.test.AR1sigTFi==1){
    lm0i<-gls(value_detrended~1,data=datai,correlation=corARMA(p=1),na.action=na.omit)
    lmi<-gls(value_detrended~season,data=datai,correlation=corARMA(p=1),na.action=na.omit)
  }
  
  #calculate AICc's, and AICc difference
  #beware that AICc function (AICcmodavg package) does unexpected things, use manual calc instead!
  AICc0i<-AICc_manual(k=1,RSS=sum(lm0i$residuals^2),n=length(lm0i$residuals))
  AICci<-AICc_manual(k=2,RSS=sum(lmi$residuals^2),n=length(lmi$residuals))
  
  #if AICc_deli > 2, there is a meaningful seasonal difference
  AICc_deli<-AICc0i-AICci
  
  #add vector for computing lag coefs (initially set to value_detrended, which is true if no seasonal difference)
  #see below: overwritten if seasonal difference
  datai$value_forlag<-datai$value_detrended
  
  # set signif.seasdiffi and sign.seasdiffi to 0
  # then re-write if there is a seasonal diff
  signif.seasdiffi<-0
  sign.seasdiffi<-0
  # if(t.test.i$p.value<=0.05){
  if(AICc_deli>=AIC_thresh){
    signif.seasdiffi<-1
    # write the sign of the seasonal diff, if diff present
    sign.seasdiffi<-sign(meandiff_detrendedi)
  }
  # if seasonal difference is present, change value_forlag to detrended+deseasoned  
  if(signif.seasdiffi==1){
    datai$value_forlag[which(datai$season=="iceon")]<-datai$value_detrended[which(datai$season=="iceon")]-mean_iceon_detrendedi
    datai$value_forlag[which(datai$season=="iceoff")]<-datai$value_detrended[which(datai$season=="iceoff")]-mean_iceoff_detrendedi
  }
  
  ############ seasonal correlation analysis ############ 
  
  #create spreaded data for residual plot, so iceon/iceoff values each get columns
  #using value_forlag (aka detrended+deseasoned when sig diff, detrended otherwise)
  datai.sub<-data.frame(year=trunc(datai$year),season=datai$season,value_forlag=datai$value_forlag)
  scatter.datai<-datai.sub %>% spread(key=season,value=value_forlag) %>% as.data.frame()
  scatter.datai.na.omit<-na.omit(scatter.datai)
  
  #using raw value
  datai.sub.raw<-data.frame(year=trunc(datai$year),season=datai$season,value=datai$value)
  scatter.datai.raw<-datai.sub.raw %>% spread(key=season,value=value) %>% as.data.frame()
  scatter.datai.raw.na.omit<-na.omit(scatter.datai.raw)
  
  #using detrended value
  datai.sub.detrend<-data.frame(year=trunc(datai$year),season=datai$season,value_detrended=datai$value_detrended)
  scatter.datai.detrend<-datai.sub.detrend %>% spread(key=season,value=value_detrended) %>% as.data.frame()
  scatter.datai.detrend.na.omit<-na.omit(scatter.datai.detrend)
  
  #create AB and BC time series for seasonal correlation analysis
  # also calculate corr coefs (relict code?)
  
  #detrended+deseasoned when sig diff
  scatter.ab.i<-scatter.datai %>% select(year, iceoff, iceon)
  scatter.bc.i<-scatter.datai %>% select(year, iceon, iceoff)
  
  scatter.ab.i$iceoff <- c(NA, scatter.ab.i$iceoff[1:(length(scatter.ab.i$iceoff)-1)])
  
  ab.cori = cor(scatter.ab.i$iceon,scatter.ab.i$iceoff,method = "spearman",use="complete.obs")
  bc.cori = cor(scatter.bc.i$iceon,scatter.bc.i$iceoff,method = "spearman", use = "complete.obs")
  
  
  #raw data
  scatter.ab.i.raw<-scatter.datai.raw %>% select(year, iceoff, iceon)
  scatter.bc.i.raw<-scatter.datai.raw %>% select(year, iceon, iceoff)
  
  scatter.ab.i.raw$iceoff <- c(NA, scatter.ab.i.raw$iceoff[1:(length(scatter.ab.i.raw$iceoff)-1)])
  
  ab.cor.rawi = cor(scatter.ab.i.raw$iceon,scatter.ab.i.raw$iceoff,method = "spearman",use="complete.obs")
  bc.cor.rawi = cor(scatter.bc.i.raw$iceon,scatter.bc.i.raw$iceoff,method = "spearman", use = "complete.obs")
  
  
  #detrended only data
  scatter.ab.i.detrend<-scatter.datai.detrend %>% select(year, iceoff, iceon)
  scatter.bc.i.detrend<-scatter.datai.detrend %>% select(year, iceon, iceoff)
  
  scatter.ab.i.detrend$iceoff <- c(NA, scatter.ab.i.detrend$iceoff[1:(length(scatter.ab.i.detrend$iceoff)-1)])
  
  ab.cor.detrendi = cor(scatter.ab.i.detrend$iceon,scatter.ab.i.detrend$iceoff,method = "spearman",use="complete.obs")
  bc.cor.detrendi = cor(scatter.bc.i.detrend$iceon,scatter.bc.i.detrend$iceoff,method = "spearman", use = "complete.obs")
  
  ######################################################
  
  # do pacf plot
  pacfi<-pacf(na.omit(datai$value_forlag))
  # collect pacf1 coef
  pacf1i<-round(as.numeric(as.character(pacfi[1])[1]),digits=3)
  # get rid of the pacf plot
  dev.off()
  
  # do AB correlation analysis
  # do linear model of detrended values
  # note this produces same slope and r2 as for detrended+deseasoned data
  lm.connect0.i<-lm(iceon~1,data=na.omit(scatter.ab.i.detrend))
  lm.connect.i<-lm(iceon~iceoff,data=na.omit(scatter.ab.i.detrend))
  
  # do box test to see if residuals are autocorrelated at AR1
  # if residuals autocorrelated at AR1, fit new model with correlation structure (corARMA(p=1))
  # note that only a minority of the time series have AR1 correlated residuals
  box.test.connect0.i<-Box.test(lm.connect0.i$residuals,lag=1)
  box.test.connect.i<-Box.test(lm.connect.i$residuals,lag=1)
  
  # set box.test.AB.AR1sigTFi=0 then change to 1 if box.test is significant
  box.test.AB.AR1sigTFi<-0
  if(box.test.connect.i$p.value<=0.05){
    box.test.AB.AR1sigTFi<-1
  }
  # run the new model if box.test is signif
  if(box.test.AB.AR1sigTFi==1){
    #lm0i<-gls(iceon~1,data=scatter.ab.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
    #lmi<-gls(iceon~iceoff,data=scatter.ab.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
    lm.connect0.i<-gls(iceon~1,data=scatter.ab.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
    lm.connect.i<-gls(iceon~iceoff,data=scatter.ab.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
  }
  
  #p-value
  pval.lm.connect.i<-NA#summary(lm.connect.i)$coefficients[2,4]
  
  #calculate AICc's, and AICc difference for above linear models
  AICc.connect.0i<-AICc_manual(k=1,RSS=sum(lm.connect0.i$residuals^2),n=length(lm.connect0.i$residuals))
  AICc.connect.i<-AICc_manual(k=2,RSS=sum(lm.connect.i$residuals^2),n=length(lm.connect.i$residuals))
  
  #if AICc_deli > 2, there is a meaningful seasonal difference
  ab.AICc_connect_deli<-AICc.connect.0i-AICc.connect.i
  lm.connect.ab.i<-lm.connect.i
  
  # now do BC correlation analysis, same procedure as AB above
  lm.connect0.i<-lm(iceoff~1,data=na.omit(scatter.bc.i.detrend))
  lm.connect.i<-lm(iceoff~iceon,data=na.omit(scatter.bc.i.detrend))
  
  box.test.connect0.i<-Box.test(lm.connect0.i$residuals,lag=1)
  box.test.connect.i<-Box.test(lm.connect.i$residuals,lag=1)
  box.test.BC.AR1sigTFi<-0
  if(box.test.connect.i$p.value<=0.05){
    box.test.BC.AR1sigTFi<-1
  }
  if(box.test.BC.AR1sigTFi==1){
    #lm0i<-gls(iceoff~1,data=scatter.bc.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
    #lmi<-gls(iceoff~iceon,data=scatter.bc.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
    lm.connect.0i<-gls(iceoff~1,data=scatter.bc.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
    lm.connect.i<-gls(iceoff~iceon,data=scatter.bc.i.detrend,correlation=corARMA(p=1),na.action=na.omit)
  }
  
  pval.lm.connect.i<-NA#summary(lm.connect.i)$coefficients[2,4]
  
  AICc.connect.0i<-AICc_manual(k=1,RSS=sum(lm.connect0.i$residuals^2),n=length(lm.connect0.i$residuals))
  AICc.connect.i<-AICc_manual(k=2,RSS=sum(lm.connect.i$residuals^2),n=length(lm.connect.i$residuals))
  
  bc.AICc_connect_deli<-AICc.connect.0i-AICc.connect.i
  lm.connect.bc.i<-lm.connect.i
  
  # create vectors indicating the presence/absence (1=present) and sign (1 or -1) of AB and BC correlations
  signif.AB.connecti<-0
  sign.AB.connecti<-0
  signif.BC.connecti<-0
  sign.BC.connecti<-0
  if(ab.AICc_connect_deli>=AIC_thresh_connect){
    signif.AB.connecti<-1
    sign.AB.connecti<-
      #sign of est. slope
      sign(as.numeric(lm.connect.ab.i$coefficients[2]))
  }
  if(bc.AICc_connect_deli>=AIC_thresh_connect){
    signif.BC.connecti<-1
    sign.BC.connecti<-
      #sign of est. slope
      sign(as.numeric(lm.connect.bc.i$coefficients[2]))
  }
  
  ############ do time series multiplot ############
  png(file=paste(c("Outputs/Figures/", i,"_",as.character(datai$variable[1]),"_",as.character(datai$lakename[1]),"_",as.character(datai$stationname[1]),".png"),collapse=""),
      width=7,height=7.5,unit="in",res=180)
  par(mfrow=c(3,3))
  plot(datai$year,datai$value,xlab="Year",ylab="value",type="b",pch=21,bg=datai$season, main = "raw values")
  plot(datai$year,datai$value_detrended,xlab="Year",ylab="value_detrended",type="b",pch=21,bg=datai$season, main = "detrended")
  plot(datai$year,datai$value_dtds,xlab="Year",ylab="value_detrend_deseas",type="b",pch=21,bg=datai$season, main = "detrended + deseasoned")
  pacf(na.omit(datai$value),main="pacf raw data")
  pacf(na.omit(datai$value_detrended),main="pacf detrended")
  pacf(na.omit(datai$value_dtds),main="pacf detrended+deseasoned")
  boxplot(datai$value_detrended~datai$season,names=c("summer","winter"),main="detrended")
  plot(scatter.ab.i.detrend$iceoff,scatter.ab.i.detrend$iceon,xlab="prev summer",ylab="winter",pch=21,main="detrended, S->W")
  plot(scatter.bc.i.detrend$iceon,scatter.bc.i.detrend$iceoff,xlab="winter",ylab="next summer",pch=21,main="detrended, W->S")
  mtext(text=paste(datai$variable[1],datai$lakename[1]),adj=0.5,side=3,outer=TRUE,padj=1.5)
  dev.off()
  
  ########### compile output vectors ############
  resultsi<-data.frame(lakei=datai$lakename[1],stationi=datai$stationname[1],variablei=datai$variable[1],
                       regionloci=datai$lakeregloc[1],stationlati=datai$stationlat[1],stationlongi=datai$stationlong[1],
                       n.icei=n.ice,n.obvsi=n.obvs,na.counti=na.count,
                       mean_iceoni, mean_iceoffi, meandiffi,meanratioi,
                       mean_iceon_detrendedi, mean_iceoff_detrendedi, meandiff_detrendedi,meanratio_detrendedi,
                       #  median_iceon_detrendedi, median_iceoff_detrendedi, mediandiff_detrendedi,medianratio_detrendedi,
                       box.test.AR1sigTFi,
                       signif.seasdiffi, sign.seasdiffi,
                       AICc_seasdiffi=AICci,AICc_noseasdiffi=AICc0i,AICc_deli,
                       ab.cor.rawi, bc.cor.rawi,
                       ab.cori, bc.cori,
                       ab.cor.detrendi, bc.cor.detrendi,
                       pacf1i, 
                       ab.AICc_connect_deli,bc.AICc_connect_deli,
                       signif.AB.connecti,sign.AB.connecti,signif.BC.connecti,sign.BC.connecti,
                       box.test.AB.AR1sigTFi,box.test.BC.AR1sigTFi)
  
  #create second results data frame for multipanel residual plots (such as multipanel AB plot for chl a, and multipanel BC plot for chl a)
  repi<-length(scatter.datai.na.omit$iceon)
  
  #using MA detrended data:
  data.abbc.i<-merge(scatter.ab.i.detrend,scatter.bc.i.detrend,by="year",suffixes=c(".abi",".bci"),all=TRUE)
  data.abbc.i<-data.abbc.i %>% select(yeari=year,iceoff.abi,iceon.abi,iceon.bci,iceoff.bci)
  repi<-length(data.abbc.i[,1])
  results.connect.pairedi<-data.frame(lakei=rep(datai$lakename[1],repi),stationi=rep(datai$stationname[1],repi),
                                      variablei=rep(datai$variable[1],repi), regionloci=rep(datai$lakeregloc[1],repi),
                                      stationlati=rep(datai$stationlat[1],repi),stationlongi=rep(datai$stationlong[1],repi),
                                      sign.seasdiffi=rep(sign.seasdiffi,repi), data.abbc.i,
                                      ab.AICc_connect_deli,bc.AICc_connect_deli,
                                      signif.AB.connecti=rep(signif.AB.connecti,repi),sign.AB.connecti=rep(sign.AB.connecti,repi),
                                      signif.BC.connecti=rep(signif.BC.connecti,repi),sign.BC.connecti=rep(sign.BC.connecti,repi))
  
  # stack results for each i, into results data frame
  if(i==1){results<-resultsi; 
  results.connect.paired<-results.connect.pairedi}
  if(i>1){results<-rbind(results,resultsi); 
  results.connect.paired<-rbind(results.connect.paired,results.connect.pairedi)}
  
}
######################### end for loop ################################


######### clean results dataframe ##########

# drop the trailing "i" from results columns
names(results)<-
  substr(names(results),1,-1+nchar(names(results)))

names(results.connect.paired)<-
  substr(names(results.connect.paired),1,-1+nchar(names(results.connect.paired)))

#convert "variable" to character
results$variable<-as.character(results$variable)

#rename variables
results$variable[which(results$variable=="chla")]<-"chl a"
results$variable[which(results$variable=="propotherphyto")]<-"other phyto"
results$variable[which(results$variable=="propcyano")]<-"cyano"
results$variable[which(results$variable=="propchloro")]<-"chloro"
results$variable[which(results$variable=="propcrypto")]<-"crypto"
results$variable[which(results$variable=="propdiatom")]<-"diatom"
results$variable[which(results$variable=="propdaphnia")]<-"daphnia"
results$variable[which(results$variable=="propcalanoid")]<-"calanoid"
results$variable[which(results$variable=="propcyclopoid")]<-"cyclopoid"
results$variable[which(results$variable=="propothercladoc")]<-"other cladoc"
results$variable[which(results$variable=="propotherzoop")]<-"other zoop"
results$variable[which(results$variable=="proprotifer")]<-"rotifer"
results$variable[which(results$variable=="totdissnitro")]<-"TDN"
results$variable[which(results$variable=="totdissphos")]<-"TDP"
results$variable[which(results$variable=="totdoc")]<-"DOC"
results$variable[which(results$variable=="totnitro")]<-"TN"
results$variable[which(results$variable=="totphos")]<-"TP"
results$variable[which(results$variable=="watertemp")]<-"water temp"

results<-results[order(results$variable,results$lake,results$station),]

results.out<-results %>% select(lake,station,variable,regionloc,stationlat,stationlong,
                                n.ice,n.obvs,
                                mean_iceon,mean_iceoff,meandiff,meanratio,
                                mean_iceon_detrended,mean_iceoff_detrended,meandiff_detrended,
                                meanratio_detrended,box.test.AR1sigTF,
                                AICc_seasdiff,AICc_noseasdiff,
                                seasdiff_aicTF=signif.seasdiff,sign_seasdiff=sign.seasdiff,AICc_del,
                                ab.AICc_connect_del,bc.AICc_connect_del,
                                sign.AB.connect,sign.BC.connect,
                                ab.cor.detrend,bc.cor.detrend,
                                box.test.AB.AR1sigTF,box.test.BC.AR1sigTF                                
                                #                             pacf1.detrend.deseas=pacf1
)

results.out<-mutate(results.out,connect.any=ifelse(sign.AB.connect==1|sign.AB.connect==-1|sign.BC.connect==1|sign.BC.connect==-1,1,0))
results.out<-data.frame(lake=results.out$lake,station=results.out$station,variable=results.out$variable,regionloc=results.out$regionloc,
                        Latitude=results.out$Latitude,Longitude=results.out$Longitude,
                        signif(results.out[,-c(1:6)],3))

# write results to file
# these are the "by lake" results
write.csv(results.out, file="Outputs/Files/timeseries_results.csv")


#convert "variable" to character
results.connect.paired$variable<-as.character(results.connect.paired$variable)
#rename variables
results.connect.paired$variable[which(results.connect.paired$variable=="chla")]<-"chl a"
results.connect.paired$variable[which(results.connect.paired$variable=="propotherphyto")]<-"other phyto"
results.connect.paired$variable[which(results.connect.paired$variable=="propcyano")]<-"cyano"
results.connect.paired$variable[which(results.connect.paired$variable=="propchloro")]<-"chloro"
results.connect.paired$variable[which(results.connect.paired$variable=="propcrypto")]<-"crypto"
results.connect.paired$variable[which(results.connect.paired$variable=="propdiatom")]<-"diatom"
results.connect.paired$variable[which(results.connect.paired$variable=="propdaphnia")]<-"daphnia"
results.connect.paired$variable[which(results.connect.paired$variable=="propcalanoid")]<-"calanoid"
results.connect.paired$variable[which(results.connect.paired$variable=="propcyclopoid")]<-"cyclopoid"
results.connect.paired$variable[which(results.connect.paired$variable=="propothercladoc")]<-"other cladoc"
results.connect.paired$variable[which(results.connect.paired$variable=="propotherzoop")]<-"other zoop"
results.connect.paired$variable[which(results.connect.paired$variable=="proprotifer")]<-"rotifer"
results.connect.paired$variable[which(results.connect.paired$variable=="totdissnitro")]<-"TDN"
results.connect.paired$variable[which(results.connect.paired$variable=="totdissphos")]<-"TDP"
results.connect.paired$variable[which(results.connect.paired$variable=="totdoc")]<-"DOC"
results.connect.paired$variable[which(results.connect.paired$variable=="totnitro")]<-"TN"
results.connect.paired$variable[which(results.connect.paired$variable=="totphos")]<-"TP"
results.connect.paired$variable[which(results.connect.paired$variable=="watertemp")]<-"water temp"

# write to file
#write.csv(results.connect.paired, file="Outputs/Files/timeseries_results.connect.paired.csv")

######################################################################################################
################ histograms of pacf1 values by time series w/ seasonal diffs #########################
######################################################################################################

# prepare data for pacf plot
results.pacf1<- results %>% select(lake,station,variable,signif.seasdiff,pacf1)
acf.long.data<-subset(results.pacf1,!results.pacf1$variable %in% c("airtemp",
                                                                   "bactcount","ciliacount","ciliamass","color","photicdepth",
                                                                   "phytomass","radiation","secchidepth","zoopmass","sim"))
dataplot<-acf.long.data

# do pacf plot
pacf1.plot<-ggplot(dataplot,aes(x=pacf1,fill=factor(signif.seasdiff)))+
  geom_histogram(alpha=0.5) +#, origin=-1) + # ylim=c(-1,1) +
  facet_wrap(~variable) + 
  scale_fill_discrete(name="seasondiff") +
  geom_vline(xintercept=0,colour="gray")+
  theme_bw()

# write pacf plot to file
png(file="Outputs/Figures/timeseries_pacf1.png",width=14,height=12,unit="in",res=180)
pacf1.plot
dev.off()


####################### aggregating results ########################
# aggregate results by variable
results.byvar <- results %>%
  group_by(variable) %>%
  dplyr::summarize(timeseriescount=length(lake),
                   signif.seasidiff.pct=100*(sum(signif.seasdiff)/timeseriescount),
                   seasdiff.pos=100*(sum(sign.seasdiff==1)/timeseriescount),
                   seasdiff.neg=100*(sum(sign.seasdiff==-1)/timeseriescount),
                   AR1.pos=100*(sum(ifelse(pacf1>0,1,0),na.rm=TRUE)/timeseriescount),
                   AR1.neg=100*(sum(ifelse(pacf1<0,1,0),na.rm=TRUE)/timeseriescount),
                   AB.pos=100*(sum(ifelse(ab.cor.raw>0,1,0),na.rm=TRUE)/timeseriescount),
                   AB.neg=100*(sum(ifelse(ab.cor.raw<0,1,0),na.rm=TRUE)/timeseriescount),
                   BC.pos=100*(sum(ifelse(bc.cor.raw>0,1,0),na.rm=TRUE)/timeseriescount),
                   BC.neg=100*(sum(ifelse(bc.cor.raw<0,1,0),na.rm=TRUE)/timeseriescount),
                   AB.dtds.pos=100*(sum(ifelse(ab.cor>0,1,0),na.rm=TRUE)/timeseriescount),
                   AB.dtds.neg=100*(sum(ifelse(ab.cor<0,1,0),na.rm=TRUE)/timeseriescount),
                   BC.dtds.pos=100*(sum(ifelse(bc.cor>0,1,0),na.rm=TRUE)/timeseriescount),
                   BC.dtds.neg=100*(sum(ifelse(bc.cor<0,1,0),na.rm=TRUE)/timeseriescount),
                   sign.AB.connect.pos=100*(sum(ifelse(sign.AB.connect==1,1,0),na.rm=TRUE)/timeseriescount),
                   sign.AB.connect.neg=100*(sum(ifelse(sign.AB.connect==-1,1,0),na.rm=TRUE)/timeseriescount),
                   sign.AB.connect.any=100*(sum(ifelse(sign.AB.connect==1|sign.AB.connect==-1,1,0),na.rm=TRUE)/timeseriescount),
                   sign.BC.connect.pos=100*(sum(ifelse(sign.BC.connect==1,1,0),na.rm=TRUE)/timeseriescount),
                   sign.BC.connect.neg=100*(sum(ifelse(sign.BC.connect==-1,1,0),na.rm=TRUE)/timeseriescount),
                   sign.BC.connect.any=100*(sum(ifelse(sign.BC.connect==1|sign.BC.connect==-1,1,0),na.rm=TRUE)/timeseriescount),
                   sign.ABBC.connect.any=100*(sum(ifelse(sign.AB.connect==1|sign.AB.connect==-1|sign.BC.connect==1|sign.BC.connect==-1,1,0),na.rm=TRUE)/timeseriescount)
  ) %>% as.data.frame()
results.byvar 

results.byvar<-data.frame(variable=results.byvar$variable,signif(results.byvar[,-1],3))

# write results to file
write.csv(file="Outputs/Files/timeseries_results_byvar.csv",results.byvar)


# create "small" results files
results.small <- results.out %>% 
  filter(variable %in% c("chl a", "phytocount", "calanoid", "cyclopoid",
                         "daphnia", "other cladoc", "rotifer", "zoopcount","zoopcount_nonrotifer", "DOC",
                         "TDN", "TDP", "TN", "TP", "water temp"))

results.byvar.small <- results.byvar %>% 
  filter(variable %in% c("chl a", "phytocount", "calanoid", "cyclopoid",
                         "daphnia", "other cladoc", "rotifer","zoopcount","zoopcount_nonrotifer", "DOC",
                         "TDN", "TDP", "TN", "TP", "water temp"))

results.byvar.small$group<-""
results.byvar.small$group[which(results.byvar.small$variable %in% c("chl a","phytocount"))]<-"1-phyto"

results.byvar.small$group[which(results.byvar.small$variable %in% c("zoopcount","zoopcount_nonrotifer","daphnia","calanoid",
                                                                    "cyclopoid", "other cladoc","rotifer"))]<-"2-zoop"
results.byvar.small$group[which(results.byvar.small$variable %in% c("TP","TN","TDP","TDN","DOC"))]<-"3-chem"
results.byvar.small$group[which(results.byvar.small$variable %in% c("water temp"))]<-"4-phys"

results.byvar.small<-results.byvar.small %>% select(group,variable,timeseriescount,seasdiff.pos,seasdiff.neg,
                                                    sign.AB.connect.pos,sign.AB.connect.neg,
                                                    sign.BC.connect.pos,sign.BC.connect.neg,
                                                    sign.ABBC.connect.any)

results.byvar.small<-results.byvar.small[order(results.byvar.small$group,results.byvar.small$variable),]

num.cols <- sapply(results.byvar.small, is.numeric)
results.byvar.small[,num.cols]<-round(results.byvar.small[,num.cols], digits = 0)


results.byvar.small

results.small<-
  data.frame(lake=results.small$lake,station=results.small$station,
             variable=results.small$variable,regionloc=results.small$regionloc,
             stationlat=results.small$stationlat,stationlong=results.small$stationlong,
             signif(results.small[,-c(1:6)],3))

# write results to file
write.csv(results.small, "Outputs/Files/timeseries_results_small.csv", row.names = FALSE) 
write.csv(results.byvar.small, "Outputs/Files/timeseries_results_byvar_small.csv", row.names = FALSE) 

#####################################################################
######################### residual plots ############################
#####################################################################

# prepare data for residual plot
results.paired<- results.connect.paired %>% select(lake,station,variable,year,
                                                   value_iceon.ab=iceon.ab,value_iceoff.ab=iceoff.ab,
                                                   value_iceon.bc=iceon.bc,value_iceoff.bc=iceoff.bc,
                                                   signif.AB.connect,signif.BC.connect)

results.paired.chla<-subset(results.paired,results.paired$variable == "chl a")
results.paired.chla$lakesta<-paste(results.paired.chla$lake, results.paired.chla$station, sep = ".")

dataplot<-results.paired.chla

dataplot$lakesta <- as.factor(dataplot$lakesta)

# do paired AB|BC plots
paired.plot.chla.AB<-ggplot(dataplot,aes(x=value_iceoff.ab,y=value_iceon.ab,colour=factor(signif.AB.connect)))+
  geom_point() + xlab("prev summer value, detrended")+ ylab("winter value, detrended")+
  facet_wrap(~lakesta, scales = "free") + 
  scale_colour_discrete(name="relatedTF")+
  theme_bw() +
  theme(strip.text.x = element_text(size = 8))

paired.plot.chla.BC<-ggplot(dataplot,aes(x=value_iceon.bc,y=value_iceoff.bc,colour=factor(signif.BC.connect)))+
  geom_point() + xlab("winter value, detrended")+ ylab("next summer value, detrended")+
  facet_wrap(~lakesta, scales = "free") + 
  scale_colour_discrete(name="relatedTF")+
  theme_bw() +
  theme(strip.text.x = element_text(size = 8))

# write paired plot to file
png(file="Outputs/Figures/timeseries_paired.plot.chla.AB.png",width=14,height=12,unit="in",res=180)
paired.plot.chla.AB
dev.off()
png(file="Outputs/Figures/timeseries_paired.plot.chla.BC.png",width=14,height=12,unit="in",res=180)
paired.plot.chla.BC
dev.off()
