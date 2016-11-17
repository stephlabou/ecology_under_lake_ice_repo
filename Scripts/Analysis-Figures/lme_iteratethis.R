##############################################################################################
# This script is called by lme_style.R This script runs and compares alternative lme models, #
# for one variable at a time (chla, DOC, etc.) and outputs variable-specific figures.        #
##############################################################################################

################################# Code by Steve Powers #######################################

#LME models
covars.lme<-lme(diff ~ lakemaxdepth+lakearea+lakeelevation+stationlat.abs,random=~1|year, data=data.lme)
depth.lme<-lme(diff ~ lakemaxdepth,random=~1|year, data=data.lme)
area.lme<-lme(diff ~ lakearea,random=~1|year, data=data.lme)
elev.lme<-lme(diff ~ lakeelevation,random=~1|year, data=data.lme)
lat.lme<-lme(diff ~ stationlat.abs,random=~1|year, data=data.lme)

depth.area.elev.lme<-lme(diff ~ lakemaxdepth+lakearea+lakeelevation,random=~1|year, data=data.lme)
depth.area.lat.lme<-lme(diff ~ lakemaxdepth+lakearea+stationlat.abs,random=~1|year, data=data.lme)
depth.elev.lat.lme<-lme(diff ~ lakemaxdepth+lakeelevation+stationlat.abs,random=~1|year, data=data.lme)
area.elev.lat.lme<-lme(diff ~ lakearea+lakeelevation+stationlat.abs,random=~1|year, data=data.lme)

depth.area.lme<-lme(diff ~ lakemaxdepth+lakearea,random=~1|year, data=data.lme)
depth.elev.lme<-lme(diff ~ lakemaxdepth+lakeelevation,random=~1|year, data=data.lme)
depth.lat.lme<-lme(diff ~ lakemaxdepth+stationlat.abs,random=~1|year, data=data.lme)

area.elev.lme<-lme(diff ~ lakearea+lakeelevation,random=~1|year, data=data.lme)
area.lat.lme<-lme(diff ~ lakearea+stationlat.abs,random=~1|year, data=data.lme)

elev.lat.lme<-lme(diff ~ lakeelevation+stationlat.abs,random=~1|year, data=data.lme)

#anova
anov<-anova(covars.lme,depth.lme,area.lme,elev.lme,lat.lme,
      depth.area.elev.lme,depth.area.lat.lme,depth.elev.lat.lme,area.elev.lat.lme,
      depth.area.lme,depth.elev.lme,depth.lat.lme,
      area.elev.lme,area.lat.lme,
      elev.lat.lme)

#make dataframe
df<-data.frame(anov)
df<-df[order(df$AIC),]

#find means
data.lme.agg<-data.lme %>% group_by(lakename,stationname) %>%
  dplyr::summarize(mean.diff=mean(diff),
                   lakemaxdepth=mean(lakemaxdepth,na.rm=TRUE),
                   lakearea=mean(lakearea,na.rm=TRUE),
                   lakeelevation=mean(lakeelevation,na.rm=TRUE),
                   stationlat.abs=mean(stationlat.abs,na.rm=TRUE)) %>% as.data.frame()

#make covariate crossplot figure
name<-paste("Outputs/Figures/covar.",data.lme$variable[1],".png",sep="")
png(name,height=8,width=8,units="in",res=300)
par(mfrow=c(2,2))

plot((data.lme.agg$lakemaxdepth),data.lme.agg$mean.diff)#,ylim=c(-1000,5000))

plot((data.lme.agg$lakearea),data.lme.agg$mean.diff)#,ylim=c(-1000,5000))

plot((data.lme.agg$lakeelevation),data.lme.agg$mean.diff)#,ylim=c(-1000,5000))

plot((data.lme.agg$stationlat.abs),data.lme.agg$mean.diff)#,ylim=c(-1000,5000))

mtext(text=data.lme$variable[1],side=3,outer=TRUE,padj=3)#,padj=-1)
dev.off()

