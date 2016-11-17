###################################################################################
# This script reformats the synthesis under ice data set to a "long" form used in # 
# subsequent analyses in Hampton et al. 2016. In the process of converting from   #
# wide to long format, this script also aggregates values across pooled stations  #
# (see Hampton et al. 2016 for more details). Subsequent references to            #
# "stationname" (in the csv outputs here, and in  later scripts using those csv   #
# outputs) are the pooled stationnames, referred to in the original data as       #
# "poolstation."                                                                  #
###################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
###################################################################################                          

########################### Code by Ryan Batt #################################

#set working directory - setwd() - to "core" folder (e.g. the name of this entire repository)

#### clear workspace ####
rm(list = ls()) 
graphics.off()


# =============================================
# = Source All Scripts from Functions Folders =
# =============================================
func.location <- "Scripts/Functions"
invisible(sapply(list.files(func.location, full.names=TRUE, recursive=TRUE), source, .GlobalEnv))


# =================
# = Load Packages =
# =================
library(data.table)
# library(reshape2)
# library(plyr)
# library(dplyr)


# =============================
# = Read in Data Set Collated =
# =============================
# read in data output from pre_edits.R

data.raw <- fread("Data/public_data_prepped.csv", stringsAsFactors = FALSE)

#Warnings are fine


# ==================
# = Make Data Long =
# ==================
raw.names <- names(data.raw) # column names

# columns that contain identifying values
id.vars <- c("lakename", "stationlat", "stationlong", "stationname", "poolstation", 
             "sampletype", "year","season", "researcher", "lakeregloc", "lakecountry", 
             "lakearea","lakemeandepth", "lakemaxdepth", "lakeelevation", "watershedarea",
             "h2oresidence", "lakefetch", "stationdistance", "stationdepth","multiplestations", 
             "startday", "startmonth", "startyear", "endday","endmonth", "endyear", "iceduration", "periodn")


# logic for which columns contain measured values vs. identifying values, numeric vs not
measure.vars <- raw.names[!raw.names%in%id.vars]
numeric.vars <- raw.names[sapply(data.raw, class)%in%c("integer","numeric")]
numeric.measure <- raw.names[raw.names%in%measure.vars & raw.names%in%numeric.vars]
not.numeric.measure <- raw.names[raw.names%in%measure.vars & !raw.names%in%numeric.vars]

# Copy 'raw' data, subset, then melt into long format
data.trim <- copy(data.raw) # copy (use this b/c of data.table)
data.trim[,c(not.numeric.measure):=NULL] # subset
data.long <- melt(data.trim, id.vars=id.vars, measure.vars=numeric.measure) # melt

#There will be a warning message - this is fine.


# ==============================================
# = Split Variables in ave, max, cv Statistics =
# ==============================================
data.long[,old.variable:=variable] # copy the original measured variable column name into a column called "old.variable"
data.long[,variable:=strip.stat(variable)] # remove the stat from the variable name (e.g., avechla goes to chla)
data.long[,stat:=guess.stat(old.variable)] # put the statistic into its own column (e.g., avechla now has stat "ave")

setkeyv(data.long, c(id.vars, "variable", "stat")) # set the "key" for the data.table

# columns that are not time-varying, thus could be used as a key
key.noTime <- c("researcher", "lakename", "lakeregloc", "lakecountry", "lakearea", "lakemeandepth", 
                "lakemaxdepth", "lakeelevation", "watershedarea", "h2oresidence", "lakefetch", 
                "stationdistance", "stationdepth", "stationname", "poolstation", "stationlat", "stationlong", 
                "multiplestations", "startday", "startmonth", "startyear", "endday", "endmonth", "endyear", 
                "iceduration", "periodn", "samplenarrat", "sampletype", "variable", "stat")


# ================
# = Ensure Order =
# ================
w.season <- which(key(data.long)=="season")
setorderv(data.long, c(key(data.long)[1:w.season]), order=c(rep(1, w.season-1), -1))


# =================================
# = Aggregate across poolstations =
# =================================

# Create data.table to contain pooled data
data.long.pool <- copy(data.long)

# Pool stations - aggregate values
 data.long.pool[,
 	c(
 	"value.pool"
 	)
 	:=
 	list(
 		mean(value, na.rm=TRUE)
 	), 
 		
 	by=c("lakename","poolstation","variable","stat","year","season")
 ]

# Pool stations - aggregate lat and long, distance and depth
data.long.pool[,
               c(
                 "stationlat",
                 "stationlong",
                 "stationdistance",
                 "stationdepth"
               )
               :=
                 list(
                   mean(stationlat, na.rm=TRUE),
                   mean(stationlong, na.rm=TRUE),
                   mean(stationdistance, na.rm=TRUE),
                   mean(stationdepth, na.rm=TRUE)
                 ), 
               
               by=c("lakename","poolstation")
               ]


# switch value.pool NaN's to NA's
data.long.pool[is.nan(value.pool), value.pool:=NA_real_]

# Set key, and drop rows that were previously only distinguished by station, but are equivalent in poolstation
setkey(data.long.pool, lakename, poolstation, variable, stat, year, season)
data.long.pool <- unique(data.long.pool)


# put value.pool into value, drop unused columns
data.long.pool[!is.na(poolstation),stationname:=poolstation]
data.long.pool[,value:=value.pool]
data.long.pool[,value.pool:=NULL]


# rename data.long.pool as data.long
data.long.unpool <- copy(data.long)
data.long <- copy(data.long.pool)


# =================================================
# = Function for reformatting before/after values =
# =================================================
get.before.after <- function(x){
	
	fix.fill.names <- c(
		"stationlat",
		"stationlong",
		"lakeregloc",
		"lakecountry",
		"lakearea",
		"lakemeandepth",
		"lakemaxdepth",
		"lakeelevation",
		"watershedarea",
		"h2oresidence",
		"lakefetch",
		"stationdistance",
		"stationdepth",
		"multiplestations",
		"old.variable"
	)
	
	fix.fill.list <- as.list(x[,list(
		mean(stationlat,na.rm=T),
		mean(stationlong,na.rm=T),
		unique(lakeregloc),
		unique(lakecountry),
		mean(lakearea,na.rm=T),
		mean(lakemeandepth,na.rm=T), # some lakes have more than 1 value entered for this
		mean(lakemaxdepth,na.rm=T),
		unique(lakeelevation),
		unique(watershedarea),
		unique(h2oresidence),
		unique(lakefetch),
		unique(stationdistance),
		unique(stationdepth),
		unique(multiplestations),
		unique(as.character(old.variable))
	)])
	
	skele <- x[,list(
		year=rep(seq(min(year),max(year)),each=2),
		season=NA
	)]
	
	skele[,season:=rep(c("iceon","iceoff"),length(unique(year)))]
	setkey(skele, year, season)
	setkey(x, year, season)
	test3 <- x[skele, allow.cartesian=TRUE]
	
	setorderv(test3, c("year", "season"), order=c(1,-1))

	
	
	test3[, c("value.before", "value", "value.after"):=list(
		value.before=c(NA, value[-length(value)]),
		value=value, 
		value.after=c(value[-1],NA)
	)]

}

#Note: this may take a few minutes
data.long2 <- data.long[, j={
	send <- copy(.SD)
	get.before.after(send)

}, by=c("lakename", "stationname", "variable", "stat")]



# ===================================
# = Add Useful Columns to data.long =
# ===================================
# number of years a particular variable was observed for a lake. Might be useful to include more values in "by" to make it site specific.
data.long2[,nYears:=length(unique(year[!is.na(value)])), by=c("lakename", "stationname","variable","stat")] 

# number of seasons that this variable was measured in this year-location-site. Should generally be 2.
data.long2[,n.season:=length(unique(season[!is.na(value)])), by=c("year","lakename", "stationname", "stat", "variable")] 

# number of years with variable measured both seasons
data.long2[,nFullYears:=length(unique(year[!is.na(value) & n.season==2])), by=c("lakename", "stationname", "stat", "variable")] 

# if the station is "NA", just assume it means that there was only 1 station (which I'll call unnamed, not that there were multiple station, but unsure of which one these data are from.
data.long2[is.na(stationname), stationname:="unnamed"] 


# ====================================
# = Aggregate data.long across years =
# ====================================

u.noNA <- function(x){unique(x[!is.na(x)])}

data.long.timeCollapse <- data.long2[n.season==2,
	list(
		stationlat=mean(stationlat,na.rm=T), 
		stationlong=mean(stationlong, na.rm=T), 
		lakecountry=u.noNA(lakecountry), 
		lakearea=mean(lakearea, na.rm=T), 
		lakemeandepth=mean(lakemeandepth, na.rm=T), 
		lakeelevation=mean(lakeelevation, na.rm=T), 
		watershedarea=mean(watershedarea, na.rm=T), 
		h2oresidence=mean(h2oresidence, na.rm=T), 
		lakefetch=mean(lakefetch, na.rm=T), 
		stationdistance=mean(stationdistance, na.rm=T), 
		iceduration=mean(iceduration, na.rm=T), 
		value=mean(value, na.rm=TRUE), 
		value.before=mean(value.before, na.rm=TRUE), 
		value.after=mean(value.after, na.rm=TRUE), 
		stationdepth=mean(stationdepth, na.rm=TRUE), 
		nFullYears=unique(nFullYears)
	), 
	by=c("season", "lakename", "stationname","variable", "stat")
]



data.long <- copy(data.long2)


setkey(data.long, lakename, stationname, year, season)
setcolorder(data.long, c("lakename", "stationname", "year", "season", "poolstation", "variable", "stat",
                         "sampletype", "stationlat","stationlong","researcher", "lakeregloc","lakecountry","lakearea",
                         "lakemeandepth","lakemaxdepth","lakeelevation","watershedarea","h2oresidence",
                         "lakefetch","stationdistance","stationdepth","multiplestations",
                         "startday","startmonth","startyear","endday","endmonth","endyear","iceduration",
                         "periodn","nYears","n.season", "nFullYears","old.variable","value.before","value","value.after"))


# ========================
# = Fix Station Matching =
# ========================
# handy expression for use in data.table "i"
upls <- expression(unique(paste(lakename,stationname2)))
pls <- expression((paste(lakename,stationname2)))

# If station name is same as lake name, rename it as "unnamed"
# If there is only one station name, rename it as "only Station"
# This will choose simpler station names that can be applied to both data sets
# thus reducing the number of chances for arbitrary mismatches in station names
data.long[,stationname2:=stationname]
data.long[stationname2==lakename, stationname2:="unnamed"]
data.long[,n.stationname2:=length(eval(upls)), by=c("lakename")]
data.long[n.stationname2==1,stationname2:="onlyStation"]


# ========
# = Save =
# ========

write.csv(data.long, "Data/data.long.csv", row.names=FALSE)  

write.csv(data.long.timeCollapse, "Data/data.long.timeCollapse.csv", row.names=FALSE)

