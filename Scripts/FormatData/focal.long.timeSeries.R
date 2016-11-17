###################################################################################
# The goal of this script is to look at lakes that have several years of data,    #
# and to focus on chlorophyll and immediately related variables (so to speak).    #
###################################################################################

###################################################################################
# Please note that results found using the publicly available dataset             #
# (available from the Knowledge Network for Biocomplexity, DOI: 10.5063/F12V2D1V) #
# will be slightly different than those published in Hampton et al. 2016, due     #
# to the exclusion of some lakes in the publicly  available dataset.              #  
###################################################################################  

########################### Code by Ryan Batt ##################################

#### clear workspace ####
rm(list = ls()) 
graphics.off()

# =============================================
# = Source All Scripts from Functions Folders =
# =============================================
func.location <- "Scripts/Functions"
invisible(sapply(list.files(func.location, full.names=TRUE, recursive=TRUE), source, .GlobalEnv))


# ==================
# = Load Libraries =
# ==================
library(data.table)


# =============
# = Load Data =
# =============

data.long <- fread("Data/data.long.csv")
data.long[,year:=as.numeric(year)]
data.long[season=="iceoff", year:=(year+0.5)] # "iceoff" comes after "iceon" of the same year; thus, add 0.5 to the "year" if the season is "iceoff"


# ===========================================
# = Subset to long time series =
# ===========================================

data.ts2 <- data.long[nYears>=10] # set 10 years of a variable's data as the cutoff

data.longTS <- data.ts2[stat=="ave"]
data.longTS <- data.longTS[, # remove years before any-non NA's were observed
	j={
		
		# next 2 lines are aimed at removing years where none of the variables of interest have had data yet (leading NA's for all), 
	  #or none of them will have data again (trailing NA's for all)
		yr <- range(year[!is.na(value)])
		out.0 <- .SD[year>=yr[1] & year<=yr[2]]
		
		# next lines are here to make the time series "regular"; even if the variable doesn't have data in a year, 
		#if it has data before and after that year, it should have an NA for the missing year, not just completely omit the year/row.
		full.year <- as.data.table(expand.grid(year=yr[1]:yr[2], variable=.SD[,unique(variable)]))
		# print(full.year)
		out <- merge(full.year, out.0, all=TRUE, by=c("year","variable"))
		out
	},
	by=c("lakename","stationname")
]
setkeyv(data.longTS, c("lakename","stationname","variable","year")) 


#write to csv
write.csv(data.longTS, "Data/data.longTS.csv", row.names=FALSE)


