#################################################################
# Function to define which row contains the corresponding data  #
# from the season before the focal row                          #
#################################################################

##################### Code by Ryan Batt #########################

# ==========================
# = Season Before Function =
# ==========================
season.before <- function(x){
	yearV <- x[,unique(year)]
	seasonV <- x[,unique(season)]
	lakenameV <- x[,unique(lakename)]
	stationnameV <- x[,unique(stationname)]
	variableV <- x[,unique(variable)]
	statV <- x[,unique(stat)]
	
	
	iceon.logic <- seasonV=="iceon"
	iceon.int <- as.integer(iceon.logic)
	
	year.out <- yearV-iceon.int
	season.out <- c("iceon","iceoff")[iceon.int+1]
	
	if(is.na(stationnameV)){
		before.ind <- data.long[,year==year.out & season==season.out & lakename==lakenameV & is.na(stationname) & variable==variableV & stat==statV]
		logic.base <- c(year.out, season.out, lakenameV, variableV, statV)
		
	}else{
		before.ind <- data.long[,year==year.out & season==season.out & lakename==lakenameV & stationname==stationnameV & variable==variableV & stat==statV]
		logic.base <- c(year.out, season.out, lakenameV, stationnameV, variableV, statV)
		
	}
	
	# print(x)
	
	
	if(!any(before.ind)|any(is.na(logic.base))){
		value.out <- NA_real_
	}else{
		value.out <- data.long[before.ind,value]
	}
	
	
	return(value.out)
}