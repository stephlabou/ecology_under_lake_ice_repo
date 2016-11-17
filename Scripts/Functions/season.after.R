##################################################################
# Function to define which row contains the corresponding data   #
# from the season after the focal row                            #
##################################################################

################ Code by Ryan Batt ########################

# =========================
# = Season After Function =
# =========================

season.after <- function(x){
	yearV <- x[,unique(year)]
	seasonV <- x[,unique(season)]
	lakenameV <- x[,unique(lakename)]
	stationnameV <- x[,unique(stationname)]
	variableV <- x[,unique(variable)]
	statV <- x[,unique(stat)]
	
	message(paste(lakenameV, yearV, seasonV, variableV, statV))

	
	iceoff.logic <- seasonV=="iceoff"
	iceoff.int <- as.integer(iceoff.logic)
	
	year.out <- yearV+iceoff.int
	season.out <- c("iceoff", "iceon")[iceoff.int+1]
	
	if(is.na(stationnameV)){
		after.ind <- data.long[,year==year.out & season==season.out & lakename==lakenameV & is.na(stationname) & variable==variableV & stat==statV]
		logic.base <- c(year.out, season.out, lakenameV, variableV, statV)
		
	}else{
		after.ind <- data.long[,year==year.out & season==season.out & lakename==lakenameV & stationname==stationnameV & variable==variableV & stat==statV]
		logic.base <- c(year.out, season.out, lakenameV, stationnameV, variableV, statV)
		
	}
	
	
	
	if(!any(after.ind)|any(is.na(logic.base))){
		value.out <- NA_real_
	}else{
		value.out <- data.long[after.ind,value]
	}
	
	
	return(value.out)	
}
