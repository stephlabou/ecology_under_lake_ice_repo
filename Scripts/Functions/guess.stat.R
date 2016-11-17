#########################################################################
# Function to figure out what statistic is assocated w/ a column name   #
# (used in create.data.long.R)                                          #
#########################################################################

#################### Code by Ryan Batt ##################################

guess.stat <- function(x){
	is.cv <- grepl("^cv", x)
	is.max <- grepl("^max", x)
	
	stat <- rep(NA, length(x))
	
	stat[is.cv] <- "cv"
	stat[is.max] <- "max"
	stat[!is.cv & !is.max] <- "ave"
	
	stat
}
strip.stat <- function(x){
	gsub("^(cv|ave|max)", "", x)
}