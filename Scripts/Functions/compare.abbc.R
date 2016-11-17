###########################################################
# Functions to compare 2 seasons, or in general, 2 values #
###########################################################

################ Code by Ryan Batt ########################

# ==========================
# = Season Index Functions =
# ==========================
dev <- function(x,y){y-x} # difference
abs.dev <- function(x,y){abs(y-x)} # absolute deviation
perc.change <- function(x, y){ # percent change
	(y-x)/x
}
abs.perc.change <- function(x, y){ # absolute percent change
	abs((y-x)/x)
}

abs.rel.dev <- function(x, y){
	mu <- (x+y)/2
	abs((y-x)/mu)
}