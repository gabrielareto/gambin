
# The random generator (r*) take n species within the octaves:

# I think it is required to set up boundaries for each octave.
# The tie method of Preston makes things more difficult. 
# This would require further thought.

# If we want "fitGambin" to work well, then we must employ
# the same octave binding as fitGambin. I.e., 1 to the first
# octave, 2-3 to the second, 4-7 to the third, etc.
rgambin <- function(n, alpha, maxoctave, round=TRUE)
{
	# Define to which octave a given species belong:
	octs <- dgambin(alpha=alpha, maxoctave=maxoctave)
	octs.n <- sample(1:length(octs), size = n, replace=TRUE, prob=octs)
	
	# Define boundaries for octaves:
	# (following the fitGambin method,
	# but this should be changed if necessary)
	u <- 2^(1:length(octs))-1		# upper limit
	l <- c(1, u+1)[-(length(u)+1)]	# lower limit
	
	# Sample uniformily within boundaries:
	out <- runif(n, min=l[octs.n], max=u[octs.n])
		
	if(round) out <- round(out)
	return(out)
}

