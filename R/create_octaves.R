create_octaves <- function(abundances, subsample = 0)  
{
	if(!subsample == 0) abundances <- .sample_abundances(abundances, subsample)
  stopifnot(is.numeric(abundances))
  abundances = abundances[abundances > 0]     # remove zeros
  octs = floor(sapply(abundances, log2))
  octs = factor(octs, levels = 0:max(octs))   # ensure that all octaves are tabled, even if no species fall in that octave
  ret <- data.frame(table(octs))
  names(ret) = c("octave","species")
  ret$octave = 0:(nrow(ret)-1)                # octaves are numbered from 0, 1, 2... etc.
  return(ret)
}


### There are different methods of binding abundances into
### abundance classes. There is only one method currently
### implemented in the gambin package. I have found useful
### other binding methods (allowing, e.g., continuous
### abundance data). These are different binding methods:

data(BCI)
abundances <- colSums(BCI)

# NOTE: for all these methods a subsampling procedure
# could be implemented. To be general, such subsampling
# must be expressed not in number of individuals, but
# as a proportion of abundances to be retained in the
# subsampling.


# Method 1: current method for gambin package:
# Changes:
	# log2(x) instead of sapply(x, log2). Faster.
	# octave names removed
	# option to return empty octaves or not
	
create_octaves1 <- function(abundances, allow.empty=TRUE)
{
	# Quality control:
	stopifnot(is.numeric(abundances))
	abundances = abundances[abundances > 0]
	
	# Binding:
	octs <- floor(log2(abundances))
	octs <- factor(octs, levels = 0:max(octs))
	freq <- as.numeric(table(octs))
	
	# Subsampling:
	#
	#
	
	# Output:
	if(!allow.empty) freq <- t[freq>0]
	return(freq)
}


# Method 2: Preston binding as in vegan, no tiesplit
# Reference code: vegan::as.preston()
create_octaves2 <- function(abundances, allow.empty=TRUE)
{
	# Quality control:
	stopifnot(is.numeric(abundances))
	abundances = abundances[abundances > 0]
	
	# Binding:
	xlog2 <- ceiling(log2(abundances))
	tmp <- table(xlog2)
	indx <- as.numeric(names(tmp)) + 1
	freq <- numeric(max(indx))
	freq[indx] <- tmp
	
	# Subsampling:
	#
	#
	
	# Output:
	if(!allow.empty) freq <- freq[freq>0]
	return(freq)
}


# Method 3: Preston binding as in vegan, with tiesplit
# Reference code: vegan::as.preston()
# This is, apparently, the original method of Preston, 
# and the default binding method employed in vegan.
create_octaves3 <- function(abundances, allow.empty=TRUE)
{
	# Quality control:
	stopifnot(is.numeric(abundances))
	abundances = abundances[abundances > 0]
	
	# Binding:
	xlog2 <- log2(abundances)
	ties <- xlog2 == ceiling(xlog2)
	tiefreq <- table(xlog2[ties])
	notiefreq <- table(ceiling(xlog2[!ties]))
	itie <- as.numeric(names(tiefreq)) + 1
	nitie <- as.numeric(names(notiefreq)) + 1
	freq <- numeric(max(itie + 1, nitie))
	freq[itie] <- tiefreq/2
	freq[itie + 1] <- freq[itie + 1] + tiefreq/2
	freq[nitie] <- freq[nitie] + notiefreq

	# Subsampling:
	#
	#

	# Output:
	if(!allow.empty) freq <- freq[freq>0]
	return(freq)
}


# Method 4: doubling classes of relative abundance in the [0, 1] interval.
# The octave with the highest abundances contains 0.5-1
# The second octave with the highest abundances contains 0.25-0.5
# Etc.
# The minimum abundance allowed is ~8.9e-16.
create_octaves4 <- function(abundances, allow.empty=TRUE)
{
	# Quality control:
	stopifnot(is.numeric(abundances))
	stopifnot(max(abundances)<=1)
	abundances = abundances[abundances > 0]
	
	# Binding:
	freq <- as.numeric(table(base::cut(abundances, breaks=1/(2^(50:0)))))
	freq <- freq[which(freq>0)[1]:length(freq)] # removes the left empty tail

	# Subsampling:
	#
	#

	# Output:
	if(!allow.empty) freq <- freq[freq>0]
	return(freq)
}


# Method 5: doubling classes of relative abundance in the [0, MAX] interval,
# being MAX the maximum observed relative abundance (never 1 if S>1).
# The octave with the highest abundances contains (0.5*MAX) to MAX
# The second octave with the highest abundances contains (0.25*MAX) to (0.5*MAX) 
# Etc.
# The minimum abundance allowed is ~8.9e-16 * MAX
create_octaves5 <- function(abundances, allow.empty=TRUE)
{
	# Quality control:
	stopifnot(is.numeric(abundances))
	stopifnot(max(abundances)<=1)
	abundances = abundances[abundances > 0]
	
	# Binding:
	MAX = max(abundances)
	freq <- as.numeric(table(base::cut(abundances, breaks=MAX/(2^(50:0)))))
	freq <- freq[which(freq>0)[1]:length(freq)] # removes the left empty tail

	# Subsampling:
	#
	#

	# Output:
	#if(!allow.empty) freq <- freq[freq>0] # never affects
	return(freq)
}











