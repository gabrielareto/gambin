.est_confint <-
function(est_alpha, est_likelihood, mydata, level)
{
  conf_logLik <- function(alpha, mydata, est_likelihood, level)
  {
    return(abs(.logLik_gamBin(alpha, mydata) - est_likelihood - exp(-qchisq(level,1)/2))) # the likelihood ratio test for confidence intervals (from "Beyond Traditional Statistical Measures")
  }
  lower <- optimise(conf_logLik, interval = c(0,est_alpha), mydata = mydata, est_likelihood = est_likelihood, level = level)$minimum
  higher <- optimise(conf_logLik, interval = c(est_alpha, 30), mydata = mydata, est_likelihood = est_likelihood, level = level)$minimum
  return(c(lower, higher))
}

.logLik_gamBin <-
function(alpha, mydata, maxoct) #the maxoct argument can be removed if we do not estimate maxoct
{
  if(missing(maxoct)) maxoct <- max(mydata$octave)
  dgamb <- dgambin(alpha, maxoct)
  exponent <- mydata$species # this line and the next can be removed if we do not estimate maxoctave
  if(length(exponent) < length(dgamb)) exponent[(length(exponent)+1):length(dgamb)] <- 0
  lik_dist <- dgamb^exponent
  logLik <- sum(-log(lik_dist))
  logLik
}


.sample_abundances <-
function(abundances, p=1)
{
  if(!is.numeric(abundances)) stop("abundances must be numeric")
  individuals = round(sum(abundances)*p)
  species <- as.factor(1:length(abundances))
  samplevector <- rep(species, abundances)
  samp <- sample(samplevector, individuals, replace = F)
  return(table(samp))
}



# This should function for continuous variables as well.
# One individual can be thought as a "cluster" of infinitesimal 
# abundances? Thus, aggregation patterns of infinitesimal abundances
# must be considered. 

# Aggregation patterns of infinitesimal abundances can be
# approximated by establishing a "cluster size". One cluster
# can be thought as one individual. The subsampling will
# affect to clusters, and not to each infinitesimal unit
# of abundance. The cluster size is a function of the
# clustering level or clustering intensity.

sample_abundances_general <- function(abundances, p=1, cl.size=0, cl.level=0, r=1, pS=1)
{
	# Quality control and workflow:
	if(!is.numeric(abundances)) stop("abundances must be numeric")
	if(is.na(p)) stop("please provide a proportion of abundances to subsample")
	abundances <- abundances[abundances>0]
	
	if(p == 1) method = "as.is"
	if(p != 1)
	{
		if(is.na(cl.size))
		{
			if(!is.na(cl.level)) method = "cl.level"
			if(is.na(cl.level))
			{
				if(!is.na(r) & !is.na(pS)) stop("please provide desired Pearson's r or proportion of species to keep, but not both")
				if(is.na(r) & !is.na(pS)) method = "find.by.pS"
				if(!is.na(r) & is.na(pS)) method = "find.by.r"
			}
		}
		
		if(!is.na(cl.size))
		{
			if(cl.size==0) method = "naive"
			if(cl.size>0 & is.na(cl.level)) method = "cl.size"
			if(cl.size>0 & !is.na(cl.level)) stop("please provide cluster size or clustering level, but not both")
			
		}
		
	}
	
	
	# method: as.is
	if(method == "as.is") out <- abundances
	
	# method: naive
	# This is based on the radical assumption of perfect
	# mixing of infinitesimal units of abundances. This
	# will keep all the species and the Pearson correlation
	# with the original abundances will be almost perfect.
	if(method == "naive") out <- abundances*p
	
	# method: cl.size
	# The user inputs directly the size of each cluster
	# and the algorithm interprets each cluster as one
	# 'individual'. I.e. it subsamples clusters of the
	# specified size.
	if(method == "cl.size")
	{
		clusters <- round(abundances/cl.size)
		n = round(sum(clusters)*p) # number of clusters to sample
		species <- as.factor(1:length(clusters))
		samplevector <- rep(species, clusters)
		samp <- sample(samplevector, n, replace = F)
		out <- as.numeric(table(samp)*cs)
	}
	
	# method: cl.level
	# This method takes the clustering level specified by
	# the user and calculates a cluster size from it. Then
	# it works as in method cl.size. If the clustering level
	# equals 0, clusters will be of size min(abundances).
	# If it equals 1, clusters will be of size max(abundances).
	if(method == "cl.level")
	{
		cs <- min(abundances) + diff(range(abundances))*cl.level
		clusters <- round(abundances/cs)
		n = round(sum(clusters)*p) # number of clusters to sample
		species <- as.factor(1:length(clusters))
		samplevector <- rep(species, clusters)
		samp <- sample(samplevector, n, replace = F)
		out <- as.numeric(table(samp)*cs)
	}
	
	# method: find.by.pS
	# This method optimizes the value of cluster size,
	# so a proportion of pS is kept.
	if(method == "find.by.pS")
	{
		S  = round(length(abundances)*pS) # desired number of species
		f <- function(cs)
		{
			clusters <- round(abundances/cs)
			n = round(sum(clusters)*p) # number of clusters to sample
			species <- as.factor(1:length(clusters))
			samplevector <- rep(species, clusters)
			samp <- sample(samplevector, n, replace = F)
			out <- abs(sum(table(samp)>0) - S)
			return(out)
		}
		
		cs <- median(unlist(replicate(100, optimise(f, interval=range(abundances)))[1,]))
		# Employ the optimized cluster size:
		clusters <- round(abundances/cs)
		n = round(sum(clusters)*p) # number of clusters to sample
		species <- as.factor(1:length(clusters))
		samplevector <- rep(species, clusters)
		samp <- sample(samplevector, n, replace = F)
		out <- as.numeric(table(samp)*cs)
	}
	
	# method: find.by.r
	# This method optimizes the value of cluster size,
	# so the subsampled abundances correlate with a
	# given Pearson's r with the original abundances.
	if(method == "find.by.r")
	{
		f <- function(cs)
		{
			clusters <- round(abundances/cs)
			n = round(sum(clusters)*p) # number of clusters to sample
			species <- as.factor(1:length(clusters))
			samplevector <- rep(species, clusters)
			samp <- sample(samplevector, n, replace = F)
			out <- abs(cor(as.numeric(table(samp)*cs), abundances) - r)
			return(out)
		}
		
		cs <- median(unlist(replicate(100, optimise(f, interval=range(abundances)))[1,]))
		# Employ the optimized cluster size:
		clusters <- round(abundances/cs)
		n = round(sum(clusters)*p) # number of clusters to sample
		species <- as.factor(1:length(clusters))
		samplevector <- rep(species, clusters)
		samp <- sample(samplevector, n, replace = F)
		out <- as.numeric(table(samp)*cs)
	}
	
	# Return output:
	cat("Method employed:", method)
	return(out)
}




# Example:
#a = colSums(BCI)
#x <- sample_abundances_general(abundances=a, p=0.15, cl.size=NA, cl.level=NA, r=NA, pS=0.50)

#plot(a, x)
#cor(a, x)
#sum(x>0)
#length(a)


