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


sample_abundances0 <-
function(abundances, p=1)
{
  if(!is.numeric(abundances)) stop("abundances must be numeric")
  individuals = round(sum(abundances)*p)
  species <- as.factor(1:length(abundances))
  samplevector <- rep(species, abundances)
  samp <- sample(samplevector, individuals, replace = F)
  return(table(samp))
}



# The following function extends the subsampling to continuous
# abundances (e.g. relative abundances, biomass, etc).
# One individual can be thought as a "cluster" of infinitesimal 
# abundances. In other cases, however, how many infinitesimal
# abundances constitute a 'natural unit' of abundances is unknown.
# Different kinds of aggregation patterns of infinitesimal abundances
# can be considered. These can be approximated by establishing a 
# "cluster size" (one cluster can be thought as one individual).
# The subsampling will affect to clusters, and not to each infinitesimal unit
# of abundance. The cluster size can be established through different methods.

# This function requires the following parameters:
# abundances: a numeric vector with abundances in any unit.
# p: proportion of abundances (e.g. of individuals) to subsample.
# method: the method that defines how to handle abundances:
	# "cl.size" (default option)
	# The user inputs directly the size of each cluster
	# and the algorithm interprets each cluster as one
	# 'individual'. I.e. it subsamples clusters of the
	# specified size in parameter cl.size.
	
	# "cl.level":
	# This method takes the clustering level specified by
	# the user and calculates a cluster size from it. Then
	# it works as in method cl.size. 

	# "find.by.pS":
	# This method optimizes the value of cluster size,
	# so a proportion of pS species is kept.
	
	# "find.by.r":
	# This method optimizes the value of cluster size,
	# so the subsampled abundances correlate with a
	# given Pearson's r with the original abundances.
	
	# "naive": 
	# This is based on the radical assumption of perfect
	# mixing of infinitesimal units of abundances. This
	# will keep all the species and the Pearson correlation
	# with the original abundances will be almost perfect.

# cl.size: Parameter for the "cl.size" method. It indicates
# how much abundance is to be contained within one cluster
# (~individual). And then performs discrete sampling of
# clusters. The default option is cl.size=min(abundances),
# so if abundances are expressed in individuals, and the
# minimum abundance is 1 individual, it will perform the
# standar sampling of individuals that one may expect.

# cl.level: Parameter for the "cl.level" method:
	# If cl.level=0, clusters will be of size min(abundances).
	# If cl.level=1, clusters will be of size max(abundances).
	# Any value within [0, 1] is allowed.

# pS: Parameter for the "find.by.pS" method. It is the
# proportion of species in [0, 1] to be kept in the sumsampling.

# r: Parameter for the "find.by.r" method. It is the Pearson's
# correlation desired between the original abundances and the
# subsampled abundances. It must be in [0, 1].

sample_abundances <- function(abundances, p=1, method=c("cl.size", "cl.level", "find.by.pS", "find.by.r", "naive"), cl.size=min(abundances), cl.level=NA, r=NA, pS=NA)
{
	# Quality control and workflow:
	if(!is.numeric(abundances)) stop("abundances must be numeric")
	if(is.na(p)) stop("please provide a proportion of abundances to subsample")
	if(length(method)>1) method = "cl.size"
	if(method == "cl.size")
	{
		if(is.na(cl.size))
		{
			warning("no cluster size provided: it will consider the minimum observed abundance as the cluster size")
			cl.size = min(abundances)
		}
		
		if(!is.na(cl.size)) if(cl.size>max(abundances)) stop("cluster size must be less than the maximum observed abundance")
	}
	
	if(method == "cl.level")
	{
		if(is.na(cl.level)) warning("no clustering level provided: it will consider the minimum observed abundance as the cluster size")
		cl.level = 0
	}
	
	if(method == "find.by.pS")
	{
		if(is.na(pS)) stop("a proportion of species to be kept must be provided with method 'find.by.pS'")
		if(!is.na(pS)){
			if(pS<=0 | pS>1) stop("a proportion of species to be kept in the interval (0,1] must be provided")
		}
	}
	
	if(method == "find.by.r")
	{
		if(is.na(r)) stop("a desired Pearson's r must be provided with method 'find.by.r'")
		if(!is.na(r)){
			if(r<0 | r>1) stop("a desired Pearson's r in the interval [0,1] must be provided")
		}
	}
	
	abundances <- abundances[abundances>0]
	if(p == 1) out <- abundances
	
	# Calculation of the cluster size:	
	if(p < 1 & method == "cl.size")
		cs <- cl.size		

	if(p < 1 & method == "cl.level")
		cs <- min(abundances) + diff(range(abundances))*cl.level

	if(p < 1 & method == "find.by.pS")
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
	}
	
	if(p < 1 & method == "find.by.r")
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
	}
	
	# Implements the calculations:
	if(p < 1 & method!="naive")
	{
		clusters <- round(abundances/cs)
		n = round(sum(clusters)*p) # number of clusters to sample
		species <- as.factor(1:length(clusters))
		samplevector <- rep(species, clusters)
		samp <- sample(samplevector, n, replace = F)
		out <- as.numeric(table(samp)*cs)
	}
	if(p < 1 & method == "naive") out <- abundances*p
	
	# Return output:
	cat("Method employed:", method, "\n")
	cat("Cluster size:", cs, "\n")
	return(out)
}




# Example:
#a = colSums(BCI)

#set.seed(165)
#x1 <- sample_abundances(a, p=0.5)
#set.seed(165)
#x2 <- sample_abundances0(a, p=0.5)
#all(x1==x2) # must be TRUE



