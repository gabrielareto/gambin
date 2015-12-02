
# The quantile function (q*) finds the 'boundary' that determines
# a given value of the distribution function (i.e. pgambin).
# This is somewhat ackward given the discrete nature of octaves.
qgambin <- function(q, alpha, maxoctave)
{
	p <- pgambin(alpha=alpha, maxoctave=maxoctave)
	qgambin0 <- function(q0) which(p>=q0)[1]
	out <- sapply(q, qgambin0)
	if(!is.null(names(out))) out = names(out) # very weird
	return(out)
}

