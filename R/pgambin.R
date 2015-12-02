# The distribution function (p*) gives the probability
# that a variable X takes on a value less than or equal to a number x.
pgambin <- function(alpha, maxoctave) cumsum(dgambin(alpha=alpha, maxoctave=maxoctave))
