\name{create_octaves}
\alias{create_octaves}

\title{
Reclassify a vector of species' abundances into abundance octaves
}
\description{
Creates abundance octaves by a log2 transform that doubles the number of abundance classes within each octave (method 3 of Gray, Bjoergesaeter & Ugland 2006). Octave 0 contains the number of species with 1 individual, octave 1 the number of species with 2 or 3 individuals, octave 2 the number of species with 4 to 7 individuals, and so forth. 
}
\usage{
create_octaves(abundances, subsample = 0)
}

\arguments{
  \item{abundances}{ A numerical vector of species abundances in a community}
  \item{subsample}{ If > 0, the community is subsampled by this number of individuals before creating octaves. This is useful for analyses where \code{alpha} is estimated from a standardized number of individuals}
}

\value{ A data.frame with two variables: \code{octave} with the name of each octave and \code{species} with the number of species in that octave
}
\references{
Gray, J.S., Bjoergesaeter, A. & Ugland, K.I. (2006) On plotting species abundance distributions. Journal of Animal Ecology, 75, 752-756.
}
\author{
Michael Krabbe Borregaard
}

\examples{

data(moths)
create_octaves(moths)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
