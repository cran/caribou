\name{GRH10}
\alias{GRH10}
\docType{data}
\title{
George River Herd 2010 Data
}
\description{
Population size of the George River (Riviere George) herd (Quebec) 
from a postcalving survey in July 2010.
}
\usage{data(GRH10)}
\format{
13 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{xi}}{number of radio-collared animals in the detected (photographed) groups}
    \item{\code{gni}}{size of the detected groups}
  }
}
\details{
During this survey, 71 collars were active.\cr
43 collars (61\%) were photographed among the 13 detected groups.
}
\source{
Couturier, S., unpubl. data
}
\examples{
petersen(GRH10, M=71)
petersen(GRH10, M=71, S=1000)
abundance(GRH10, n=71, model="H")
abundance(GRH10, n=71, model="I")
abundance(GRH10, n=71, model="T", B=2)
abundance(GRH10, n=71, model="T", B=4)
abundance(GRH10, n=71, model="T", B=6)
}
\keyword{datasets}
