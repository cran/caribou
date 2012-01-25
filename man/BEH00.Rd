\name{BEH00}
\alias{BEH00}
\docType{data}
\title{
Bluenose-East Caribou Herd 2000 Data
}
\description{
Population size of the Bluenose-East caribou herd (Northwest Territories, Canada) 
from a postcalving survey in June-July 2000.
}
\usage{data(BEH00)}
\format{
16 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{xi}}{number of radio-collared animals in the detected (photographed) groups}
    \item{\code{gni}}{size of the detected groups}
  }
}
\details{
During this survey, 33 collars were active (Patterson et al. also made calculation considering that
30 collars were active because 3 caribou were never located after collaring).\cr
23 collars (70\%) were photographed among the 16 detected groups.
}
\note{
This data set excludes the six groups without radio-collared animal from Table 1 of Patterson et al. 
(2004) since they cannot be included in the estimates of the total population size.\cr

For the last three groups, the size \code{gni} has been modified according to what Patterson et al. 
(2004) did (38\% of what appears in Table 1).  
}
\source{
Patterson, B. R., Olsen, B. T. and Joly, D. O. (2004). Population estimate for the Bluenose-East 
caribou herd using post-calving photography. \emph{Arctic}, \bold{57}, 47-58. 
}
\examples{
data(BEH00)
petersen(BEH00, M=33)
abundance(BEH00, n=33, model="H")
}
\keyword{datasets}
