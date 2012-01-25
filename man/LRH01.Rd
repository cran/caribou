\name{LRH01}
\alias{LRH01}
\docType{data}
\title{
Leaf River Herd 2001 Data
}
\description{
Population size of the Leaf River (Riviere aux Feuilles) herd (Quebec) 
from a postcalving survey in July 2001.
}
\usage{data(LRH01)}
\format{
17 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{xi}}{number of radio-collared animals in the detected (photographed) groups}
    \item{\code{gni}}{size of the detected groups}
  }
}
\details{
During this survey, 120 collars were active.\cr
23 collars (19\%) were photographed among the 17 detected groups.
}
\note{
The small sample size of this census was caused by technical and weather 
related problems in July 2001. This provided an opportunity 
to see the behaviour of the different models under low sampling regime. 
}
\source{
Couturier, S., Jean, D., Otto, R. and Rivard, S. (2004). \emph{Demography of the 
migratory tundra caribou (Rangifer tarandus) of the Nord-du-Quebec 
region and Labrador}. Min. Ressources naturelles, Faune et Parcs, Quebec. 68 p. 
ISBN: 2-550-43725-X
}
\examples{
data(LRH01)
petersen(LRH01, M=120)
abundance(LRH01, n=120, model="H")
abundance(LRH01, n=120, model="I")
abundance(LRH01, n=120, model="T", B=2)
# The threshold model with B >= 3 is equivalent
# to the homogeneity model for this data set
# because max(LRH01$xi)=2
}
\keyword{datasets}
