\name{GRH01}
\alias{GRH01}
\docType{data}
\title{
George River Herd 2001 Data
}
\description{
Population size of the George River (Riviere George) herd (Quebec) 
from a postcalving survey in July 2001.
}
\usage{data(GRH01)}
\format{
27 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{xi}}{number of radio-collared animals in the detected (photographed) groups}
    \item{\code{gni}}{size of the detected groups}
  }
}
\details{
During this survey, 109 collars were active.\cr
56 collars (57\%) were photographed among the 27 detected groups.
}
\source{
Couturier, S., Jean, D., Otto, R. and Rivard, S. (2004). \emph{Demography of the 
migratory tundra caribou (Rangifer tarandus) of the Nord-du-Quebec 
region and Labrador}. Min. Ressources naturelles, Faune et Parcs, Quebec. 68 p. 
ISBN: 2-550-43725-X
}
\examples{
data(GRH01)
petersen(GRH01, M=109)
abundance(GRH01, n=109, model="H")
abundance(GRH01, n=109, model="I")
abundance(GRH01, n=109, model="T", B=2)
abundance(GRH01, n=109, model="T", B=4)
abundance(GRH01, n=109, model="T", B=6)
}
\keyword{datasets}
