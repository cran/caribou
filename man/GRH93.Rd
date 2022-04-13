\name{GRH93}
\alias{GRH93}
\docType{data}
\title{
George River Herd 1993 Data
}
\description{
Population size of the George River (Riviere George) herd (Quebec and Labrador) 
from a postcalving survey in July 1993.
}
\usage{data(GRH93)}
\format{
28 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{xi}}{number of radio-collared animals in the detected (photographed) groups}
    \item{\code{gni}}{size of the detected groups}
  }
}
\details{
During this survey, 92 collars were active.\cr
73 collars (79\%) were photographed among the 28 detected groups.
}
\source{
Rivest, L.-P., Couturier, S. and Crepeau, H. (1998). Statistical Methods for estimating caribou abundance 
using postcalving aggregations detected by radio telemetry. \emph{Biometrics}, \bold{54}(3), 865-876.
}
\examples{
petersen(GRH93, M=92)
petersen(GRH93, M=92, S=4000)
abundance(GRH93, n=92, model="H")
abundance(GRH93, n=92, model="I")
abundance(GRH93, n=92, model="T", B=2)
abundance(GRH93, n=92, model="T", B=4)
abundance(GRH93, n=92, model="T", B=6)
}
\keyword{datasets}
