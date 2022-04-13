\name{WAH11}
\alias{WAH11}
\docType{data}
\title{
Western Arctic Herd 2011 Data
}
\description{
Population size of the Western Arctic Herd (Alaska)  
from a postcalving survey in 2011.
}
\usage{data(WAH11)}
\format{
10 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{xi}}{number of radio-collared animals in the detected (photographed) groups}
    \item{\code{gni}}{size of the detected groups}
  }
}
\details{
During this survey, 97 collars were active.\cr
96 collars (99\%) were photographed among the 10 detected groups.
}
\note{
This data set provide the opportunity to test the models under very high sampling regime. 
}
\source{
Dau, J., unpubl. data
}
\examples{
petersen(WAH11, M=97)
abundance(WAH11, n=97, model="H")
abundance(WAH11, n=97, model="I")
abundance(WAH11, n=97, model="T", B=2)
}
\keyword{datasets}

