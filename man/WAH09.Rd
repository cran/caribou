\name{WAH09}
\alias{WAH09}
\docType{data}
\title{
Western Arctic Herd 2009 Data
}
\description{
Population size of the Western Arctic Herd (Alaska)  
from a postcalving survey in July 2009.
}
\usage{data(WAH09)}
\format{
9 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{xi}}{number of radio-collared animals in the detected (photographed) groups}
    \item{\code{gni}}{size of the detected groups}
  }
}
\details{
During this survey, 89 collars were active.\cr
89 collars (100\%) were photographed among the 9 detected groups.
}
\note{
This data set provide the opportunity to test the models under very high sampling regime. 
}
\source{
Parrett, L., unpubl. data
}
\examples{
data(WAH09)
petersen(WAH09, M=89)
abundance(WAH09, n=89, model="H")
abundance(WAH09, n=89, model="I")
abundance(WAH09, n=89, model="T", B=2)
}
\keyword{datasets}
