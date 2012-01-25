\name{petersen}
\alias{petersen}
\alias{print.petersen}
\title{ Lincoln-Petersen estimator for caribou abundance }
\description{
The function \code{petersen} estimates the total population size based on the Lincoln-Petersen Index 
as applied to radio telemetry data by White and Garrott (1990). It uses the Lincoln-Petersen estimator 
with Chapman's (1951) bias correction and the bias corrected standard error estimator of 
Seber (1970) and Wittes (1972).
}
\usage{
petersen(mat, M, S = 0)

\method{print}{petersen}(x,\dots)
}
\arguments{
  \item{mat}{ A matrix containing in the first column the number of radio-collared animals in the 
              detected (photographed) groups and in the second column the corresponding size of the 
              detected  groups. }
  \item{M}{ A numeric: the total number of active collars during the census 
           (equivalent to the argument \code{n} in the function \code{\link{abundance}}).}  
  \item{S}{ A numeric: the minimum size that define well aggregated groups. Only observations
            from well aggregated groups (containing at least \code{S} animals) 
            are taken into account in the computation. The default is zero (every observation is 
            included in the computation).}
  \item{x}{ An object, produced by the \code{petersen} function, to print.}
  \item{\dots}{ Further arguments to be passed to methods (see \code{print.default}).}
}
\value{
  \item{G }{The number of well aggregated groups.}
  \item{R }{The total number of radio-collared animal observed in the well aggregated groups.}
  \item{C }{The total number of animals observed in the well aggregated groups containing at least one radio-collared animal during the survey.}    
  \item{T.hat }{The modified lincoln-Petersen estimator for the total number of animals in a herd.} 
  \item{se_T.hat }{The estimated standard error of \code{LP_T.hat}.}
  \item{mat_aggre }{A matrix containing a subset of the input matrix \code{mat}: the data for the well aggregated groups only.}
   \item{call }{The function call (object of class "call").}
}
\references{
Chapman, D. G. (1951). Some properties of the hypergeometric distribution with applications to 
zoological sample censuses. \emph{University of California Publications in Statistics}, 
\bold{1}(7), 131-160.

Seber, G.A.F. (1970). The effects of trap response on tag recapture estimates. \emph{Biometrics}, 
\bold{26}, 13-22.

White, G.C. and Garrott, R.A. (1990). \emph{Analysis of wildlife radio-tracking data}. San Diego: 
Academic Press.

Wittes, J.T. (1972). On the bias and estimated variance of Chapman's two-sample capture-recapture 
population estimate. \emph{Biometrics}, \bold{28}, 592-597.
}
\author{ Helene Crepeau \email{Helene.Crepeau@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and\cr
Serge Couturier \email{tuttu@videotron.ca} and\cr 
Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} }
\seealso{ \code{\link{abundance}} }
\examples{
data(GRH93)
petersen(GRH93, M=92) # default S=0
petersen(GRH93, M=92, S=4000)
}
