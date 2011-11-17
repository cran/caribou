\name{abundance}
\alias{abundance}
\alias{print.abundance}
\title{ Post-calving method for caribou abundance estimation }
\description{
The function \code{abundance} applies the methodology found in Rivest, L.-P., Couturier, S. and Crepeau, H. (1998) for estimating caribou abundance using postcalving aggregations detected by radio telemetry.}
\usage{
abundance(mat, n, model = c("H", "I", "T"), B)

\method{print}{abundance}(x,\dots)
}
\arguments{
  \item{mat}{ A matrix containning in the first column the number of radio-collared animals in the detected (photographed) groups and in the second column the corresponding size of the detected  groups. }
  \item{n}{ A numeric: the total number of active collars during the census.}  
  \item{model}{ A character string indicating the model to determine the probability that a group with collared animals is detected \eqn{p_i}{pi}. It can be either "H" = homogeneity model, "I" = independence model or "T" = threshold model (see \bold{Details}). The default is "H". }
  \item{B}{ A numeric: a bound for the threshold model. }
  \item{x}{ An object, produced by the \code{abundance} function, to print.}
  \item{\dots}{ Further arguments to be passed to methods (see \code{print.default}).}
}
\details{
DETECTION MODELS\cr\cr

\bold{- homogeneity model} (\code{model="H"}):
\deqn{p_i = r \mbox{    provided} X_i>1}{pi = r provided xi>1}\cr

\bold{- independence model} (\code{model="I"}):
\deqn{p_i = 1-(r^{X_i})}{pi = 1-(r^xi)}\cr
	
\bold{- threshold model} (\code{model="T"}):
\deqn{
p_i = \left\{ \begin{array}{ll}
1 & \mbox{if } X_i \geq \mbox{\code{B}} \\
r & \mbox{if } 1 \leq X_i < \mbox{\code{B}} \end{array} \right.
}{pi = 1 if xi>=B, r if 1<=xi<B}

where \eqn{p_i}{pi} is the probability that a group with collared animals is detected, \eqn{X_i}{xi} is the number of radio-collared in the detected (photographed) groups and \eqn{r} is a parameter related to the probability of detection. For the threshold model, \code{B} is a bound given as a function's argument.
}
\value{
  \item{mp }{The number of detected groups having radio-collared animals.}		
  \item{xt }{The total number of radio-collared animals found in the detected groups.}
  \item{gnt }{The total number of animals counted in the detected groups.}
  \item{rr }{The estimated parameter related to the probability of detection.}
  \item{se_rr }{The estimated standard error of \code{rr}.}    
  \item{mat_pi }{A matrix containing a sorted copy of the input matrix \code{mat} in the first two columns, the detection probabilities \eqn{p_i}{pi} in the third column and the probabilities that the group has at least one collared animal \eqn{\pi_i}{pi.i} in the last column.}
  \item{T.hat }{The estimator for the total number of animals in a herd.} 
  \item{se_T.hat }{The estimated standard error of \code{T.hat}.}
  \item{loglikelihood }{The maximum value of the loglikelihood function for the detected model.}
  \item{randomness_test }{A vector with the statistic and the p-value of a score test for the randomness assumption available only for the homogeneity, independance and threshold model with B=2 or 3.}
  \item{call }{The function call (object of class "call").}
}
\references{
Rivest, L.-P., Couturier, S. and Crepeau, H. (1998). Statistical Methods for estimating caribou abundance using postcalving aggregations detected by radio telemetry. \emph{Biometrics}, \bold{54}(3), 865-876.
}
\author{ Helene Crepeau \email{Helene.Crepeau@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and\cr
Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} }
\seealso{ \code{\link{petersen}} }
\examples{
data(GRH93)
abundance(GRH93, n=92, model="H")
abundance(GRH93, n=92, model="I")
abundance(GRH93, n=92, model="T", B=2)
abundance(GRH93, n=92, model="T", B=4)
abundance(GRH93, n=92, model="T", B=6)
}
