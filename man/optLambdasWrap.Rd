\name{optLambdasWrap}
\alias{optLambdasWrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal ridge penalties with sequential optimization.
}
\description{
Sequentially optimizes a cross-validated score w.r.t. ridge penalties for multiple data blocks.
Also implements preferential ridge, which allows to first optimize for the preferential data types.
}
\usage{
optLambdasWrap(penaltiesinit = NULL, XXblocks, Y, X1 = NULL, pairing = NULL,
  folds, intercept = ifelse(class(Y) == "Surv", FALSE, TRUE), frac1 = NULL,
  score = "loglik", model = NULL, epsIWLS = 0.001, maxItrIWLS = 25,
  traceCV = TRUE, reltol = 1e-04, optmethod1 = "SANN", optmethod2 =
  ifelse(length(penaltiesinit) == 1, "Brent", "Nelder-Mead"), maxItropt1 = 10,
  maxItropt2 = 25, save = FALSE, parallel = FALSE, pref = NULL, fixedpen = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{penaltiesinit}{
Numeric vector. Initial values for penaltyparameters. May be obtained from \code{\link{fastCV2}}.
}
  \item{XXblocks}{
List of \code{nxn} matrices. Usually output of \code{\link{createXXblocks}}.
}
  \item{Y}{
Response vector: numeric, binary, factor or \code{survival}.
}
  \item{X1}{
Matrix. Dimension \code{n x p_0, p_0 < n}, representing unpenalized covariates }
  \item{pairing}{
Numerical vector of length 3 or \code{NULL} when pairs are absent. Represents
the indices (in \code{XXblocks}) of the two data blocks involved in pairing,
plus the index of the paired block.
}
    \item{folds}{
List, containing the splits of the samples. Usually obtained by
\code{\link{CVfolds}}
}
 \item{intercept}{
Boolean. Should an intercept be included?
}
  \item{frac1}{
Scalar. Prior fraction of cases. Only relevant for \code{model=" logistic"}.
}
  \item{score}{
Character. See Details.
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
}
  \item{epsIWLS}{
Scalar. Numerical bound for IWLS convergence.
}
  \item{maxItrIWLS}{
Integer. Maximum number of iterations used in IWLS.
}
  \item{traceCV}{
Boolean. Should the output of the IWLS algorithm be traced?
}
  \item{reltol}{
Scalar. Relative tolerance for optimization methods.
}
 \item{optmethod1}{
Character. First, global search method. Any of the methods \code{c("Brent",
"Nelder-Mead", "Sann")} may be used, but
simulated annealing by \code{"Sann"} is recommended to search a wide landscape.
Other unconstrained methods
offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{optmethod2}{
Character. Second, local search method. Any of the methods \code{c("Brent",
"Nelder-Mead", "Sann")} may be used, but
\code{"Nelder-Mead"} is generally recommended. Other unconstrained methods
offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{maxItropt1}{
Integer. Maximum number of iterations for \code{optmethod1}.
}
  \item{maxItropt2}{
Integer. Maximum number of iterations for \code{optmethod2}.
}

\item{save}{
Boolean. If TRUE appends the penalties and resulting CVscore to global variable \code{allscores}
}
  \item{parallel}{
Boolean. Should computation be done in parallel? If \code{TRUE}, requires to run \code{\link{setupParallel}} first.
}
  \item{pref}{
Integer vector or \code{NULL}. Contains indices of data types in \code{XXblocks}
that are preferential.
}
  \item{fixedpen}{
Integer vector or \code{NULL}. Contains indices of data types of which penalty
is fixed to the corresponding value in \code{penaltiesinit}.
}
}
\details{
As opposed to \code{\link{optLambdas}} this function first searches globally,
then locally.
Hence, more time-consuming, but better guarded against multiple local optima.

See \code{\link{Scoring}} for details on \code{score}. We highly recommend to
use smooth scoring functions, in particular \code{"loglik"}.
For ranking-based criteria like \code{"auc"} and \code{"cindex"} we advise to
use repeated CV (see \code{\link{CVfolds}}) to avoid ending up in any of the
many local optima.
}
\value{
List, with components:
\item{res}{Outputs of all optimizers used}
\item{lambdas}{List of penalties found by the  optimizers}
\item{optpen}{Numerical vector with final, optimal penalties}
}

\examples{
data(dataXXmirmeth)
resp <- dataXXmirmeth[[1]]
XXmirmeth <- dataXXmirmeth[[2]]

# Find initial lambdas: fast CV per data block separately.
cvperblock2 <- fastCV2(XXblocks=XXmirmeth,Y=resp,kfold=10,fixedfolds = TRUE)
lambdas <- cvperblock2$lambdas

# Create (repeated) CV-splits of the data.
leftout <- CVfolds(Y=resp,kfold=10,nrepeat=3,fixedfolds = TRUE)

# Optimizes cross-validate criterion (default: log-lik)
# Increase the number of iterations for optimal results
jointlambdas <- optLambdasWrap(penaltiesinit=lambdas, XXblocks=XXmirmeth,Y=resp,
folds=leftout,score="loglik",save=T,maxItropt1=5, maxItropt2=5)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optLambdas}} for one-pass optimization. \code{\link{fastCV2}} for initialization of penalties.A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}
