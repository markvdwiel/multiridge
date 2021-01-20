\name{optLambdas}
\alias{optLambdas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal ridge penalties.
}
\description{
Optimizes a cross-validated score w.r.t. ridge penalties for multiple data blocks.
}
\usage{
optLambdas(penaltiesinit = NULL, XXblocks, Y, X1 = NULL, pairing = NULL, folds,
  intercept = ifelse(class(Y) == "Surv", FALSE, TRUE), frac1 = NULL, score = "loglik",
  model = NULL, epsIWLS = 0.001, maxItrIWLS = 25, traceCV = TRUE, reltol = 1e-04,
  optmethod = ifelse(length(penaltiesinit) == 1, "Brent", "Nelder-Mead"), maxItropt = 500,
  save = FALSE, parallel = FALSE, fixedpen = NULL, fixedseed = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{penaltiesinit}{
Numeric vector. Initial values for penaltyparameters. May be obtained from \code{\link{fastCV}}.
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
Numerical vector of length 3 or \code{NULL} when pairs are absent. Represents the indices (in \code{XXblocks}) of the two data blocks involved in pairing,
plus the index of the paired block.
}
  \item{folds}{
List, containing the splits of the samples. Usually obtained by \code{\link{CVfolds}}
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
  \item{optmethod}{
Character. Optimization method. Any of the methods \code{c("Brent", "Nelder-Mead", "Sann")} may be used, but
\code{"Nelder-Mead"} is generally recommended. Other unconstrained methods offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{maxItropt}{
Integer. Maximum number of iterations for \code{optmethod}.
}
\item{save}{
Boolean. If TRUE appends the penalties and resulting CVscore to global variable \code{allscores}
}
  \item{parallel}{
Boolean. Should computation be done in parallel? If \code{TRUE}, requires to run \code{\link{setupParallel}} first.
}
  \item{fixedpen}{
Integer vector or \code{NULL}. Contains indices of data types of which penalty is fixed to the corresponding value in \code{penaltiesinit}.
}
  \item{fixedseed}{
Boolean. Should the initialization be fixed? For reproducibility.
}
}
\details{
See \code{\link{Scoring}} for details on \code{score}.
 We highly recommend to use smooth scoring functions, in particular \code{"loglik"}.
For ranking-based criteria like \code{auc} and \code{cindex} we advise to use repeated CV (see \code{\link{CVfolds}}) to avoid ending up in any of the many local optima.
}
\value{
List, with components:
\item{optres}{Output of the optimizer}
\item{optpen}{Vector with determined optimal penalties}
\item{allsc}{Matrix with CV scores for all penalty parameter configurations used by the optimizer}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optLambdasWrap}} for i) (recommended) optimization in two steps: first global, then local; and ii) sequential optimization
when some data types are preferred over others. \code{\link{fastCV}} for initialization of penalties.
}
\examples{

}
