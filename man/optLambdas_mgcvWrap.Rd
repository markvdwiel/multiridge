\name{optLambdas_mgcvWrap}
\alias{optLambdas_mgcvWrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal ridge penalties with sequential optimization.
}
\description{
Sequentially optimizes a marginal likelihood score w.r.t. ridge penalties for multiple data blocks.
}
\usage{
optLambdas_mgcvWrap(penaltiesinit=NULL, XXblocks,Y, pairing=NULL, model=NULL, reltol=1e-4,
  optmethod1= "SANN", optmethod2 =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),
  maxItropt1=10,maxItropt2=25,tracescore=TRUE,fixedseed =TRUE, pref=NULL, fixedpen=NULL,
  sigmasq = 1, opt.sigma=ifelse(model=="linear",TRUE, FALSE))
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
\item{pairing}{Numerical vector of length 3 or \code{NULL} when pairs are absent. Represents the indices (in \code{XXblocks}) of the two data blocks involved in pairing, plus the index of the paired block.
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
}
  \item{reltol}{
Scalar. Relative tolerance for optimization methods.
}
 \item{optmethod1}{
Character. First, global search method. Any of the methods \code{c("Brent", "Nelder-Mead", "Sann")} may be used, but
simulated annealing by \code{"Sann"} is recommended to search a wide landscape. Other unconstrained methods
offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{optmethod2}{
Character. Second, local search method. Any of the methods \code{c("Brent", "Nelder-Mead", "Sann")} may be used, but
\code{"Nelder-Mead"} is generally recommended. Other unconstrained methods
offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{maxItropt1}{
Integer. Maximum number of iterations for \code{optmethod1}.
}
  \item{maxItropt2}{
Integer. Maximum number of iterations for \code{optmethod2}.
}
 \item{tracescore}{
Boolean. Should the output of the scores be traced?
}
  \item{fixedseed}{
Boolean. Should the initialization be fixed? For reproducibility.
}
  \item{pref}{
Integer vector or \code{NULL}. Contains indices of data types in \code{XXblocks}
that are preferential.
}
\item{fixedpen}{
Integer vector or \code{NULL}. Contains indices of data types of which penalty
is fixed to the corresponding value in \code{penaltiesinit}.
}
\item{sigmasq}{Default error variance.
}

\item{opt.sigma}{Boolean. Should the error variance be optimized as well? Only relevant for \code{model="linear"}.
}


}
\details{
As opposed to \code{\link{optLambdas_mgcv}} this function first searches globally, then locally.
Hence, more time-consuming, but better guarded against multiple local optima.
See \code{\link[mgcv]{gam}} for details on how the marginal likelihood is computed.
}
\value{
List, with components:
\item{res}{Outputs of all optimizers used}
\item{lambdas}{List of penalties found by the  optimizers}
\item{optpen}{Numerical vector with final, optimal penalties}
}
\references{
Wood, S. N. (2011), Fast stable restricted maximum likelihood and marginal likelihood
estimation of semiparametric generalized linear models, J. Roy. Statist. Soc., B 73(1), 3-36.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optLambdas_mgcv}} for one-pass optimization. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}
