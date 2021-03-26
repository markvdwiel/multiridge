\name{fastCV2}
\alias{fastCV2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast cross-validation per data block
}
\description{
Fast cross-validation for high-dimensional data. Finds optimal penalties separately per data block. Useful for initialization.
}
\usage{
fastCV2(XXblocks, Y, X1 = NULL, kfold = 10, intercept = ifelse(class(Y) == "Surv", FALSE, TRUE), parallel   = FALSE, fixedfolds = TRUE, model = NULL, eps = 1e-10, lambdamax= 10^6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{XXblocks}{
List of data frames or matrices, representing \code{b=1,...,B} data blocks of dimensions \code{n x p_b}.
}
\item{Y}{
Response vector: numeric, binary, factor or \code{survival}.
}
   \item{X1}{
Matrix. Dimension \code{n x p_0, p_0 < n}, representing unpenalized covariates.
}
  \item{kfold}{
Integer. Desired fold.}

   \item{intercept}{
Boolean. Should an intercept be included?
}
  \item{parallel}{
Boolean. Should computation be done in parallel? If \code{TRUE}, requires to run \code{\link{setupParallel}} first.
}
 \item{fixedfolds}{
Boolean. Should fixed splits be used for reproducibility?
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
}
  \item{eps}{
Scalar. Numerical bound for IWLS convergence.
}

 \item{lambdamax}{
Numeric. Upperbound for lambda.
}

}
\details{
This function is basically a wrapper for applying \code{\link{optLambdas}} per data block separately using Brent optimization.
Unlike \code{\link{fastCV}} it does not require loading the (orphaned) penalized package.
}
\value{
Numerical vector containing penalties optimized separately per data block. Useful for initialization.
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
\code{\link{fastCV}}, \code{\link{optLambdas}}, \code{\link{optLambdasWrap}} which optimize the penalties jointly.}
\examples{

}
