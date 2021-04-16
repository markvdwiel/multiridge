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
fastCV2(XXblocks, Y, X1 = NULL, kfold = 10, intercept = ifelse(class(Y) ==
"Surv", FALSE, TRUE), parallel   = FALSE, fixedfolds = TRUE, model = NULL,
eps = 1e-10, reltol = 0.5, lambdamax= 10^6, traceCV=TRUE)
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
 \item{reltol}{
Scalar. Relative tolerance for optimization method.
}
 \item{lambdamax}{
Numeric. Upperbound for lambda.
}

 \item{traceCV}{
Boolean. Should the CV results be traced and printed?
}

}
\details{
This function is basically a wrapper for applying \code{\link{optLambdas}} per data block separately using Brent optimization.
}
\value{
Numerical vector containing penalties optimized separately per data block. Useful for initialization.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optLambdas}}, \code{\link{optLambdasWrap}} which optimize the penalties jointly.
A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}}

\examples{
data(dataXXmirmeth)
resp <- dataXXmirmeth[[1]]
XXmirmeth <- dataXXmirmeth[[2]]

cvperblock2 <- fastCV2(XXblocks=XXmirmeth,Y=resp,kfold=10,fixedfolds = TRUE)
lambdas <- cvperblock2$lambdas
}
