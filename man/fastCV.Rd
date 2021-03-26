\name{fastCV}
\alias{fastCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast, SVD-based cross-validation per data block
}
\description{
Fast, SVD-based cross-validation for high-dimensional data. Finds optimal penalties separately per data block. Useful for initialization.
}
\usage{
fastCV(Xblocks, Y, X1 = NULL, kfold = 10, intercept = ifelse(class(Y) == "Surv", FALSE, TRUE), parallel   = FALSE, fixedfolds = TRUE, model = NULL, eps = 1e-10, lambdamax= 10^6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Xblocks}{
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
Note that this function requires loading the orphaned \code{penalized} package. 
This function is kept for reproducibility of results. 
Use \code{\link{fastCV2}} if you do not have access to \code{penalized}. 
}
\value{
Numerical vector containing penalties optimized separately per data block. Useful for initialization.
This function uses SVD for fast cross-validation.
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
\code{\link{optLambdas}}, \code{\link{optLambdasWrap}} which optimize the penalties jointly. \code{\link{fastCV2}}, which does not depend on the penalized package}
\examples{

}
