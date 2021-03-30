\name{predictIWLS}
\alias{predictIWLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions from ridge fits
}
\description{
Produces predictions from ridge fits for new data.
}
\usage{
predictIWLS(IWLSfit, X1new = NULL, Sigmanew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IWLSfit}{
List, containing fits from either \code{\link{IWLSridge}} (linear, logistic
ridge) or \code{\link{IWLSCoxridge}}
}
  \item{X1new}{
Matrix. Dimension \code{nnew x p_0}, representing unpenalized covariates for
new data.
}
  \item{Sigmanew}{
Matrix. Dimensions \code{nnew x n}. Sample cross-product from penalized
variables, usually computed by first applying \code{\link{createXXblocks}}
and then \code{\link{SigmaFromBlocks}}.
}
}
\details{
Predictions rely purely on the linear predictors, and do not require producing
the parameter vector.
}
\value{
Numerical vector of linear predictor for the test samples.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{IWLSridge}} (\code{\link{IWLSCoxridge}}) for fitting linear and
logistic ridge (Cox ridge). \code{\link{betasout}} for obtaining parameter
estimates.
\code{\link{Scoring}} to evaluate the predictions. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}
\examples{
#Example below shows how to create the input argument Sigmanew (for simulated data)
#Simulate
Xbl1 <- matrix(rnorm(1000),nrow=10)
Xbl2 <- matrix(rnorm(2000),nrow=10)
Xbl1new <- matrix(rnorm(200),nrow=2)
Xbl2new <- matrix(rnorm(400),nrow=2)

#check whether dimensions are correct
nrow(Xbl1)==nrow(Xbl1new)
nrow(Xbl2)==nrow(Xbl2new)
ncol(Xbl1)==nrow(Xbl2)
ncol(Xbl1new)==ncol(Xbl2new)

#create cross-product
XXbl <- createXXblocks(list(Xbl1,Xbl2),list(Xbl1new,Xbl2new))

#suppose penalties for two data types equal 5,10, respectively
Sigmanew <- SigmaFromBlocks(XXbl,c(5,10))

#check dimensions (should be nnew x n)
dim(Sigmanew)
}
