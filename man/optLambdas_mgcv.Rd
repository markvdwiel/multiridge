\name{optLambdas_mgcv}
\alias{optLambdas_mgcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal ridge penalties with maximimum marginal likelihood
}
\description{
Optimizes a marginal likelihood score w.r.t. ridge penalties for multiple data blocks.
}
\usage{
optLambdas_mgcv(penaltiesinit=NULL, XXblocks,Y, pairing=NULL, model=NULL, reltol=1e-4,
  optmethod=ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),maxItropt=500,
  tracescore=TRUE, fixedpen=NULL, fixedseed =TRUE, sigmasq = 1,
  opt.sigma=ifelse(model=="linear",TRUE, FALSE))
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
  \item{optmethod}{
Character. Optimization method. Any of the methods \code{c("Brent", "Nelder-Mead", "Sann")} may be used, but
\code{"Nelder-Mead"} is generally recommended. Other unconstrained methods offered by \code{\link{optim}} may also be used, but have not been tested.
}
  \item{maxItropt}{
Integer. Maximum number of iterations for \code{optmethod}.
}
  \item{tracescore}{
Boolean. Should the output of the scores be traced?
}
\item{fixedpen}{
Integer vector or \code{NULL}. Contains indices of data types of which penalty
is fixed to the corresponding value in \code{penaltiesinit}.
}
  \item{fixedseed}{
Boolean. Should the initialization be fixed? For reproducibility.
}

\item{sigmasq}{Default error variance.
}

\item{opt.sigma}{Boolean. Should the error variance be optimized as well? Only relevant for \code{model="linear"}.
}

}
\details{
See \code{\link[mgcv]{gam}} for details on how the marginal likelihood is computed.

}
\value{
List, with components:
\item{optres}{Output of the optimizer}
\item{optpen}{Vector with determined optimal penalties}
\item{allsc}{Matrix with marginal likelihood scores for all penalty parameter configurations used by the optimizer}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optLambdas_mgcvWrap}} for i) (recommended) optimization in two steps: first global, then local; and ii) sequential optimization
when some data types are preferred over others. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
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

# Compute cross-validated score for initial lambdas
CVscore(penalties=lambdas, XXblocks=XXmirmeth,Y=resp,folds=leftout,
score="loglik")

# Optimize by using marginal likelihood criterion
jointlambdas2 <- optLambdas_mgcvWrap(penaltiesinit=lambdas, XXblocks=XXmirmeth,
Y=resp)

# Optimal lambdas
optlambdas <- jointlambdas2$optpen

}

