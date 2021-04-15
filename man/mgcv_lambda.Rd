\name{mgcv_lambda}
\alias{mgcv_lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum marginal likelihood score
}
\description{
Computed maximum marginal likelihood score for given penalty parameters using \code{mgcv}.
}
\usage{
mgcv_lambda(penalties, XXblocks,Y, model=NULL, printscore=TRUE, pairing=NULL, sigmasq = 1,
  opt.sigma=ifelse(model=="linear",TRUE, FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{penalties}{
Numeric vector.
}
  \item{XXblocks}{
List of \code{nxn} matrices. Usually output of \code{\link{createXXblocks}}.
}
  \item{Y}{
Response vector: numeric, binary, factor or \code{survival}.
}
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}. Is inferred from
\code{Y} when \code{NULL}.
}

\item{printscore}{Boolean. Should the score be printed?
}

\item{pairing}{Numerical vector of length 3 or \code{NULL} when pairs are absent. Represents the indices (in \code{XXblocks}) of the two data blocks involved in pairing, plus the index of the paired block.
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
Numeric, marginal likelihood score for given penalties
}
\references{
Wood, S. N. (2011), Fast stable restricted maximum likelihood and marginal likelihood
estimation of semiparametric generalized linear models, J. Roy. Statist. Soc., B 73(1), 3-36.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CVscore}} for cross-validation alternative. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}
