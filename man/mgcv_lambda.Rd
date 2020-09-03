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
mgcv_lambda(penalties, XXblocks,Y, model, printscore=TRUE))
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
Character. Any of \code{c("linear", "logistic", "cox")}.
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
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CVscore}} for cross-validation alternative
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
