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
optLambdas_mgcv(penaltiesinit=NULL, XXblocks,Y, model="logistic", reltol=1e-4, optmethod=ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"), maxItropt=500,tracescore=TRUE,fixedseed =TRUE)
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
  \item{model}{
Character. Any of \code{c("linear", "logistic", "cox")}.
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
  \item{fixedseed}{
Boolean. Should the initialization be fixed? For reproducibility.
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
\code{\link{optLambdas_mgcvWrap}} for i) (recommended) optimization in two steps: first global, then local; and ii) sequential optimization
when some data types are preferred over others.
}
\examples{

}
