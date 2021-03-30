\name{betasout}
\alias{betasout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coefficient estimates from (converged) IWLS fit
}
\description{
Extracts estimated regression coefficients from the final Iterative Weighted Least Squares fit, as obtained from linear, logistic, or Cox ridge regression.
}
\usage{
betasout(IWLSfit, Xblocks, penalties, pairing = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IWLSfit}{
List object, see details
}
  \item{Xblocks}{
List of data frames or matrices, representing \code{b=1,...,B} data blocks of dimensions \code{n x p_b}.
}
  \item{penalties}{
Numerical vector.
}
  \item{pairing}{
Numerical vector of length 3 or \code{NULL}.
}
}
\details{
\code{IWLSfit} should be the output of either \code{\link{IWLSridge}} or \code{\link{IWLSCoxridge}}. \code{Xblocks} may be created by \code{\link{createXblocks}}.

}
\value{
List. Number of components equals number of components of \code{Xblocks} plus one, as the output is augmented with an intercept estimate (first component, \code{NULL} if absent).
Each component is a numerical vector representing regression parameter estimates. Lengths of vectors match column dimensions of \code{Xblocks} (nr of variables for given data type)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createXblocks}}
}

