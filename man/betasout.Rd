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
betasout(IWLSfit, Xblocks, X1=NULL, penalties, pairing = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IWLSfit}{
List object, see details
}
  \item{Xblocks}{
List of data frames or matrices, representing \code{b=1,...,B} data blocks of dimensions \code{n x p_b}.
}
   \item{X1}{
Matrix. Dimension \code{n x p_0, p_0 < n}, representing unpenalized covariates.
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
\code{\link{createXblocks}}. A full demo and data are available from:\cr
\url{https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4}
}

\examples{
data(dataXXmirmeth)
resp <- dataXXmirmeth[[1]]
XXmirmeth <- dataXXmirmeth[[2]]
lambdas <- c(100,1000)

# Prepare fitting for the specified penalties.
XXT <- SigmaFromBlocks(XXmirmeth,penalties=lambdas)

# Fit. fit$etas contains the n linear predictors
fit <- IWLSridge(XXT,Y=resp)

# Computation of the regression coefficients requires the original
# (large!) nxp data sets, available from link above
\dontrun{
Xbl <- createXblocks(list(datamir,datameth))
betas <- betasout(fit, Xblocks=Xbl, penalties=lambdas)
}
}
